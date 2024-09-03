library(MAST)
library(Seurat)
library(Matrix)
library(SingleCellExperiment)


#' @importMethodsFrom SummarizedExperiment assay
#' @importMethodsFrom SummarizedExperiment assayNames
#' @importMethodsFrom SummarizedExperiment 'assay<-'
NULL




#' Perform Differential Gene Expression using MAST 
#' Adapted from https://www.sc-best-practices.org/conditions/differential_gene_expression.html#single-cell-specific
#'
#' @param sca SinglecellAssay object over which to perform the DGE
#' @param formulaToUse DE Formula to apply
#' @param conditionOfInterest Contrast for which the p-value and fold change should be computed
#' @param minProportionOfCell Minimum global proportion of the cells expression a given gene.  Genes expressed in less cells than that are not tested for DGE.
#' @param libSizeVar Variable name representing the number of genes expressed per cells (as used in the formulaToUse)
#' @param filterFDRThr Filter output results based on this FDR threshold
#' @param mixed Set to true if the formula represents a Mixed model.  Set to false otherwise.
#' @returns A data frame with DGE results.
#' @import MAST
#' @import Seurat
#' @import Matrix
#' @import SingleCellExperiment
#' @export
find_de_MAST_RE <- function(sca,formulaToUse,conditionOfInterest,minProportionOfCell = 0.1,libSizeVar="ngeneson", filterFDRThr=0.05,mixed = T){

    print("Dimensions before subsetting:")
    print(dim(sca))
    print("")
    # keep genes that are expressed in more than 10% of all cells
    sca <- sca[freq(sca)>minProportionOfCell,]
    print("Dimensions after subsetting:")
    print(dim(sca))
    print("")
    # add a column to the data which contains scaled number of genes that are expressed in each cell
    cdr2 <- Matrix::colSums(c(sca)>0)
    colData(sca)[[libSizeVar]] <- scale(cdr2)

    # define and fit the model (if mixed or only fixed)
    if (mixed)
    {
        zlmCond <- MAST::zlm(formula = as.formula(formulaToUse), 
                       sca=sca, 
                       method='glmer', 
                       ebayes=F, 
                       strictConvergence=F,
                       fitArgsD=list(nAGQ = 0)) # to speed up calculations
    
    } else {
        zlmCond <- MAST::zlm(formula = as.formula(formulaToUse), 
                       sca=sca, 
                       ebayes=T, silent=T) 
    }
    
    # perform likelihood-ratio test for the condition that we are interested in    
    summaryCond <- summary(zlmCond, doLRT=conditionOfInterest,parallel = TRUE)
    # get the table with log-fold changes and p-values
    summaryDt <- summaryCond$datatable
    result <- merge(summaryDt[contrast==conditionOfInterest & component=='H',.(primerid, `Pr(>Chisq)`)], # p-values
                     summaryDt[contrast==conditionOfInterest & component=='logFC', .(primerid, coef)],
                     by='primerid') # logFC coefficients
    # MAST uses natural logarithm so we convert the coefficients to log2 base to be comparable to edgeR
    result[,coef:=result[,coef]/log(2)]
    # do multiple testing correction
    result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    result = result[result$FDR<filterFDRThr,, drop=F]

    result <- stats::na.omit(as.data.frame(result))
    return(result)
}


#' Compute Fold Changes  between 2 conditions and multiple samples using a pseudobulk approach.  First, the assay is aggregated by sample and the data is normalized using a scaling factor corresponding to the average library depth of all the samples.  Then the FoldChange is compute over the data slot of the pseudobulk aggregated assay using Seurat standard "FoldChange" function.
#'
#' @param seurat_object Seurat Object, with the condition of interest setup as "Ident"
#' @param ident.1 Ident level of the First group.  Fold change is compute as (ident.1/ident.2)
#' @param ident.2 Ident level of the Second group.  Fold change is compute as (ident.1/ident.2)
#' @param assay Assay name over which to compute fold changes.
#' @param sampleId Meta data field name containing the Id of the individual samples.
#' @returns A Fold change dataframe, with the assay feature names as rownames.
#' @import MAST
#' @import Seurat
#' @import Matrix
#' @import SingleCellExperiment
#' @examples
#'\dontrun{
#' fcdf = ComputeFCPseudobulk(seurat_object,'condition1','condition2',assay='SCT',sampleId='orig.ident')
#' head(fcdf)
#' }
#' @export
ComputeFCPseudobulk <- function(seurat_object, 
                                ident.1,
                                ident.2,
                                assay = "RNA",
                                sampleId = "SampleId"
                                )
{
    
    pb <- AggregateExpression(seurat_object, 
                             assays = assay, 
                             return.seurat = T, 
                             group.by = unique(c(sampleId,"ident")))
    #Here, a scale factor representing the mean sequencing depth of all samples is used
    # Using a proper scalefactor is important because of the +1 pseudocount added to fold change computation
    pb <- NormalizeData(pb, scale.factor = mean(Matrix::colSums(GetAssayData(pb[[assay]],slot = "counts") )))
    Idents(pb) <- "orig.ident"
    fcpb<- FoldChange(pb,ident.1 = ident.1,ident.2 = ident.2)
    fc <- FoldChange(seurat_object,ident.1 = ident.1,ident.2 = ident.2)
    fcpb$pct.1 <- fc[rownames(fcpb),"pct.1"]
    fcpb$pct.2 <- fc[rownames(fcpb),"pct.2"]
    return (fcpb)

}

#' Perform Differential Gene Expression using MAST using a mixed effect model containing a random effect for individual biological replicates.  The fold changes are however computed using a pseudo-bulk approach (see ComputeFCPseudobulk)
#' 
#' @param seurat_object Seurat Object, with the condition of interest setup as "Ident"
#' @param ident.1 Ident level of the First group.  Fold change is compute as (ident.1/ident.2)
#' @param ident.2 Ident level of the Second group.  Fold change is compute as (ident.1/ident.2)
#' @param formula DE Formula to apply
#' @param sampleId Meta data field name containing the Id of the individual samples.  This variable will be added to the formula as a random effect.
#' @param assay Assay name over which to compute fold changes.
#' @param identVar Name of the Ident var within the provided formula
#' @param libSizeVar Variable name representing the number of genes expressed per cells (as used in the formulaToUse)
#' @param minProportionOfCell Minimum global proportion of the cells expression a given gene.  Genes expressed in less cells than that in one of the 2 groups are not tested for DGE.
#' @param filterFDRThr Filter output results based on this FDR threshold
#' @returns A data frame with DGE results.
#' @import MAST
#' @import Seurat
#' @import Matrix
#' @import SingleCellExperiment
#' @export
RunDGE_MAST_MixedModel <- function(seurat_object,
                                   ident.1,
                                   ident.2,
                                   formula = "~ ngeneson + group",
                                   sampleId="sampleId",
                                   assay="RNA",
                                   identVar = "group",
                                   librarysize_var = "ngeneson",
                                   minProportionOfCell = 0.05,
                                   filterFDRThr=2)
{
    
        cells <- which(Idents(seurat_object) %in% c(ident.1,ident.2))
        
        label <- factor(Idents(seurat_object)[cells])

        pbFC <- ComputeFCPseudobulk(seurat_object,
                                   ident.1,
                                   ident.2,
                                   assay = assay,
                                   sampleId = sampleId)
        
        sca = MAST::SceToSingleCellAssay(as.SingleCellExperiment(seurat_object[rownames(pbFC)[pbFC$pct.1 >= minProportionOfCell |
                                                                                        pbFC$pct.2 >= minProportionOfCell],
                                                                         cells],assay =assay), 
                                   class = "SingleCellAssay",check_sanity=FALSE)
 
        

        #Set the base group
        label <- relevel(label,ident.2)
    
        colData(sca)[[sampleId]] <- factor(colData(sca)[[sampleId]])

        colData(sca)[[identVar]] <- label
        #Our condition of interest (contrast) is on the group variable 

        
        contrast = paste0(gsub(" ", "_",identVar),ident.1)

        mixedFormula = paste0(formula," + (1|",sampleId,")")
       

        #We control for the library size (ngeneson),
        #First run with the biological replicate as a fixed effect
        res = find_de_MAST_RE(sca,
                              formula = mixedFormula,
                              conditionOfInterest = contrast, 
                              mixed=TRUE,
                              minProportionOfCell = 0,
                              libSizeVar=librarysize_var, 
                              filterFDRThr=filterFDRThr)

        

        colnames(res) <- c('gene','p.value','log2FC','fdr')
        res$log2FC <- pbFC[res$gene,"avg_log2FC"]
        res$pct.1 <- pbFC[res$gene,"pct.1"]
        res$pct.2 <- pbFC[res$gene,"pct.2"]
        return(res)
}

#' Perform Differential Gene Expression using MAST using a GLM model c The fold changes are however computed using a pseudo-bulk approach (see ComputeFCPseudobulk)
#' 
#' @param seurat_object Seurat Object, with the condition of interest setup as "Ident"
#' @param ident.1 Ident level of the First group.  Fold change is compute as (ident.1/ident.2)
#' @param ident.2 Ident level of the Second group.  Fold change is compute as (ident.1/ident.2)
#' @param formula DE Formula to apply
#' @param sampleId Meta data field name containing the Id of the individual samples.  
#' @param assay Assay name over which to compute fold changes.
#' @param identVar Name of the Ident var within the provided formula
#' @param libSizeVar Variable name representing the number of genes expressed per cells (as used in the formulaToUse)
#' @param minProportionOfCell Minimum global proportion of the cells expression a given gene.  Genes expressed in less cells than that in one of the 2 groups are not tested for DGE.
#' @param filterFDRThr Filter output results based on this FDR threshold
#' @returns A data frame with DGE results.
#' @import MAST
#' @import Seurat
#' @import Matrix
#' @import SingleCellExperiment
#' @export
RunDGE_MAST_GLMModel <- function(seurat_object,
                                   ident.1,
                                   ident.2,
                                   formula = "~ ngeneson + group",
                                   sampleId="sampleId",
                                   assay="RNA",
                                   identVar = "group",
                                   librarysize_var = "ngeneson",
                                   minProportionOfCell = 0.05,
                                   filterFDRThr=2)
{
    
        cells <- which(Idents(seurat_object) %in% c(ident.1,ident.2))
        
        label <- factor(Idents(seurat_object)[cells])

        pbFC <- ComputeFCPseudobulk(seurat_object,
                                   ident.1,
                                   ident.2,
                                   assay = assay,
                                   sampleId = sampleId)
        
        sca = MAST::SceToSingleCellAssay(as.SingleCellExperiment(seurat_object[rownames(pbFC)[pbFC$pct.1 >= minProportionOfCell |
                                                                                        pbFC$pct.2 >= minProportionOfCell],
                                                                         cells],assay =assay), 
                                   class = "SingleCellAssay",check_sanity=FALSE)
 
        

        #Set the base group
        label <- relevel(label,ident.2)
    
       colData(sca)[[identVar]] <- label
        #Our condition of interest (contrast) is on the group variable 

        
        contrast = paste0(gsub(" ", "_",identVar),ident.1)


        #We control for the library size (ngeneson),
        #First run with the biological replicate as a fixed effect
        res = find_de_MAST_RE(sca,
                              formula = formula,
                              conditionOfInterest = contrast, 
                              mixed=False,
                              minProportionOfCell = 0,
                              libSizeVar=librarysize_var, 
                              filterFDRThr=filterFDRThr)

        

        colnames(res) <- c('gene','p.value','log2FC','fdr')
        res$log2FC <- pbFC[res$gene,"avg_log2FC"]
        res$pct.1 <- pbFC[res$gene,"pct.1"]
        res$pct.2 <- pbFC[res$gene,"pct.2"]
        return(res)
}