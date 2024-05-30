require(fgsea)
require(Seurat)
require(future)
require(BiocParallel)
require(AUCell)
require(Matrix)
require(qvalue)


#' Load Genesets files in GMT format from a folder
#'
#' @param folder Path where the GMT files are stored (directory)
#' @returns A list of list of geneset, the top level list is a named list of GMT files and the lower level list is a named list of the genesets within this file.
#' @examples
#'\dontrun{
#' genesets = LoadGeneSetsGMTs("../genesets")
#' names(genesets)
#' genesets[[1]]
#' }
#' @export
LoadGeneSetsGMTs <- function(folder = ".")
{
    ## Loading Genesets
    genesetFiles <- list.files(folder, pattern='gmt')
    genesets = list()
    for (f in genesetFiles)
    {
        print(paste0("Loading ",f,"..."))
        g = fgsea::gmtPathways(file.path(folder,f))
        print(paste0("Loaded ", length(g), " genesets"))
        genesetType <- gsub(pattern = "\\.gmt$", "", f)
        genesets[[genesetType]] <- g
    }
    print(paste0("Total number of genesets files loaded: ",length(genesets)))

    return(genesets)

}

#' Collapse a list of list of genesets into a single level list
#'
#' @param genesetsLists List of list of genesets, as returned by LoadGeneSetsGMTs
#' @param usePrefix If TRUE, the top level list named a prefixed to the collapsed list names.  Default to FALSE to preserve the original geneset names.
#' @returns A single-level list genesets, ready to be used for pathway activity scoring.
#' @examples
#'\dontrun{
#' genesetsLists = LoadGeneSetsGMTs("../genesets")
#' genesets <- CollapseGenesetList(genesetsLists)
#' genesets[[1]]
#' }
#' @export
CollapseGenesetList <- function(genesetsLists,usePrefix=FALSE)
{
    collapsed_gs <-  unlist(genesetsLists, recursive=FALSE,use.names = TRUE)
    if(!usePrefix)
    {
        for(g in names(genesetsLists))
        {
            names(collapsed_gs) <- gsub(paste0('^',g,'.'),'',names(collapsed_gs))
        }
    }
    return(collapsed_gs)

}

#' Create a Ranking of Genes from DE results.  The ranking is a based on a multiplication of the log2FC and the -log10(pvalue)
#'
#' @param deResultDF Data Frame containing the Differential expression results
#' @param pvalColumn Column name containing the p.value
#' @param log2FCColumn Column name containing the log2FC
#' @param geneColumn Column name containing the gene name
#' @param twoSided Is the ranking 2 sided (default) or not.  If not 2-sided, the absoluted value of the 2-sided ranking is returned.
#' @param groupByColumn If the input deResultDF contains multiple sets of results, the group by column should be specified to identify each unique resultset in the dataframe.  A ranking will be defined for each resultset.
#' @returns A list of rankings, one per resultset defined by the groupByColumn.  If the groupByColumn column is NULL, a single ranking is in the returned list and is named "All"
#' @examples
#'\dontrun{
#' rankings = CreateRanking(myDEResultsDF,pvalColumn='p_val',log2FCColumn="avg_Log2FC",geneColumn='gene',groupByColumn="cluster")
#' names(rankings)
#' rankings[[1]]
#' }
#' @export
CreateRanking <- function(deResultDF, 
              pvalColumn="p.value", 
              log2FCColumn="log2FC",
              geneColumn="gene",
              twoSided=TRUE,
              groupByColumn=NULL)
{

    CreateRankingInternal <- function(df)
    {
        pval = df[,pvalColumn]
        log2FC = df[,log2FCColumn]
        pval[pval==0] = 1e-300
        FGSEA_Ranking = log2FC * -log10(pval)
        if(!twoSided)
        {
            FGSEA_Ranking = abs(FGSEA_Ranking)
        }
        if(anyDuplicated(df[,geneColumn]))
        {
            stop("Genes must be unique within each group.  Perhaps you should be specifying a groupByColumn???")
        }

        names(FGSEA_Ranking) = df[,geneColumn]
        return(FGSEA_Ranking)
    }
    results = list()
    if(is.null(groupByColumn))
    {
        results['All'] = CreateRankingInternal(deResultDF)
    } else {
        groups = unique(deResultDF[,groupByColumn])
        for(g in groups)
        {
            results[[g]] <- CreateRankingInternal(deResultDF[deResultDF$groupByColumn == g,])
        }
    }
    return(results)
}


#' Run FGSEAMultilevel in batch mode, for each rankedGeneLists and for each type of Genesets
#'
#' @param genesetLists A lists of list of genesets, as returned by LoadGeneSetsGMTs
#' @param rankedGeneLists A list of ranked gene lists, as returned by CreateRanking
#' @param removeGenePattern Optional pattern to remove unwanted genes from the ranking before running FGSEA.  Could be used to remove Ribosomal or mitochondrial genes, for instance.
#' @param geneListsGroupedBy Column name in the output data frame used to identify the different ranked gene lists
#' @param ... Additional arguments passed on to fgseaMultilevel
#' @returns A dataframe containing the FGSEA results.  The leadingEdge column is collapsed into a string delimited by ; 
#' @examples
#'\dontrun{
#' genesets = LoadGeneSetsGMTs("../genesets")
#' rankings = CreateRanking(myDEResultsDF,pvalColumn='p_val',log2FCColumn="avg_Log2FC",geneColumn='gene',groupByColumn="cluster")
#' fgseaResults <- RunFGSEAMultilevelBatch(genesets,rankings,removeGenePattern="RPS|RPL|MT-")
#' head(fgseaResults)
#' }
#' @export
RunFGSEAMultilevelBatch<- function(genesetLists,
              rankedGeneLists,
              removeGenePattern = NULL,
              geneListsGroupedBy="RankedGeneList",
              ...)
{
    results = NULL
    for(g in names(genesetLists))
    {
        genesets = genesetLists[[g]]
        for(group in names(rankedGeneLists))
        {
            print(paste0("Running ",group,"; on ", g))
            ranks = rankedGeneLists[[group]]

            ranks = na.omit(ranks) # Just in case

            if(!is.null(removeGenePattern))
            {
                #Remove genes that match the pattern (ribosomal or mitochondrial genes, for instance)
                ranks <- ranks[grep(removeGenePattern,names(ranks),invert=TRUE,value = TRUE)]
            }
            fgseaRes <- fgsea::fgseaMultilevel(genesets,ranks, ...)
            fgseaRes[,geneListsGroupedBy]  = group
            fgseaRes$genesetType = g
            if(is.null(results))
            {
                results = fgseaRes
            } else {
                results = rbind(results,fgseaRes)
            }            
        }
    }

    #Collapse the leadingEdge column into a string
    results$leadingEdge <- vapply(results$leadingEdge, paste, collapse = ";", character(1L))

    return(results)
}


#' Split data matrix into smaller sub-matrices ('chunks')
#' Note : Code copied from https://rdrr.io/github/carmonalab/UCell/man/split_data.matrix.html
#' 
#'
#' @param   matrix      Input data matrix 
#' @param   chunk.size  How many cells to include in each sub-matrix
#' 
#' @return  A list of sub-matrices, each with size {n_features x chunk_size}
split_data.matrix <- function(matrix, chunk.size=100) {
    ncols <- dim(matrix)[2]
    nchunks <- (ncols-1) %/% chunk.size + 1
    
    split.data <- list()
    min <- 1
    for (i in seq_len(nchunks)) {
        if (i == nchunks-1) {  #make last two chunks of equal size
            left <- ncols-(i-1)*chunk.size
            max <- min+round(left/2)-1
        } else {
            max <- min(i*chunk.size, ncols)
        }
        split.data[[i]] <- matrix[,min:max,drop=FALSE]
        min <- max+1    #for next chunk
    }
    return(split.data)
}

#' Run AUCell Pathway Activity scoring and threshold the results based on a null distribution obtained from a permuted gene matrix.
#'
#' @param dataMat A sparse scRNAseq gene expression matrix  (genes * cells)
#' @param genesets A list of genesets to compute the activity score
#' @param p_val_thr Pvalue threshold to use to threshold the results.  Pathway activity score values less significant than this p-value will be set to 0.
#' @param minSize Only use genesets having at least that many detected genes in the dataMatrix
#' @param ncores Use that number of cores when running AUCell
#' @param ... Additional arguments passed on to UCell::AUCell_run
#' @returns A sparse pathway activity matrix (pathway * cells), which has been thresholded so that non significant scores are 0.
#' @examples
#'\dontrun{
#' genesets = CollapseGenesetList(LoadGeneSetsGMTs("../genesets"))
#' geneExpMat = as(GetAssayData(object = seurat_object, assay = 'RNA', layer = "data"),"dgCMatrix")
#' pasMat = ComputeAUCellPathwayActivity(geneExpMat,genesets,cores=5)
#' seurat_object[['PAS']] <- CreateAssayObject(data = pasMat)
#' }
#' @export
ComputeAUCellPathwayActivity <- function(dataMat,
                              genesets,
                              p_val_thr = 0.05,
                              minSize = 10,                               
                              ncores = 1,
                              ...)
{


    detectedGeneCounts <- Matrix::rowSums(as(dataMat,"dgCMatrix") > 0)
    detectedGeneCounts <- detectedGeneCounts[detectedGeneCounts >0]
    detectedGeneCounts <- sapply(genesets,FUN=function(x){sum(x %in% names(detectedGeneCounts))})
    #Filter genesets that have too few detected genes in our dataset
    genesets <- genesets[detectedGeneCounts > minSize]
    
    
    
    #Adjust the thrshold if we do not have enough cells.
    if (p_val_thr <= 1/NCOL(dataMat))
    {
        p_val_thr = 1/NCOL(dataMat)
    }
    
    split.data <- split_data.matrix(matrix=dataMat, chunk.size=2000)
    

    #work in chunks of data, so that the dense matrix created by the apply is not too large.
    shuffledMat <- NULL
    for(i in split.data)
    {
         #Shuffle the gene expression independently for each cell
        mat <- as(apply(i,MARGIN=2,FUN=sample),"dgCMatrix")
        if(is.null(shuffledMat))
        {
            shuffledMat <- mat
        } else {
            shuffledMat <- cbind(shuffledMat,mat)
        }
    }

    #Put back the row names
    rownames(shuffledMat) <- rownames(seurat_object[[assay_name]])

    # Run AUCells on both the shuffled matrix, and the normal matrix
    shuffledMat_AUC <- AUCell::AUCell_run(shuffledMat,
                        genesets, BPPARAM = MulticoreParam(ncores), ...)
    AUC <- AUCell::AUCell_run(dataMat,
                        genesets, BPPARAM = MulticoreParam(ncores), ...)

    
    split.data <- split(rownames(AUC@assays@data$AUC),2000)
    res <- NULL
    
    #Apply Threshold in batches, because a dense matrix is created at each batch

    for(i in split.data)
    {
        mat <- as(t(sapply(i,
                           FUN=function(x){
                               #Compute empirical p-values based on the 
                               pvals <- qvalue::empPvals(AUC@assays@data$AUC[x,],
                                                         shuffledMat_AUC@assays@data$AUC[x,], 
                                                         pool = TRUE)
                               pas <- AUC@assays@data$AUC[x,]
                               qs <- quantile(shuffledMat_AUC@assays@data$AUC[x,],prob=c(1-p_val_thr))
                               pas[pvals > p_val_thr] <- 0

                               pas
                               })),"dgCMatrix")

        if(is.null(res))
        {
            res <- mat
        } else {
            res <- rbind(res,mat)
        }
    }    
    
    colnames(res) <- colnames(shuffledMat_AUC@assays@data$AUC)
    rownames(res) <- rownames(shuffledMat_AUC@assays@data$AUC)

    return(res)

}