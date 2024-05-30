require(fgsea)

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