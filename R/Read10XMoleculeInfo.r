require(Matrix)
require(rhdf5)

#' Read the 10X Molecule_info.h5 file and outputs separated EXONS/INTRONS count matrices along with basic cell and feature metadata
#'
#' @param moleculeFile Path to the molecule_info.h5 file.  Defaults to the molecule_info.h5 file in the current directory.
#' @param featureType Feature type to extract.  Defaults to 'Gene Expression'
#' @param filtered True to return the Filtered matrix (only cells).  False to return the full matrix of all barcodes.  Defaults to True.
#' @returns A list with the following 4 elements : EXONS , INTRONS , cells.meta.data, features.meta.data
#' @examples
#'\dontrun{
#' counts = Read10XMoleculeInfo(moleculeFile="cellranger/sample/outs/molecule_info.h5")
#' counts$EXONS # The exon count matrix
#' counts$INTRONS # The intron count matrix
#' counts$cells.meta.data # data.frame of cells meta data
#' counts$features.meta.data # data.frame of features meta data
#' seurat_object <- CreateSeuratObject(counts$EXONS+counts$INTRONS,
#' meta.data=counts$cells.meta.data,
#' assay = "RNA")
#' seurat_object[['INTRONS']] <- CreateAssayObject(counts = counts$INTRONS )
#' seurat_object[['EXONS']] <- CreateAssayObject(counts = counts$EXONS )
#' seurat_object[['RNA']] <- AddMetaData(seurat_object[['INTRONS']], counts$features.meta.data)
#' seurat_object[['INTRONS']] <- AddMetaData(seurat_object[['INTRONS']], counts$features.meta.data)
#' seurat_object[['EXONS']] <- AddMetaData(seurat_object[['EXONS']], counts$features.meta.data)
#' }
#' @export
Read10XMoleculeInfo <- function(moleculeFile="molecule_info.h5",
                                featureType = 'Gene Expression',
                                filtered = TRUE)
{


    #Read the feature types, and determines which features to keep
    featureTypes = rhdf5::h5read(moleculeFile,"features/feature_type")
    featuresToKeep = which(featureTypes == featureType)
    filt = rhdf5::h5read(moleculeFile, "feature_idx") %in% (featuresToKeep-1)


    #Determines which barcodes are called cells
    barcodes = rhdf5::h5read(moleculeFile,"barcodes")
    barcode_info <- rhdf5::h5read(moleculeFile,"barcode_info/pass_filter")
    passFilterBarcodes =  barcode_info[1,barcode_info[2,]==0]

    #Reads Feature Id and Features Name
    featureName = rhdf5::h5read(moleculeFile,"features/name")
    featureId = rhdf5::h5read(moleculeFile,"features/id")


    name =  make.unique(featureName)
    if(filtered)
    {
        filt = filt & rhdf5::h5read(moleculeFile, "barcode_idx") %in% passFilterBarcodes
    } else {
        #keep all barcodes
        passFilterBarcodes = barcode_info[1,]
    }
    #Filter on Exon UMI
    filt_exon = rhdf5::h5read(moleculeFile,"umi_type") == 1

    #Size of output matrix
    ncells = NROW( barcodes)
    ngenes = NROW(name)

    #Extract Exons count matrix into a sparse matrix
    countExons <- Matrix::sparseMatrix(i = rhdf5::h5read(moleculeFile,"feature_idx")[filt & filt_exon]+1,
                 j = rhdf5::h5read(moleculeFile,"barcode_idx")[filt& filt_exon]+1,
                 x = rep(1, sum(filt& filt_exon)),
                dims = c(ngenes,ncells),
                 use.last.ij = FALSE,
                giveCsparse = TRUE)
    #Extract Introns count matrix matrix
    countIntrons <- Matrix::sparseMatrix(i = rhdf5::h5read(moleculeFile,"feature_idx")[filt & !filt_exon]+1,
                 j = rhdf5::h5read(moleculeFile,"barcode_idx")[filt& !filt_exon]+1,
                 x = rep(1, sum(filt& !filt_exon)),
                dims = c(ngenes,ncells),
                 use.last.ij = FALSE,
                giveCsparse = TRUE)

    #Filter sparse matrix to only the features to keep and the barcodes to keep.
    countExons <- countExons[featuresToKeep,passFilterBarcodes+1]
    countIntrons <- countIntrons[featuresToKeep,passFilterBarcodes+1]
    barcodesNames <-barcodes[passFilterBarcodes+1]
    if(nchar(barcodesNames[1]) == 16)
    {
        barcodesNames <- paste0(barcodesNames,"-1")
    }
    colnames(countExons) <- colnames(countIntrons) <- barcodesNames
    rownames(countExons) <- rownames(countIntrons) <- name[featuresToKeep]
    cells.meta.data <- data.frame(nCount_EXONS =Matrix::colSums(countExons),
                                  nCount_INTRONS = Matrix::colSums(countIntrons))
    cells.meta.data$pct_INTRONS = 100*cells.meta.data$nCount_INTRONS / pmax(1,cells.meta.data$nCount_INTRONS + cells.meta.data$nCount_EXONS)

    features.meta.data <- data.frame(feature_id = featureId[featuresToKeep],
                                  feature_name = featureName[featuresToKeep],
                                  genome = rhdf5::h5read(moleculeFile,"features/genome")[featuresToKeep],
                                  feature_type = featureTypes[featuresToKeep],
                                  nCount_EXONS =Matrix::rowSums(countExons),
                                  nCount_INTRONS = Matrix::rowSums(countIntrons))
    features.meta.data$pct_INTRONS = 100*features.meta.data$nCount_INTRONS / pmax(1,features.meta.data$nCount_INTRONS + features.meta.data$nCount_EXONS)

    return (list('EXONS'=countExons,
            'INTRONS'=countIntrons,
            'cells.meta.data'=cells.meta.data,
            'features.meta.data'=features.meta.data))
}
