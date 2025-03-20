library(openxlsx)
library(Seurat)
library(dplyr)
library(rols)

#' Export a list of tables to 1 excel file
#'
#' @param tableList Named List of tables to export.  The name will be used as the sheet name in the output Excel file
#' @param file Excel file name (and path)
#' @param colNames Export column names (TRUE / FALSE)
#' @param rowNames Export row names (TRUE / FALSE)
#' @param freezeFirstRow Set to TRUE to freeze the first row.
#' @param freezeFirstCol Set to TRUE to freeze the first column
#' @param tableStyle Style of the data table in the excel workbook
#' @examples
#'\dontrun{
#' ExportTables(tableList = List("table1" = df1,
#'                               "table2" = df2),
#'              file = "./path/to/export/file.xlsx")
#' }
#' @export
ExportTables <- function(tableList,
                         file,
                         colNames = TRUE, 
                         rowNames = FALSE,
                         freezeFirstRow = TRUE,
                         freezeFirstCol = TRUE,
                         tableStyle = "TableStyleLight9")
{
    wb <- createWorkbook()
    #Create one sheet for each table
    for(sheetName in names(tableList))
    {
        addWorksheet(wb, sheetName = sheetName, gridLines = FALSE)
        freezePane(wb, sheet = sheetName, firstRow = freezeFirstRow, firstCol = freezeFirstCol) ## freeze first row and column
        writeDataTable(wb, sheet =sheetName, x = tableList[[sheetName]],colNames = TRUE, rowNames = TRUE,tableStyle = tableStyle)
    }

    saveWorkbook(wb, file, overwrite = TRUE) 
}

#' Export a pseudobulk assay to an Excel workbook.  2 sheets will be created, one with the assay counts data (PseudoBulk_RawCounts) and another one with the normalized log2CPM data (PseudoBulk_log2CPMNorm)
#'
#' @param seurat_object Seurat Object, with the condition of interest setup as "Ident"
#' @param file Excel file name (and path)
#' @param assay Assay name to export
#' @param features Subset of features to export (or NULL to export all features)
#' @param group.by Metadata fields over which to compute the pseudobulk assay.  Set to "ident" by default.
#' @param tableStyle Style of the data table in the excel workbook
#' @returns a seurat object representing the pseudo_bulk assay
#' @examples
#'\dontrun{
#' pb_seurat <- ExportPseudoBulkAssay(seurat_object = seurat_object,
#'                                      assay = "RNA",
#'                                      group.by = c("SampleId", "Condition"),
#'                                      file = "./path/to/export/file.xlsx")
#' }
#' @export
ExportPseudoBulkAssay <- function(seurat_object,
                                  file,
                                  assay = "RNA",
                                  features = NULL,
                                  group.by= "ident",
                                  tableStyle = "TableStyleLight9")
{
    wb <- createWorkbook()
    agg <- AggregateExpression(
    seurat_object,
    assays = assay,
    features = features,
    return.seurat = TRUE,
    group.by = group.by,
    normalization.method = "LogNormalize",
    scale.factor = 1000000,
    margin = 1,
    verbose = TRUE)

    #Sccale to log2(CPM)
    agg <- NormalizeData(agg,scale.factor=1000000)

    addWorksheet(wb, sheetName = "PseudoBulk_RawCounts", gridLines = FALSE)
    freezePane(wb, sheet = 1, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
    writeDataTable(wb, sheet = 1, x = as.data.frame(agg[[assay]]$counts),colNames = TRUE, rowNames = TRUE,tableStyle = tableStyle)

    addWorksheet(wb, sheetName = "PseudoBulk_log2CPMNorm", gridLines = FALSE)
    freezePane(wb, sheet = 2, firstRow = TRUE, firstCol = TRUE) ## freeze first row and column
    writeDataTable(wb, sheet = 2, x = as.data.frame(agg[[assay]]$data),colNames = TRUE, rowNames = TRUE,tableStyle = tableStyle)

    saveWorkbook(wb, file, overwrite = TRUE) 
    return(agg)

}

#' Export a matrix to the mtx format accepted by the Broad single cell Portal
#'
#' @param seurat_object Seurat Object
#' @param folder Folder where to export the matrix
#' @param assay Assay name to export
#' @param layer Layer to export (data) by default
#' @returns nothing
#' @examples
#'\dontrun{
#' ExportMatrixSingleCellPortal(seurat_object = seurat_object,
#'                                      folder = "./processed"
#'                                      assay = "RNA",
#'                                      layer = "data")
#' }
#' @export
ExportMatrixSingleCellPortal <- function(seurat_object,
                                  folder = "./processed",
                                  assay = "RNA",
                                  layer = "data")
{

    DropletUtils::write10xCounts(folder,
                                 GetAssayData(seurat_object[[assay]],
                                              layer=layer),
                                overwrite=T,
                                version="3")
    #Broad portal doesnt like matrix.mtx.gz name
    file.rename(file.path(folder, "matrix.mtx.gz"), 
                file.path(folder, "matrix.txt.gz"))

}

#' Export cell metadata to a format accepted by the Broad single cell Portal
#'
#' @param seurat_object Seurat Object,
#' @param file meta data file name
#' @param columns Columns to export (default to all metadata columns)
#' @returns nothing
#' @examples
#'\dontrun{
#'  ExportMetadataSingleCellPortal(seurat_object = seurat_object,
#'                                      file = "./metadata.txt")
#' }
#' @export
ExportMetadataSingleCellPortal <- function(seurat_object,
                                  file = "./metadata.txt",
                                  columns = colnames(seurat_object@meta.data))
{

    #Export Metadata file
    types <- sapply(columns, FUN=function(x){
        if(typeof(seurat_object[[x]][[1]]) == "character")
        {
            return("group")
        } else {
            return("numeric")
        }
    })
    types <- append(list('NAME' = 'TYPE'),types)
    file1 <- tempfile(pattern = "metadata1", tmpdir = tempdir(), fileext = "")
    file2 <- tempfile(pattern = "metadata2", tmpdir = tempdir(), fileext = "")

    write.table(data.frame(types),file1,quote=F, sep = "\t",row.names=F)
    write.table(seurat_object@meta.data[,columns],file2,col.names=F,quote=F, sep = "\t")
    system(paste0("cat ",file1," ",file2," > ", file))

}

#' Export cell clustering to a format accepted by the Broad single cell Portal
#'
#' @param seurat_object Seurat Object,
#' @param file meta data file name
#' @param reduction the data reduction name to export (UMAP by default)
#' @param extraColumns Extra Columns to export alongside the UMAP
#' @returns nothing
#' @examples
#'\dontrun{
#'  ExportClusteringSingleCellPortal(seurat_object = seurat_object,
#'                                      file = "./metadata.txt",
#'                                      reduction = "UMAP")
#' }
#' @export
ExportClusteringSingleCellPortal <- function(seurat_object,
                                  file = "./clustering.txt",
                                  reduction = "UMAP",
                                  extraColumns = c("seurat_clusters"))
{

    reduction = as.data.frame(Embeddings(seurat_object, reduction = reduction))
    #Only keep first 2 dimensions
    reduction = reduction[,1:2]

    columns = extraColumns[extraColumns %in% colnames(seurat_object@meta.data)]
    types <- sapply(columns, FUN=function(x){
        if(typeof(seurat_object[[x]][[1]]) == "character")
        {
            return("group")
        } else {
            return("numeric")
        }
    })

    if(NROW(columns) > 0)
    {
        types <- append(list("NAME" = "TYPE","X"="numeric","Y"="numeric"),types)
        for(c in columns)
        {
            reduction[,c] = seurat_object@meta.data[,c]
        }
    } else {
        types <- list("NAME" = "TYPE","X"="numeric","Y"="numeric")
    }

    
    file1 <- tempfile(pattern = "metadata1", tmpdir = tempdir(), fileext = "")
    file2 <- tempfile(pattern = "metadata2", tmpdir = tempdir(), fileext = "")

    write.table(data.frame(types),file1,quote=F, sep = "\t",row.names=F)
    write.table(reduction,file2,col.names=F,quote=F, sep = "\t")
    system(paste0("cat ",file1," ",file2," > ", file))

}

#' Export coordinate label files to a format accepted by the Broad single cell Portal
#'
#' @param seurat_object Seurat Object,
#' @param file meta data file name
#' @param reduction The reduction to base the coordinate on
#' @param label The Group label column to use
#' @returns nothing
#' @examples
#'\dontrun{
#'  ExportLabelCoordinatesSingleCellPortal(seurat_object = seurat_object,
#'                                      file = "./labels.txt",
#'                                      reduction = "UMAP",
#'                                      label = "cellType"
#' )
#' }
#' @export
ExportLabelCoordinatesSingleCellPortal <- function(seurat_object,
                                  file = "./labels.txt",
                                  reduction = "UMAP",
                                  label = "seurat_clusters")
{

    reduction = as.data.frame(Embeddings(seurat_object, reduction = reduction))
    #Only keep first 2 dimensions
    reduction = reduction[,1:2]

    reduction[,"LABELS"] = seurat_object@meta.data[,label]
    colnames(reduction) <- c("X","Y","LABELS")
    labels <- reduction %>% group_by(LABELS) %>% summarize(X=mean(X),Y=mean(Y))
    labels <- labels[,c("X","Y","LABELS")]
    colnames(labels) <- c("X","Y","LABELS")
    write.table(labels,file,row.names=F,quote=F, sep = "\t")

}


#' Add the Required meta data field for uploading the dataset to the Broad Single Cell Portal.  See (see https://singlecell.zendesk.com/hc/en-us/articles/360060609852-Required-metadata).  The values will be added in a "best-effort" way, but should be manually reviewed afterwards.
#'
#' @param seurat_object Seurat Object,
#' @param maleCells id of male cells 
#' @param femaleCells id of female cells
#' @param mixedCells id of mixed cells
#' @param biosampleField Existing field containing the biosample identifier
#' @param donorField Existing field containing the donor identifier
#' @param species Text field identifying the specie.  This will be matched against the taxonomy ontology
#' @param normalCells id of normal (non-diseased) cells
#' @param diseaseName Name of the disease affecting the non-normal cells. This will be matched against the MONDO ontology
#' @param organName Name of the organ under study. This will be matched against the UBERON ontology
#' @param library_preparation_protocol Name of library preparation protocol. This will be matched against the EFO ontology 
#' @param cellTypeField Cell meta data field from the seurat_object containing the cell type names.  This will be matched against the cell ontology.
#' @returns a seurat object with the required metadata fields added (see https://singlecell.zendesk.com/hc/en-us/articles/360060609852-Required-metadata).  The use is expected to validate the results before uploading the data.
#' @examples
#'\dontrun{
#'  seurat_object <- AnnotateOntologiesSingleCellPortal(seurat_object,femaleCells=colnames(seurat_object),
#'                                   biosampleField="sample",
#'                                   donorField="condition",
#'                                  species = "mouse",
#'                                  normalCells = seurat_object$condition == "NORMAL",
#'                                  diseaseName = "AML",
#'                                  organName = "Blood",
#'                                  library_preparation_protocol = "Flex",
#'                                  cellTypeField = "ClusterNameLevel1"  )  
#' }
#' @export
AnnotateOntologiesSingleCellPortal <- function(seurat_object,
                                      maleCells = c(),
                                      femaleCells = c(),
                                      mixedCells = c(),
                                      biosampleField = "sample",
                                      donorField = "donor",
                                      species = "human",
                                      normalCells = c(),
                                      diseaseName = "AML",
                                      organName = "Blood",
                                      library_preparation_protocol = "10X 3' v3 Gene Expression",
                                      cellTypeField = "cellType"
                                      )
{

    OLSearch <- function(x, ontology = "cl",default='CL_0000255')
    {
        qry <- OlsSearch(q = x, ontology = ontology, exact = F, rows = 1)
        (qry <- olsSearch(qry))
        qdrf <- as(qry, "data.frame")
        if(NROW(qdrf) > 0)
        {
            return(qdrf$short_form[1])
        } else {

                return(default)

        }
        
    }

    
    # First map cell names
    cellTypeNames = unique(seurat_object@meta.data[,cellTypeField,T])
    CellOntologyMapping <- sapply(cellTypeNames, FUN=function(x){
        return(OLSearch(x,ontology = "cl",default='CL_0000255'))   
    })

    cl <- Ontology("cl")
    CellOntologyLabelMapping = sapply(names(CellOntologyMapping),FUN=function(x) {termLabel(Term(cl, gsub("_",":",CellOntologyMapping[[x]],fixed=T)))})

    ct <-seurat_object@meta.data[,cellTypeField,T]
    names(ct) <- colnames(seurat_object)
    seurat_object$cell_type =   sapply(ct,FUN=function(x) {CellOntologyMapping[[as.character(x)]]})                       
    seurat_object$cell_type__ontology_label =   sapply(ct,FUN=function(x) {CellOntologyLabelMapping[[as.character(x)]]})  
        
    # Assign Sex
    seurat_object$sex = "unknown"
    seurat_object$sex[maleCells] = "male"
    seurat_object$sex[femaleCells] = "female"
    seurat_object$sex[mixedCells] = "mixed"

    # Assign biosample and donor
    seurat_object$biosample_id = seurat_object[[biosampleField]]
    seurat_object$donor_id = seurat_object[[donorField]]

    taxon <- Ontology("ncbitaxon")
    #Assign Taxon
    taxon = OLSearch(species,ontology = "ncbitaxon",default='NCBITaxon_131567')
    seurat_object$species  = taxon
    seurat_object$species__ontology_label = termLabel(Term("ncbitaxon", gsub("_",":",taxon,fixed=T)))

    #Assign Disease
    disease = OLSearch(diseaseName,ontology = "mondo",default='MONDO_0700096')
    seurat_object$disease  = disease
    seurat_object$disease__ontology_label = termLabel(Term("mondo", gsub("_",":",disease,fixed=T)))    
    seurat_object$disease[normalCells]  = "PATO_0000461"
    seurat_object$disease__ontology_label[normalCells] = termLabel(Term("pato", gsub("_",":","PATO_0000461",fixed=T)))

    #Assign organ
    organ = OLSearch(organName,ontology = "uberon",default='MONDO_0700096')
    seurat_object$organ  = organ
    seurat_object$organ__ontology_label = termLabel(Term("uberon", gsub("_",":",organ,fixed=T)))    

    #Assign library_preparation_protocol
    protocol = OLSearch(library_preparation_protocol,ontology = "efo",default='EFO_0010183')
    seurat_object$library_preparation_protocol  = protocol
    seurat_object$library_preparation_protocol__ontology_label = termLabel(Term("efo", gsub("_",":",protocol,fixed=T))) 
    return(seurat_object)
}