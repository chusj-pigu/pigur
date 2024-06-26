library(openxlsx)
library(Seurat)

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
#'                                      file = "./path/to/export/file.xlsx"
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