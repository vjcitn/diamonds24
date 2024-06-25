
#' app demo for SingleR application
#' @import shiny
#' @import shinytoastr
#' @import celldex
#' @importFrom SingleR SingleR
#' @importFrom scRNAseq MuraroPancreasData
#' @import scran
#' @import scater
#' @import scuttle
#' @import BiocParallel
#' @note We work just with the MuraroPancreasData, leaving
#' more general upload capability as a project.
#' @export
scapp = function() {
# derived from ls("package:celldex")
 optfuns = c("BlueprintEncodeData", "DatabaseImmuneCellExpressionData", 
"HumanPrimaryCellAtlasData", "ImmGenData", "MonacoImmuneData", 
"MouseRNAseqData", "NovershternHematopoieticData")
 ui = fluidPage(
  shinytoastr::useToastr(),
  sidebarLayout(
   sidebarPanel(
    helpText("app for labeling single cells with selected references"),
    radioButtons("ref", "refs", optfuns),
# consider option for label.main, label.fine
    numericInput("ncomp", "npcs", min=2, max=5, value=2),
    helpText("provide an upload function here"), width=2
    ),
    mainPanel(
     plotOutput("view")
    )
   )
  )
 server = function(input, output) {
  # needs help here -- 1) use upload method, 2) verify gene symbols present, get from rowData if not
  build_sce = reactive({
   given = scRNAseq::MuraroPancreasData() 
   rownames(given) = rowData(given)$symbol
   dups = which(duplicated(rownames(given)))
   if (length(dups)>0) given = given[-dups,]
   given = scuttle::logNormCounts(given)
   ref2use = get(input$ref)()
   myb = BiocParallel::MulticoreParam(4)
   shinytoastr::toastr_info("starting SingleR")
   sing = SingleR::SingleR(given, ref2use, ref2use$label.main, BPPARAM=myb)
   shinytoastr::toastr_info("done")
   given$celltype = sing$labels
   scater::runPCA(given)
   })
  output$view = renderPlot({
   given = build_sce()
   scater::plotPCA(given, colour_by = "celltype", 
        ncomponents=input$ncomp, theme_size=14)
   })
 }
 runApp(list(ui=ui, server=server))
}

   
 
 
