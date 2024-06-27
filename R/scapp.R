
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
 optfuns = c(ENCOblue="BlueprintEncodeData", DICE="DatabaseImmuneCellExpressionData", 
HPCA="HumanPrimaryCellAtlasData", ImmGen="ImmGenData", Monaco="MonacoImmuneData", 
mouse="MouseRNAseqData", hematopo="NovershternHematopoieticData")
 ui = fluidPage(
  shinytoastr::useToastr(),
  sidebarLayout(
   sidebarPanel(
    helpText("app for labeling single cells with selected references"),
    helpText(sprintf("version %s", packageVersion("diamonds24"))),
    radioButtons("ref", "refs", optfuns, selected="HumanPrimaryCellAtlasData"),
# consider option for label.main, label.fine
    numericInput("ncomp", "npcs", min=2, max=5, value=2),
    helpText("provide an upload function here"), width=2
    ),
    mainPanel(
     tabsetPanel(
      tabPanel("main", plotOutput("view")),
      tabPanel("interact", plotly::plotlyOutput("viewly")),
      tabPanel("author", plotOutput("auth")),
      tabPanel("ref comp", verbatimTextOutput("called")),
      tabPanel("ref comp2", verbatimTextOutput("called2"))
     )
    )
   )
  )
 server = function(input, output) {
  # needs help here -- 1) use upload method, 2) verify gene symbols present, get from rowData if not
  # setup, runPCA once, defines "given"
   given = scRNAseq::MuraroPancreasData() 
   rownames(given) = rowData(given)$symbol
   dups = which(duplicated(rownames(given)))
   if (length(dups)>0) given = given[-dups,]
   given = scuttle::logNormCounts(given)
   given = scater::runPCA(given)
  #
  run_SingleR = reactive({
   ref2use = get(input$ref)()
   myb = BiocParallel::MulticoreParam(4)
   shinytoastr::toastr_info("starting SingleR")
   sing = SingleR::SingleR(given, ref2use, ref2use$label.main, BPPARAM=myb)
   shinytoastr::toastr_info("done")
   given$celltype = sing$labels
   #scater::runPCA(given)
   given
   })
  output$view = renderPlot({
   given = run_SingleR()
   scater::plotPCA(given, colour_by = "celltype", 
        ncomponents=input$ncomp, theme_size=14)
   })
  output$viewly = plotly::renderPlotly({
   given = run_SingleR()
   dfr = SingleCellExperiment::reducedDim(given)
   mydf = data.frame(PC1 = dfr[,1], PC2=dfr[,2], type=given$celltype)
   gg = ggplot2::ggplot(mydf, aes(x=PC1, y=PC2, text=type,
      colour=type)) +
     ggplot2::geom_point()
   plotly::ggplotly(gg)
   })
  output$auth = renderPlot({
   scater::plotPCA(given, colour_by = "label",
        ncomponents=input$ncomp, theme_size=14)
   })
  output$called = renderPrint({
   ref2use = get(input$ref)()
   sort(table(ref2use$label.main))
   })
  output$called2 = renderPrint({
   ref2use = get(input$ref)()
   sort(table(ref2use$label.fine))
   })
 }
 runApp(list(ui=ui, server=server))
}

   
 
 
