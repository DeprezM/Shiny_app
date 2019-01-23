library(shiny)
library(Seurat)
library(dplyr)
library(Matrix)


colors.BBInt <- c("Cycling Basal" = "#F3766E",
                  "Basal" = "#A31E22",
                  "Basal 2" = "#A31E22", "Basal 1" = "#d7282d",
                  "Suprabasal" = "#FCCC0A",
                  "Suprabasal 3" = "#ecc113",
                  "Suprabasal 1" = "#FCCC0A",
                  "Suprabasal 2" = "#e0c96c",
                  "Secretory" = "#53c653", 
                  "Secretory 1" = "#84C553", "Secretory 2" = "#a0d279", 
                  "Secretory 3" = "#53c653",
                  "Multiciliated" = "#466cb9",  "MCCs" = "#466cb9",
                  "Deuterosomal" = "#2da9d2",
                  "PNEC" = "#904fd1",
                  "Precursor" = "#e067b5",
                  "Ionocyte" = "#c63563",
                  "Surfactant +" = "#a7d6a7","BPIFB2 +" = "#75ff75",
                  "Goblet" = "#006600",
                  "Glandular"="#3e937f","Endothelial"="#281c1c","Fibroblast"="#aea9ce",
                  "Immune"="#e2bae2","Mesenchymal"="#7c5959","Smooth muscle"="#7c5959",
                  "Macrophage 2"="#ca81ca", "Myofibroblast" = "#a98989",
                  "Macrophage"="#e2bae2","Lymphocyte"="#ffd3ff","Plasmocyte"="#7e6384","Mastocyte"="#eea5ff")


genes_list <- list(
  "Cycling Basal" = c("MKI67", "CDC20", "TOP2A", "HIST1H4C", "KRT5", "TP63"),
  "Basal" = c("KRT5", "KRT15", "MIR205HG", "TP63", "DLK2", "SERPINF1"),
  "Suprabasal" = c("SERPINB3", "SERPINB4", "SPRR3", "KRT13", "MT1X", "NTS"),
  "Secretory" = c("SCGB1A1", "VMO1", "MSMB", "BPIFA1", "TFF3", "LYZ"),
  "Goblet" = c("MUC5AC", "CXCL3", "FCGBP", "SCGB1A1", "VMO1", "TFF3"),
  "Deuterosomal" = c("CDC20Bshort", "CCNO", "HES6", "CCDC67", "FOXJ1", "NEK2"),
  "Multiciliated" = c("FOXJ1", "TPPP3", "PIFO", "IFT57", "OMG", "LRRC23"),
  "Ionocyte" = c("CFTR", "ASCL3", "ATP6V1B1", "ATP6V1G3", "FOXI1", "PDE1C"),
  "PNEC" = c("GRP", "PCSK1N", "SCG5", "SCG2", "CHGA", 'ASCL1'),
  "Mucous" = c("BPIFB2", "MUC5B", "TFF1", "PLA2G2A", "MSLN", "CLDN22"),
  "Glandular" = c("LTF", "LYZ", "TCN1", "PIP", "DMBT1", "FOLR1"),
  "Fibroblast" = c("FBLN1", "DPT", "PTGDS", "LUM", "COL1A1", "DCN"),
  "Smooth muscle" = c("ACTA2", "ACTG2", "TAGLN", "MYL9", "CALD1", "RERGL"),
  "Endothelial" = c("DARC", "AQP1", "GNG11", "VIM", "GIMAP7", "VWF"),
  "Macrophage" = c("CD68", "TYROBP", "LST1", "FCER1G", "LAPTM5", "AIF1"),
  "Plasmocyte" = c("CD79A", "PTPN7", "IGLL5", "IGJ", "CD79B", "MZB1"),
  "Mastocyte" = c("TPSAB1", "SLC18A2", "CPA3"),
  "Lymphocyte" = c("NKG7", "GZMA", "GZMB", "CD2", "CCL5", "CD3D"),
  "Monocyte" = c("MS4A6A", "AIF1", "LST1", "COTL1")
)

## Create App layout
ui <- fluidPage(
    # Create a full page for single dataset analysis or full aggregate analysis
    fileInput(inputId="file_name", label="Select Saved Seurat Object File"),
    navbarPage(title = textOutput("file_name"),
               #### INDIVIDUAL ####
               tabPanel("Individual",
                        fluidRow( # split page layout
                          column(3,  # sidebar
                                 wellPanel(
                                   radioButtons(inputId = "plotType", label="Plot type?", choices = c("t-SNE" = "tsne","Umap" = "umap"), selected = "tsne"),
                                   tags$hr(),
                                   sliderInput(inputId = "dotSize", label = "Set point size", value=1, min=0.01, max=10),
                                   sliderInput(inputId = "labelSize", label = "Set label size", value=4, min=0, max=10)
                                 ),
                                 wellPanel(
                                   conditionalPanel(
                                     condition = "output.fileUploaded == true",
                                     selectizeInput(inputId = "geneName", label = "Enter Gene Name", choices = NULL, multiple = T),
                                     #textInput(inputId = "geneName", label="Enter Gene Name", value="MALAT1"),
                                     actionButton("goButton", "Plot Gene")
                                   )
                                 ),
                                 wellPanel(
                                   selectizeInput(inputId = "markerGene", label = "Select Cell Type",
                                                  choices = names(genes_list)),
                                   actionButton("goButtons", "Plot Genes")
                                 )
                          ), # end sidebar
                          column(9, # main panel
                                 plotOutput("TSNE"),
                                 tags$hr(),
                                 tags$h4("Gene Expression"),
                                 tabsetPanel(
                                   tabPanel("TSNE", plotOutput("featureTSNE")),
                                   tabPanel("UMAP", plotOutput("featureUmap")),
                                   tabPanel("Violins", plotOutput("vlnPlot"))),
                                 tags$hr(),
                                 tags$h4("Marker Genes Expression"),
                                 tabsetPanel(
                                   tabPanel("TSNE", plotOutput("featureTSNEs")),
                                   tabPanel("UMAP", plotOutput("featureUmaps")),
                                   tabPanel("Violins", plotOutput("vlnPlots")))
                          ) # end main panel
                        ) # end fluid row 
                        
                        
               ), # end tabPanel Individual
               
               #### AGGREGATE ####
               tabPanel("Aggregate",
                        fluidRow( # split page layout
                          column(3,  # sidebar
                                 wellPanel(
                                   radioButtons(inputId = "plotType_aggr", label="Plot type?", choices = c("t-SNE" = "tsne","Umap" = "umap"), selected="tsne"),
                                   tags$hr(),
                                   sliderInput(inputId = "dotSize_aggr", label = "Set point size", value=1, min=0.01, max=10)
                                 ),
                                 wellPanel(
                                   selectizeInput(inputId = "geneName_aggr", label = "Enter Gene Name", choices = NULL, multiple = T),
                                   uiOutput("Group_1"),
                                   # uiOutput("Group_2"),
                                   actionButton("goButtons_aggr", "Plot Gene")
                                 )
                          ), # end sidebar
                          column(9, # main panel
                                 plotOutput("TSNE_aggr"),
                                 plotOutput("TSNE_aggr_manip"),
                                 tags$hr(),
                                 tags$h4("Gene Expression per Group"),
                                 tabsetPanel(
                                   tabPanel("Violins", plotOutput("vlnPlot_aggr")))
                          ) # end main panel
                        ) # end fluid row 
                      
               ) # end tabPanel Aggregate
               
    ) # end navbarPage
  ) # end fluidPage

# )


#### SERVER ####

server <- function(input, output, session) {
  
  # increase the max upload file size
  options(shiny.maxRequestSize=10000000*1024^2)
  
  SeuratObject <- reactive({
    if(is.null(input$file_name)) return(NULL)
    # Load user-defined dataset
    load(input$file_name$datapath)
  })
  
  output$fileUploaded <- reactive({
    return(!is.null(SeuratObject()))
  })
  outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)
  
  # Retrieve user-defined dataset name
  output$file_name <- renderText({
    unlist(strsplit(input$file_name$name, ".", fixed= T))[[1]]
  })
  
  
  #### INDIVIDUAL ####
  # Display TSNE Plot as overview toggled options for display legends/labels and setting size. Default displays demo dataset
  
  output$TSNE <- renderPlot({
    load(input$file_name$datapath)
      DimPlot(get(SeuratObject()), reduction = input$plotType, label = T, 
                                       pt.size = input$dotSize, label.size = input$labelSize, cols = colors.BBInt)
  })
  
  observe({
    req(input$file_name)
    gene_name_list <- rownames(GetAssayData(object = get(SeuratObject()), assay = "RNA"))
    updateSelectizeInput(session, "geneName", "Enter Gene Name", choices = gene_name_list)

  })

  fPlots <- eventReactive(input$goButtons, {
    load(input$file_name$datapath)
    FeaturePlot(get(SeuratObject()), 
                genes_list[[input$markerGene]][genes_list[[input$markerGene]] %in% rownames(GetAssayData(object = get(SeuratObject())))], 
                ncol = 3, reduction = "tsne", pt.size = input$dotSize, cols = c("grey70", "red4"))
  })
  
  uPlots <- eventReactive(input$goButtons, {
    load(input$file_name$datapath)
    FeaturePlot(get(SeuratObject()), 
                genes_list[[input$markerGene]][genes_list[[input$markerGene]] %in% rownames(GetAssayData(object = get(SeuratObject())))], 
                ncol = 3, reduction = "umap", pt.size = input$dotSize, cols = c("grey70", "red4"))
  })
  
  vPlots <- eventReactive(input$goButtons, {
    load(input$file_name$datapath)
    VlnPlot(get(SeuratObject()), features = 
              genes_list[[input$markerGene]][genes_list[[input$markerGene]] %in% rownames(GetAssayData(object = get(SeuratObject())))], ncol = 3) +
      ggplot2::scale_fill_manual(values = colors.BBInt)
  })
  
  output$featureTSNEs <- renderPlot(fPlots())
  output$featureUmaps <- renderPlot(uPlots())
  output$vlnPlots <- renderPlot(vPlots())
  
  
  
  fPlot <- eventReactive(input$goButton, {
    load(input$file_name$datapath)
    FeaturePlot(get(SeuratObject()), input$geneName, reduction = "tsne",  pt.size = input$dotSize, 
                cols = c("grey70", "red4")) + NoLegend()
  })
  
  uPlot <- eventReactive(input$goButton, {
    load(input$file_name$datapath)
    FeaturePlot(get(SeuratObject()), input$geneName, reduction = "umap", 
                pt.size = input$dotSize, cols = c("grey70", "red4")) + NoLegend()
  })
  
  vPlot <- eventReactive(input$goButton, {
    load(input$file_name$datapath)
    VlnPlot(get(SeuratObject()), features = input$geneName) +
      ggplot2::scale_fill_manual(values = colors.BBInt) +
      ggplot2::xlab(paste0("Expr cells : ",sum(GetAssayData(object = get(SeuratObject()))[input$geneName,] > 0),
                  " / ", ncol(GetAssayData(object = get(SeuratObject())))))
  })
  
  output$featureTSNE <- renderPlot(fPlot())
  output$featureUmap <- renderPlot(uPlot())
  output$vlnPlot <- renderPlot(vPlot())
  
  #### AGGREGATE ####
  
  output$TSNE_aggr <- renderPlot({
    load(input$file_name$datapath)
    DimPlot(get(SeuratObject()), reduction = input$plotType_aggr,
                        pt.size = input$dotSize_aggr, cols = colors.BBInt)
  })
  
  output$TSNE_aggr_manip <- renderPlot({
    load(input$file_name$datapath)
    DimPlot(get(SeuratObject()), reduction = input$plotType_aggr, group.by = "manip",
            pt.size = input$dotSize_aggr)
  })
  
  observe({
    req(input$file_name)
    gene_name_list <- rownames(GetAssayData(object = get(SeuratObject()), assay = "RNA"))
    updateSelectizeInput(session, "geneName_aggr", "Enter Gene Name", choices = gene_name_list)
    
  })
  
  # create group list to check
  output$Group_1 <- renderUI({
    celltypes <- unique(get(SeuratObject())$celltype)
    checkboxGroupInput("group1_ckd", "Choose Cell Type", celltypes)
  })
  
  # output$Group_2 <- renderUI({
  #   manips <- unique(get(SeuratObject())$manip)
  #   checkboxGroupInput("group2_ckd", "Choose Dataset", manips)
  # })
  
  vPlot_aggr <- eventReactive(input$goButtons_aggr, {
    load(input$file_name$datapath)
    VlnPlot(get(SeuratObject()), features = input$geneName_aggr, idents = input$group1_ckd, 
            group.by = "manip", assay = "RNA")# +
      # ggplot2::scale_fill_manual(values = colors.BBInt) +
      # ggplot2::xlab(paste0("Expr cells : ",sum(GetAssayData(object = get(SeuratObject()))[input$geneName,] > 0),
      #                      " / ", ncol(GetAssayData(object = get(SeuratObject())))))
  })
  
  output$vlnPlot_aggr <- renderPlot(vPlot_aggr())
  
  
}

shinyApp(ui = ui, server = server)
