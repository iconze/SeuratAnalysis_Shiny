library(Matrix)
library(DT)
library(shiny)
library(shinybusy)
library(shinyWidgets)
library(dplyr)
library(plyr)
library(Seurat)
library(ggplot2)
library(scater)
library(BiocManager)
library(glmGamPoi)
library(sctransform)
library(gotop)
library(stringr)
library(viridis)

options(shiny.maxRequestSize=100*1024^2)
### Start server code 
{
ui <- shinyUI(fluidPage(
  use_gotop(),
  navbarPage(
    title = "scRNA analysis: Seurat workflow",
    
    tabPanel(
      "Seurat Workflow",
      fixedRow(
        titlePanel("Seurat Workflow"),
        p("Process your data in this first tab, before having a look at its Differential Expression and Cluster formation."),
        hr()
      ),
      fixedRow(
        h2("Upload 10x genomics data for processing!"),
        p("To upload your data, you need one .out folder for each dataset. In that folder, you need a folder called", em("GeneFull"), " whithin which you need at least one folder with either", em("filtered")," or", em("raw"), 
          "data (barcodes.tsv, genes.tsv and matrix.mtx) with the folders named accoridngly. Put all .out folders in the same location and specify the directory in the", 
          em("Directory"),
          "field. End with a", em("/"),"."), 
        p("Next, type the initial folder name in the", em("Foldername"), "field. In the", em("Condition"), "field, 
        specify which of your two conditions this dataset belongs to and select the appropriate replicate in the", em("Replicate"), "field.
          Lastly, you can chose if you want to use the",
          em("Filtered"), "or", em("Raw"), "data available in your .out folder."),
        p("Click", em("Use this dataset!"), "to use the specified dataset. Check if the directory and 
          characteristics are correct in the given table. You can modify them by addressing them using the", em("Dataset Number"), "selection.
          Add as many datasets as you like. Once all datasets are listed in the table, click", em("Load and merge these datasets"), "to continue by merging the
          datasets.")
        
      ),
      fixedRow(
        textInput("inputdir", "Directory", "C:/Users/isabe/Elective_BeyerLab/Module_Code/"),
        verbatimTextOutput("Directory")
      ),
      
      fixedRow(
        column(2, uiOutput("DatasetNo")),
        #column(2, numericInput("DatasetNo_I", "Dataset Number...", 1)),
        column(3, textInput("inputfile", "Foldername of the dataset","Y1K1Solo.out" )),
        column(3, textInput("type_I", "Type/Condition of the dataset", "Young")),
        column(2, numericInput("replicate_I", "Replicate", 1)),
        column(2, selectInput("filter_status_I", "Filtered or Raw?",
                              choices = c("Filtered", "Raw"), selected = "Filtered"))
      ),
      
      fixedRow(
        actionButton("Input_button", "Use this dataset!")
      ),
      
      fixedRow(
        tableOutput("filetable")
      ),
      
      fixedRow(
        actionButton("load_button", "Load and merge these datasets")
      ),
      
        ## Output for initial violin plots to compare the datasets

        #Input for Quality Control 
        fixedRow(
          h2("Quality Control"),
          p("To see the data with the currently set threshholds click 'Show Plots'. To filter the data using the threshholds click 'Apply' - 
            the removal of data cannot be reversed prior to subsequent analysis." )
        ),
        fixedRow(
          column(2, offset = 1,
                 br(),
                 uiOutput("nFeat_selection")),
          column(2,
                 br(),
                 br(),
                 uiOutput("nCount_selection")),
          column(2, 
                 br(),
                 br(),
                 uiOutput("nCount_doublet")),
          column(2, 
                 uiOutput("mt_selection")),
          column(2,
                 br(),
                 uiOutput("rbs_selection"))
        ),
        #Button for starting quality control
        fixedRow(
          hr(),
          column(3, offset = 3,
                 helpText("See your data with the current threshholds without modifying the data itself."),
                 actionButton("show_button", "Show plots"),
                 add_busy_spinner(spin = "fading-circle")#progress indicator
                 
          ),
          column(3,
                 helpText("SCTransform and PCA is started automatically"),
                 br(),
                 actionButton("start_button", "Apply")
                 
          )
          
        ),
        #Violin Plots of data  
        fixedRow(
          hr(),
          h3("Violinplots:"),
          column(4, 
                 uiOutput("selection_meta"),
                 uiOutput("selection_split")),
          column(4, plotOutput("vplot"))
          
        ),
        
        ## Output for adjustable ggplots 
        fixedRow(
          hr(),
          h3("Adjustable ggplots:"),
          column(3,
                 uiOutput("plot_selection"),
                 #for bar-plot: type/rep
                 uiOutput("plot_data"),
                 #for density-plot & point: metadata & intercepts
                 uiOutput("plot_data2"),
                 uiOutput("color")
          ),
          column(7,
                 plotOutput("ggplot")
          ),
          column(2,
                 uiOutput("intercept2"),
                 uiOutput("intercept1")
          )
        ),
        
        
        #Plot after SCTransform: variable features
        fixedRow(
          hr(),
          h3("Variable Feature Plot"),
          helpText("Is shown after SCTransform."),
          plotOutput("VarFeat_plot")
        ),
        
        #Elbow Plot once PCA is finished
        fixedRow(
          h3("Elbow Plot"),
          helpText("Is shown after PCA and can be used to assess how many dimensions might be 
                   suitable for downstream analysis."),
          plotOutput("Elbow")
        ),
        fixedRow(
          hr(),
          h1("Clustering"),
          p("Can be used to find appropriate clusters - settings can be changed and clusters are calculated.
            Once you are happy with the shown clusters, continue with finding Cluster markergenes
            and differentially expressed genes.")
        ),
        #Clustering settings + button to start clustering
        fixedRow(
          column(3,offset = 1,
                 uiOutput("npcs")),
          column(3, uiOutput("dims")),
          column(3, 
                 br(),
                 uiOutput("res"))

        ),
        
        #Settings for Dim-/Featureplots
        fixedRow(
          column(3, offset = 1,
          uiOutput("clusfeat")
          ),
          column(3, offset = 3,
             br(),
             br(),
             actionButton("Cluster", "Start Clustering"))
        ),
      fixedRow(
        
      ),
      
      fixedRow(
          column(6,plotOutput("Dimplot")),
          column(6, plotOutput("Dimplot1"))
        ),
      
      fixedRow(
        h2("Rename your clusters!"),
        p("If you don't want to use simple numbers as cluster names, you can rename them here. (If you want to keep the numbers, please click",
        em("Use these names"), "and continue with your analysis. If you want to change the names,
          select the cluster that you want to rename in the", em("Cluster to be renamed"),
          "box. Next, you can type the name for that cluster into the", em("Clustername"), "box. Names have to be unique within your dataset,
          but not all clusters have to be renamed.
          To use the specified name, click", em("Give name"), ". Once you are happy with all given names click", 
          em("Use these names"), "and continue with your analysis.")
      ),
      fixedRow(
        column(3, uiOutput("Clus_select")),
        column(6, uiOutput("Clus_name"),
               actionButton("name_button", "Give name"),
        )
      ),
      fixedRow(
        uiOutput("name_error"),
        uiOutput("name_table")
      ),
      fixedRow(
        uiOutput("Names_given"),
        actionButton("apply_names", "Use these names")
      ),
      fixedRow(
        plotOutput("Dimplot_newnames")
      ),
      
      fixedRow(
          h2("Find differentially expressed genes and cluster markers"),
          p("Happy with the clusters? Let's find cluster markers and genes that are differentially expressed between the clusters.")
        ),
      fixedRow(
          #marker button starts findallmarkers()
          column(4, offset = 4,
          actionButton("marker_button", "Find DE and markers")
          )
        ),
      fixedRow(
          br(),
          p("Once analysis is done you can have a closer look at your data in the other Tabs!")
      )
    ), #close panel
    
    tabPanel(
      "Differential Expression",
      fixedRow(
        titlePanel("Differentially expressed genes and cluster marker genes")
      ),
      fixedRow(
        #select which cluster to show from markerselection
        uiOutput("selectmclus")
      ),
      
      fixedRow(
        #table for markergenes
        h4("Table with cluster marker genes"),
        tableOutput("mytable")
      ), 
      
      fixedRow(
        h4("Table with genes that are differentially expressed"),
        uiOutput("DEselection")
      ),
      
      fixedRow(
        DT::dataTableOutput("DEtable")
      ),
      
      fixedRow(
        h4("Plot gene expression"),
        p("Plot the expression of any gene you want in a Violin Plot or a Feature Plot, or look at the most variable cluster marker genes in a heatmap."),
        column(3, uiOutput("plotsel"))
      ),
      
      fixedRow(
        column(3, uiOutput("genesel")),
        column(3, uiOutput("splitby"))
      ),
      fixedRow(
        plotOutput("vplot2")
      )
      #use_gotop(),
    ),
    
    tabPanel(
      "Clustering",
      fixedRow(
        h4("Cellinfo/Genes vs. Cellinfo/Genes in a Featureplot"),
        p("Select if you want to display cellinfo or genes as features in each plot!")
      ),
      fixedRow(
        column(6, uiOutput("sel1_left"),
               uiOutput("sel2_left"),
               plotOutput("plot1")),
        column(6, uiOutput("sel1_right"),
               uiOutput("sel2_right"),
               plotOutput("plot2"))
      ),
    
    fixedRow(
      h4("Proportion plot"),
      column(3,
             uiOutput("cellinfosel"),
             uiOutput("groupsel")),
      column(9, 
             plotOutput("propplot"))
    )
    #use_gotop()
    ),
    tabPanel(
      "Type Comparison",
      
      fixedRow(
        h4("Find genes differentially expressed between conditions.")),
      add_busy_spinner(spin = "fading-circle"),#progress indicator
      fixedRow(
        column(10, textOutput("deexpl"))
      ),
      fixedRow(
               DT::dataTableOutput("DEtable_type")
        ),
      
      
      fixedRow(
        column(8, plotOutput("scattertype", click = "scattertype_click")),
        column(4, 
               uiOutput("genesel1"),
               htmlOutput("x_value")
        )
      ),
     
        #Dim-/Featureplots
        fixedRow(
          plotOutput("Dimplot2")
        )
       
      #use_gotop(),
    )
    
   
  ) #Close: navbarPage 

))} #Close: fluidPage


# Define server logic ----
#sce$combined <- readRDS("/data/user/apapada1/Downloads/Seuratobject_bfQC_merged.rda")
#sce$combined$log10GenesPerUMI <- log10(sce$combined$nFeature_RNA) / log10(sce$combined$nCount_RNA)

findpercentmt <- function(seurat){
  seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
  return(seurat)
}
findpercentrbs <- function(seurat){
  seurat[["percent.rbs"]] <- PercentageFeatureSet(seurat, pattern = "Rp[sl]")
  return(seurat)
}

metas1 <- list("No. of UMIs" = "nCount_RNA", "No. of features" = "nFeature_RNA", "Mitochondrial percentage" = "percent.mt", "Ribosomal percentage" = "percent.rbs")

server <- function(input, output) {
  
  output$Directory <- renderText({ input$inputdir })
  
  vecreps <- reactiveValues(reps = c())
  vectypes <- reactiveValues(types = c())
  levtypes <- reactiveValues(levels = c())
  vecdirectories <- reactiveValues(directories = c())
  vecfiles <- reactiveValues(files = c())
  filter <- reactiveValues(status = c())
  dataf <- reactiveValues(frame = data.frame(c(), c()))
  seurat <- reactiveValues(list = c())
  cell <- reactiveValues(ID = c())
  sce <- reactiveValues(combined = seurat)
  new <- reactiveValues(seurats = c())
  veclev <- reactiveValues(levels=c())
  number <- reactiveValues(inp = 1)
  height <- reactiveValues(n=800)
  y <- reactive(
  height$n
  )
  
  output$DatasetNo <- renderUI({
    numericInput("DatasetNo_I", "Dataset Number...", number$inp)
    
  })
  
  observeEvent(input$Input_button,{
    vecreps$reps[input$DatasetNo_I] <- as.character(input$replicate_I)
    vectypes$types[input$DatasetNo_I] <- input$type_I
    levtypes$levels <- levels(vectypes$types)
    veclev$levels[input$DatasetNo_I] <- levtypes$levels
    vecdirectories$directories[input$DatasetNo_I] <- paste(input$inputdir, input$inputfile, "/GeneFull/", input$filter_status_I, sep ="")
    vecfiles$files[input$DatasetNo_I] <- input$inputfile
    filter$state[input$DatasetNo_I] <- input$filter_status_I
    cell$ID[input$DatasetNo_I] <- paste(input$type_I, input$replicate_I, sep = "")
    dataf$frame <- data.frame(vecfiles$files, vectypes$types, vecreps$reps, filter$state, vecdirectories$directories)
    number$inp <- number$inp + 1
  })
  
  output$filetable <- renderTable({
    validate(
      need(length(dataf$frame) != 0, "Please chose your first data set"),
      need(str_sub(input$inputdir,-1) == "/", "Make sure the directory is correct and the typed directory ends with a /"),
      need(str_sub(input$inputfile, -4) == ".out", "Make sure your inputfolder is an .out folder")
    )
    colnames(dataf$frame) <- c("Folder", "Type", "Replicate", "Filtered or Raw?", "Directory")
    dataf$frame
  })
  
  addtype <- function(seuratcolumn, chars){
    seuratcolumn$type <- chars
    return(seuratcolumn)
  }
  
  addrep <- function(seurat, chars){
    seurat$rep <- chars
    return(seurat)
  }
  
  
  observeEvent(input$load_button,{
    seurat$list <- lapply(vecdirectories$directories, Read10X)
    seurat$list <- lapply(seurat$list, CreateSeuratObject, project = "kimmel.kidney", min.cells = 6, min.features = 200)
    seurat$list <- lapply(seurat$list, findpercentmt)
    seurat$list <- lapply(seurat$list, findpercentrbs)
    seurat$list <- mapply(addtype, seurat$list, vectypes$types)
    seurat$list <- mapply(addrep, seurat$list, vecreps$reps)
    
    print(length(seurat$list))
    for (i in 2:length(seurat$list)){
      new$seurats[[i-1]] <- seurat$list[[i]]
    }
    
    sce$combined <- merge(seurat$list[[1]], y=new$seurats, add.cell.ids = cell$ID)
  })
  
  #Quality Control
  
    
    #Input: QC - Features
    output$nFeat_selection <- renderUI(
      numericInput("nFeat", label = paste("Remove cells with less than ... features"), value = 0)
    )
    #Input: QC - UMIs
    output$nCount_selection <- renderUI(
      numericInput("nCount", label = "Remove cells with less than ... UMIs", value = 0)
    )
    
    output$nCount_doublet <- renderUI(
      numericInput("nCount_doub", label = "Remove cells with more than ... UMIs", value = 1000000)
    )
    
    #Input QC - mitochondrial %
    output$mt_selection <- renderUI(
      numericInput("mt_cut", label = "Remove cells with more than ...% mitochondrial genes", value = 100)
    )
    
    #Input QC - ribosomal %
    output$rbs_selection <- renderUI(
      numericInput("rbs_cut", label = "Remove cells with more than ...% ribosomal genes", value = 100)
    )
    
    #Start with "Show plots"/show_button:
    observeEvent(input$show_button,{
      sce.combined1 <-sce$combined
      sce.combined1$log10GenesPerUMI <- log10(sce.combined1$nFeature_RNA) / log10(sce.combined1$nCount_RNA)
      sce.combined1 <- subset(sce$combined, subset =
                               (nFeature_RNA > input$nFeat) &
                               (nCount_RNA > input$nCount) &
                               (nCount_RNA < input$nCount_doub) &
                               (percent.mt < input$mt_cut) &
                               (percent.rbs < input$rbs_cut))
      metadata <- sce.combined1@meta.data

      #Generate input for Violin plots
      output$selection_meta <- renderUI({
        #factormeta1 <<- c("type", "rep", "orig.ident")
        selectInput("metas", "Select metadata to display", 
                    choices = metas1, 
                    selected=metas1[1])
      })
      #Decide split-variable of Violin plots
      output$selection_split <- renderUI({

          selectInput("split", "Select split variable",
                             choices = list("type", "rep"),
                             selected = "type")
        })
      

      #Print violin plots
      output$vplot <-  
        renderPlot({
          VlnPlot(object = sce.combined1, features = input$metas, split.by = input$split)
        })
      
      #Generate modifiable ggplots
      {
        #Input: What kind of plot do you want?
        output$plot_selection <- 
          renderUI(selectInput("plotselected", "Select plot to display", 
                               choices = list("bar" = "bar1", "density" = "density1", "point" = "point1"), 
                               selected="bar")
          )
        
        #Render second selection based on what kind of plot you want
        
        output$plot_data <- renderUI({
          req(input$plotselected)
          #     If you want a barplot: which bars?
          if(input$plotselected=="bar1"){
            selectInput("bar_input",label="Number of cells per...", choices=list("type" = "type1", "rep" = "rep1"),
                        selected="type", multiple = F)
            #If you want a density plot: Counts, Features or %mt?
          }else if (input$plotselected=="density1"){
            selectInput("density_input",label="Density of... (split by type)", choices=metas1, 
                        selected=metas1[1])
          }else if (input$plotselected=="point1"){
            selectInput("density_input",label="Plot...", choices=metas1, 
                        selected=metas1[1])
          }else{return(NULL)}
          
          #render 2nd selection
        })
        
        #Selection of second variable for point-plot
        
        output$plot_data2 <- renderUI({
          if(input$plotselected == "point1"){
            selectInput("plot_input2",label="against...", choices=metas1, 
                        selected=metas1[2])
          }
        })
        
        output$color <- renderUI({
          if(input$plotselected == "point1"){
            selectInput("plot_color", label="Gradient by...", choices=metas1, 
                        selected=metas1[3])
          }
        })
        
        #Position of the intercepts (to assess values for QC)
        output$intercept1 <- renderUI({
          req(input$plotselected)
          req(input$density_input)
          if(input$plotselected=="density1"){
            numericInput("interceptA", label="Show Intercept at...", value = 0)
          }else if (input$plotselected=="point1"){
            numericInput("intercept_hor", label="Show horizontal Intercept at...", value = 0)
          }else{return(NULL)}
        })
        
        #Second intercept for point plot
        output$intercept2 <- renderUI({
          req(input$plotselected)
          req(input$density_input)
          if(input$plotselected=="point1"){
            numericInput("intercept_ver", label="Show vertical Intercept at...", value = 0)
          }else{return(NULL)}
        })
        
        output$ggplot<-renderPlot({
          req(input$plotselected)

          #Plot for barplot
          if(input$plotselected=="bar1"){
            #VlnPlot(object = sce$combined, features = "nCount_RNA")
            req(input$bar_input)
            
            if(input$bar_input == "type1"){
              ggplot(metadata, aes(x=type, fill = type)) +
                geom_bar() +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                theme(plot.title = element_text(hjust=0.5, face="bold")) +
                ggtitle("Number of Cells")
              
            }else if (input$bar_input == "rep1"){
              ggplot(metadata, (x=rep)) +
                geom_bar() +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
                theme(plot.title = element_text(hjust=0.5, face="bold")) +
                ggtitle("Number of Cells")
            }
            #Plot for density plot
          }else if (input$plotselected=="density1"){
            req(input$density_input)
            ggplot(metadata, aes_string(color="type", x=input$density_input, fill="type")) + 
              geom_density(alpha = 0.2) + 
              theme_classic() +
              scale_x_log10() + 
              geom_vline(xintercept = input$interceptA)
            
          }else if (input$plotselected=="point1"){
            req(input$plot_input2)
            ggplot(metadata, aes_string(x=input$density_input, y=input$plot_input2, color=input$plot_color)) + 
              geom_point() + 
              scale_colour_gradient(low = "gray90", high = "black") +
              stat_smooth(method=lm) +
              scale_x_log10() + 
              scale_y_log10() + 
              theme_classic() +
              geom_hline(yintercept = input$intercept_hor) +
              geom_vline(xintercept = input$intercept_ver) +
              facet_wrap(~type)
          }else{return(NULL)}
          
        })
        
        #ggplots
        }              
      
      
      #button show
    })
    
    observeEvent(input$start_button,{
      #sce$combined <-sce$combined
      sce$combined <<- subset(sce$combined, subset =
                                (nFeature_RNA > input$nFeat) &
                                (nCount_RNA > input$nCount) &
                                (nCount_RNA < input$nCount_doub) &
                                (percent.mt < input$mt_cut) &
                                (percent.rbs < input$rbs_cut))
      

      #Run SCTransform
      sce$combined <<- SCTransform(sce$combined, vars.to.regress = c("percent.mt", "percent.rbs"), method = "glmGamPoi", verbose = FALSE)
      
      #Visualize variable features (detected by SCTransform):
      top10 <- head(VariableFeatures(sce$combined), 10)
      #plot variable features:
      output$VarFeat_plot <- renderPlot({
        plot1 <- VariableFeaturePlot(sce$combined) +
          NoLegend()
        LabelPoints(plot = plot1, points = top10, repel = TRUE)
      })
      
      #Run PCA
      sce$combined <<- RunPCA(sce$combined, npcs = 50, verbose = FALSE)
      #Render Elbowplot
      output$Elbow <- renderPlot({
        ElbowPlot(sce$combined, ndims = 50)
      })
    
      #start button  
    })
    
    #SC Transformation
    # observeEvent(input$Transform_button,{
    
    #transformbutton
    # })
    
    #run PCA -> started by button
    # observeEvent(input$PCA_button,{
    
    #PCA button
    # })
    
    #Settings for clustering: npcs, dimensions and resolution
    {
      output$npcs <- renderUI({
        numericInput(
          "npcs_i", label = "Total number of PCs to compute and store", value = 50
        )})
      output$dims <- renderUI({
        numericInput(
          "dims_i", label = "Number of dimensions to consider for clustering", value = 20
        )})
      output$res <- 
        renderUI({numericInput(
          "res_i", label = "Resolution for clustering", value = 0.5
        )})
      output$clusfeat <- renderUI({
        selectInput("clus_feature", label = "Select gradient for feature plot", 
                    choices = metas1, 
                    selected=metas1[1])})
    }
    #Start clustering with the settings - started by button:
    observeEvent(input$Cluster, {
      sce$combined <<- RunPCA(sce$combined, npcs = input$npcs_i, verbose = FALSE)
      sce$combined <<- RunUMAP(sce$combined, reduction = "pca", dims = 1:input$dims_i)
      sce$combined <<- FindNeighbors(sce$combined, reduction = "pca", dims = 1:input$dims_i)
      sce$combined <<- FindClusters(sce$combined, resolution = input$res_i)
      no.clusters <<- c(levels(sce$combined$seurat_clusters))
      #Render dimplot/featureplot: ability to adjust color gradient
      output$Dimplot <- renderPlot({
        FeaturePlot(sce$combined, pt.size=1, label=F, reduction="umap",
                    features= input$clus_feature, shape.by = "type", 
                    cols=c("lightgrey", "darkblue"))
      })
      
      output$Dimplot1 <- renderPlot({
        DimPlot(sce$combined, pt.size=1, reduction="umap")
      })

    #Ability to add Cluster Names
    output$Clus_select <- renderUI({
      name_vec <<- no.clusters
      numericInput("Clus_select_I", "Cluster to be renamed:", 0, min = 0, max = length(name_vec)-1)
    })
    
    output$Clus_name <- renderUI({
      textInput("Clus_name_I", "Name for the cluster")
    })
    
    observeEvent(input$name_button,{
      name_vec[sum(input$Clus_select_I, 1)] <- input$Clus_name_I
      name_vec <<- name_vec
      repl <- data.frame("Cluster" = levels(sce$combined$seurat_clusters), "Given name" = name_vec)
      output$name_error <- renderUI({
        validate(
          need(any(any(input$Clus_name_I == name_vec)==FALSE), "Make sure each name is unique")
        )
      })
        output$name_table <- renderTable({
        repl
      })
    })
    
    observeEvent(input$apply_names,{
      cell.labels <- sce$combined$seurat_clusters
      cell.labels <- mapvalues(cell.labels,
                               from=0:(length(name_vec)-1),
                               to=name_vec)
      sce$combined <<- AddMetaData(sce$combined, cell.labels, 'Cluster.Names')
      print(sce$combined$Cluster.Names)
      name_vec <- levels(sce$combined$Cluster.Names)
      output$Dimplot_newnames <- renderPlot({
        Idents(sce$combined) <- "Cluster.Names"
        DimPlot(sce$combined, label = T, repel = T)
      })
    })
    
    
    #cluster button
    })
    
    #Find all markers (=clustermarkers) - started by button 
    observeEvent(input$marker_button,{
      Idents(sce$combined) <- "Cluster.Names"
      combined.markers <<- FindAllMarkers(sce$combined, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
      combined.markers %>%
        group_by(cluster) %>%
        slice_max(n = 2, order_by = avg_log2FC)
      ClusterNames <- combined.markers$cluster
      ClusterNames <- mapvalues(ClusterNames,
                                from=0:(length(name_vec)-1),
                                to=name_vec)
      combined.markers$NamedCluster <- ClusterNames
      ClusterNames <- as.list(levels(ClusterNames))
      saveRDS(combined.markers, file = "clustermarkers_t.rda")
      saveRDS(sce$combined, file = "Seuratobject_t.rda")
      #Determine cluster for which you want to see the markers
      output$selectmclus <- renderUI({
        selectInput("clusterselect", "Top markers for cluster...", 
                    choices = ClusterNames, selected = ClusterNames[1])
      })
      
      # render table of cluster markers
      output$mytable <- renderTable({
        print("hello3")
        combined.markers <- subset(combined.markers, combined.markers$NamedCluster == input$clusterselect)
        combined.markers <- slice_max(combined.markers, n = 15, order_by = avg_log2FC)
        combined.markers})

      no.clusters <<- c(levels(sce$combined$seurat_clusters))
      Idents(sce$combined) <- "Cluster.Names"
      ClusterNames <- as.list(levels(sce$combined$Cluster.Names))
      DEexpr <- rep(c("1"), length(no.clusters))
      for (i in 1:length(no.clusters)){
        x <- ClusterNames[i]
        markers <- FindMarkers(sce$combined, ident.1 = x, min.pct = 0.25)
        DEexpr[i] <- list(markers)
        #saveRDS(DEexpr, file = "DEAnalysis.rda")
      }
      
      
      # Find differentially expressed genes in each cluster 
      #- chose which cluster to look at 
      output$DEselection <- renderUI({
        selectInput("DEclus", "Differentially expressed genes for cluster...", 
                    choices = ClusterNames, selected = ClusterNames[1])
      })
      
      output$DEtable <- DT::renderDataTable({
        allTab <- Reduce(rbind, lapply(seq(DEexpr), function(ii) {
          DEexpr[[ii]]$cluster <- levels(sce$combined$Cluster.Names)[ii]
          return(DEexpr[[ii]])
        }) )
        
        saveRDS(allTab, "typemarkers.rda")
        allTab <- subset(allTab, cluster == input$DEclus)
        names(DEexpr) <- levels(sce$combined$seurat_clusters)
        
        datatable( allTab, rownames = T)  } ,  extensions= 'Buttons', options = list( scrollY = 25 , dom = 'Bfrtip',
                                                                                      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),paging=T))
            #DE button/#marker button
    })

    output$plotsel <- renderUI ({
      selectInput("plotsel_I", label = "Which plot do you want to see?", choices = list("Violin Plot","Feature Plot", "Heatmap"),
                  selected = "Violin Plot")
    })
    
    output$genesel <- renderUI ({
      textInput("genesel_I", "Which gene do you want to see? (not applicable for heatmap)")
    })
    
    output$splitby <- renderUI ({
      radioButtons("splitby_I", "Split by... (not applicable for heatmap)", choices = list("type", "rep", "none"),
                   selected = "type")
    })
    
    output$vplot2 <- renderPlot({
      Idents(sce$combined) <- "Cluster.Names"
      if(input$plotsel_I == "Violin Plot"){
        validate(
          need(input$genesel_I, "Type in your gene of interest"),
        )
        validate(
          need(any(input$genesel_I == rownames(sce$combined)), "Gene not found") #continue here
        )
        height$n <- 500
        if (input$splitby_I == "none"){
          VlnPlot(object = sce$combined, features = input$genesel_I) +
            scale_x_discrete(labels=sce$combined$Cluster.Names)
        }else{
          VlnPlot(object = sce$combined, features = input$genesel_I, split.by = input$splitby_I)
        }
      }else if (input$plotsel_I == "Feature Plot"){
        validate(
          need(input$genesel_I, "Type in your gene of interest"),
        )
        validate(
          need(any(input$genesel_I == rownames(sce$combined)), "Gene not found") #continue here
        )
        height$n <- 500
        y <<- height$n
        if (input$splitby_I == "none"){
          sce.combined <- sce$combined
          #sce.combined$type <- 
          FeaturePlot(sce$combined, features = input$genesel_I)
        }else{
          FeaturePlot(sce$combined, features = input$genesel_I, split.by = input$splitby_I, reduction="umap")}
      }else if(input$plotsel_I == "Heatmap"){
        top10 <- combined.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
        sce.combined_subs <- subset( sce$combined ,  downsample = 200)
        sce.combined_subs <- ScaleData( sce.combined_subs, features = rownames( sce.combined_subs) )
        height$n <- length(levels(sce$combined$seurat_clusters))*110
        
        # make a plot
        DoHeatmap(sce.combined_subs, features = top10$gene) + 
        NoLegend() 
      }
    }, height = y)
    
    output$sel1_left <- renderUI({
      selectInput("sel1_left_I", label = "Cell info or Genes?",
                  choices = list("Cell info", "Genes"))
    })
    
    output$sel2_left <- renderUI({
      if(input$sel1_left_I == "Genes"){
        textInput("sel2_left_I", label = "Gene to display", value = "Kap")
      }else if (input$sel1_left_I =="Cell info"){
        selectInput("sel2_left_I", label = "What cellinfo to display?",
                    choices = metas1,
                    selected = metas1[1])
      }
    })
    
    output$plot1 <- renderPlot({
      if(input$sel1_left_I == "Genes"){
      validate(
        need(input$sel2_left_I, "Type in your gene of interest"),
      )
      validate(
        need(any(input$sel2_left_I == rownames(sce$combined)), "Gene not found") #continue here
      )
      }
      Idents(sce$combined) <- "Cluster.Names"
      FeaturePlot(sce$combined, features = input$sel2_left_I, reduction = "umap")
    })
    
    output$sel1_right <- renderUI({
      selectInput("sel1_right_I", label = "Cell info or Genes?",
                  choices = list("Cell info", "Genes"),
                  selected = "Genes")
    })
    output$sel2_right <- renderUI({
      if(input$sel1_right_I == "Genes"){
        textInput("sel2_right_I", label = "Gene to display", value = "Kap")
      }else if (input$sel1_right_I =="Cell info"){
        selectInput("sel2_right_I", label = "What cellinfo to display?",
                    choices = metas1,
                    selected = metas1[1])
      }
    })
    
    output$plot2 <- renderPlot({
      if(input$sel1_right_I == "Genes"){
      validate(
        need(input$sel2_right_I, "Type in your gene of interest"),
      )
      validate(
        need(any(input$sel2_right_I== rownames(sce$combined)), "Gene not found") #continue here
      )
      }
      Idents(sce$combined) <- "Cluster.Names"
      FeaturePlot(sce$combined, features = input$sel2_right_I, reduction = "umap")
    })
    
    output$cellinfosel <- renderUI({
      selectInput("cellinfosel_I", label = "Cell info to plot (x axis)",
                  choices = list("type", "rep", "Cluster.Names"),
                  selected = "type")
    })
    
    output$groupsel <- renderUI({
      selectInput("groupsel_I", label = "Cell info to plot (x axis)",
                  choices = list("type", "rep", "Cluster.Names"),
                  selected = "Cluster.Names")
    })
    
    output$propplot <- renderPlot({
      data <- as.data.frame(sce$combined@meta.data)
      Idents(sce$combined) <- "Cluster.Names"
      ggplot(data, aes_string(x = input$cellinfosel_I, y = "nCount_RNA", fill = input$groupsel_I)) +
        geom_col(position = "fill") +
        scale_y_continuous(labels = scales::percent) +
        theme_classic()+
        scale_fill_viridis(discrete = TRUE)
    })
    
    output$deexpl <- renderText({
      types <- levels(factor(sce$combined$type))
      x1 <- types[1]
      x2 <- types[2]
      return(paste("Genes differentially expressed when comparing ", x1, " with ", x2,
            ". Positive SCT indicate stronger expression in ", x1, ".", sep=""))
    })
    
    output$DEtable_type <- DT::renderDataTable({
      #Determine differential expression (young vs. old): use "corrected counts" that are stored in the data slot of the SCT assay
      Idents(sce$combined) <- "type"
      types <- levels(factor(sce$combined$type))
      x1 <- types[1]
      x2 <- types[2]
      agemarkers <- FindMarkers(sce$combined, assay = "SCT", ident.1 = x1, ident.2 = x2, verbose = FALSE, future.seed=TRUE)
      agemarkers <<- subset(agemarkers, p_val < 0.05)
      agemarkers
    })
    
    output$genesel1 <- renderUI({
      textInput("genesel1_I", label = "Gene to highlight in scatterplot:")
    })
    
    output$scattertype <- renderPlot({
      Idents(sce$combined) <- "type"

      avexpr1 <<- as.data.frame(AverageExpression(sce$combined, assays = "SCT", verbose = FALSE))
      colnames(avexpr1) <- c("Condition.1", "Condition.2") #Condition1: young, Condition2: old
      avexpr <- data.frame(log1p(avexpr1$Condition.1), log1p(avexpr1$Condition.2))
      rownames(avexpr)<-rownames(avexpr1)
      colnames(avexpr)<-colnames(avexpr1)
      types <- levels(factor(sce$combined$type)) #c(old, young, ordered by alphabet)
      x1 <- types[1] #old
      x2 <- types[2] #young
      avexpr$genenames <- rownames(avexpr)
      avexpr <<- avexpr
      markerfeatures <- rownames(agemarkers)[1:15]
      rm(avexpr1)
      highlightdf <- subset.data.frame(avexpr, genenames == input$genesel1_I)
      plot1 <- ggplot(avexpr, aes(Condition.1, Condition.2)) + #young, old
        geom_point(color = "lightblue") +
        geom_point(data=highlightdf, aes(x=Condition.1,y=Condition.2), color='red', size=3) +
        ggtitle("DE between the two conditions")+
        theme_classic() +
        labs(x=x2, y=x1) #young, old
      LabelPoints(plot = plot1, points = markerfeatures, repel = TRUE)
      
    })
    
    output$x_value <- renderText({
      if (is.null(input$scattertype_click$x)) return("")
      else {
        
        #avexpr$genenames <- rownames(avexpr)
        genename <- nearPoints(avexpr, input$scattertype_click) #does not really work so far
        HTML("The selected gene is", genename)
      }
    })
    
    output$Dimplot2 <- renderPlot({
      DimPlot(sce$combined, pt.size=1, reduction="umap",
              split.by = "type")
    })
    
  } #QC  
 #Server 
  

# Run the app ----
shinyApp(ui = ui, server = server)
 