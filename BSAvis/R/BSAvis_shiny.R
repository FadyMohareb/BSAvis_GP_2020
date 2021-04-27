#Run BSAvis - interactive BSA visualisation tool

#' @title Run BSAvis 
#' @description Run the interactive BSA visualisation tool.
#'
#' @param vcf.list object containing meta information and vcf data frame
#' 
#' @importFrom dplyr %>%
#' @export BSAvis_shiny
#'
#' @examples 
#' BSAvis_shiny(vcf.list=vcf_list)


BSAvis_shiny <- function(vcf.list){

  #User Interface
  ui <-
    fluidPage(useShinyalert(),
      
      # ==============================================================================
      #                                     HOME
      # ==============================================================================
      
      fluidRow(
      id = "body",
      column(width = 12, hr(),
        h1("BSAvis", align="center"),
        br(),
        #HTML('<center><img src="www/BSAvis_logo.png" width="200"></center>'),
        #shiny::img(src="www/BSAvis_small_logo.png"),
        h5("Welcome to BSAvis, an interactive Bulk Segregant Analysis (BSA) tool.", 
        align ="center")
      ),
      column(width = 12, hr())
    ),
      
      # ==============================================================================
      # FIRST METHOD
      #                               SNP-index Method
      # ==============================================================================
      
      #Title
      titlePanel(strong("BSAvis")),
      tags$h5("SNP-index Method"),
      br(),
      
      sidebarLayout(
        sidebarPanel(
          
          #Select Chromosome
          selectInput(
            inputId = "chrMethod1",
            label = "Chromosome Number:",
            choices = (BSAvis::extract_chrIDs(vcf.list$meta))[!grepl("0$", BSAvis::extract_chrIDs(vcf.list$meta))] #remove chromosome ID if ends with 0 (chromosome 0)
          ),
          
          hr(),
          
          #Select pool names
          selectInput(
            inputId = "wtMethod1",
            label = "Wild-type bulk:",
            choices = unique((vcf.list$df)$Indiv),
            selected = unique((vcf.list$df)$Indiv)[grep("*minus", unique((vcf.list$df)$Indiv))]
          ),
          
          selectInput(
            inputId = "mMethod1",
            label = "Mutant bulk:",
            choices = unique((vcf.list$df)$Indiv),
            selected = unique((vcf.list$df)$Indiv)[grep("*plus", unique((vcf.list$df)$Indiv))]
          ),
          
          hr(),
          
          radioButtons(
            inputId = "bulk",
            label = "Bulk:",
            choices = list("Both" = 0,
                           "Wild-type bulk" = 1,
                           "Mutant bulk" = 2),
            selected = 0
          ),
          
          radioButtons(
            inputId = "variants_snpindex",
            label = "Variants:",
            choices = list("SNPs" = 1,
                           "SNPs+InDels" = 2),
            
            selected = 1
          ),
          
          hr(),
          
          #Set Window Size
          sliderInput(
            inputId = "window",
            label = "Window Size:",
            min = 500000,
            max = 3000000,
            value = 1000000, 
            step = 500000
          ),
          
          #Set Step Size
          sliderInput(
            inputId = "step",
            label = "Step Size:",
            min = 5000,
            max = 10000,
            value = 10000,
            step = 1000
          ),
          
          hr(),
          
          #Filtering
          
          #min SNP-index
          numericInput(
            inputId = "min_snpindex",
            label = "Min SNP-index:",
            min = 0,
            max = 1,
            step = 0.1,
            value = 0.3
          ),
          
          #max SNP-index
          numericInput(
            inputId = "max_snpindex",
            label = "Max SNP-index:",
            min = 0,
            max = 1,
            step = 0.1,
            value = 0.9
          ),
          
          #min DP
          numericInput(
            inputId = "min_DP1",
            label = "Min DP:",
            min = 0,
            max = 80,
            step = 1,
            value = 50
          ),
          
          #max DP
          numericInput(
            inputId = "max_DP1",
            label = "Max DP:",
            min = 0,
            max = 500,
            step = 1,
            value = 200
          ),
          
          #min GQ
          numericInput(
            inputId = "min_GQ",
            label = "Min GQ:",
            min = 0,
            max = 100,
            step = 1,
            value = 99
          )
        ),
        
        #SNP-index - plot panels
        mainPanel(tabsetPanel(type = "tabs",
                              tabPanel("SNP-index",
                                       #loading spinner before generating the plot
                                       withSpinner(plotOutput("snp_index",
                                                              dblclick = "plot1_dblclick",  #add double-click event
                                                              brush = brushOpts(            #enable brush
                                                                id = "plot1_brush",
                                                                resetOnNew = TRUE))),
                                       #add save button
                                       actionButton("save1", "Save Plot")),
                              
                              tabPanel("delta(SNP-index)",
                                       #loading spinner before generating the plot
                                       withSpinner(plotOutput("delta_snp_index",
                                                              dblclick = "plot2_dblclick",  #add double-click event
                                                              brush = brushOpts(            #enable brush
                                                                id = "plot2_brush",
                                                                resetOnNew = TRUE))),
                                       #add save button
                                       actionButton("save2", "Save Plot"))
          ))
      ),
      
      # ==============================================================================
      # SECOND METHOD
      #                               SNP-ratio Method
      # ==============================================================================
      
      titlePanel(strong("BSAvis")),
      tags$h5("SNP-ratio Method"),
      br(),
      
      sidebarLayout(
        sidebarPanel(
          
          #Select Chromosome
          selectInput(
            inputId = "chrMethod2",
            label = "Chromosome Number:",
            choices = (BSAvis::extract_chrIDs(vcf.list$meta))[!grepl("0$", BSAvis::extract_chrIDs(vcf.list$meta))]
          ),
          
          hr(),
          
          #Select pool names
          selectInput(
            inputId = "wtMethod2",
            label = "Wild-type bulk:",
            choices = unique((vcf.list$df)$Indiv),
            selected = unique((vcf.list$df)$Indiv)[grep("*minus", unique((vcf.list$df)$Indiv))]
          ),
          
          selectInput(
            inputId = "mMethod2",
            label = "Mutant bulk:",
            choices = unique((vcf.list$df)$Indiv),
            selected = unique((vcf.list$df)$Indiv)[grep("*plus", unique((vcf.list$df)$Indiv))]
          ),
          
          hr(),
          
          radioButtons(
            inputId = "variants_snpratio",
            label = "Variants:",
            choices = list("SNPs" = 1,
                           "SNPs+InDels" = 2),
            
            selected = 1
          ),
          
          hr(),
          
          #Min SNP-ratio value to consider
          numericInput(
            inputId = "min_snpratio",
            label = "Min SNP-ratio:",
            min = 0,
            step = 0.1,
            value = 0.1
          ),
          
          #min DP
          numericInput(
            inputId = "min_DP2",
            label = "Min DP:",
            min = 0,
            max = 80,
            step = 1,
            value = 50
          ),
          
          #max DP
          numericInput(
            inputId = "max_DP2",
            label = "Max DP:",
            min = 0,
            max = 500,
            step = 1,
            value = 200
          ),
          
          hr(),
          
          #LOESS Smoothing - Degree
          radioButtons(
            inputId = "degree",
            label = "Degree:",
            choices = list("0" = 0,
                           "1" = 1,
                           "2" = 2),
            selected = 2
          ),
          
          #LOESS Smoothing - Span
          numericInput(
            inputId = "span",
            label = "Span:",
            min = 0,
            value = 0.07,
            step = 0.01
          )
        ),
        
        #SNP-ratio - plot panel
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("SNP-ratio",
                     #loading spinner before generating the plot
                     withSpinner(plotOutput("snp_ratio",
                                 dblclick = "plot3_dblclick",   #add double-click event
                                 brush = brushOpts(             #enable brush
                                   id = "plot3_brush",
                                   resetOnNew = TRUE))),
                     #add action button
                     actionButton("save3", "Save Plot"))
        ))
      )
    )
  
  
  #Server
  server <- function(input, output) {
    
    #Create objects to store brush bounds for each of the plots
    ranges1 <- reactiveValues(x = NULL, y = NULL)
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    ranges3 <- reactiveValues(x = NULL, y = NULL)
    
    # ==============================================================================
    # FIRST METHOD
    #                               SNP-index Method
    # ==============================================================================
    
    # -------------------------------- SNP-Index ----------------------------------- 
    #Create reactive expression to plot SNP-index
    SNPindexPlot <- reactive({

      #Assign value to chosenVariants based on user input regarding the type of variants to show
      chosenVariants <- ifelse(input$variants_snpindex == 1, "SNP", "all")
      
      #Assign value to chosenBulk based on user input regarding which bulk's SNP-index plot to visualise
      if (input$bulk==0) {
        chosenBulk <- 0
      }
      
      else if (input$bulk==1) {
        chosenBulk <- 1
      }
      
      else if (input$bulk==2) {
        chosenBulk <- 2
      }
      
      #Call SNP-index plot function
      BSAvis::shiny_SNPindex(vcf.list = vcf.list, 
                             wtBulk = input$wtMethod1, 
                             mBulk = input$mMethod1, 
                             variants = chosenVariants,
                             min.SNPindex = input$min_snpindex, 
                             max.SNPindex = input$max_snpindex, 
                             min.DP = input$min_DP1, 
                             max.DP = input$max_DP1, 
                             min.GQ = input$min_GQ,
                             chrID = input$chrMethod1, 
                             chr = input$chrMethod1, 
                             windowSize = as.numeric(input$window), 
                             windowStep = as.numeric(input$step), 
                             bulk = chosenBulk, 
                             ranges = ranges1)
    })
    
    #Add SNP-index plot to its corresponding place in the UI
    output$snp_index <- renderPlot({
      SNPindexPlot()
    })
    
    # -------------------------------- delta(SNP-Index) -----------------------------------
    #Create reactive expression to plot delta(SNP-index)
    deltaSNPindexPlot <- reactive({

 

      chosenVariants <- ifelse(input$variants_snpindex == 1, "SNP", "all")
      
      BSAvis::shiny_deltaSNPindex(vcf.list = vcf.list, 
                                  wtBulk = input$wtMethod1, 
                                  mBulk = input$mMethod1, 
                                  variants = chosenVariants,
                                  min.SNPindex = input$min_snpindex, 
                                  max.SNPindex = input$max_snpindex, 
                                  min.DP = input$min_DP1, 
                                  max.DP = input$max_DP1, 
                                  min.GQ = input$min_GQ,
                                  chrID = input$chrMethod1, 
                                  chr = input$chrMethod1, 
                                  windowSize = as.numeric(input$window), 
                                  windowStep = as.numeric(input$step), 
                                  ranges = ranges2)
    })
    
    #Add delta(SNP-index) plot to its corresponding place in the UI
    output$delta_snp_index <- renderPlot({
      deltaSNPindexPlot()
    })
    
    # ==============================================================================
    # SECOND METHOD
    #                               SNP-ratio PLOT
    # ==============================================================================  
    
    #Create reactive expression to plot SNP-ratio
    SNPratioPlot <- reactive({

      chosenVariants <- ifelse(input$variants_snpratio == 1, "SNP", "all")
      
      BSAvis::shiny_SNPratio(vcf.list = vcf.list, 
                             wtBulk = input$wtMethod2, 
                             mBulk = input$mMethod2, 
                             variants = chosenVariants,
                             min.SNPratio = input$min_snpratio, 
                             min.DP = input$min_DP2, 
                             max.DP = input$max_DP2,
                             chrID = input$chrMethod2, 
                             chr = input$chrMethod2, 
                             degree = input$degree, 
                             span = input$span, 
                             ranges = ranges3)
    })
    
    #Add SNP-ratio plot to its corresponding place in the UI
    output$snp_ratio <- renderPlot({
      SNPratioPlot()
    })
    
    # ==============================================================================
    # 
    #                                    Save plots
    # ==============================================================================  
    #Save SNP-index plot when the corresponding button is clicked. 
    #First, show pop-up window allowing the user to set saving parameters
    observeEvent(input$save1, {
      shinyalert(inputId="savingInfoAlert1", html = TRUE, text = tagList(
        helpText(h6("Default values are being shown. Please, type to customise the parameters.")),
        br(),
        textInput("directory1", "Directory to save the plot", getwd()),
        textInput("filename1", "File name (without file extesion)", paste0("plot_SNPindex_ch", input$chrMethod1)),
        numericInput("dpi1", "DPI", 1200),
        numericInput("width1", "Width (inches)", 7.5),
        numericInput("height1", "Height (inches)", 5)
      ))
    })
    
    observeEvent(input$savingInfoAlert1, {
      ggplot2::ggsave(filename=paste0(input$filename1, ".tiff"), path=input$directory1, plot=SNPindexPlot(), device = "tiff", dpi = input$dpi1, width = input$width1, height = input$height1)
    })
    
    
    #Save delta(SNP-index) plot when the corresponding button is clicked. 
    #First, show pop-up window allowing the user to set saving parameters
    observeEvent(input$save2, {
      shinyalert(inputId="savingInfoAlert2", html = TRUE, text = tagList(
        helpText(h6("Default values are being shown. Please, type to customise the parameters.")),
        br(),
        textInput("directory2", "Directory to save the plot", getwd()),
        textInput("filename2", "File name (without file extesion)", paste0("plot_deltaSNPindex_ch", input$chrMethod1)),
        numericInput("dpi2", "DPI", 1200),
        numericInput("width2", "Width (inches)", 7.5),
        numericInput("height2", "Height (inches)", 5)
      ))
    })
    
    observeEvent(input$savingInfoAlert2, {
      ggplot2::ggsave(filename=paste0(input$filename2, ".tiff"), path=input$directory2, plot=deltaSNPindexPlot(), device = "tiff", dpi = input$dpi2, width = input$width2, height = input$height2)
    })
    
    
    #Save SNP-ratio plot when the corresponding button is clicked. 
    #First, show pop-up window allowing the user to set saving parameters
    observeEvent(input$save3, {
      shinyalert(inputId="savingInfoAlert3", html = TRUE, text = tagList(
        helpText(h6("Default values are being shown. Please, type to customise the parameters.")),
        br(),
        textInput("directory3", "Directory to save the plot", getwd()),
        textInput("filename3", "File name (without file extesion)", paste0("plot_SNPratio_ch", input$chrMethod2)),
        numericInput("dpi3", "DPI", 1200),
        numericInput("width3", "Width (inches)", 7.5),
        numericInput("height3", "Height (inches)", 5)
      ))
    })
    
    observeEvent(input$savingInfoAlert3, {
      ggplot2::ggsave(filename=paste0(input$filename3, ".tiff"), path=input$directory3, plot=SNPratioPlot(), device = "tiff", dpi = input$dpi3, width = input$width3, height = input$height3)
    })
    
    # ==============================================================================
    # 
    #                                 Zoom functionality
    # ==============================================================================  
    #For each of the plots, check if there is a brush when a double-click happens. If so, zoom to the brush bounds; if not, reset the zoom.
    
    #Implement double-click to zoom in SNP-index plot
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      if (!is.null(brush)) {
        ranges1$x <- c(brush$xmin, brush$xmax)
        ranges1$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges1$x <- NULL
        ranges1$y <- NULL
      }
    })
    
    #Implement double-click to zoom in delta(SNP-index) plot
    observeEvent(input$plot2_dblclick, {
      brush <- input$plot2_brush
      if (!is.null(brush)) {
        ranges2$x <- c(brush$xmin, brush$xmax)
        ranges2$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges2$x <- NULL
        ranges2$y <- NULL
      }
    })
    
    #Implement double-click to zoom in SNP-ratio plot
    observeEvent(input$plot3_dblclick, {
      brush <- input$plot3_brush
      if (!is.null(brush)) {
        ranges3$x <- c(brush$xmin, brush$xmax)
        ranges3$y <- c(brush$ymin, brush$ymax)
        
      } else {
        ranges3$x <- NULL
        ranges3$y <- NULL
      }
    })
    
  }
  
  shinyApp(ui, server)
}
