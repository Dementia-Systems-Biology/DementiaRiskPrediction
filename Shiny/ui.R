ui <- tagList(
  
  tags$head(tags$style(HTML("
                           .navbar-nav {
                           float: none !important;
                           }
                           .navbar-nav > li:nth-child(4) {
                           float: right;
                           }
                            .navbar-nav > li:nth-child(3) {
                           float: right;
                           }
                           .my_style_1{ 
                             background-image: url(Background1.jpg);
                           }
                           .my_style_1 { margin-top: -20px; }
                           
                           .my_style_1 { width: 100%; }
                           
                           .container-fluid { padding-left: 0; padding-right: 0; }
                           
                           .my_style_1 { position: absolute; top: 0; }
                           
                           "))),
  
  
  fluidPage(
    
    navbarPage(title = "", id = "navbar",
               
               
               ###################################################################
               
               #  Calculate risk score: Home tab    
               
               ###################################################################
               tabPanel("Home", 
                        value = "input_panel", 
                        icon = icon("fas fa-home"), class = "my_style_1",
                        
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        
                        #******************************************************#
                        # Title and welcome message
                        #******************************************************#
                        fluidRow(
                          column(6, offset = 3,
                                 align = "center", 
                                 style = "background-color:#F5F5F5;",
                                 
                                 br(),
                                 br(),
                                 h1(strong(span(style = "color:#000000", "Multivariate methylation risk score for cognitive impairment and dementia"))),
                                 hr(),
                                 h5(span(style = "color:#000000", "Get started by uploading a .csv file with beta- or M-values.")),
                                 
                                 br()
                          )
                        ),
                        
                        #******************************************************#
                        # File upload
                        #******************************************************#
                        fluidRow(
                          column(6, offset = 3, 
                                 align = "center", 
                                 style = "background-color:#F5F5F5;",
                                 
                                 fileInput(inputId = "uploadcsv",
                                           label = NULL,
                                           accept = ".csv",
                                           placeholder = "Select .csv data file")
                                 
                          )
                          
                        ),
                        
                        #******************************************************#
                        # Start analysis
                        #******************************************************#
                        fluidRow(
                          
                          column(6, offset = 3, 
                                 align = "center", 
                                 style = "background-color:#F5F5F5;",
                                 
                                 
                                 #Use example data
                                 actionBttn(inputId = "example", 
                                            label = "Example",
                                            style = "jelly",
                                            color = "royal",
                                            size = "md",
                                            icon = icon("refresh")),
                                 
                                 #Start the analysis
                                 actionBttn(inputId = "startAnalysis",
                                            label = "Start",
                                            style = "jelly",
                                            color = "primary",
                                            size = "md",
                                            icon = icon("arrow-right")),
                                 
                                 
                                 
                                 br(),
                                 br(),
                                 br(),
                                 br(),
                                 
                          )
                        ),
                        
                      
                        
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br(),
                        br()
                        
                        
                        
               ), # End of home tab
               
               
               
               ###################################################################
               
               #  Outputs
               
               ###################################################################
               
               tabPanel("Output", value = "output_panel", 
                        icon = icon("fas fa-layer-group"),
                        
                        navlistPanel(id = "tabs_output",
                                     
                                     # Output table
                                     tabPanel("Methylation Profile Scores (MPSs)", value = "MRS_table",
                                              h1(strong("Methylation Profile Scores")),
                                              hr(),
                                              dataTableOutput("predictedScore_factors"),
                                              downloadButton(outputId ="download_predictedScore_factors", 
                                                             label = "Download MPS Table")),
                                     
                                     # Output table
                                     tabPanel("Multivariate Methylation Risk Score (MMRS)", value = "epiMCI_table",
                                              h1(strong("Multivariate Methylation Risk Score")),
                                              hr(),
                                              dataTableOutput("predictDF"),
                                              downloadButton(outputId = "download_predictDF", 
                                                             label = "Download MMRS Table"))

                                     
                                     
                                    
                                     
                        ) # navlist panel
               ) # outputs
               

               
    ) #navbarpage
  ) # fluidpage
) # taglist
