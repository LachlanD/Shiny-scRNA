#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(DT)
library(shinycssloaders)



# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 500*1024^2)
shinyUI(
    dashboardPage(
        dashboardHeader(title = "scRNA GUI"),

        dashboardSidebar(
            sidebarMenu( id = "tabs",
                
                         menuItem("Files", tabName = "Files", icon = icon("file-upload"), selected = TRUE),
                         menuItem("Filter",tabName = "Filter", icon = icon("filter")),
                         menuItem("Normalization", tabName = "Normalization", icon = icon("chart-line")),
                         menuItem("Variance", tabName = "Variance", icon = icon("chart-bar")),
                         menuItem("Expression", tabName = "Expression", icon = icon("feather")),
                         menuItem("Principal Component Analysis", tabName = "PCA", icon = icon("ruler-combined")),
                         menuItem("t-SNE Analysis", tabName = "TSNE", icon = icon("project-diagram")),
                         menuItem("Cluster", tabName = "Cluster", icon = icon("code-branch")),
                         menuItem("Markers", tabName = "Markers", icon = icon("flag")),
                         menuItem("Export", tabName = "Export", icon = icon("file-download")),
                         menuItem("Random Seed", tabName = "Seed", icon = icon("dice")),
                         menuItem("About", tabName = "About", icon = icon("info-circle"))
            )
        ),
        dashboardBody(
            tabItems(
                tabItem(tabName = "Files",
                        fluidRow(
                            column(3,
                                   fileInput("mtx", "Upload mtx file", multiple = FALSE, accept = c(".mtx")),
                                   uiOutput("mtxUI")
                            ),
                            column(3,
                                   fileInput("genes", "Upload genes file", multiple = FALSE, accept = c(".tsv", ".gz")),
                                   uiOutput("genesUI")
                            ),
                            column(3,
                                   fileInput("barcodes", "Upload barcodes file", multiple = FALSE, accept = c(".tsv", ".gz")),
                                   uiOutput("barcodesUI")
                            ),
                            column(3,
                                   fileInput("gene_info", "Upload gene annotation file", multiple = FALSE, accept = c(".gz")),
                                   uiOutput("annUI")
                            )
                        ),
                        fluidRow(
                             withSpinner(plotOutput("dgePlot"))
                        ),
                        fluidRow(
                            column(2,
                                actionButton("load", "Load"),
                            offset = 10)
                        ),
                        fluidRow(
                            column(2,
                                   actionButton("next_1", "Next", icon = icon("forward")), 
                            offset = 10)
                        )
                ),
                tabItem(tabName = "Filter",
                        checkboxInput("f1", label = "Remove genes with all zero counts", value = FALSE),
                        checkboxInput("f2", label = "Remove genes with non-zero counts in fewer than x percent of cells", value = FALSE),
                        numericInput("f2_value", "Minimum percentage of cells with non-zero counts", value = 0.01, min = 0.0, max = 1.0, step = 0.001),
                        checkboxInput("f3", label = "Remove unAnnotated genes", value = TRUE),
                        checkboxInput("f4", label = "Remove duplicates", value = TRUE),
                        uiOutput("filterUI"),
                        plotOutput("filterPlot"),
                        fluidRow(
                            column(2,actionButton("back_2", "Back", icon = icon("backward")), offset = 8),
                            column(2,actionButton("next_2", "Next", icon = icon("forward")))
                        )
                ),
                tabItem(tabName = "Normalization",
                        fluidRow(
                            column(2,
                                radioButtons("clus", "cluster method", choices = c("igraph", "hclust")),
                                actionButton("cluster_1", "Cluster")
                            ),
                            column(10,
                                withSpinner(plotOutput("normUI"))
                            )
                        ),
                    
                        withSpinner(plotOutput("libsize")),
                        actionButton("normalize", "Compute Normalization Factors"),
                        fluidRow(
                            column(2,actionButton("back_3", "Back", icon = icon("backward")), offset = 8),
                            column(2,actionButton("next_3", "Next", icon = icon("forward")))
                        )
                ),
                tabItem(tabName = "Variance",
                        fluidRow(
                            column(width =4,
                                checkboxInput("select", "select/unselect", width = '150px')
                            ),
                            column(width =2, tags$b("top")),
                            column(width = 6,
                                numericInput("top", label = NULL, value = 20, min = 0, step =1, width = '100px')
                            )
                        ),
                        DTOutput("varTable"),
                        plotOutput("means", click = "plot_click", brush = "brush"),
                        fluidRow(
                            column(2,actionButton("back_4", "Back", icon = icon("backward")), offset = 8),
                            column(2,actionButton("next_4", "Next", icon = icon("forward")))
                        )
                ),
                tabItem(tabName = "Expression",
                    plotOutput("expression"),
                    DTOutput("expTable"),
                    fluidRow(
                        column(2,actionButton("back_5", "Back", icon = icon("backward")), offset = 8),
                        column(2,actionButton("next_5", "Next", icon = icon("forward")))
                    )
                ),
                tabItem(tabName = "PCA",
                    fluidRow(    
                        column(4,
                            numericInput("pcas", label = "Number of components", value =3, min = 1, max = 5, step = 1)
                        ),
                        column(4,
                            radioButtons("subset_type", label = "use genes", choiceValues = c(1,2,3), choiceNames = c("all", "selected", "highest variance"), selected = 3)
                        ),
                        column(4,
                            numericInput("hvg", "no. genes", value = 10, min = 5, step = 1, width = '150px')
                        )
                    ),
                    
                    withSpinner(plotOutput("pca")),
                    fluidRow(column(2, actionButton("pca", "Calculate Principal Components"), offset =8)),
                    fluidRow(
                        column(2,actionButton("back_6", "Back", icon = icon("backward")), offset = 8),
                        column(2,actionButton("next_6", "Next", icon = icon("forward")))
                    )        
                ),
                tabItem(tabName = "TSNE",
                    withSpinner(plotOutput("tsne")),
                    fluidRow(column(2, actionButton("tsne", "Plot t-SNE"), offset =8)),
                    DTOutput("colTable"),
                    fluidRow(
                        column(2,actionButton("back_7", "Back", icon = icon("backward")), offset = 8),
                        column(2,actionButton("next_7", "Next", icon = icon("forward")))
                    )
                ),
                tabItem("Cluster",
                    fluidRow(
                        column(4,
                            radioButtons("cluster_type", label = "Cluster", choices = c("Hierarchical", "k-Nearest Neighbour")),
                        ),
                        column(4,
                            numericInput("knn", "k", value = 60, min = 10, step = 1, width = "150px"),
                            numericInput("min_c", "min. heirarchial cluster size", value = 50, min = 10, step = 1, width = "150px") 
                        )
                    ),
                    withSpinner(plotOutput("clusterTSNE")),
                    fluidRow(column(4, actionButton("cluster_2", "Cluster"), offset = 8)),
                    fluidRow(
                        column(2,actionButton("back_8", "Back", icon = icon("backward")), offset = 8),
                        column(2,actionButton("next_8", "Next", icon = icon("forward")))
                    )
                ),
                tabItem(tabName = "Markers",
                    htmlOutput("heatmap_output"),
                    actionButton("show_heatmap", "Show/Refresh Heatmap"),
                    fluidRow(
                        column(2,actionButton("back_9", "Back", icon = icon("backward")), offset = 8),
                        column(2,actionButton("next_9", "Next", icon = icon("forward")))
                    )    
                ),
                tabItem(tabName = "Export",
                    tags$div("export")      
                ),
                tabItem(tabName = "Seed",
                        numericInput("seed", label = "random seed", value = 100)
                )
            )
    )
))
