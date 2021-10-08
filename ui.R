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

# Define UI for application that draws a histogram
options(shiny.maxRequestSize = 500*1024^2)
shinyUI(
    dashboardPage(
        dashboardHeader(title = "scRNA GUI"),

        dashboardSidebar(
            sidebarMenu( id = "tabs",
                
                         menuItem("Files", tabName = "Files", icon = icon("file-upload"), selected = TRUE),
                         menuItem("Filter",tabName = "Filter", icon = icon("chart-bar")),
                         menuItem("Normalization", tabName = "Normalization", icon = icon("search")),
                         menuItem("Export", tabName = "Export", icon = icon("file-download")),
                         menuItem("About", tabName = "About", icon = icon("info-circle"))
            )
        ),
        dashboardBody(
            tabItems(
                tabItem(tabName = "Files",
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
                               uiOutput("annUI"),
                               actionButton("load", "Load",)
                        ),
                        plotOutput("dgePlot")
                ),
                tabItem(tabName = "Filter",
                        checkboxInput("f1", label = "Remove genes with all zero counts", value = FALSE),
                        checkboxInput("f2", label = "Remove genes with non-zero counts in fewer than x percent of cells", value = FALSE),
                        numericInput("f2_value", "Minimum percentage of cells with non-zero counts", value = 0.01, min = 0.0, max = 1.0, step = 0.001),
                        checkboxInput("f3", label = "Remove unAnnotated genes", value = TRUE),
                        checkboxInput("f4", label = "Remove duplicates", value = TRUE),
                        uiOutput("filterUI"),
                        plotOutput("filterPlot")
                ),
                tabItem(tabName = "Normalization",
                        radioButtons("clus", "cluster method", choices = c("igraph", "hclust")),
                        tableOutput("normUI"),
                        plotOutput("libsize")
                        
                )
                        
            )
                               
    )
))
