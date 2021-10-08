#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(scater)
library(scran)
library(edgeR)
library(Matrix)


shinyServer(function(input, output) {

    mtx <- reactive({
        validate(need(input$mtx, message = FALSE))
        
        input$mtx$datapath
    })
    
    genes <- reactive({
        validate(need(input$genes, message = FALSE))
        
        
        input$genes$datapath
    })
    
    barcodes <- reactive({
        validate(need(input$barcodes, message = FALSE))
        
        input$barcodes$datapath
    })
    
    gene_info <- reactive({
        validate(need(input$gene_info, message = FALSE))
        
        input$gene_info$datapath    
    })
    
    output$genesUI <- renderUI({
        validate(need(genes(), message = FALSE))
        
        g <- read.delim(genes(), header = FALSE)
        
        tags$div(nrow(g), " genes")
    })
    
    output$barcodesUI <- renderUI({
        validate(need(barcodes(), message = FALSE))
        
        b <- read.delim(barcodes(), header = FALSE)
        
        tags$div(nrow(b), " cells")
    })
    
    output$mtxUI <- renderUI({
        validate(need(mtx(), message = FALSE))
        
        m <- readMM(mtx())
        
        tags$div(nrow(m), "X", ncol(m), " Matrix")
    })
    
    output$annUI <- renderUI({
        validate(need(gene_info(), message = FALSE))
        
        a <- read.delim(gene_info())
        
        tags$div(nrow(a), " gene annotations")
    })

    dge <- eventReactive(input$load,{
        validate(need(mtx(), message = FALSE),
                 need(genes(), message = FALSE),
                 need(barcodes(), message = FALSE))
        
        N <- scan(mtx(),skip=2,what=0L,sep=" ",nmax=3,quiet=TRUE)
        ngenes <- N[1]
        ncells <- N[2]
        nmtx <- N[3]
        
        # Read gene Ids
        Genes <- read.table(genes(),header=FALSE,comment.char="",sep="\t",row.names=1,colClasses="character")
        if(nrow(Genes) != ngenes) stop("Number of feature IDs doesn't agree with header information in mtx file")
        names(Genes)[1] <- "Symbol"
        if(ncol(Genes) > 1L) names(Genes)[2] <- "Type"
        
        
        #	Read mtx file of counts
        m <- read.table(mtx(),skip=3,header=FALSE,comment.char="",sep=" ",colClasses="integer",nrows=nmtx)
        
        #	Convert Market Exchange Format to ordinary matrix
        y <- matrix(0L,ngenes,ncells)
        i <- m[,1]+(m[,2]-1L)*ngenes
        y[i] <- m[,3]
        dimnames(y) <- list(Gene=row.names(Genes),Cell=1:ncells)
        
        #	Optionally read barcodes
        if(is.null(barcodes)) {
            Samples <- NULL
        } else {
            Barcodes <- scan(barcodes(),what="",quiet=TRUE)
            if(length(Barcodes) != ncells) stop("Number of barcodes doesn't agree with header information in mtx file")
            Samples <- data.frame(Barcode=Barcodes)
        }
        
        
        dge <- DGEList(count=y,genes=Genes,samples=Samples)
        
        

        ann <- alias2SymbolUsingNCBI(dge$genes$Symbol, required.columns = c("GeneID", "Symbol"), gene.info.file = gene_info()) 
        
        dge$genes <- cbind(dge$genes, Official=ann$Symbol, GeneID=ann$GeneID)
        
        dimnames(dge$counts)<-list(dge$genes$Official, colnames(dge$counts))
        
        mito<- grep("^mt-", dge$genes$Symbol)
        
        percent.mito <- colSums(dge$counts[mito, ])/dge$samples$lib.size
        
        nGenes <- colSums(dge$counts != 0)
        
        dge$samples <- cbind(dge$samples, percent.mito=percent.mito, nGenes=nGenes)
        
        o <- order(rowSums(dge$counts), decreasing = TRUE)
        
        dge <- dge[o,]
        
        dge
    })
    
    dge.filtered <- reactive({
        validate(need(dge(), message = FALSE))
        
        dge.filtered<-dge()
        
        
        if(input$f1){
            f1<- rowSums(dge.filtered$counts > 0) > 0
            dge.filtered <- dge.filtered[f1,]
        }
        if(input$f2){
            f2<- rowSums(dge.filtered$counts > 0) >= ncol(dge.filtered)*input$f2_value
            dge.filtered <- dge.filtered[f2,]    
        }
        if(input$f3){
            f3 <- !is.na(dge.filtered$genes$Official)
            dge.filtered <- dge.filtered[f3,]
        }
        if(input$f4){
            f4 <- !duplicated(dge.filtered$genes$Official)
            dge.filtered <- dge.filtered[f4,]
        }
        
        dge.filtered
    })
    
    output$dgePlot <- renderPlot({
        validate(need(dge(), message = FALSE))
        
        d<-dge()
        
        
        par(mfrow=c(1,2))
        plot(d$samples[,c("lib.size","nGenes")], pch=16, cex=0.7)
        plot(d$samples[,c("lib.size","percent.mito")], pch=16, cex=0.7)
        
    })
    
    output$filterUI <- renderUI({
        validate(need(dge(), message = FALSE),
                 need(dge.filtered(), message = FALSE))
        
        tags$div(nrow(dge.filtered()), " out of ", nrow(dge()))
        
    })
    
    output$filterPlot <- renderPlot({
        validate(need(dge.filtered(), message = FALSE))
        
        hist(log(rowSums(dge.filtered()$counts)))
        
    })
    
    sce <- reactive({
        validate(need(dge.filtered(), message = FALSE))
                 
        SingleCellExperiment(list(counts=dge.filtered()$counts))
    })
    
    clus <- reactive({
        validate(need(input$clus, message = FALSE),
                 need(sce(), message = FALSE))
        
        quickCluster(sce(), method=input$clus, min.mean=0.5)
    })
    
    output$normUI <- renderTable({
        table(clus())
    }) 
    
    sce.c <- reactive({
        sce.c <- computeSumFactors(sce(), cluster=clst)
    })
    
    output$libsize <- renderPlot({
        libSize <- dge()$samples$lib.size
        plot(libSize/1000, sizeFactors(sce.c()), log="xy", pch=16, cex=0.7, xlab="library size (x1000)", ylab="Size fact")
    })
})
