#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

library(BiocManager)
options(repos = BiocManager::repositories())

library(scater)
library(scran)
library(edgeR)
#library(pheatmap)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)


shinyServer(function(input, output, session) {

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
        
        m <- read.table(mtx(), skip =2, nrows = 1)
        
        tags$div(m$V1, "X", m$V2, " Matrix")
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
        
        ###########################################################
        ## Following Code modifed from read10X in edgeR package
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
        ## End edgeR code
        ########################################################
        
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
        
        hist(log(rowSums(dge.filtered()$counts)), xlab = "Log counts", ylab = "Frequency", main = "Gene Log Counts Frequency")
        
    })
    
    sce <- reactive({
        validate(need(dge.filtered(), message = FALSE))
                 
        SingleCellExperiment(list(counts=dge.filtered()$counts))
    })
    
    clus <- eventReactive(input$cluster_1, {
        validate(need(input$clus, message = FALSE),
                 need(sce(), message = FALSE))
        
        set.seed(input$seed)
        
        quickCluster(sce(), method=input$clus, min.mean=0.5)
    })
    
    output$normUI <- renderPlot({
        validate(need(clus(), message = FALSE))
        
        plot(clus(), xlab = "clusters", ylab = "Frequency")
    }) 
    
    sce.n <- eventReactive(input$normalize,{
        s <- computeSumFactors(sce(), cluster=clus())
        s <-logNormCounts(s)
        
        colData(s) <- cbind(colData(s), log10LibSize=log10(libSize()))
        
        s
    })
    
    libSize <- reactive({
        dge()$samples$lib.size
    })
    
    output$libsize <- renderPlot({
        validate(need(sce.n(), message = FALSE),
                 need(libSize(), message = FALSE))
        
        plot(libSize()/1000, sizeFactors(sce.n()), log="xy", pch=16, cex=0.7, xlab="library size (x1000)", ylab="Size fact")
    })
    
    var <- reactive({
        validate(need(sce.n(), message = FALSE))
        var <- modelGeneVar(sce.n())
        
        var<-var[order(var$bio, decreasing = TRUE),]
        var<- as.data.frame(var)
        var$ID <- seq.int(nrow(var))
        var
    })
        
    
    output$varTable <- renderDT({
        validate(need(var(), message = FALSE))
        
        var()[1:6]
    })
    
    output$means <- renderPlot({
        validate(need(var(), message = FALSE))
        v <- var()
        
        s <- input$varTable_rows_selected
        
        fit <- fitTrendVar(v$mean, v$total)
        
        plot(v$mean, v$total,pch=16,cex=0.7)
        curve(fit$trend(x), col="blue", add=TRUE, lwd=2)
        if (length(s)) points(v[s,]$mean, v[s,]$total, col = "red", pch = 19, cex = 2)
    })
    
    observeEvent(input$plot_click,{
        validate(need(var(), message = FALSE))
        v <- var()
        
        np <- nearPoints(v, input$plot_click, xvar = "mean", yvar = "total")
        
        if(nrow(np)){
            dtp <- dataTableProxy(
                "varTable",
                session = shiny::getDefaultReactiveDomain(),
                deferUntilFlush = TRUE
            )
            s <- input$varTable_rows_selected
            if(length(s)){
                if(!(np$ID %in% s))
                    s <- c(s, np$ID)
                else
                    s <- s[s != np$ID]
            } else {
                s <- c(np$ID)
            }
            
            selectRows(dtp, s)
        }
    })
    
    observeEvent(input$select,{
        t <- input$top
        
        dtp <- dataTableProxy(
            "varTable",
            session = shiny::getDefaultReactiveDomain(),
            deferUntilFlush = TRUE
        )

        
        if(input$select){
            selectRows(dtp, 1:t)
        } else {
            selectRows(dtp, NULL)    
        }
    })
    
    observeEvent(input$top,{
        t <- input$top
        
        dtp <- dataTableProxy(
            "varTable",
            session = shiny::getDefaultReactiveDomain(),
            deferUntilFlush = TRUE
        )
        
        if(input$select){
            selectRows(dtp, 1:t)
        }
    })
    
    var.e <-reactive({
        validate(need(var(), message = FALSE))
        
        s <- input$varTable_rows_selected
        
        if(length(s))
            v<-var()[s,]
        v
    })
    
    output$expression <- renderPlot({
        validate(need(sce.n(), message = FALSE))
        sce <- sce.n()
        
        s <- input$expTable_rows_selected
        if(length(s))
        {
            plotExpression(sce, features=rownames(var.e())[s])
        }
    })
    
    output$expTable <- renderDT({
        validate(need(var.e(), message = FALSE))
        
        var.e()[1:6]
    })
    
    output$colTable <- renderDT({
        validate(need(var.e(), message = FALSE))
        
        var.e()[1:6]
    }, selection = 'single')
    
    output$pca <- renderPlot({
        validate(
            need(sce.p(), message = FALSE),
            need(libSize(), message = FALSE)
        )
        
        sce <- sce.p()

        ncol(reducedDim(sce, "PCA"))
        
        n <- input$pcas
        
        colData(sce) <- cbind(colData(sce), log10LibSize=log10(libSize()))
        plotPCASCE(sce, colour_by = "log10LibSize", ncomponents = n)    
    })
    
    sce.p <- eventReactive(input$pca, {
        validate(
            need(sce.n(), message = FALSE),
            need(var(), message = FALSE),
            need(hvg(), message = FALSE)
        )
        
        sce <- denoisePCA(sce.n(), subset.row = hvg() ,technical=var()$tech, min.rank=10, max.rank=30)
        
        sce
    })
    
    sce.t <- eventReactive(input$tsne,{
        validate(need(sce.n(), message = FALSE))
        
        set.seed(input$seed)
        sce <- runTSNE(sce.n(), use_dimred="PCA", perplexity=100, theta=0)
        
        sce
    })
    
    output$tsne <- renderPlot({
        validate(need(sce.t(), message = FALSE))
        
        s <- input$colTable_rows_selected
        
        if(length(s)){
            plotTSNE(sce.t(), colour_by=rownames(var.e())[s])
        } else {
            plotTSNE(sce.t(), colour_by="log10LibSize")
        }
    })
    
    
    hvg <- reactive({
        validate(need(var(), message = FALSE))
        
        v<-var()
        
        
        if(input$subset_type == 1){
            s <- NULL
        } else if(input$subset_type == 2){
            s <- v[input$varTable_rows_selected,]
        } else {
            s <- v[1:input$hvg,]
        }
        
        rownames(s)
    })
    
    cluster <- eventReactive(input$cluster_2, {
        validate(need(sce.p(), message = FALSE),
                 need(hvg(), message = FALSE))
        

        if(input$cluster_type == "Hierarchical"){
            h <- quickCluster(sce.p(), subset.row=hvg(), assay.type="logcounts", method="hclust", min.size=input$min_c, min.mean=0.5)
            c <- factor(h)
        } else {
            snn <- buildSNNGraph(sce.p(), k=input$knn, use.dimred="PCA")
            w <- igraph::cluster_walktrap(snn)
            c <- factor(w$membership)
        }
        
        c
    })
    
    sce.c <- eventReactive(input$cluster_2, {
        validate(need(sce.t(), message = FALSE),
            need(cluster(), message = FALSE)         
        )
        
        s <- sce.t()
        
        s$Cluster <- cluster()
        
        s
    })
    
    output$clusterPlot <- renderPlot({
        validate(need(cluster(), message = FALSE))
        
        plot(cluster(), xlab = "clusters", ylab = "Frequency")
    })
    
    output$clusterTSNE <- renderPlot({
        validate(need(sce.c(), message = FALSE))
        
        plotTSNE(sce.c(), colour_by="Cluster")  
    })
    
    markers <- reactive({
        validate(need(sce.c(), message = FALSE))
        
        markers <- findMarkers(sce.c(), group=sce.c()$Cluster, direction="up")
        
        markers
    })
    
    gset <- reactive({
        validate(need(markers(), message = FALSE))
        
        top10 <- lapply(markers(), "[", 1:10, )
        g <- as.vector(sapply(top10, "rownames"))
        g <- g[!duplicated(g)]
        
        
        g
    })
    
    heatmap <- reactive({
        validate(need(sce.c(), message = FALSE),
                 need(gset(), message = FALSE))
        
        s <- sce.c()[gset()]

        Heatmap(logcounts(s), column_split = s$Cluster, show_row_dend = FALSE, show_column_dend = FALSE, cluster_columns = FALSE)
    })
    
    observeEvent(input$show_heatmap,{
        InteractiveComplexHeatmapWidget(input, output, session, heatmap(),
                                        output_id = "heatmap_output")   
    })
    
    observeEvent(input$next_1, {
        validate(need(dge(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Filter")
    })
    
    observeEvent(input$next_2, {
        validate(need(dge.filtered(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Normalization")
    })
    
    observeEvent(input$back_2, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Files")
    })
    
    observeEvent(input$next_3, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Variance")
    })
    
    observeEvent(input$back_3, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Filter")
    })
    
    observeEvent(input$next_4, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Expression")
    })
    
    observeEvent(input$back_4, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Normalization")
    })
    
    observeEvent(input$next_5, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "PCA")
    })
    
    observeEvent(input$back_5, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Variance")
    })
    
    observeEvent(input$next_6, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "TSNE")
    })
    
    observeEvent(input$back_6, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Expression")
    })
    
    observeEvent(input$next_7, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Cluster")
    })
    
    observeEvent(input$back_7, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "PCA")
    })
    
    observeEvent(input$next_8, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Markers")
    })
    
    observeEvent(input$back_8, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "TSNE")
    })
    
    observeEvent(input$next_9, {
        validate(need(sce.n(), message = FALSE))
        
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Export")
    })
    
    observeEvent(input$back_9, {
        updateTabItems(shiny::getDefaultReactiveDomain(), "tabs", selected = "Cluster")
    })
})
