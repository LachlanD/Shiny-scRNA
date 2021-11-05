# Shiny-scRNA
Siny/R app which provides a graphical interface for exploratory data analysis of single cell RNA-Seq data using R tools such as Bioconductor, Scater, Scran, EdgeR, Interactive Complex Heatamps & DataTables.

[EdgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
[interactive Complex Heatmaps](https://bioconductor.org/packages/release/bioc/html/InteractiveComplexHeatmap.html)
[Bioconductor](https://bioconductor.org/)
[Scran](https://bioconductor.org/packages/release/bioc/html/scran.html)
[Scater](https://bioconductor.org/packages/release/bioc/html/scater.html)
[DataTables for R](https://rstudio.github.io/DT/)


## Launching App

To start the app download or clone the ui.R and server.R files and open them with RStudio the click the run app button.

![Launch](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/run.PNG?raw=true)

The app can also be hosted on a server and allowing user to run it through a browser based client.  (Free servers on shinyapps.io have instufficent memory however) 

## Files
Upload the three files from the cell ranger director (Matrix file, genes file, barcodes file) and a gene annotations file for the organsism. 

![Files](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/files.PNG?raw=true)

After uploading the file inspecting the library size distributions should show if there was a major problem in the data.

## Filter
The next stage we filter the data. Removing duplicated and unAnnotated genes is necessary for following data analysis.  Removing genes where only a small percentage of cells have non-zero counts is recommended.  Adjusting the minimum percentage should yield an appropriate distribution of counts.

![Filter](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/filter.PNG?raw=true)

## Normalization
The next step is to normailze counts with respect each cells library size.  Since we expect biological as well as techinical differences in library size we first perform a basic clustering of cells of similar expression profiles and only normalize with respect to cells in the same cluster.


![cluster](https://user-images.githubusercontent.com/5520490/140042369-72bb47fe-b108-4e13-8491-2e4e30bbbbf0.png)

## Variance

Scran computes the technical and biological variance of each gene.  Using the table or the graph we can select the genes we wish to include in further analysis.

![varaince1](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/variance1.gif?raw=true)


![variance2](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/variance2.gif?raw=true)

## Expression

Plot the exprssion of the genes by selecting them from the table.

![expression](https://github.com/LachlanD/Shiny-scRNA/blob/main/img/expression.gif?raw=true)




