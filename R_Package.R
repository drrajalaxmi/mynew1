
library(gplots)
library(ggplot2)
library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyBS)
library(DT)
library(clusterProfiler)
library(DESeq2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(edgeR)
library(RColorBrewer)
library(XML)
library(annotate)
library(stringr)
library(reshape2)
library(data.table)
library(plotly)
library(heatmaply)
library(DOSE)
library(Harman)
library(dplyr)
library(tidyverse)
library(tidyr)

#' debrowserdeanalysis
#'
#' Module to perform and visualize DE results.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @param columns, columns
#' @param conds, conditions
#' @param params, de parameters
#' @return DE panel
#' @export
#'
#' @examples
#'     x <- debrowserdeanalysis()
#'
debrowserdeanalysis <- function(input = NULL, output = NULL, session = NULL,
                                data = NULL, columns = NULL, conds = NULL, params = NULL) {
  if(is.null(data)) return(NULL)
  deres <- reactive({
    runDE(data, columns, conds, params)
  })
  prepDat <- reactive({
    applyFiltersNew(addDataCols(data, deres(), columns, conds), input)
  })
  observe({
    dat <-  prepDat()[prepDat()$Legend == input$legendradio,]
    dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat)
    getTableDetails(output, session, "DEResults", dat2, modal=FALSE)
  })
  list(dat = prepDat)
}
#' getDEResultsUI
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getDEResultsUI("batcheffect")
#'
#' @export
#'
getDEResultsUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "DE Results",
                          solidHeader = T, status = "info",  width = 12,
                          fluidRow(
                            column(12,
                                   uiOutput(ns("DEResults"))
                            ),
                            actionButtonDE("goMain", "Go to Main Plots", styleclass = "primary")
                          )
      )
    )
  )
}

#' cutOffSelectionUI
#'
#' Gathers the cut off selection for DE analysis
#'
#' @param id, namespace id
#' @note \code{cutOffSelectionUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- cutOffSelectionUI("cutoff")
#' @export
#'
cutOffSelectionUI <- function(id){
  ns <- NS(id)
  list(
    getLegendRadio(id),
    textInput(ns("padj"), "padj value cut off", value = "0.01" ),
    textInput(ns("foldChange"), "or foldChange", value = "2" )
  )
}

#' applyFiltersNew
#'
#' Apply filters based on foldChange cutoff and padj value.
#' This function adds a "Legend" column with "Up", "Down" or
#' "NS" values for visualization.
#'
#' @param data, loaded dataset
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- applyFiltersNew()
#'
applyFiltersNew <- function(data = NULL, input = NULL) {
  if (is.null(data)) return(NULL)
  padj_cutoff <- as.numeric(input$padj)
  foldChange_cutoff <- as.numeric(input$foldChange)
  m <- data
  if (!("Legend" %in% names(m))) {
    m$Legend <- character(nrow(m))
    m$Legend <- "NS"
  }
  m$Legend[m$foldChange >= foldChange_cutoff &
             m$padj <= padj_cutoff] <- "Up"
  m$Legend[m$foldChange <= (1 / foldChange_cutoff) &
             m$padj <= padj_cutoff] <- "Down"
  return(m)
}


#' runDE
#'
#' Run DE algorithms on the selected parameters.  Output is
#' to be used for the interactive display.
#'
#' @param data, A matrix that includes all the expression raw counts,
#'     rownames has to be the gene, isoform or region names/IDs
#' @param columns, is a vector that includes the columns that are going
#'     to be analyzed. These columns has to match with the given data.
#' @param conds, experimental conditions. The order has to match
#'     with the column order
#' @param params, all params for the DE methods
#' @return de results
#'
#' @export
#'
#' @examples
#'     x <- runDE()
#'
runDE <- function(data = NULL, columns = NULL, conds = NULL, params = NULL) {
  if (is.null(data)) return(NULL)
  de_res <- NULL

  if (startsWith(params[1], "DESeq2"))
    de_res <- runDESeq2(data, columns, conds, params)
  else if (startsWith(params[1], "EdgeR"))
    de_res <- runEdgeR(data, columns, conds, params)
  else if (startsWith(params[1], "Limma"))
    de_res <- runLimma(data, columns, conds, params)
  data.frame(de_res)
}

#' runDESeq2
#'
#' Run DESeq2 algorithm on the selected conditions.  Output is
#' to be used for the interactive display.
#'
#' @param data, A matrix that includes all the expression raw counts,
#'     rownames has to be the gene, isoform or region names/IDs
#' @param columns, is a vector that includes the columns that are going
#'     to be analyzed. These columns has to match with the given data.
#' @param conds, experimental conditions. The order has to match
#'     with the column order
#' @param params, fitType: either "parametric", "local", or "mean" for the type
#'     of fitting of dispersions to the mean intensity.
#'     See estimateDispersions for description.
#'  betaPrior: whether or not to put a zero-mean normal prior
#'     on the non-intercept coefficients See nbinomWaldTest for
#'     description of the calculation of the beta prior. By default,
#'     the beta prior is used only for the Wald test, but can also be
#'     specified for the likelihood ratio test.
#' testType: either "Wald" or "LRT", which will then use either
#'     Wald significance tests (defined by nbinomWaldTest), or the
#'     likelihood ratio test on the difference in deviance between a
#'     full and reduced model formula (defined by nbinomLRT)
#' rowsum.filter: regions/genes/isoforms with total count
#'      (across all samples) below this value will be filtered out
#' @return deseq2 results
#'
#' @export
#'
#' @examples
#'     x <- runDESeq2()
#'
runDESeq2 <- function(data = NULL, columns = NULL, conds = NULL, params = NULL) {
  if (is.null(data)) return(NULL)
  if (length(params)<3)
    params <- strsplit(params, ",")[[1]]

  fitType <- if (!is.null(params[2])) params[2]
  betaPrior <-  if (!is.null(params[3])) params[3]
  testType <- if (!is.null(params[4])) params[4]
  rowsum.filter <-  if (!is.null(params[5])) as.integer(params[5])

  data <- data[, columns]

  data[, columns] <- apply(data[, columns], 2,
                           function(x) as.integer(x))

  coldata <- prepGroup(conds, columns)
  # Filtering non expressed genes
  filtd <- data
  if (is.numeric(rowsum.filter) && !is.na(rowsum.filter))
    filtd <- subset(data, rowSums(data) > rowsum.filter)

  # DESeq data structure is going to be prepared
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(filtd),
                                colData = coldata, design = ~group)
  # Running DESeq
  if (testType == "LRT")
    dds <- DESeq(dds, fitType = fitType, betaPrior = as.logical(betaPrior), test=testType, reduced= ~ 1)
  else
    dds <- DESeq(dds, fitType = fitType, betaPrior = as.logical(betaPrior), test=testType)

  res <- results(dds)
  return(res)
}

#' runEdgeR
#'
#' Run EdgeR algorithm on the selected conditions.  Output is
#' to be used for the interactive display.
#'
#' @param data, A matrix that includes all the expression raw counts,
#'     rownames has to be the gene, isoform or region names/IDs
#' @param columns, is a vector that includes the columns that are going
#'     to be analyzed. These columns has to match with the given data.
#' @param conds, experimental conditions. The order has to match
#'     with the column order
#' @param params, normfact: Calculate normalization factors to scale the raw
#'     library sizes. Values can be "TMM","RLE","upperquartile","none".
#' dispersion: either a numeric vector of dispersions or a character
#'     string indicating that dispersions should be taken from the data
#'     object. If a numeric vector, then can be either of length one or
#'     of length equal to the number of genes. Allowable character
#'     values are "common", "trended", "tagwise" or "auto".
#'     Default behavior ("auto" is to use most complex dispersions
#'     found in data object.
#' testType: exactTest or glmLRT. exactTest: Computes p-values for differential
#'     abundance for each gene between two digital libraries, conditioning
#'     on the total count for each gene. The counts in each group as a
#'     proportion of the whole are assumed to follow a binomial distribution.
#'     glmLRT: Fit a negative binomial generalized log-linear model to the read
#'     counts for each gene. Conduct genewise statistical tests for a given
#'     coefficient or coefficient contrast.
#' rowsum.filter: regions/genes/isoforms with total count
#'      (across all samples) below this value will be filtered out
#' @return edgeR results
#'
#' @export
#'
#' @examples
#' x <- runEdgeR()
#'
runEdgeR<- function(data = NULL, columns = NULL, conds = NULL, params = NULL){
  if (is.null(data)) return(NULL)
  if (length(params)<3)
    params <- strsplit(params, ",")[[1]]
  normfact <- if (!is.null(params[2])) params[2]
  dispersion <- if (!is.null(params[3])) params[3]
  testType <- if (!is.null(params[4])) params[4]
  rowsum.filter <-  if (!is.null(params[5])) as.integer(params[5])

  data <- data[, columns]
  data[, columns] <- apply(data[, columns], 2,
                           function(x) as.integer(x))
  dispersion <- as.numeric(dispersion)
  conds <- factor(conds)
  filtd <- data
  if (is.numeric(rowsum.filter) && !is.na(rowsum.filter))
    filtd <- subset(data, rowSums(data) > rowsum.filter)

  d<- edgeR::DGEList(counts = filtd, group = conds)
  d <- edgeR::calcNormFactors(d, method = normfact)
  # If dispersion is 0, it will estimate the dispersions.
  de.com <- c()
  cnum = summary(conds)[levels(conds)[1]]
  tnum = summary(conds)[levels(conds)[2]]
  des <- c(rep(1, cnum),rep(2, tnum))
  design <- model.matrix(~des)
  if (testType == "exactTest"){
    if (dispersion == 0){
      d <- edgeR::estimateDisp(d, design)
      de.com <- edgeR::exactTest(d)
    }else{
      de.com <- edgeR::exactTest(d, dispersion=dispersion)
    }
  }else if (testType == "glmLRT"){
    if (dispersion == 0){
      d <- edgeR::estimateDisp(d, design)
      fit <- edgeR::glmFit(d, design)
    }else{
      fit <- edgeR::glmFit(d, design, dispersion=dispersion)
    }
    de.com <- edgeR::glmLRT(fit,coef=2)
  }

  options(digits=4)

  padj<- p.adjust(de.com$table$PValue, method="bonferroni")
  res <-data.frame(cbind(de.com$table$logFC/log(2),de.com$table$PValue, padj))
  colnames(res) <- c("log2FoldChange", "pvalue", "padj")
  rownames(res) <- rownames(filtd)
  return(res)
}

#' runLimma
#'
#' Run Limma algorithm on the selected conditions.  Output is
#' to be used for the interactive display.
#'
#' @param data, A matrix that includes all the expression raw counts,
#'     rownames has to be the gene, isoform or region names/IDs
#' @param columns, is a vector that includes the columns that are going
#'     to be analyzed. These columns has to match with the given data.
#' @param conds, experimental conditions. The order has to match
#'     with the column order
#' @param params, normfact: Calculate normalization factors to scale the raw
#'     library sizes. Values can be "TMM","RLE","upperquartile","none".
#' fitType, fitting method; "ls" for least squares or "robust"
#'     for robust regression
#' normBet: Normalizes expression intensities so that the
#'     intensities or log-ratios have similar distributions across a set of arrays.
#' rowsum.filter: regions/genes/isoforms with total count
#'     (across all samples) below this value will be filtered out
#' @return Limma results
#'
#' @export
#'
#' @examples
#'     x <- runLimma()
#'
runLimma<- function(data = NULL, columns = NULL, conds = NULL, params = NULL){
  if (is.null(data)) return(NULL)
  if (length(params)<3)
    params <- strsplit(params, ",")[[1]]
  normfact = if (!is.null(params[2])) params[2]
  fitType = if (!is.null(params[3])) params[3]
  normBet = if (!is.null(params[4])) params[4]
  rowsum.filter <-  if (!is.null(params[5])) as.integer(params[5])

  data <- data[, columns]
  data[, columns] <- apply(data[, columns], 2,
                           function(x) as.integer(x))
  conds <- factor(conds)

  cnum = summary(conds)[levels(conds)[1]]
  tnum = summary(conds)[levels(conds)[2]]
  filtd <- data
  if (is.numeric(rowsum.filter) && !is.na(rowsum.filter))
    filtd <- as.matrix(subset(data, rowSums(data) > rowsum.filter))

  des <- factor(c(rep(levels(conds)[1], cnum),rep(levels(conds)[2], tnum)))
  names(filtd) <- des
  design <- cbind(Grp1=1,Grp2vs1=des)

  dge <- DGEList(counts=filtd, group = des)

  dge <- calcNormFactors(dge, method=normfact, samples=columns)

  v <- voom(dge, design=design, normalize.method = normBet, plot=FALSE)

  fit <- lmFit(v, design=design)
  fit <- eBayes(fit)

  options(digits=4)
  tab <- topTable(fit,coef=2, number=dim(fit)[1],genelist=fit$genes$NAME)
  res <-data.frame(cbind(tab$logFC, tab$P.Value, tab$adj.P.Val))

  colnames(res) <- c("log2FoldChange", "pvalue", "padj")
  rownames(res) <- rownames(tab)
  return(res)
}

#' prepGroup
#'
#' prepare group table
#'
#' @param cols, columns
#' @param conds, inputconds
#' @return data
#' @export
#'
#' @examples
#'     x <- prepGroup()
#'
prepGroup <- function(conds = NULL, cols = NULL) {
  if (is.null(conds) || is.null(cols)) return (NULL)
  coldata <- data.frame(cbind(cols, conds))
  coldata$conds <- factor(coldata$conds)
  colnames(coldata) <- c("libname", "group")
  coldata
}

#' addDataCols
#'
#' add aditional data columns to de results
#'
#' @param data, loaded dataset
#' @param de_res, de results
#' @param cols, columns
#' @param conds, inputconds
#' @return data
#' @export
#'
#' @examples
#'     x <- addDataCols()
#'
addDataCols <- function(data = NULL, de_res = NULL, cols = NULL, conds = NULL) {
  if (is.null(data) || is.null(de_res)) return (NULL)
  norm_data <- data[, cols]

  coldata <- prepGroup(conds, cols)

  mean_cond_first <- getMean(norm_data, as.vector(coldata[coldata$group==levels(coldata$group)[1], "libname"]))
  mean_cond_second <- getMean(norm_data, as.vector(coldata[coldata$group==levels(coldata$group)[2], "libname"]))

  m <- cbind(rownames(de_res), norm_data[rownames(de_res), cols],
             log10(unlist(mean_cond_second) + 1),
             log10(unlist(mean_cond_first) + 1),
             de_res[rownames(de_res),
                    c("padj", "log2FoldChange", "pvalue")],
             2 ^ de_res[rownames(de_res),
                        "log2FoldChange"],
             -1 * log10(de_res[rownames(de_res), "padj"]))
  colnames(m) <- c("ID", cols, "x", "y",
                   "padj", "log2FoldChange", "pvalue",
                   "foldChange", "log10padj")
  m <- as.data.frame(m)
  m$padj[is.na(m[paste0("padj")])] <- 1
  m$pvalue[is.na(m[paste0("pvalue")])] <- 1
  m
}

#' getMean
#'
#' Gathers the mean for selected condition.
#'
#' @param data, dataset
#' @param selcols, input cols
#' @return data
#' @export
#'
#' @examples
#'     x <- getMean()
#'
getMean<-function(data = NULL, selcols=NULL) {
  if (is.null(data)) return (NULL)
  mean_cond<-NULL
  if (length(selcols) > 1)
    mean_cond <-list(rowMeans( data[, selcols]))
  else
    mean_cond <-list(norm_data[selcols])
  mean_cond
}

#' getLegendRadio
#'
#' Radio buttons for the types in the legend
#' @param id, namespace id
#' @note \code{getLegendRadio}
#' @return radio control
#'
#' @examples
#'
#'     x <- getLegendRadio("deprog")
#'
#' @export
#'
getLegendRadio <- function(id) {
  ns <- NS(id)
  types <- c("Up", "Down", "NS")
  radioButtons(inputId=ns("legendradio"),
               label="Data Type:",
               choices=types
  )
}




#' debrowserall2all
#'
#' Module for a bar plot that can be used in data prep, main plots
#' low count removal modules or any desired module
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @param cex, the size of the dots
#' @return all2all plot
#' @export
#'
#' @examples
#'     x <- debrowserall2all()
#'
debrowserall2all <- function(input, output, session, data = NULL,
                             cex=2) {
  if(is.null(data)) return(NULL)
  output$all2allplot <- renderPlot({
    all2all(data, cex)
  })
  output$all2allUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = "All2all plot", status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotOutput(session$ns("all2allplot"),
                                    width = input$width, height=input$height))
  })
}

#' getAll2AllPlotUI
#'
#' all2all plots UI.
#'
#' @note \code{getAll2AllPlotUI}
#' @param id, namespace id
#' @return the panel for all2all plots;
#'
#' @examples
#'     x <- getAll2AllPlotUI("bar")
#'
#' @export
#'
getAll2AllPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("all2allUI"))
}

#' all2allControlsUI
#'
#' Generates the controls in the left menu for an all2all plot
#'
#' @note \code{all2allControlsUI}
#' @param id, namespace id
#' @return returns the controls for left menu
#' @examples
#'     x <- all2allControlsUI("bar")
#' @export
#'
all2allControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(id, " - Options"),
                           sliderInput("cex", "corr font size",
                                       min = 0.1, max = 10,
                                       step = 0.1, value = 2)
  )
}

#' all2all
#'
#' Prepares all2all scatter plots for given datasets.
#'
#' @param data, data that have the sample names in the header.
#' @param cex text size
#' @return all2all scatter plots
#' @examples
#'     plot<-all2all(mtcars)
#'
#' @export
#'
all2all <- function(data, cex=2) {
  pcor <- function(x, y, ...) panel.cor(x, y, cex.cor = cex)
  nr <- nrow(data)
  if (nr > 1000)
    nr <- 1000
  pairs(log10(data[1:nr, ]), cex = 0.25,
        diag.panel = panel.hist, lower.panel = pcor)
}

#' panel.hist
#'
#' Prepares the historgram for the all2all plot.
#'
#' @param x, a vector of values for which the histogram is desired
#' @param ..., any additional params
#' @return all2all histogram plots
#' @examples
#'     panel.hist(1)
#'
#' @export
#'
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nb <- length(breaks)
  y <- h$counts
  y <- y / max(y)
  rect(breaks[-nb], 0, breaks[-1], y, col = "red", ...)
}

#' panel.cor
#'
#' Prepares the correlations for the all2all plot.
#'
#' @param x, numeric vector x
#' @param y, numeric vector y
#' @param prefix, prefix for the text
#' @param cex.cor, correlation font size
#' @param ..., additional parameters
#' @return all2all correlation plots
#' @examples
#'     panel.cor(c(1,2,3), c(4,5,6))
#'
#' @export
#'
panel.cor <- function(x, y, prefix = "rho=", cex.cor=2, ...){
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor.test(x, y, method = "spearman",
                na.rm = TRUE, exact = FALSE)$estimate
  txt <- round(r, digits = 2)
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}


#' debrowserbatcheffect
#'
#' Module to correct batch effect
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param ldata, loaded data
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowserbatcheffect()
#'
debrowserbatcheffect <- function(input, output, session, ldata = NULL) {
  if(is.null(ldata)) return(NULL)
  batchdata <- reactiveValues(count=NULL, meta = NULL)
  observeEvent(input$submitBatchEffect, {
    if (is.null(ldata$count)) return (NULL)

    countData <- ldata$count
    withProgress(message = 'Normalization', detail = "Normalization", value = NULL, {
      if (input$norm_method != "none"){
        countData <- getNormalizedMatrix(ldata$count, method=input$norm_method)
      }
    })
    withProgress(message = 'Batch Effect Correction', detail = "Adjusting the Data", value = NULL, {
      if (input$batchmethod == "Combat"){
        batchdata$count <- correctCombat(input, countData, ldata$meta)
      }
      else if (input$batchmethod == "Harman"){
        batchdata$count <- correctHarman(input, countData, ldata$meta)
      }
      else{
        batchdata$count <-  countData
      }
    })

    batchdata$meta <- ldata$meta
  })

  output$batchfields <- renderUI({
    if (!is.null(ldata$meta))
      list( conditionalPanel(condition = paste0("input['", session$ns("batchmethod"),"']!='none'"),
                             selectGroupInfo( ldata$meta, input, session$ns("treatment"), "Treatment"),
                             selectGroupInfo( ldata$meta, input, session$ns("batch"), "Batch")))
  })

  batcheffectdata <- reactive({
    ret <- NULL
    if(!is.null(batchdata$count)){
      ret <- batchdata
    }
    return(ret)
  })

  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", ldata)
    getSampleDetails(output, "filteredSummary", "filteredDetails", batcheffectdata())
    getTableDetails(output, session, "beforebatchtable", ldata$count, modal = TRUE)
    callModule(debrowserpcaplot, "beforeCorrectionPCA", ldata$count, ldata$meta)
    callModule(debrowserIQRplot, "beforeCorrectionIQR",  ldata$count)
    callModule(debrowserdensityplot, "beforeCorrectionDensity", ldata$count)
    if ( !is.null(batcheffectdata()$count ) && nrow(batcheffectdata()$count)>2 ){
      withProgress(message = 'Drawing the plot', detail = "Preparing!", value = NULL, {
        getTableDetails(output, session, "afterbatchtable", batcheffectdata()$count, modal = TRUE)
        callModule(debrowserpcaplot, "afterCorrectionPCA",  batcheffectdata()$count, batcheffectdata()$meta)
        callModule(debrowserIQRplot, "afterCorrectionIQR",  batcheffectdata()$count)
        callModule(debrowserdensityplot, "afterCorrectionDensity", batcheffectdata()$count)
      })
    }
  })

  list(BatchEffect=batcheffectdata)
}


#' batchEffectUI
#' Creates a panel to coorect batch effect
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- batchEffectUI("batcheffect")
#'
#' @export
#'
batchEffectUI <- function (id) {
  ns <- NS(id)

  list(
    fluidRow(
      shinydashboard::box(title = "Batch Effect Correction and Normalization",
                          solidHeader = TRUE, status = "info",  width = 12,
                          fluidRow(
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("uploadSummary")),
                                         DT::dataTableOutput(ns("sampleDetails"))),
                                   uiOutput(ns("beforebatchtable"))
                            ),
                            column(2,
                                   shinydashboard::box(title = "Options",
                                                       solidHeader = TRUE, status = "info",
                                                       width = 12,
                                                       normalizationMethods(id),
                                                       batchMethod(id),
                                                       uiOutput(ns("batchfields")),
                                                       actionButtonDE(ns("submitBatchEffect"), label = "Submit", styleclass = "primary")
                                   )
                            ),
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("filteredSummary")),
                                         DT::dataTableOutput(ns("filteredDetails"))),
                                   uiOutput(ns("afterbatchtable"))
                            )
                          ),
                          conditionalPanel(condition = paste0("input['", ns("submitBatchEffect"),"']"),
                                           actionButtonDE("goDE", "Go to DE Analysis", styleclass = "primary"),
                                           actionButtonDE("goQCplots", "Go to QC plots", styleclass = "primary"))),
      shinydashboard::box(title = "Plots",
                          solidHeader = TRUE, status = "info",  width = 12,
                          fluidRow(column(1, div()),
                                   tabsetPanel( id = ns("batchTabs"),
                                                tabPanel(id = ns("PCA"), "PCA",
                                                         column(5,
                                                                getPCAPlotUI(ns("beforeCorrectionPCA"))),
                                                         column(2,
                                                                shinydashboard::box(title = "PCA Controls",
                                                                                    solidHeader = T, status = "info",  width = 12,
                                                                                    tabsetPanel( id = ns("pcacontrols"),
                                                                                                 tabPanel ("Before",
                                                                                                           pcaPlotControlsUI(ns("beforeCorrectionPCA"))),
                                                                                                 tabPanel ( "After",
                                                                                                            pcaPlotControlsUI(ns("afterCorrectionPCA")))))),
                                                         column(5,
                                                                getPCAPlotUI(ns("afterCorrectionPCA")))
                                                ),
                                                tabPanel(id = ns("IQR"), "IQR",
                                                         column(5,
                                                                getIQRPlotUI(ns("beforeCorrectionIQR"))),
                                                         column(2, div()),
                                                         column(5,
                                                                getIQRPlotUI(ns("afterCorrectionIQR")))
                                                ),
                                                tabPanel(id = ns("Density"), "Density",
                                                         column(5,
                                                                getDensityPlotUI(ns("beforeCorrectionDensity"))),
                                                         column(2, div()),
                                                         column(5,
                                                                getDensityPlotUI(ns("afterCorrectionDensity")))
                                                )
                                   )
                          )
      )
    ), getPCAcontolUpdatesJS())
}
#' normalizationMethods
#'
#' Select box to select normalization method prior to batch effect correction
#'
#' @note \code{normalizationMethods}
#' @param id, namespace id
#' @return radio control
#'
#' @examples
#'
#'     x <- normalizationMethods("batch")
#'
#' @export
#'
normalizationMethods <- function(id) {
  ns <- NS(id)
  selectInput(ns("norm_method"), "Normalization Method:",
              choices = c("none", "MRN", "TMM", "RLE", "upperquartile"))
}

#' batchMethod
#'
#' select batch effect method
#' @param id, namespace id
#' @note \code{batchMethod}
#' @return radio control
#'
#' @examples
#'
#'     x <- batchMethod("batch")
#'
#' @export
#'
batchMethod <- function(id) {
  ns <- NS(id)
  selectInput(ns("batchmethod"), "Correction Method:",
              choices = c("none", "Combat", "Harman"),
              selected='none'
  )
}

#' Correct Batch Effect using Combat in sva package
#'
#' Batch effect correction
#' @param input, input values
#' @param idata, data
#' @param metadata, metadata
#' @return data
#' @export
#'
#' @examples
#'     x<-correctCombat ()
correctCombat <- function (input = NULL, idata = NULL, metadata = NULL) {
  if (is.null(idata) || input$batch == "None") return(NULL)
  batch <- metadata[, input$batch]
  treatment <- metadata[, input$treatment]
  columns <- colnames(idata)
  meta <- data.frame(cbind(columns, treatment, batch))
  datacor <- data.frame(idata[, columns])
  datacor[, columns] <- apply(datacor[, columns], 2,
                              function(x) as.integer(x))

  datacor[, columns] <- apply(datacor[, columns], 2,
                              function(x) return(x + runif(1, 0, 0.01)))

  modcombat = model.matrix(~1, data = meta)

  combat_blind = sva::ComBat(dat=as.matrix(datacor), batch=batch)

  a <- cbind(idata[rownames(combat_blind), 2], combat_blind)

  a[, columns] <- apply(a[, columns], 2, function(x) ifelse(x<0, 0, x))
  colnames(a[, 1]) <- colnames(idata[, 1])
  a[,columns]
}

#' Correct Batch Effect using Harman
#'
#' Batch effect correction
#' @param input, input values
#' @param idata, data
#' @param metadata, metadata
#' @return data
#' @export
#'
#' @examples
#'     x<-correctHarman ()
correctHarman <- function (input = NULL, idata = NULL, metadata = NULL) {
  if (is.null(idata)) return(NULL)
  batch.info <- data.frame(metadata[, c(input$treatment, input$batch)])
  rownames(batch.info) <- rownames(metadata)
  colnames(batch.info) <- c("treatment", "batch")

  harman.res <- harman(idata, expt= batch.info$treatment, batch= batch.info$batch, limit=0.95)
  harman.corrected <- reconstructData(harman.res)
  harman.corrected[harman.corrected<0] <- 0
  harman.corrected
}


#' getGeneList
#'
#' Gathers the gene list to use for GOTerm analysis.
#"
#' @note \code{GOTerm}
#'
#' @export
#'
#' @note \code{getGeneList}
#' symobol to ENTREZ ID conversion
#' @param genes, gene list
#' @param org, orgranism for gene symbol entrez ID conversion
#' @return ENTREZ ID list
#'
#' @examples
#'     x <- getGeneList(c('OCLN', 'ABCC2'))
#'
getGeneList <- function(genes = NULL, org = "org.Hs.eg.db") {
  # Get the entrez gene identifiers that are mapped to a gene symbol
  if (!installpack(org)) return(NULL)
  allkeys <- AnnotationDbi::keys(eval(parse(text = org)),
                                 keytype="SYMBOL")
  existinggenes <- unique(as.vector(unlist(lapply(toupper(genes),
                                                  function(x){ allkeys[x == toupper(allkeys)] }))))
  if (length(existinggenes) < 1){
    txt <- paste0("Please check the gene names! DEBrowser only accepts gene symbols. Ex:",
                  paste(head(allkeys), sep=","))
    print(txt)
    showNotification(txt)
    return(NULL)
  }

  mapped_genes <- mapIds(eval(parse(text = org)), keys = existinggenes,
                         column="ENTREZID", keytype="SYMBOL",
                         multiVals = "first")
  genelist <- unique(as.vector(unlist(mapped_genes)))
  genelist
}

#' getEntrezTable
#'
#' Gathers the entrezIds of the genes in given list and their data
#"
#' @note \code{GOTerm}
#'
#' @export
#'
#' @note \code{getEntrezTable}
#' symobol to ENTREZ ID conversion
#' @param genes, gene list
#' @param dat, data matrix
#' @param org, orgranism for gene symbol entrez ID conversion
#' @return table with the entrez IDs in the rownames
#'
#' @examples
#'     x <- getEntrezTable()
#'
getEntrezTable <- function(genes = NULL, dat = NULL, org = "org.Hs.eg.db") {
  if (is.null(genes)) return(NULL)
  if (!installpack(org)) return(NULL)
  allkeys <- AnnotationDbi::keys(eval(parse(text = org)),
                                 keytype="SYMBOL")
  entrezIDs <- unlist(strsplit(genes, "/"))

  mapped_genes <- mapIds(eval(parse(text = org)), keys = rownames(dat),
                         column="ENTREZID", keytype="SYMBOL",
                         multiVals = "first")
  mapped_genes <- mapped_genes[mapped_genes %in% entrezIDs]
  genelist <- cbind(mapped_genes, dat[names(mapped_genes), ])

  genelist <- data.frame(genelist)
}

#' getEntrezIds
#'
#' Gathers the gene list to use for GOTerm analysis.
#"
#' @note \code{GOTerm}
#'
#' @export
#'
#' @note \code{getEntrezIds}
#' symobol to ENTREZ ID conversion
#' @param genes, gene list with fold changes
#' @param org, orgranism for gene symbol entrez ID conversion
#' @return ENTREZ ID list
#'
#' @examples
#'     x <- getEntrezIds()
#'
getEntrezIds <- function(genes = NULL, org = "org.Hs.eg.db") {
  if (is.null(genes)) return(NULL)
  if (!installpack(org)) return(NULL)
  allkeys <- AnnotationDbi::keys(eval(parse(text = org)),
                                 keytype="SYMBOL")

  mapped_genes <- mapIds(eval(parse(text = org)), keys = rownames(genes),
                         column="ENTREZID", keytype="SYMBOL",
                         multiVals = "first")
  mapped_genes <- mapped_genes[!is.na(mapped_genes)]
  genelist <- cbind(mapped_genes, names(mapped_genes), genes[names(mapped_genes), "log2FoldChange"])

  colnames(genelist) <- c("ENTREZID", "SYMBOL", "log2FoldChange")
  genelist <- data.frame(genelist)
  genelist$log2FoldChange <- as.numeric(as.character(genelist$log2FoldChange))
  rownames(genelist) <- genelist$ENTREZID
  genelist
}

#' getEnrichGO
#'
#' Gathers the Enriched GO Term analysis data to be used within the
#' GO Term plots.
#'
#' @note \code{getEnrichGO}
#' @param genelist, gene list
#' @param pvalueCutoff, p value cutoff
#' @param org, the organism used
#' @param ont, the ontology used
#' @return Enriched GO
#' @examples
#'     x <- getEnrichGO()
#'
#' @export
#'
getEnrichGO <- function(genelist = NULL, pvalueCutoff = 0.01,
                        org = "org.Hs.eg.db", ont="CC") {
  if (is.null(genelist)) return(NULL)
  if (!installpack(org)) return(NULL)
  res <- c()
  res$enrich_p <- clusterProfiler::enrichGO(gene = genelist, OrgDb = org,
                                            # res$enrich_p <- enrichGO(gene = genelist, organism = "human",
                                            ont = ont, pvalueCutoff = pvalueCutoff)

  res$p <- barplot(res$enrich_p, title = paste("Enrich GO", ont),
                   font.size = 12)
  res$table <- NULL
  if (!is.null(nrow(res$enrich_p@result)) )
    res$table <- res$enrich_p@result[,c("ID", "Description",
                                        "GeneRatio", "pvalue", "p.adjust", "qvalue")]
  return(res)
}

#' getEnrichKEGG
#'
#' Gathers the Enriched KEGG analysis data to be used within the
#' GO Term plots.
#'
#' @note \code{getEnrichKEGG}
#' @param genelist, gene list
#' @param org, the organism used
#' @param pvalueCutoff, the p value cutoff
#' @return Enriched KEGG
#' @examples
#'     x <- getEnrichKEGG()
#' @export
#'
getEnrichKEGG <- function(genelist = NULL, pvalueCutoff = 0.01,
                          org = "org.Hs.eg.db") {
  if (is.null(genelist)) return(NULL)
  res <- c()
  res$enrich_p <- enrichKEGG(gene = genelist, organism = getOrganism(org),
                             pvalueCutoff = pvalueCutoff)
  res$p <- barplot(res$enrich_p, title =
                     paste("KEGG Enrichment: p-value=", pvalueCutoff))
  res$table <- NULL
  if (!is.null(nrow(res$enrich_p@result)) )
    res$table <- res$enrich_p@result[,c("ID", "Description",
                                        "GeneRatio", "pvalue", "p.adjust", "qvalue")]
  return(res)
}

#' clusterData
#'
#' Gathers the Cluster analysis data to be used within the
#' GO Term plots.
#'
#' @note \code{clusterData}
#' @param dat, the data to cluster
#' @return clustered data
#' @examples
#'     mycluster <- clusterData()
#'
#' @export
#'
clusterData <- function(dat = NULL) {
  if (is.null(dat)) return(NULL)
  ret <- list()
  itemlabels <- rownames(dat)
  norm_data <- getNormalizedMatrix(dat)
  mydata <- na.omit(norm_data)  # listwise deletion of missing
  mydata <- scale(mydata)  # standardize variables

  wss <- (nrow(mydata) - 1) * sum(apply(mydata, 2, var))
  for (i in 2:15) wss[i] <- sum(kmeans(mydata, centers = i)$withinss)
  plot(1:15, wss, type = "b",
       xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
  k <- 0
  for (i in 1:14) {
    if ( ( wss[i] / wss[i + 1] ) > 1.2 ) {
      k <- k + 1
    }
  }
  # K-Means Cluster Analysis
  fit <- kmeans(mydata, k)  # 5 cluster solution
  # get cluster means
  aggregate(mydata, by = list(fit$cluster), FUN = mean)
  # append cluster assignment
  mydata_cluster <- data.frame(mydata, fit$cluster)

  # distance <- dist(mydata, method = 'euclidean')
  # distance matrix fit <- hclust(distance,
  # method='ward.D2') plot(fit, cex = 0.1)
  # display dendogram groups <- cutree(fit, k=k) rect.hclust(fit,
  # k=k, border='red')
  return(mydata_cluster)
}

#' compareClust
#'
#' Compares the clustered data to be displayed within the GO Term
#' plots.
#'
#' @note \code{compareClust}
#' @param dat, data to compare clusters
#' @param ont, the ontology to use
#' @param org, the organism used
#' @param fun, fun
#' @param title, title of the comparison
#' @param pvalueCutoff, pvalueCutoff
#' @return compared cluster
#' @examples
#'     x <- compareClust()
#'
#' @export
#'
compareClust <- function(dat = NULL, ont = "CC", org = "org.Hs.eg.db",
                         fun = "enrichGO", title = "Ontology Distribution Comparison",
                         pvalueCutoff = 0.01) {
  if (is.null(dat)) return(NULL)
  if (!installpack(org)) return(NULL)
  res <- c()
  genecluster <- list()
  k <- max(dat$fit.cluster)
  for (i in 1:k) {
    clgenes <- rownames(dat[dat$fit.cluster == i, ])
    genelist <- getGeneList(clgenes, org)
    genecl <- list()
    genecl <- push(genecl, genelist)
    genecluster[c(paste("X", i, sep = ""))] <- genecl
  }
  res$table <- NULL
  p <- tryCatch({
    title <- paste(fun, title)
    xx <- c()
    if (fun == "enrichKEGG"){
      xx <- compareCluster(genecluster, fun = fun,
                           organism = getOrganism(org),
                           pvalueCutoff = pvalueCutoff)
    } else if (fun == "enrichDO") {
      if (!installpack("DOSE")) return(NULL)
      xx <- compareCluster(genecluster, fun = fun,
                           pvalueCutoff = pvalueCutoff)
    } else {
      title <- paste(ont, title)
      xx <- compareCluster(genecluster, fun = fun,
                           ont = ont, OrgDb = org, pvalueCutoff = pvalueCutoff)
      #ont = ont, organism = "human", pvalueCutoff = pvalueCutoff)
    }
    if (!is.null(xx@compareClusterResult) )
      res$table <- xx@compareClusterResult[,
                                           c("Cluster", "ID", "Description", "GeneRatio", "BgRatio",
                                             "pvalue", "p.adjust", "qvalue")]
    res$p <- dotplot(xx, title = title)
  })
  res
}
#' getEnrichDO
#'
#' Gathers the Enriched DO Term analysis data to be used within the
#' GO Term plots.
#'
#' @note \code{getEnrichDO}
#' @param genelist, gene list
#' @param pvalueCutoff, the p value cutoff
#' @return enriched DO
#' @examples
#'     x <- getEnrichDO()
#' @export
#'
getEnrichDO <- function(genelist = NULL, pvalueCutoff = 0.01) {
  if (is.null(genelist)) return(NULL)
  res <- c()
  res$enrich_p <- enrichDO(gene = genelist, ont = "DO",
                           pvalueCutoff = pvalueCutoff)

  res$p <- barplot(res$enrich_p, title = "Enrich DO", font.size = 12)
  res$table <- NULL
  if (!is.null(nrow(res$enrich_p@result)) )
    res$table <- res$enrich_p@result[,c("ID", "Description",
                                        "GeneRatio", "pvalue", "p.adjust", "qvalue")]
  res
}

#' drawKEGG
#'
#' draw KEGG patwhay with expression values
#'
#' @note \code{drawKEGG}
#' @param input, input
#' @param dat, expression matrix
#' @param pid, pathway id
#' @return enriched DO
#' @examples
#'     x <- drawKEGG()
#' @importFrom pathview pathview
#' @export
#'
drawKEGG <- function(input = NULL, dat = NULL, pid = NULL) {
  if (is.null(dat) && is.null(pid)) return(NULL)
  tryCatch({
    if (installpack("pathview")){
      org <- input$organism
      genedata <- getEntrezIds(dat[[1]], org)
      foldChangeData <- data.frame(genedata$log2FoldChange)
      rownames(foldChangeData) <- rownames(genedata)
      pathview(gene.data = foldChangeData,
               pathway.id = pid,
               species = substr(pid,0,3),
               gene.idtype="entrez",
               out.suffix = "b.2layer", kegg.native = TRUE)

      unlink(paste0(pid,".png"))
      unlink(paste0(pid,".xml"))
    }
  })
}


#' debrowserlowcountfilter
#'
#' Module to filter low count genes/regions
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param ldata, loaded data
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowserlowcountfilter()
#'
debrowserlowcountfilter <- function(input = NULL, output = NULL, session = NULL, ldata = NULL) {
  if (is.null(ldata)) return(NULL)
  fdata <- reactiveValues(count=NULL, meta = NULL)
  observeEvent(input$submitLCF, {
    if (is.null(ldata$count)) return (NULL)
    filtd <- ldata$count
    filtd[, colnames(filtd)] <- apply(filtd[, colnames(filtd)], 2, function(x) as.integer(x))

    if (input$lcfmethod == "Max"){
      filtd <- subset(filtd, apply(filtd, 1, max, na.rm = TRUE)  >=  as.numeric(input$maxCutoff))
    } else if (input$lcfmethod == "Mean") {
      filtd <- subset(filtd, rowMeans(filtd, na.rm = TRUE) >= as.numeric(input$meanCutoff))
    }
    else if (input$lcfmethod == "CPM") {
      cpmcount <- edgeR::cpm(filtd)
      filtd <- subset(filtd, rowSums(cpmcount > as.numeric(input$CPMCutoff),
                                     na.rm = TRUE) >= as.numeric(input$numSample))
    }
    fdata$count <- filtd
    fdata$meta <- ldata$meta
  })

  output$cutoffLCFMet <- renderUI({
    ret<-textInput(session$ns("maxCutoff"), "Filter features where Max Value <", value = "10" )
    if (input$lcfmethod == "Mean"){
      ret<-textInput(session$ns("meanCutoff"), "Filter features where Row Means <", value = "10" )
    }
    else if (input$lcfmethod == "CPM"){
      ret <- list(textInput(session$ns("CPMCutoff"), "Filter features where CPM <", value = "1" ),
                  textInput(session$ns("numSample"), "at least # of samples", value = toString(ncol(ldata$count)-1) ))
    }
    ret
  })

  filtereddata <- reactive({
    ret <- NULL
    if(!is.null(fdata$count)){
      ret <- fdata
    }
    return(ret)
  })

  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", ldata)
    getSampleDetails(output, "filteredSummary", "filteredDetails", filtereddata())
    getTableDetails(output, session, "loadedtable",  data = ldata$count,  modal = TRUE)
    callModule(debrowserhistogram, "beforeFiltering", ldata$count)

    if ( !is.null(filtereddata()$count ) && nrow(filtereddata()$count)>2 ) {
      getTableDetails(output, session, "filteredtable",  data = filtereddata()$count, modal = TRUE)
      callModule(debrowserhistogram, "afterFiltering", filtereddata()$count)
    }
  })

  list(filter=filtereddata)
}

#' dataLCFUI
#' Creates a panel to filter low count genes and regions
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- dataLCFUI("lcf")
#'
#' @export
#'
dataLCFUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "Low Count Filtering",
                          solidHeader = TRUE, status = "info",  width = 12,
                          fluidRow(
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("uploadSummary")),
                                         DT::dataTableOutput(ns("sampleDetails"))),
                                   uiOutput(ns("loadedtable"))
                            ),
                            column(2,
                                   shinydashboard::box(title = "Filtering Methods",
                                                       solidHeader = TRUE, status = "info",
                                                       width = 12,
                                                       lcfMetRadio(id),
                                                       uiOutput(ns("cutoffLCFMet")),
                                                       actionButtonDE(ns("submitLCF"), label = "Filter", styleclass = "primary")
                                   )
                            ),
                            column(5,div(style = 'overflow: scroll',
                                         tableOutput(ns("filteredSummary")),
                                         DT::dataTableOutput(ns("filteredDetails"))),
                                   uiOutput(ns("filteredtable"))
                            )
                          ),
                          conditionalPanel(condition = paste0("input['", ns("submitLCF"),"']"),
                                           actionButtonDE("Batch", label = "Batch Effect Correction", styleclass = "primary"),
                                           conditionalPanel(condition = "!(input.Batch)",
                                                            actionButtonDE("goDEFromFilter", "Go to DE Analysis", styleclass = "primary"),
                                                            actionButtonDE("goQCplotsFromFilter", "Go to QC plots", styleclass = "primary")))
      ),
      shinydashboard::box(title = "Histograms",
                          solidHeader = TRUE, status = "info",  width = 12,
                          fluidRow(
                            column(6,histogramControlsUI(ns("beforeFiltering")),
                                   getHistogramUI(ns("beforeFiltering"))),
                            column(6,histogramControlsUI(ns("afterFiltering")),
                                   getHistogramUI(ns("afterFiltering")))
                          ))
    ))
}

#' lcfMetRadio
#'
#' Radio buttons for low count removal methods
#'
#' @param id, namespace id
#' @note \code{lcfMetRadio}
#' @return radio control
#'
#' @examples
#'
#'     x <- lcfMetRadio("lcf")
#'
#' @export
#'
lcfMetRadio <- function(id) {
  ns <- NS(id)
  radioButtons(inputId=ns("lcfmethod"),
               label="Low count filtering method:",
               choices=c(Max='Max',
                         Mean='Mean',
                         CPM='CPM'
               ),
               selected='Max'
  )
}


#' debrowserheatmap
#'
#' Heatmap module to create interactive heatmaps and get selected list from
#' a heatmap
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param expdata, a matrix that includes expression values
#' @return heatmapply plot
#'
#' @examples
#'     x <- debrowserheatmap()
#'
#' @export
#'
#'
debrowserheatmap <- function( input, output, session, expdata = NULL){
  if(is.null(expdata)) return(NULL)
  output$heatmap <- renderPlotly({
    shinyjs::onevent("mousemove", "heatmap", js$getHoverName(session$ns("hoveredgenename")))
    shinyjs::onevent("click", "heatmap", js$getHoverName(session$ns("hoveredgenenameclick")))

    withProgress(message = 'Drawing Heatmap', detail = "interactive", value = 0, {
      runHeatmap(input, session, orderData())
    })
  })
  output$heatmap2 <- renderPlot({
    withProgress(message = 'Drawing Heatmap', detail = "non-interactive", value = 0, {
      runHeatmap2(input, session, orderData())
    })
  })
  heatdata <- reactive({
    cld <- prepHeatData(expdata, input)
    if (input$kmeansControl)
    {
      res <- niceKmeans(cld, input)
      cld <- res$clustered
    }
    cld
  })

  button <- reactiveVal(FALSE)
  orderData <- reactive({
    newclus <- heatdata()
    if (input$changeOrder && isolate(button()) && !is.null(input$clusterorder)){
      newclus <- changeClusterOrder(isolate(input$clusterorder), newclus)
    }
    button(FALSE)
    newclus
  })
  observeEvent(input$changeOrder,{
    button(TRUE)
  })
  output$heatmapUI <- renderUI({
    if (is.null(input$interactive)) return(NULL)
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("Heatmap"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,   getPlotArea(input, session))
  })

  hselGenes <- reactive({
    if (is.null(input$selgenenames)) return("")
    unlist(strsplit(input$selgenenames, split=","))
  })
  shg <- reactive({
    if (is.null(input$hoveredgenename)) return("")
    js$getSelectedGenes(session$ns("heatmap"), session$ns("selgenenames"))
    input$hoveredgenename
  })
  observe({
    if(!input$changeOrder)
      updateTextInput(session, "clusterorder", value = paste(seq(1:input$knum), collapse=","))

    if (is.null(shg()))
      js$getSelectedGenes()
  })
  shgClicked <- reactive({
    if (is.null(input$hoveredgenenameclick) || input$hoveredgenenameclick == "") return(input$hoveredgenename)
    input$hoveredgenenameclick
  })

  list( shg = (shg), shgClicked=(shgClicked), selGenes=(hselGenes), getSelected = (orderData))
}

#' getPlotArea
#'
#' returns plot area either for heatmaply or heatmap.2
#' @param input, input variables
#' @param session, session
#' @return heatmapply/heatmap.2 plot area
#'
#' @examples
#'     x <- getPlotArea()
#'
#' @export
#'
#'
getPlotArea <- function(input = NULL, session = NULL){
  if (is.null(input)) return(NULL)
  ret <- c()

  if (input$interactive){
    ret <- plotlyOutput(session$ns("heatmap"),
                        height=input$height, width=input$width)
  }
  else{
    ret <- plotOutput(session$ns("heatmap2"),
                      height = input$height, input$width)
  }
  ret
}

#' runHeatmap
#'
#' Creates a heatmap based on the user selected parameters within shiny
#' @param input, input variables
#' @param session, session
#' @param expdata, a matrix that includes expression values
#' @return heatmapply plot
#'
#' @examples
#'     x <- runHeatmap()
#'
#' @export
#'
#'
runHeatmap <- function(input = NULL, session = NULL, expdata = NULL){
  if (is.null(expdata)) return(NULL)
  cld <-expdata
  hclustfun_row <- function(x, ...) hclust(x, method = input$hclustFun_Row)
  hclustfun_col <- function(x, ...) hclust(x, method = input$hclustFun_Col)
  distfun_row <- function(x, ...) {
    if (input$distFun_Row != "cor") {
      return(dist(x, method = input$distFun_Row))
    } else {
      return(as.dist(1 - cor(t(x))))
    }
  }
  distfun_col <- function(x, ...) {
    if (input$distFun_Col != "cor") {
      return(dist(x, method = input$distFun_Col))
    } else {
      return(as.dist(1 - cor(t(x))))
    }
  }

  if (!input$customColors ) {
    heatmapColors <- eval(parse(text=paste0(input$pal,
                                            '(',input$ncol,')')))
  }
  else{
    if (!is.null(input$color1))
      heatmapColors <- colorRampPalette(c(input$color1,
                                          input$color2, input$color3))(n = 1000)
  }

  if (!input$kmeansControl){
    p <- heatmaply(cld,
                   main = input$main,
                   xlab = input$xlab,
                   ylab = input$ylab,
                   row_text_angle = input$row_text_angle,
                   column_text_angle = input$column_text_angle,
                   dendrogram = input$dendrogram,
                   branches_lwd = input$branches_lwd,
                   seriate = input$seriation,
                   colors = heatmapColors,
                   distfun_row =  distfun_row,
                   hclustfun_row = hclustfun_row,
                   distfun_col = distfun_col,
                   hclustfun_col = hclustfun_col,
                   showticklabels = c(input$labCol, input$labRow),
                   k_col = input$k_Col,
                   k_row = input$k_Row
    )
  }else {
    if (!input$showClasses){
      cld <- data.frame(cld)
      cld <- as.matrix(cld [, -match("class",names(cld))])
    }
    rhcr <- hclust(dist(cld))
    chrc <- hclust(dist(t(cld)))
    p <- heatmaply(cld,
                   main = input$main,
                   xlab = input$xlab,
                   ylab = input$ylab,
                   row_text_angle = input$row_text_angle,
                   column_text_angle = input$column_text_angle,
                   #dendrogram = input$dendrogram,
                   dendrogram = "none",
                   branches_lwd = input$branches_lwd,
                   seriate = input$seriation,
                   colors = heatmapColors,
                   showticklabels = c(input$labCol, input$labRow),
                   Rowv = as.dendrogram(rhcr),
                   Colv = as.dendrogram(chrc),
                   k_col = input$k_Col,
                   k_row = input$knum
    )
  }
  p <- p %>%
    plotly::layout(
      height=input$height, width=input$width,
      margin = list(l = input$left,
                    b = input$bottom,
                    t = input$top,
                    r = input$right
      ))
  p$elementId <- NULL
  p
}

#' runHeatmap2
#'
#' Creates a heatmap based on the user selected parameters within shiny
#' @param input, input variables
#' @param session, session
#' @param expdata, a matrix that includes expression values
#' @return heatmap.2
#'
#' @examples
#'     x <- runHeatmap2()
#'
#' @export
#'
#'
runHeatmap2 <- function(input = NULL, session = NULL, expdata = NULL){
  if(is.null(expdata)) return(NULL)
  if (nrow(expdata)>5000)
    expdata <- expdata[1:5000, ]

  if (!input$customColors ) {
    heatmapColors <- eval(parse(text=paste0(input$pal,
                                            '(',input$ncol,')')))
  }
  else{
    if (!is.null(input$color1))
      heatmapColors <- colorRampPalette(c(input$color1,
                                          input$color2, input$color3))(n = 1000)
    #heatmapColors <- colorRampPalette(c("red", "white", "blue"))(n = 1000)
  }

  hclustfun_row <- function(x, ...) hclust(x, method = input$hclustFun_Row)
  distfun_row <- function(x, ...) {
    if (input$distFun_Row != "cor") {
      return(dist(x, method = input$distFun_Row))
    } else {
      return(as.dist(1 - cor(t(x))))
    }
  }
  if (!input$showClasses && "class" %in% names(expdata) ){
    expdata <- data.frame(expdata)
    expdata <- as.matrix(expdata [, -match("class",names(expdata))])
  }
  if (input$kmeansControl){
    m <- heatmap.2(as.matrix(expdata), Rowv = FALSE, main = input$main, dendrogram = input$dendrogram,
                   Colv = FALSE, col = heatmapColors, labRow = input$labRow,
                   distfun = distfun_row, hclustfun = hclustfun_row, density.info = "none",
                   trace = "none", margins = c(10,10))
  }else{
    m <- heatmap.2(as.matrix(expdata), main = input$main, dendrogram = input$dendrogram,
                   col = heatmapColors, labRow = input$labRow,
                   distfun = distfun_row, hclustfun = hclustfun_row, density.info = "none",
                   trace = "none", margins = c(10,10))
  }
  m
}


#' changeClusterOrder
#'
#' change order of K-means clusters
#'
#' @note \code{changeClusterOrder}
#' @param order, order
#' @param cld, data
#' @return heatmap plot area
#' @examples
#'     x <- changeClusterOrder()
#' @export
#'
changeClusterOrder <- function(order = NULL, cld = NULL){
  if (is.null(order) || is.null(cld) ) return(NULL)
  newcluster <- c()
  idx <- as.integer(as.vector(unlist(strsplit(order, ","))))
  da <- data.frame(cld)
  for (i in 1:length(idx)) {
    newcluster <- rbind(newcluster, da[da$class == idx[i], ])
  }
  newcluster
}


#' niceKmeans
#'
#' Generates hierarchially clustered K-means clusters
#'
#' @note \code{niceKmeans}
#' @param df, data
#' @param input, user inputs
#' @param iter.max, max iteration for kmeans clustering
#' @param nstart, n for kmeans clustering
#' @return heatmap plot area
#' @examples
#'     x <- niceKmeans()
#' @export
#'
niceKmeans <-function (df = NULL, input = NULL, iter.max = 1000, nstart=100) {
  if(is.null(df)) return(NULL)
  source <-df
  kmeans <- kmeans(source, centers = input$knum, iter.max = iter.max, algorithm=input$kmeansalgo, nstart=nstart)
  clustered <- data.frame()
  distfun_row <- function(x, ...) {
    if (input$distFun_Row != "cor") {
      return(dist(x, method = input$distFun_Row))
    } else {
      return(as.dist(1 - cor(t(x))))
    }
  }
  breaks <- c();
  for (i in 1:input$knum) {
    cluster <- source[kmeans$cluster==i,]
    rows <- row.names(cluster)
    clust <- hclust(distfun_row(as.matrix(cluster)), method = input$hclustFun_Row)
    clust$rowInd <- clust[[3]]
    cluster.ordered <- cluster[clust$rowInd,]
    cluster.ordered.genes <- rows[clust$rowInd]
    row.names(cluster.ordered) <- cluster.ordered.genes
    class <- data.frame(row.names = cluster.ordered.genes)
    class[,"class"] <- i
    cluster.ordered <- cbind(cluster.ordered, class)
    clustered <- rbind(clustered, cluster.ordered)
    if(i > 1 & i < input$knum) {
      breaks[i] <- as.numeric(breaks[i-1]) + length(rows)
    } else if(i==1) {
      breaks[i] <- length(rows);
    }
  }

  result <- list();
  result$clustered <- clustered;
  result$breaks <- breaks;
  return(result);
}


#' getHeatmapUI
#'
#' Generates the left menu to be used for heatmap plots
#'
#' @note \code{getHeatmapUI}
#' @param id, module ID
#' @return heatmap plot area
#' @examples
#'     x <- getHeatmapUI("heatmap")
#' @export
#'
getHeatmapUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("heatmapUI"))
}

#' heatmapControlsUI
#'
#' Generates the left menu to be used for heatmap plots
#'
#' @note \code{heatmapControlsUI}
#' @param id, module ID
#' @return HeatmapControls
#' @examples
#'     x <- heatmapControlsUI("heatmap")
#' @export
#'
heatmapControlsUI <- function(id) {
  ns <- NS(id)
  list(
    checkboxInput(ns('interactive'), 'Interactive', value = FALSE),
    kmeansControlsUI(id),
    shinydashboard::menuItem("Scale Options",
                             checkboxInput(ns('scale'), 'Scale', value = TRUE),
                             checkboxInput(ns('center'), 'Center', value = TRUE),
                             checkboxInput(ns('log'), 'Log', value = TRUE),
                             textInput(ns('pseudo'),'Pseudo Count','0.1')
    ),
    dendControlsUI(id, "Row"),
    dendControlsUI(id, "Col"),
    shinydashboard::menuItem("Heatmap Colors",
                             conditionalPanel(paste0("!input['", ns("customColors"), "']"),
                                              palUI(id),
                                              sliderInput(ns("ncol"), "# of Colors",
                                                          min = 1, max = 256, value = 256)),
                             customColorsUI(id)
    ),
    shinydashboard::menuItem("Heatmap Dendrogram",
                             selectInput(ns('dendrogram'),'Type',
                                         choices = c("both", "row", "column", "none"),selected = 'both'),
                             selectizeInput(ns("seriation"), "Seriation",
                                            c(OLO="OLO", GW="GW", Mean="mean", None="none"),selected = 'OLO'),
                             sliderInput(ns('branches_lwd'),'Branch Width',
                                         value = 0.6,min=0,max=5,step = 0.1)
    ),
    shinydashboard::menuItem("Heatmap Layout",
                             textInput(ns('main'),'Title',''),
                             textInput(ns('xlab'),'Sample label',''),
                             sliderInput(ns('row_text_angle'),'Sample Text Angle',
                                         value = 0,min=0,max=180),
                             textInput(ns('ylab'), 'Gene/Region label',''),
                             sliderInput(ns('column_text_angle'),'Gene/Region Text Angle',
                                         value = 45,min=0,max=180)
    ))
}

#' kmeansControlsUI
#'
#' get kmeans controls
#'
#' @note \code{kmeansControlsUI}
#' @param id, module ID
#' @return controls
#' @examples
#'     x <- kmeansControlsUI("heatmap")
#' @export
#'
kmeansControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem("kmeans",
                           checkboxInput(ns('kmeansControl'), 'kmeans clustering', value = FALSE),
                           conditionalPanel(paste0("input['", ns("kmeansControl"), "']"),
                                            sliderInput(ns("knum"), "k: # of Clusters",
                                                        min = 2, max = 20, value = 2),
                                            selectizeInput(ns("kmeansalgo"), "kmeans.algorithm",
                                                           c("Hartigan-Wong", "Lloyd", "Forgy",
                                                             "MacQueen"), selected = 'Lloyd'),
                                            textInput(ns('clusterorder'),
                                                      'The order of the clusters', ""),
                                            actionButtonDE(ns("changeOrder"), label = "Change Order", styleclass = "primary"),
                                            checkboxInput(ns('showClasses'), 'Show Classes', value = FALSE)))
}

#' dendControlsUI
#'
#' get distance metric parameters
#'
#' @note \code{dendControlsUI}
#' @param id, module ID
#' @param dendtype, Row or Col
#' @return controls
#' @examples
#'     x <- dendControlsUI("heatmap")
#' @export
#'
dendControlsUI <- function(id, dendtype = "Row") {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(dendtype, " dendrogram"),
                           selectizeInput(ns(paste0("distFun_", dendtype)), "Dist. method",
                                          distFunParamsUI(),
                                          selected = 'euclidean'),
                           selectizeInput(ns(paste0("hclustFun_", dendtype)), "Clustering linkage",
                                          clustFunParamsUI(),
                                          selected = 'complete'),
                           sliderInput(ns(paste0("k_", dendtype)), "# of Clusters",
                                       min = 1, max = 10, value = 2),
                           checkboxInput(ns(paste0('lab',dendtype)), paste0(dendtype, ' Labels'), value = TRUE))
}

#' clustFunParamsUI
#'
#' get cluster function parameter control
#'
#' @note \code{clustFunParamsUI}
#' @return cluster params
#' @examples
#'     x <- clustFunParamsUI()
#' @export
#'
clustFunParamsUI <- function() {
  c(Complete= "complete",Single= "single",Average= "average",
    Mcquitty= "mcquitty",Median= "median",Centroid= "centroid",
    Ward.D= "ward.D",Ward.D2= "ward.D2")
}

#' distFunParamsUI
#'
#' get distance metric parameters
#'
#' @note \code{distFunParamsUI}
#' @return funParams
#' @examples
#'     x <- distFunParamsUI()
#' @export
#'
distFunParamsUI <- function() {
  c(Cor="cor", Euclidean="euclidean",Maximum='maximum',
    Manhattan='manhattan',Canberra='canberra',
    Binary='binary',Minkowski='minkowski')
}

#' palUI
#'
#' get pallete
#'
#' @note \code{palUI}
#' @param id, namespace ID
#' @return pals
#' @examples
#'     x <- palUI("heatmap")
#' @export
#'
palUI <- function(id) {
  ns <- NS(id)
  colSel='RdBu'
  selectizeInput(inputId = ns("pal"),
                 label ="Select Color Palette",
                 choices = c('RdBu' = 'RdBu',
                             'BlueRed' = 'bluered',
                             'RedBlue' = 'redblue',
                             'RdYlBu' = 'RdYlBu',
                             'RdYlGn' = 'RdYlGn',
                             'BrBG' = 'BrBG',
                             'Spectral' = 'Spectral',
                             'BuGn' = 'BuGn',
                             'PuBuGn' = 'PuBuGn',
                             'YlOrRd' = 'YlOrRd',
                             'Heat' = 'heat.colors',
                             'Grey' = 'grey.colors'),
                 selected=colSel)
}

#' customColorsUI
#'
#' get Custom Color controls
#'
#' @note \code{getColRng}
#' @param id, namespace ID
#' @return color range
#' @examples
#'     x <- customColorsUI("heatmap")
#' @export
#'
customColorsUI <- function(id) {
  ns <- NS(id)
  list(
    checkboxInput(ns('customColors'), 'Custom Colors', value = FALSE),
    conditionalPanel(paste0("input['", ns("customColors"), "']"),
                     colourpicker::colourInput(ns("color1"), "Choose min colour", "blue"),
                     colourpicker::colourInput(ns("color2"), "Choose median colour", "white"),
                     colourpicker::colourInput(ns("color3"), "Choose max colour", "red")))
}

#' prepHeatData
#'
#' scales the data
#'
#' @param expdata, a matrixthat includes expression values
#' @param input, input variables
#' @return heatdata
#'
#' @examples
#'     x <- prepHeatData()
#'
#' @export
#'
prepHeatData <- function(expdata = NULL, input = NULL)
{
  if(is.null(expdata)) return(NULL)
  ld <- expdata
  if (!is.null(input$pseudo))
    ld <- ld + as.numeric(input$pseudo)
  if (!is.null(input$log) && input$log)
    ld <- log2(ld)
  cldt <- scale(t(ld), center = input$center, scale = input$scale)
  cld <- t(cldt)
  return(cld)
}

#' getSelHeat
#'
#' heatmap selection functionality
#'
#' @param expdata, selected genes
#' @param input, input params
#' @return plot
#' @export
#'
#' @examples
#'     x <- getSelHeat()
#'
getSelHeat <- function(expdata = NULL, input = NULL) {
  if (is.null(input)) return(NULL)
  getSelected <- reactive({
    expdata[unlist(strsplit(input, ",")), ]
  })
  list( getSelected = isolate(getSelected) )
}


#' heatmapJScode
#'
#' heatmap JS code for selection functionality
#'
#' @return JS Code
#' @export
#'
#' @examples
#'     x <- heatmapJScode()
#'
heatmapJScode <- function() {
  'shinyjs.getHoverName = function(params){

    var defaultParams = {
    controlname : "hoveredgenename"
    };
    params = shinyjs.getParams(params, defaultParams);
    var out = ""

    if (typeof  document.getElementsByClassName("nums")[0] != "undefined"){
    if (typeof  document.getElementsByClassName("nums")[0].querySelectorAll("tspan.line")[0] != "undefined"){
    out = document.getElementsByClassName("nums")[0].querySelectorAll("tspan.line")[0].innerHTML.match("row: (.*)")[1]
    $("#heatmap-heatmap").attr("gname", out)
    }
    }
    Shiny.onInputChange(params.controlname, $("#heatmap-heatmap").attr("gname"));
    }
    shinyjs.resetInputParam = function(params){
        var defaultParams = {
                controlname : "hoveredgenename"
        };
        params = shinyjs.getParams(params, defaultParams);
        console.log(params.controlname)
        Shiny.onInputChange(params.controlname, "");
    }

    shinyjs.getSelectedGenes = function(params){
    var defaultParams = {
    plotId : "heatmap",
    controlname : "selgenenames"
    };
    params = shinyjs.getParams(params, defaultParams);
    var count = document.getElementById(params.plotId).querySelectorAll("g.y2tick").length
    var start = 0
    var out = ""

    for (i = start; i < count; i++)
    {
        if (typeof document.getElementById(params.plotId).querySelectorAll("g.y2tick")[i] != "undefined"){
        out += document.getElementById(params.plotId).querySelectorAll("g.y2tick")[i].innerHTML.match(">(.*)</text>")[1]  + ","
        }
    }
    Shiny.onInputChange(params.controlname, out);
    }'
}

#' getJSLine
#'
#' heatmap JS code for selection functionality
#'
#' @return JS Code
#' @export
#'
#' @examples
#'     x <- getJSLine()
#'
getJSLine <-function()
{
  list(shinyjs::useShinyjs(),
       shinyjs::extendShinyjs(text = heatmapJScode(), functions = c("getHoverName", "getSelectedGenes", "resetInputParam")))
}


#' heatmapServer
#'
#' Sets up shinyServer to be able to run heatmapServer interactively.
#'
#' @note \code{heatmapServer}
#' @param input, input params from UI
#' @param output, output params to UI
#' @param session, session variable
#' @return the panel for main plots;
#'
#' @examples
#'     heatmapServer
#'
#' @export
#'
heatmapServer <- function(input, output, session) {
  updata <- reactiveVal()
  selected <- reactiveVal()
  expdata <- reactiveVal()
  observe({
    updata(callModule(debrowserdataload, "load", "Submit"))
  })
  observe({
    if(!is.null(updata()$load()$count))
      if (nrow(updata()$load()$count) > 1000){
        updateCheckboxInput(session, "mostvaried", value = TRUE)
        expdata(getMostVariedList(updata()$load()$count,
                                  colnames(updata()$load()$count), input))
      }
    else
      expdata(updata()$load()$count)
  })

  observeEvent (input$Submit, {
    updateTabItems(session, "DEBrowserHeatmap", "Heatmap")
  })
  observe({
    if (!is.null(expdata())){
      withProgress(message = 'Creating plot', style = "notification", value = 0.1, {
        selected(callModule(debrowserheatmap, "heatmap", expdata()))
      })
    }
  })
  output$heatmap_hover <- renderPrint({
    if (!is.null(selected()) && !is.null(selected()$shgClicked()) &&
        selected()$shgClicked() != "")
      return(paste0("Clicked: ",selected()$shgClicked()))
    else
      return(paste0("Hovered:", selected()$shg()))
  })
  output$heatmap_selected <- renderPrint({
    if (!is.null(selected()))
      selected()$selGenes()
  })
  output$topn <- renderPrint({
    if (!is.null(input$topn))
      input$topn
  })
  output$mincount <- renderPrint({
    if (!is.null(input$mincount))
      input$mincount
  })
}

#' heatmapUI
#'
#' Creates a shinyUI to be able to run DEBrowser interactively.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#'
#' @note \code{heatmapUI}
#' @return the panel for heatmapUI;
#'
#' @examples
#'     x<-heatmapUI()
#'
#' @export
#'
heatmapUI <- function(input, output, session) {
  header <- dashboardHeader(
    title = "DEBrowser Heatmap"
  )
  sidebar <- dashboardSidebar(  getJSLine(),
                                sidebarMenu(id="DEBrowserHeatmap",
                                            menuItem("Upload", tabName = "Upload"),
                                            menuItem("Heatmap", tabName = "Heatmap"),
                                            menuItem("Options", tabName = "Heatmap",
                                                     checkboxInput('mostvaried', 'Most Varied Set', value = FALSE),
                                                     conditionalPanel( (condition <- "input.mostvaried"),
                                                                       textInput("topn", "top-n", value = "500" ),
                                                                       textInput("mincount", "total min count", value = "10" )),
                                                     plotSizeMarginsUI("heatmap"),
                                                     heatmapControlsUI("heatmap"))))

  body <- dashboardBody(
    tabItems(
      tabItem(tabName="Upload", dataLoadUI("load")),
      tabItem(tabName="Heatmap",  getHeatmapUI("heatmap"),
              column(4,
                     verbatimTextOutput("heatmap_hover"),
                     verbatimTextOutput("heatmap_selected"),
                     verbatimTextOutput("topn"),
                     verbatimTextOutput("mincount")
              ))
    ))

  dashboardPage(header, sidebar, body, skin = "blue")
}



#' getPCAPlotUI
#'
#' PCA plots UI.
#' @param id, namespace id
#' @note \code{getPCAPlotUI}
#' @return the panel for PCA plots;
#'
#' @examples
#'     x <- getPCAPlotUI("pca")
#'
#' @export
#'
getPCAPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("pcaplot"))
}

#' debrowserpcaplot
#'
#' Module for a pca plot with its loadings
#' as a mainplot in debrowser
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param pcadata, a matrix that includes expression values
#' @param metadata, metadata to color the plots
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowserpcaplot()
#'
debrowserpcaplot <- function(input = NULL, output = NULL, session = NULL, pcadata = NULL, metadata = NULL) {
  if(is.null(pcadata)) return(NULL)
  qcplots <-  reactive({
    sc <- getShapeColor(input)
    plot_pca(pcadata, input$pcselx, input$pcsely,
             metadata = metadata, color = sc$color,
             size = 5, shape = sc$shape,
             textonoff = sc$textonoff,
             legendSelect = sc$legendSelect, input = input )
  })
  output$pcaplot <- renderUI({
    list(fluidRow(
      column(12,
             shinydashboard::box(
               collapsible = TRUE, title = "PCA Plot", status = "primary",
               solidHeader = TRUE, width = NULL,
               draggable = TRUE, plotlyOutput(session$ns("pca1"),
                                              height= input$height, width=input$width, inline=TRUE)
             ),
             shinydashboard::box(
               collapsible = TRUE, title = "Loadings", status = "primary",
               solidHeader = TRUE, width = NULL,
               draggable = TRUE,  plotlyOutput(session$ns("pca2"),
                                               height= input$height, width=input$width, inline=TRUE) )) )
    )
  })
  output$pca1 <- renderPlotly({
    p <- qcplots()$plot1
    p$elementId <- NULL
    p
  })
  output$pca2 <- renderPlotly({
    p <- qcplots()$plot2
    p$elementId <- NULL
    p
  })
  output$colorShapeSelect <- renderUI({
    getColorShapeSelection(metadata, input, session)
  })
}

#' pcaPlotControlsUI
#'
#' Generates the PCA PLots Left menu to be displayed within the DEBrowser.
#'
#' @param id, namespace id
#' @note \code{pcaPlotControlsUI}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- pcaPlotControlsUI("pca")
#' @export
#'
pcaPlotControlsUI <- function(id  = "pca") {
  ns <- NS(id)
  list(fluidRow(column(12, getPCselection(id, 1, "x")),
                column(12, getPCselection(id, 2, "y"))),
       fluidRow(
         column(12, getHideLegendOnOff(id)),
         column(12, getTextOnOff(id)),
         column(12, getLegendSelect(id))),
       uiOutput(ns("colorShapeSelect")))
}

#' run_pca
#'
#' Runs PCA on the selected dataset.
#'
#' @param x, dataframe with experiment data
#' @param retx, specifies if the data should be returned
#' @param center, center the PCA (Boolean)
#' @param scale, scale the PCA (Boolean)
#' @return pca list
#' @examples
#'     load(system.file("extdata", "demo", "demodata.Rda",
#'         package="debrowser"))
#'     pca_data<-run_pca(getNormalizedMatrix(
#'         demodata[rowSums(demodata[,1:6])>10,1:6]))
#'
#' @export
#'
run_pca <- function(x=NULL, retx = TRUE,
                    center = TRUE, scale = TRUE) {
  if ( is.null(x) || ncol(x) < 2) return (NULL)
  x <- subset(x, apply(x, 1, var, na.rm = TRUE) >  0)
  pca <- prcomp(t(x), retx = retx,
                center = center, scale. = scale)
  variances <- pca$sdev ^ 2
  explained <- variances / sum(variances)

  return(list(PCs = pca$x, explained = explained, pca = pca))
}

#' plot_pca
#'
#' Plots the PCA results for the selected dataset.
#'
#' @param dat, data
#' @param pcx, x axis label
#' @param pcy, y axis label
#' @param metadata, additional data
#' @param color, color for plot
#' @param shape, shape for plot
#' @param size, size of the plot
#' @param textonoff, text on off
#' @param legendSelect, select legend
#' @param input, input param
#' @return pca list
#' @examples
#'     load(system.file("extdata", "demo", "demodata.Rda",
#'             package="debrowser"))
#'     metadata<-cbind(colnames(demodata[,1:6]),
#'             colnames(demodata[,1:6]),
#'             c(rep("Cond1",3), rep("Cond2",3)))
#'     colnames(metadata)<-c("samples", "color", "shape")
#'
#'     a <- plot_pca(getNormalizedMatrix(
#'             demodata[rowSums(demodata[,1:6])>10,1:6]),
#'             metadata = metadata, color = "samples",
#'             size = 5, shape = "shape")
#'
#' @export
#'
plot_pca <- function(dat = NULL, pcx = 1, pcy = 2,
                     metadata = NULL, color = NULL, shape = NULL,
                     size = NULL, textonoff = "Off", legendSelect = "samples", input = NULL) {
  if ( is.null(dat) || ncol(dat) < 2) return(NULL)

  pca_data <- run_pca(dat)
  p_data <- prepPCADat(pca_data, metadata, input, pcx, pcy)

  # Prepare axis labels
  xaxis <- sprintf("PC%d (%.2f%%)", pcx,
                   round(pca_data$explained[pcx] * 100, 2))
  yaxis <- sprintf("PC%d (%.2f%%)", pcy,
                   round(pca_data$explained[pcy] * 100, 2))

  plot1 <- ggplot(data=p_data, aes(x=x, y=y))

  if (legendSelect == "color") {
    plot1 <-  plot1 + geom_point(mapping=aes(shape=shape, color=color), size=3 )
  }else{
    plot1 <-  plot1 + geom_point(mapping=aes(shape=shape, color=shape), size=3 )
  }
  if (textonoff == "On")
    plot1 <- plot1 + geom_text(aes(label=samples), vjust = 0, nudge_y = 1)
  plot1 <- plot1 + theme(legend.title = element_blank())
  plot1 <- plot1 +  labs(x = xaxis, y = yaxis)
  if (!is.null(input$top))
    plot1 <- plot1 + theme( plot.margin = margin(t = input$top, r =input$right, b =input$bottom, l = input$left, "pt"))
  plot1 <- ggplotly(plot1, width = input$width, height = input$height)
  if (!is.null(input$legendonoff) && input$legendonoff=="Off")
    plot1 <- plotly::hide_legend(plot1)

  plot1$elementId <- NULL

  pcaExp <- getPCAexplained(dat, pca_data, input)
  plot2 <- drawPCAExplained(pcaExp$plotdata)
  plot2$elementId <- NULL
  return (list(plot1 =  plot1, plot2 =  plot2, pcaset = pcaExp$pcaset))
}

#' prepPCADat
#'
#' prepares pca data with metadata. If metadata doesn't exists
#' it puts all the sampels into a signlge group; "Conds".
#'
#' @param pca_data, pca run results
#' @param metadata, additional meta data
#' @param input, input
#' @param pcx, x axis label
#' @param pcy, y axis label
#' @return Color and shape from selection boxes or defaults
#' @examples
#'     x <- prepPCADat()
#' @export
#'
prepPCADat <- function(pca_data = NULL, metadata = NULL, input = NULL, pcx = 1, pcy = 2){
  if (is.null(pca_data)) return (NULL)
  rownames(metadata) <- metadata[,1]

  x <- pca_data$PCs
  plot_data <- data.frame(x)
  # Prepare data frame to pass to ggplot
  xaxis <- paste0("PC", pcx)
  yaxis <- paste0("PC", pcy)
  if (!is.null(metadata)) {
    samples <- rownames(plot_data)
    color  <- rownames(plot_data)
    shape <- "Conds"
    if (!is.null(input$color_pca) && input$color_pca != "None")
      color <- as.character(metadata[samples, input$color_pca])
    if (!is.null(input$shape_pca) && input$shape_pca != "None")
      shape <- as.character(metadata[samples, input$shape_pca])

    metadata <- cbind(samples, color, shape)
    plot_data <- cbind(plot_data, metadata)
    p_data <- plot_data[,c(xaxis, yaxis, "samples", "color", "shape")]
  } else {
    samples <- rownames(plot_data)
    color  <- rownames(plot_data)
    shape <- "Conds"
    p_data <- cbind( plot_data[,c(xaxis, yaxis)], samples, color, shape)
  }
  colnames(p_data) <- c("x", "y", "samples", "color", "shape")
  p_data
}

#' getPCAexplained
#'
#' Creates a more detailed plot using the PCA results from
#' the selected dataset.
#'
#' @param datasetInput, selected data
#' @param pca_data, from user
#' @param input, input params
#' @return explained plot
#' @examples
#' load(system.file("extdata", "demo", "demodata.Rda", package="debrowser"))
#' input<-c()
#' input$qcplot<-"pca"
#' input$col_list<-colnames(demodata[,1:6])
#' dat <- getNormalizedMatrix(demodata[,1:6])
#' pca_data <- run_pca(dat)
#' x <- getPCAexplained(dat, pca_data, input)
#'
#' @export
#'
getPCAexplained <- function(datasetInput = NULL,
                            pca_data = NULL, input = NULL) {
  if (is.null(datasetInput)) return(NULL)
  datexp <- NULL
  pcaset <- NULL
  size <- length(pca_data$explained)
  if (size>9)
    size <- 9
  datexp <- data.frame(cbind(unlist(lapply(
    c(1:size),
    function(x){paste0("PC", x)})),
    round(pca_data$explained * 100, 2)))
  colnames(datexp) <- c("PCs", "explained")
  datexp$explained <- as.numeric( as.character(datexp$explained) )
  datexp <- datexp[1:size,]
  var <- pca_data$pca$sdev^2/sum(pca_data$pca$sdev^2)

  ## Select the genes for PCA, removing the least variable

  dThresh.pctile <- 1 - as.numeric(input$pctile)     # distance threshold
  gList.dThresh <- c()

  d <- pca_data$pca$rotation[,c(input$pcselx)]
  dThresh<-quantile(d, dThresh.pctile)
  gList.dThresh <- names(which(d>dThresh))
  pcaset <-  datasetInput[gList.dThresh, ]
  return (list(plotdata =  datexp, pcaset = pcaset))
}

#' getShapeColor
#'
#' Generates the fill and shape selection boxes for PCA plots.
#' metadata file has to be loaded in this case
#'
#' @param input, input values
#' @return Color and shape from selection boxes or defaults
#' @examples
#'     x <- getShapeColor()
#' @export
#'
getShapeColor <- function(input = NULL) {
  if (is.null(input)) return (NULL)
  sc <-  c()
  if (!is.null(input$color_pca))
    sc$color <- input$color_pca
  if (!is.null(input$shape_pca))
    sc$shape <- input$shape_pca

  sc$textonoff <- input$textonoff
  sc$legendSelect <- input$legendSelect
  return(sc)
}

#' Creates a more detailed plot using the PCA results from
#' the selected dataset.
#'
#' @param explainedData, selected data
#' @return explained plot
#' @examples
#'     x <- drawPCAExplained()
#'
#' @export
#'
drawPCAExplained <- function(explainedData = NULL){
  p <- NULL
  if (is.null(explainedData)) return(NULL)

  p<- plot_ly(data=explainedData, x=~PCs, y=~explained,
              type = 'bar')

  p$elementId <- NULL
  p
}


#' getPCselection
#'
#' Generates the PC selection number to be used within DEBrowser.
#' @param id, namespace id
#' @param num, PC selection number
#' @param xy, x or y coordinate
#' @note \code{getPCselection}
#' @return PC selection for PCA analysis
#' @examples
#'     x <- getPCselection("pca")
#' @export
#'
getPCselection <- function(id, num = 1, xy = "x" ) {
  ns <- NS(id)
  numericInput(ns(paste0("pcsel", xy)),
               paste0("PC selection[", xy, "]"), num, 1, 6)
}

#' getColorShapeSelection
#'
#' Generates the fill and shape selection boxes for PCA plots.
#' metadata file has to be loaded in this case
#'
#' @param metadata, metadata table
#' @param input, input
#' @param session, session
#' @return Color and shape selection boxes
#' @examples
#'     x <- getColorShapeSelection()
#' @export
#'
getColorShapeSelection <- function(metadata = NULL, input = NULL, session = NULL) {
  if (is.null(metadata) ||  is.null(session)) return (NULL)
  list(fluidRow(column(12, selectGroupInfo(metadata, input, session$ns("color_pca"), "Color field")),
                column(12, selectGroupInfo(metadata, input, session$ns("shape_pca"), "Shape field"))))
}

#' getLegendSelect
#'
#' select legend
#' @param id, namespace id
#' @note \code{getLegendSelect}
#' @examples
#'     x <- getLegendSelect("pca")
#' @export
#'
getLegendSelect <- function(id = "pca") {
  ns <- NS(id)
  lst.choices <- as.list(c("color", "shape"))
  selectInput(ns("legendSelect"), label = "Select legend",
              choices = lst.choices,
              selected = "color")
}

#' getTextOnOff
#'
#' text on PCA plot on and off
#' @param id, namespace id
#' @note \code{getTextOnOff}
#' @examples
#'     x <- getTextOnOff("pca")
#' @export
#'
getTextOnOff <- function(id = "pca") {
  ns <- NS(id)
  lst.choices <- as.list(c("On", "Off"))
  selectInput(ns("textonoff"), label = "Text On/Off",
              choices = lst.choices,
              selected = "Off")
}

#' getHideLegendOnOff
#'
#' hide legend
#' @param id, namespace id
#' @examples
#'     x <- getHideLegendOnOff("pca")
#' @export
#'
getHideLegendOnOff <- function(id = "pca") {
  ns <- NS(id)
  lst.choices <- as.list(c("On", "Off"))
  selectInput(ns("legendonoff"), label = "Legend On/Off",
              choices = lst.choices,
              selected = "On")
}



#' installpack
#'
#' install packages if they don't exist
#' display.
#'
#' @param package_name, package name to be installed
#' @note \code{installpack}
#'
#' @examples
#'     x <- installpack()
#'
#' @export
#'
installpack <- function(package_name = NULL) {
  if (is.null(package_name)) return(NULL)

  if (!loadpack(package_name))
  {
    txt <- paste0("Please install ",package_name, " to use this function.")
    print(txt)
    showNotification(txt, type="error")
    return(FALSE)
  }else{
    loadpack(package_name)
  }
  return(TRUE)
}

#' loadpack
#'
#' load packages
#'
#' @param package_name, package name to be loaded
#' @note \code{loadpack}
#'
#' @examples
#'     x <- loadpack()
#'
#' @export
#'
loadpack <- function (package_name = NULL){
  if(isTRUE(package_name %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", package_name)))
    return (TRUE)
  }
  return(FALSE)
}



#' getIntroText
#' Intro text
#'
#' @return the JS for tab updates
#'
#' @examples
#'     x<- getIntroText()
#'
#' @export
getIntroText<-function(){

  list(
    h2("Quick Start Guide"),
    p("Differential expression (DE) analysis has become an increasingly popular tool
in determining and viewing up and/or down experssed genes between two sets of
samples."),
    p("The goal of DEBrowser is to provide an easy way to perform and visualize DE analysis."),
    #####################################################
    ##Data input section
    h3("1. Data input"),
    p("The input data for DEBrowser is 2 files in '.txt', '.csv' or '.tsv' format,
  they are named as 'Count Data File' and 'Metadata File'."),
    h4("1.1 Count Data File"),
    p("This input file should contain summarized count results of all samples in the experiment,
  an example of the expected input data format is presented as below:"),
    tableOutput("countFile"),
    p("Where columns are samples, rows are the mapped genomic features
  (e.g. genes, isoforms, miRNAs, ATAC or Chip regions etc.)."),
    p("If you do not have a dataset to upload, you can use the built in demo data file by
  clicking on the 'Load Demo' buttons from two different publications. To view the entire demo data file, you can
  download the ", a("demo set 1 (Vernia et. al).", href="https://bioinfo.umassmed.edu/pub/debrowser/simple_demo.tsv"),
      "For another example, try our ", a("full dataset for demo set 1 (Vernia et. al)", href="https://bioinfo.umassmed.edu/pub/debrowser/advanced_demo.tsv")),
    h4("1.2 Metadata File"),
    p("In addition to the count data file; you may also upload metadata file to correct for
batch effects or any other normalizing conditions you might want to address that might be
within your results. To handle for these conditions, simply create a metadata file by
using the example table at below "),
    tableOutput("metaFile"),
    p("or download sample metadata file used for:",
      a("metadata file for set 1 (Vernia et. al).", href="https://bioinfo.umassmed.edu/pub/debrowser/simple_demo_meta.txt")),
    p("To be able to upload the data please press",
      actionButton("UploadBut", label = "Upload", styleclass = "primary"), "button in the data upload page."),
    p("After sucessfull upload, you should see the summary of your data in the 'Upload Summary' section.
  To move the filtering section please click",
      actionButton("FilterBut", label = "Filter", styleclass = "primary"), "button in the upload page."),
    p("If you are ready to use DEBrowser, please click 'Upload' menu on the left to start using DEBrowser"))
}

#' getDataAssesmentText
#' DataAssesment text
#'
#' @return help text for data assesment
#'
#' @examples
#'     x<- getDataAssesmentText()
#'
#' @export
getDataAssesmentText<-function(){
  list(
    h3("2. Data Assesment"),
    h4("2.1 Low Count Filtering"),
    p("In this section, you can simultaneously visualise the changes of your dataset
    while filtering out the low count genes. Choose your filtration criteria from
    Filtering Methods box which is located just center of the screen.
    Three methods are available to be used:"),
    p(strong("Max:"), "Filters out genes where maximum count for each gene across all samples are less
   than defined threshold."),
    p(strong("Mean:"), "Filters out genes where mean count for each gene are less than defined threshold."),
    p(strong("CPM:"), "First, counts per million (CPM) is calculated as the raw counts divided by the
    library sizes and multiplied by one million. Then it filters out genes where at least
    defined number of samples is less than defined CPM threshold."),
    withMathJax(),
    p("The expression cutoff value is determined according to the library size
    and normalization factors with formula $$\\text{CPM} = \\frac{\\text{raw counts}}{\\text{library size} * \\text{normalization factors} * 10^{-6}}$$
    For example, if the cutoff CPM value is 10,
    the library size and normalization factors are estimated approximately equal to \\(\\ 3 \\text{ x} 10 ^ 6\\) and 1 for at least 4 samples,
    then 10 CPM expression cutoff corresponds to about 30 read counts.
    Therefore, in this example features in more than 4 samples have less than
    30 read counts (10 CPM) is going to be low expression features
    and will be removed for batch effect correction and DE analysis."),
    p("To be able to filter out the low expression counts please press",
      actionButton("FilterBut", label = "Filter", styleclass = "primary"), "button in data filtering page."),
    h4("2.2 Quality control(QC)"),
    p("After filtering low count features, you may continue your analysis with Batch Effect
    Detection & Correction or directly jump to differential expression analysis or view
    quality control (QC) information of your dataset."),
    p("If specified metadata file containing your treatment and batch fields, by clicking ",
      actionButton("BatchBut", label = "Batch effect correction", styleclass = "primary"),
      "button, you have the option to conduct PCA, interquartile range (IQR) and density
    plots to asses if the data requires batch effect correction or not."),
    p("If user wants to skip batch effect assesment and correction step, they can either click",
      actionButton("goDEFromFilterBut", "Go to DE Analysis", styleclass = "primary"),
      " button to perform DE Analysis or ",
      actionButton("goQCplotsFromFilterBut", "Go to QC plots", styleclass = "primary"),
      " button for QC plots to draw PCA, all2all scatter, heatmaps, IQR and density plots.")
  )
}

#' getDataPreparationText
#' DataPreparation text
#'
#' @return help text for data preparation
#'
#' @examples
#'     x<- getDataPreparationText()
#'
#' @export
getDataPreparationText<-function(){
  list(
    h3("3. Data Preparation"),
    p("With metadata file containing your batch correction fields
    then you have the option to conduct ",strong("batch effect correction"), " prior to
    your analysis. By adjusting parameters of Options box, you can investigate
    your character of your dataset. These parameters of the options box are
    explained as following:"),
    p(strong("Normalization Method:"), "DEBrowser allows performing normalization
    prior the batch effect correction. You may choose your normalization method
    (among MRE, TMM, RLE, upperquartile), or if you don't want to normalize your
    data you can select none for this item."),
    p(strong("Correction Method:"), "DEBrowser uses ComBat (part of the SVA
    bioconductor package) or Harman to adjust for possible batch effect or conditional
    biases. For more information, you can visit following links for
    documentation: ComBat, Harman"),
    p(strong("Treatment:"), "Please select the column that is specified in
    metadata file for comparision, such as cancer vs control. It is named
    condition for our sample metadata."),
    p(strong("Batch:"), "Please select the column name in metadata file
    which differentiate the batches. For example in our metadata, it is called batch.
    Upon clicking submit button, comparison tables and plots will be created on the right
    part of the screen as shown below."),
    p("You can investigate the changes on the data by comparing following features:"),
    p("1. Read counts for each sample.",br(),
      "2. PCA, IQR and Density plot of the dataset.", br(),
      "3. Gene/region vs samples data"),
    p("After batch effect correction, user can click ",
      actionButton("goDEFromFilterBut", "Go to DE Analysis", styleclass = "primary"),
      " button to perform DE Analysis or ",
      actionButton("goQCplotsFromFilterBut", "Go to QC plots", styleclass = "primary"),
      " button for QC plots to draw PCA, all2all scatter, heatmaps, IQR and density plots.")
  )
}

#' getDEAnalysisText
#' DEAnalysis text
#'
#' @return help text for DE Analysis
#'
#' @examples
#'     x<- getDEAnalysisText()
#'
#' @export
getDEAnalysisText<-function(){
  list(
    h3("4. DE analysis"),
    p("The goal of differential gene expression analysis is to find genes
or transcripts whose difference in expression, when accounting for the
variance within condition, is higher than expected by chance."),
    p(a("DESeq2", href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
      "is an R package available via Bioconductor and is designed to normalize count
data from high-throughput sequencing assays such as RNA-Seq and test for
differential expression (Love et al. 2014).  With multiple parameters such as
padjust values, log fold changes, plot styles, and so on, altering plots
created with your DE data can be a hassle as well as time consuming. The
Differential Expression Browser uses DESeq2 (Love et al., 2014)",
      a("EdgeR",href="https://bioconductor.org/packages/release/bioc/html/edgeR.html"),
      "(Robinson et al., 2010), and ",
      a("Limma", href="https://bioconductor.org/packages/release/bioc/html/limma.html"),
      "(Ritchie et al., 2015) coupled with shiny (Chang, W. et al., 2016)
to produce real-time changes within your
plot queries and allows for interactive browsing of your DE results.
In addition to DE analysis, DEBrowser also offers a variety of other plots
and analysis tools to help visualize your data even further."),
    p("If you are ready to discover and visualize your data, please click ",
      actionButton("mainPlotsBut", "Go to Main Plots", styleclass = "primary"),
      "button in DE Results section."),
    h4("4.1 Used parameters for DESeq2"),
    p(strong("fitType:")),
    p("Either 'parametric', 'local', or 'mean' for the type of fitting of
dispersions to the mean intensity. See estimateDispersions for description."),
    p(strong("betaPrior:")),
    p("Whether or not to put a zero-mean normal prior on the non-intercept
coefficients See nbinomWaldTest for description of the calculation of
the beta prior. By default, the beta prior is used only for the Wald test,
but can also be specified for the likelihood ratio test."),
    p(strong("testType:")),
    p("Either 'Wald' or 'LRT', which will then use either Wald significance tests
(defined by nbinomWaldTest), or the likelihood ratio test on the difference in
deviance between a full and reduced model formula (defined by nbinomLRT)"),
    h4("4.2 Used parameters for EdgeR"),
    p(strong("Normalization:")),
    p("Calculate normalization factors to scale the raw library sizes. Values
can be 'TMM','RLE','upperquartile','none'."),
    p(strong("Dispersion:")),
    p("Either a numeric vector of dispersions or a character string indicating
that dispersions should be taken from the data object."),
    p(strong("testType:")),
    p("ExactTest or glmLRT. ",strong("exactTest:")," Computes p-values for differential
abundance for each gene between two samples, conditioning on the total
count for each gene. The counts in each group are assumed to follow a
binomial distribution. ",strong("glmLRT:")," Fits a negative binomial generalized
log-linear model to the read counts for each gene and conducts
genewise statistical tests."),
    h4("4.3 Used parameters for Limma"),
    p(strong("Normalization:")),
    p("Calculate normalization factors to scale the raw library sizes.
Values can be 'TMM','RLE','upperquartile','none'."),
    p(strong("Fit Type:")),
    p("fitting method; 'ls' for least squares or 'robust' for robust regression"),
    p(strong("Norm. Bet. Arrays:")),
    p("Normalization Between Arrays; Normalizes expression intensities so that the
intensities or log-ratios have similar distributions across a set of arrays.")
  )
}
#' getQAText
#' Some questions and answers
#'
#' @return help text for QA
#'
#' @examples
#'     x<- getQAText()
#'
#' @export
getQAText<-function(){
  list(
    h3("5. Frequently asked questions (FAQ)"),
    h4("5.1 Why un-normalized counts?"),
    p("DESeq2 requires count data as input obtained from
          RNA-Seq or another high-thorughput sequencing experiment
          in the form of matrix values. Here we convert un-integer
          values to integer to be able to run DESeq2. The matrix values
          should be un-normalized, since DESeq2 model internally corrects for
          library size. So, transformed or normalized values such as counts
          scaled by library size should not be used as input. Please use edgeR
          or limma for normalized counts."),
    h4("5.2 Why am I getting error while uploading files?"),
    p("* DEBrowser supports tab, comma or semi-colon separated files. However spaces or characters in numeric regions not supported and causes an error while uploading files. It is crutial to remove these kind of instances from the files before uploading files."),
    p("* Another reason of getting an error is using same gene name multiple times. This may occurs after opening files in programs such as Excel, which tends to automatically convert some gene names to dates (eg. SEP9 to SEP.09.2018). This leads numerous problems therefore you need to disable these kind of automatic conversion before opening files in these kind of programs."),
    p("* Some files contain both tab and space as an delimiter which lead to error. It is required to be cleaned from these kind of files before loading."),
    h4("5.3 Why some columns not showed up after upload?"),
    p("If a character in numeric area or space is exist in one of your column, either column will be eliminated or you will get an error. Therefore it is crutial to remove for these kind of instances from your files before uploading."),
    h4("5.4 Why am I getting error while uploading CSV/TSV files exported from Excel?"),
    p("* You might getting an error, because of using same gene name multiple times. This may occurs after opening files in programs such as Excel, which tends to automatically convert some gene names to dates (eg. SEP9 to SEP.09.2018). Therefore you need to disable these kind of automatic conversion before opening files in these kind of programs."),
    h4("5.5 Why can't I see all the background data in Main Plots?"),
    p("In order to increase the performance, by default 10% of non-significant(NS) genes are used to generate plots. We strongly suggest you to use all of the NS genes in your plots while publishing your results. You can easily change this parameter by clicking **Main Options** button and change Background Data(%) to 100% on the left sidebar."),
    h4("5.6 Why am I getting error when I click on DE Genes in Go Term Analysis?"),
    p("To start ", strong("Go Term"), " analysis, it is important to select correct organism from ", strong("Choose an organism"), " field. After selecting other desired parameters, you can click ", strong("Submit")," button to run Go Term analysis. After this stage, you will able to see", strong(" categories")," regarding to your selected gene list in the ", strong("Table")," Tab. Once you select this category, you can click DE Genes button to see gene list regarding to selected category."),
    h4("5.7 How to download selected data from Main plots/QC Plots/Heatmaps?"),
    p("First, you need to choose ", strong("Choose dataset"), " field as ",strong("selected")," under ",strong("Data Options")," in the left sidebar. When you select this option, new field: ",strong("The plot used in selection")," will appear under ", strong("Choose dataset")," field. You need to specify the plot you are interested from following options: Main plot, Main Heatmap, QC Heatmap. Finally you can click ", strong("Download Data"), " button to download data, or if you wish to see the selected data, you can click ",strong("Tables")," tab.")
  )
}


#' getHistogramUI
#'
#' Histogram plots UI.
#'
#' @note \code{getHistogramUI}
#' @param id, namespace id
#' @return the panel for PCA plots;
#'
#' @examples
#'     x <- getHistogramUI("histogram")
#'
#' @export
#'
getHistogramUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("histogramUI"))
}

#' debrowserhistogram
#'
#' Module for a histogram that can be used in data prep and
#' low count removal modules
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @return histogram
#' @export
#'
#' @examples
#'     x <- debrowserhistogram()
#'
debrowserhistogram <- function(input = NULL, output = NULL, session = NULL, data = NULL) {
  if(is.null(data)) return(NULL)
  output$histogram <- renderPlotly({

    h <- hist(log10(rowSums(data)), breaks = as.numeric(input$breaks), plot = FALSE)

    p <- plot_ly(x = h$mids, y = h$counts,
                 width = input$width, height=input$height) %>%
      add_bars() %>%
      plotly::layout(
        margin = list(l = input$left,
                      b = input$bottom,
                      t = input$top,
                      r = input$right
        ))
    p$elementId <- NULL
    p
  })
  output$histogramUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("plot"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotlyOutput(session$ns("histogram"),
                                      width = input$width, height=input$height))
  })
}

#' histogramControlsUI
#'
#' Generates the controls in the left menu for a histogram
#'
#' @note \code{histogramControlsUI}
#' @param id, namespace id
#' @return returns the left menu
#' @examples
#'     x <- histogramControlsUI("histogram")
#' @export
#'
histogramControlsUI <- function(id) {
  ns <- NS(id)
  textInput(ns("breaks"), "Breaks", value = "100" )
}




#' getGoPanel
#'
#' Creates go term analysis panel within the shiny
#' display.
#'
#' @note \code{getGoPanel}
#' @return the panel for go term analysis;
#'
#' @examples
#'     x <- getGoPanel()
#'
#' @export
#'
getGoPanel <- function(){
  gopanel <- list(
    wellPanel(helpText( "Please select parameters and press the
                    submit button in the left menu for the plots"),
              getHelpButton("method",
                            "http://debrowser.readthedocs.io/en/master/examples/examples.html#go-term-plots")),
    tabsetPanel(id = "gotabs", type = "tabs",
                tabPanel(title = "Plot", value = "gopanel1", id = "gopanel1",
                         column(12, wellPanel( plotOutput("GOPlots1")))),
                tabPanel(title = "Table", value = "gopanel2", id = "gopanel2",
                         column(12, wellPanel( DT::dataTableOutput("gotable"))))
    ),
    getKEGGModal(),
    getTableModal()
  )
  return(gopanel)
}



#' getGOPlots
#'
#' Go term analysis panel.  Generates appropriate GO plot
#' based on user selection.
#'
#' @param dataset, the dataset used
#' @param input, input params
#' @note \code{getGOPlots}
#' @return the panel for go plots;
#'
#' @examples
#'     x<- getGOPlots()
#' @export
#'
getGOPlots <- function(dataset = NULL, input = NULL){
  if (is.null(dataset)) return(NULL)
  goplots <- NULL
  org <- input$organism
  if (input$goplot == "disease")
    org <- "org.Hs.eg.db"
  genelist <- getGeneList(rownames(dataset), org)
  gopval <- as.numeric(input$gopvalue)
  if (input$goplot == "enrichGO"){
    res <- getEnrichGO(genelist, ont = input$ontology,
                       pvalueCutoff = gopval, org = input$organism)
    goplots<-res
    if (input$goextplot == "Dotplot")
      goplots$p <- dotplot(res$enrich_p, showCategory=30)
  }
  else if (input$goplot == "enrichKEGG"){
    res <- getEnrichKEGG(genelist, pvalueCutoff=
                           gopval, org = input$organism)
    goplots<-res
    if (input$goextplot == "Dotplot")
      goplots$p <- dotplot(res$enrich_p, showCategory=30)
  }
  else if (input$goplot == "compare"){
    cl <- clusterData(dataset)
    res <- compareClust(cl, fun=input$gofunc, input$ontology,
                        org = input$organism)
    goplots<-res
  }
  else if (input$goplot == "disease"){
    res <- getEnrichDO(genelist, pvalueCutoff = gopval )
    goplots<-res
    if (input$goextplot == "Dotplot")
      goplots$p <- dotplot(res$enrich_p, showCategory=30)
  }
  return(goplots)
}

#' getOrganismBox
#'
#' Get the organism Box.
#"
#' @note \code{getOrganismBox}
#'
#' @export
#'
#' @note \code{getOrganismBox}
#' makes the organism box
#' @return selectInput
#'
#' @examples
#'     x <- getOrganismBox()
#'
getOrganismBox <- function(){
  organismBox <- list(
    conditionalPanel( ( condition <- "input.goplot!='disease' &&
                            input.gofunc != 'enrichDO'"),
                      selectInput("organism", "Choose an organism:",
                                  choices =  c( "Human" = "org.Hs.eg.db",
                                                "Mouse" = "org.Mm.eg.db",
                                                "Rat" = "org.Rn.eg.db",
                                                "Zebrafish" = "org.Dr.eg.db",
                                                "Fly" = "org.Dm.eg.db",
                                                "Worm" = "org.Ce.eg.db",
                                                "Yeast" = "org.Sc.sgd.db",
                                                "Arabidopsis" = "org.At.tair.db"
                                  ))))
  return(organismBox)
}

#' getOrganism
#'
#' @param org, organism
#' @note \code{getOrganism}
#'
#' @export
#' @return organism name for keg
#'
#' @examples
#'     x <- getOrganism()
#'
getOrganism <- function(org){
  organisms <-  list("hsa", "mmu", "rno",
                     "dre", "dme", "cel", "sce", "At")
  names(organisms) <- c("org.Hs.eg.db",
                        "org.Mm.eg.db",
                        "org.Rn.eg.db",
                        "org.Dr.eg.db",
                        "org.Dm.eg.db",
                        "org.Ce.eg.db",
                        "org.Sc.sgd.db",
                        "org.At.tair.db")
  organisms[org][[1]]
}

#' getOrganismPathway
#'
#' @param org, organism
#' @note \code{getOrganismPathway}
#'
#' @export
#' @return organism name for pathway
#'
#' @examples
#'     x <- getOrganismPathway()
#'
getOrganismPathway <- function(org){
  organisms <- list("human", "mouse", "rat",
                    "zebrafish", "fly", "celegans", "yeast", "arabidopsis")
  names(organisms) <- c("org.Hs.eg.db",
                        "org.Mm.eg.db",
                        "org.Rn.eg.db",
                        "org.Dr.eg.db",
                        "org.Dm.eg.db",
                        "org.Ce.eg.db",
                        "org.Sc.sgd.db",
                        "org.At.tair.db")
  organisms[org][[1]]
}




#' debrowsercondselect
#'
#' Condition selection
#' This is not a module. Module construction didn't used here, just use it
#' as functions not in a module.
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, count data
#' @param metadata, metadata
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowsercondselect()
#'
debrowsercondselect <- function(input = NULL, output = NULL, session = NULL, data = NULL, metadata = NULL) {
  if (is.null(data)) return(NULL)
  choicecounter <- reactiveVal(0)
  output$conditionSelector <- renderUI({
    selectConditions(data, metadata, choicecounter(), session, input)
  })

  observeEvent(input$add_btn, {
    choicecounter(choicecounter() + 1)
  })
  observeEvent(input$rm_btn, {
    if (choicecounter() > 0)
      choicecounter(choicecounter() - 1)
  })
  list(cc = choicecounter)
}

#' condSelectUI
#' Creates a panel to select samples for each condition
#'
#' @return panel
#' @examples
#'     x <- condSelectUI()
#'
#' @export
#'
condSelectUI<- function () {
  list(
    shinydashboard::box(title = "Comparison Selection",
                        solidHeader = TRUE, status = "info",  width = NULL, height = NULL, collapsible = TRUE,
                        fluidRow(
                          uiOutput("conditionSelector"),
                          column(12,actionButtonDE("add_btn", "Add New Comparison",styleclass = "primary"),
                                 actionButtonDE("rm_btn", "Remove", styleclass = "primary"),
                                 getHelpButton("method", "http://debrowser.readthedocs.io/en/master/deseq/deseq.html"),
                                 conditionalPanel(condition = ("output.condReady>0"),
                                                  actionButtonDE("startDE", "Start DE", styleclass = "primary")))
                        ))
  )
}
#' getMethodDetails
#'
#' get the detail boxes after DE method selected
#'
#' @param num, panel that is going to be shown
#' @param input, user input
#' @examples
#'     x <- getMethodDetails()
#'
#' @export
#'
#'
getMethodDetails <- function(num = 0, input = NULL) {
  if (num > 0)
    list(
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'DESeq2'")),
        getSelectInputBox("fitType", "Fit Type", num,
                          c("parametric", "local", "mean"),
                          selectedInput("testType", num, "parametric",
                                        input), 3),
        getSelectInputBox("betaPrior", "betaPrior", num,
                          c(FALSE, TRUE),
                          selectedInput("betaPrior", num,
                                        FALSE, input),2),
        getSelectInputBox("testType", "Test Type", num,
                          c("LRT", "Wald"),
                          selectedInput("testType", num, "LRT", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," == 'EdgeR'")),
        getSelectInputBox("edgeR_normfact", "Normalization", num,
                          c("TMM","RLE","upperquartile","none"),
                          selectedInput("edgeR_normfact", num, "TMM", input), 3),
        column(2,textInput(paste0("dispersion", num), "Dispersion",
                           value = isolate(selectedInput("dispersion",
                                                         num, "0", input) ))),
        getSelectInputBox("edgeR_testType", "Test Type", num,
                          c("exactTest", "glmLRT"),
                          selectedInput("edgeR_testType", num,
                                        "exactTest", input))),
      conditionalPanel(
        (condition <- paste0("input.demethod",num," ==  'Limma'")),
        getSelectInputBox("limma_normfact", "Normalization", num,
                          c("TMM","RLE","upperquartile","none"),
                          selectedInput("limma_normfact", num, "TMM", input), 3),
        getSelectInputBox("limma_fitType", "Fit Type", num,
                          c("ls", "robust"),
                          selectedInput("limma_fitType", num, "ls", input)),
        getSelectInputBox("normBetween", "Norm. Bet. Arrays", num,
                          c("none", "scale", "quantile", "cyclicloess",
                            "Aquantile", "Gquantile", "Rquantile","Tquantile"),
                          selectedInput("normBetween", num, "none", input))),
      br())
}

#' getConditionSelector
#'
#' Selects user input conditions to run in DESeq.
#'
#' @param num, panel that is going to be shown
#' @param choices, sample list
#' @param selected, selected smaple list
#' @examples
#'     x <- getConditionSelector()
#'
#' @export
#'
getConditionSelector<- function(num=0, choices = NULL, selected = NULL) {
  if (!is.null(choices))
    list(column(6, selectInput(paste0("condition", num),
                               label = paste0("Condition ", num),
                               choices = choices, multiple = TRUE,
                               selected = selected),
                textInput(paste0("conditionname", num),
                          label=paste0("Name the Condition ",num),
                          value=paste0("Condition ", num))))
}

#' getConditionSelectorFromMeta
#'
#' Selects user input conditions to run in DESeq from metadata
#'
#' @param metadata, meta data table
#' @param input, input
#' @param index, index
#' @param num, num
#' @param choices, choices
#' @param selected, selected
#'
#' @examples
#'     x <- getConditionSelectorFromMeta()
#'
#' @export
#'





getConditionSelectorFromMeta <- function(metadata = NULL, input = NULL, index = 1, num=0,
                                         choices = NULL, selected = NULL) {
  a <- list(column(6, selectInput(paste0("condition", num),
                                  label = paste0("Condition ", num),
                                  choices = choices, multiple = TRUE,
                                  selected = selected),
                   textInput(paste0("conditionname", num),
                             label=paste0("Name the Condition ",num),
                             value=paste0("Condition ", num))))

  if (!is.null(metadata)){
    selected_meta <- selectedInput("conditions_from_meta",
                                   index, NULL, input)

    if (is.null(selected_meta)) selected_meta <- "No Selection"

    if (selected_meta != "No Selection"){
      old_selection <- ""

      if(!is.null(input[[paste0("condition", num)]])){
        selected <- input[[paste0("condition", num)]]
      }
      meta_choices_all <- NULL
      if (!is.null(selected_meta))
        meta_choices_all <- get_conditions_given_selection(metadata,
                                                           selected_meta)
      if(old_selection != selected_meta){
        if(typeof(meta_choices_all) == "character"){
          meta_choices <- list("There must be exactly 2 groups.")
        } else{
          meta1 <- meta_choices_all[[2 - (num %% 2)]]
          meta_choices <- unlist(meta1, recursive=FALSE)
        }
        selected <- meta_choices
        condname<- isolate(input[[paste0("conditionname",num)]])
      }

      a <- list(column(6, selectInput(paste0("condition", num),
                                      label = paste0("Condition ", num),
                                      choices = choices, multiple = TRUE,
                                      selected = selected),

                       textInput(paste0("conditionname", num),
                                 label=paste0("Name the Condition ",num),
                                 value= condname)))

    }
  }
  return(a)

}


#' selectedInput
#'
#' Selects user input conditions to run in DESeq.
#'
#' @param id, input id
#' @param num, panel that is going to be shown
#' @param default, default text
#' @param input, input params
#' @examples
#'     x <- selectedInput()
#'
#' @export
#'
selectedInput <- function(id = NULL, num = 0, default = NULL,
                          input = NULL) {
  if (is.null(id)) return(NULL)
  m <- NULL
  if (is.null(input[[paste0(id, num)]]))
    m <- default
  else
    m <- input[[paste0(id, num)]]
  m
}

#' getSelectInputBox
#'
#' Selects user input conditions to run in DESeq.
#'
#' @param id, input id
#' @param name, label of the box
#' @param num, panel that is going to be shown
#' @param choices, sample list
#' @param selected, selected smaple list
#' @param cw, column width
#' @examples
#'     x <- getSelectInputBox()
#'
#' @export
#'
getSelectInputBox <- function(id = NULL, name = NULL,
                              num = 0, choices = NULL, selected = NULL,
                              cw = 2) {
  if (is.null(id)) return(NULL)
  if (!is.null(choices))
    list(column(cw, selectInput(paste0(id, num),
                                label = name,
                                choices = choices, multiple = FALSE,
                                selected = selected)))
}


#' selectConditions
#'
#' Selects user input conditions, multiple if present, to be
#' used in DESeq.
#'
#' @param Dataset, used dataset
#' @param metadata, metadatatable to select from metadata
#' @param choicecounter, choicecounter to add multiple comparisons
#' @param session, session
#' @param input, input params
#' @note \code{selectConditions}
#' @return the panel for go plots;
#'
#' @examples
#'     x<- selectConditions()
#'
#' @export
#'
selectConditions<-function(Dataset = NULL,
                           metadata = NULL,
                           choicecounter = NULL,
                           session = NULL,
                           input = NULL) {
  if (is.null(Dataset)) return(NULL)

  selectedSamples <- function(num){
    if (is.null(input[[paste0("condition", num)]]))
      getSampleNames(colnames(Dataset), num %% 2 )
    else
      input[[paste0("condition", num)]]
  }
  nc <- choicecounter

  if (nc >= 0) {
    allsamples <- getSampleNames( colnames(Dataset), "all" )

    lapply(seq_len(nc), function(i) {

      selected1 <- selectedSamples(2 * i - 1)
      selected2 <- selectedSamples( 2 * i )
      to_return <- list(column(12, getMetaSelector(metadata = metadata, input=input, n = i),
                               getConditionSelectorFromMeta(metadata, input, i,
                                                            (2 * i - 1), allsamples, selected1),
                               getConditionSelectorFromMeta(metadata, input, i,
                                                            (2 * i), allsamples, selected2)
      ),
      column(12,
             column(1, helpText(" ")),
             getSelectInputBox("demethod", "DE Method", i,
                               c("DESeq2", "EdgeR", "Limma"),
                               selectedInput("demethod", i, "DESeq2", input)),
             getMethodDetails(i, input)))
      if (!is.null(selectedInput("conditions_from_meta",
                                 i, NULL, input)) && selectedInput("conditions_from_meta",
                                                                   i, NULL, input) != "No Selection"){
        facts <- levels(factor(metadata[,selectedInput("conditions_from_meta",
                                                       i, NULL, input)]))
        facts <- facts[facts != "" & facts != "NA"]
        if (length(facts) != 2) {
          showNotification("There must be exactly 2 groups in the selected condition.
                         Please use NA or space to remove extra sample groups from metadata selection.",
                           type = "error")
          updateSelectInput(session, paste0("conditions_from_meta", i), selected="No Selection" )
        }
      }
      return(to_return)
    })
  }
}

#' getMetaSelector
#'
#' Return the sample selection box using meta data table
#'
#' @param metadata, meta data table
#' @param input, input params
#' @param n, the box number
#' @return meta select box
#'
#' @examples
#'     x<-getMetaSelector()
#' @export
#'
getMetaSelector <- function(metadata = NULL, input = NULL, n = 0){
  if(!is.null(metadata)){
    df <- metadata
    col_count <- length(colnames(df))

    list(HTML('<hr style="color: white; border:solid 1px white;">'),
         br(), column(10, selectInput(paste0("conditions_from_meta",
                                             n), label = "Select Meta",
                                      choices = as.list(c("No Selection",
                                                          colnames(df)[2:col_count])),
                                      multiple = FALSE,
                                      selected =  selectedInput("conditions_from_meta",
                                                                n, "Selection 2", input))))
  }
}

#' get_conditions_given_selection
#'
#' Return the two set of conditions given the selection of meta select box
#'
#' @param metadata, meta data table
#' @param selection, selection
#' @return meta select box
#'
#' @examples
#'     x<-get_conditions_given_selection()
#' @export
#'
get_conditions_given_selection <- function(metadata = NULL, selection = NULL){
  if(is.null(metadata)) return(NULL)
  df <- metadata
  if(selection == "No Selection"){
    return(NULL)
  }
  facts <- levels(factor(df[,selection]))
  facts <- facts[facts != "" & facts != "NA"]
  if(length(facts) != 2){
    return(NULL)
  } else {
    # Assuming the first column has samples
    sample_col_name <- colnames(df)[1]

    condition1 <- facts[1]
    condition2 <- facts[2]

    # In case the conditions are integers
    if(is.null(condition2)){
      condition1 <- factor(condition1)
      condition2 <- factor(condition2)
    }

    condition1_filtered <- df[df[,selection] == condition1, ]
    a <- condition1_filtered[,sample_col_name]

    condition2_filtered <- df[df[,selection] == condition2, ]
    b <- condition2_filtered[,sample_col_name]

    both_groups <- list(a, b)
    return(both_groups)
  }
}

#' getSampleNames
#'
#' Prepares initial samples to fill condition boxes.
#' it reads the sample names from the data and splits into two.
#'
#' @param cnames, sample names in the header of a dataset
#' @param part, c(1,2). 1=first half and 2= second half
#' @return sample names.
#'
#' @examples
#'     x<-getSampleNames()
#' @export
#'
getSampleNames <- function(cnames = NULL, part = 1) {
  if (is.null(cnames)) return(NULL)

  startpos <- 1
  endpos <- length(cnames)
  if (part == 1)
    endpos <- floor(endpos / 2)
  else if (part == 0)
    startpos <- floor(endpos / 2) + 1

  cn <- cnames[startpos:endpos]
  m <- as.list(NULL)
  for (i in seq(cn)) {
    m[[i]] <- cn[i]
  }
  m
}

#' prepDataContainer
#'
#' Prepares the data container that stores values used within DESeq.
#'
#' @param data, loaded dataset
#' @param counter, the number of comparisons
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataContainer()
#'
prepDataContainer <- function(data = NULL, counter=NULL,
                              input = NULL) {
  if (is.null(data)) return(NULL)

  inputconds <- reactiveValues(demethod_params = list(), conds_1= list(),cond_name= list(),conds = list(), dclist = list())

  inputconds$conds_1 <- list()

  inputconds$cond_name<-list()
  for (cnt in seq(1:(2*counter))){
    inputconds$conds_1[cnt] <- list(isolate(input[[paste0("condition",cnt)]]))
    inputconds$cond_name[cnt]<- list(isolate(input[[paste0("conditionname",cnt)]]))
  }
  #Get parameters for each method
  inputconds$demethod_params <- NULL
  for (cnt in seq(1:counter)){
    if (isolate(input[[paste0("demethod",cnt)]]) == "DESeq2"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("fitType",cnt)]]),
        isolate(input[[paste0("betaPrior",cnt)]]),
        isolate(input[[paste0("testType",cnt)]]), sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "EdgeR"){
      inputconds$demethod_params[cnt]<- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("edgeR_normfact",cnt)]]),
        isolate(input[[paste0("dispersion",cnt)]]),
        isolate(input[[paste0("edgeR_testType",cnt)]]), sep=",")
    }
    else if (isolate(input[[paste0("demethod",cnt)]]) == "Limma"){
      inputconds$demethod_params[cnt] <- paste(
        isolate(input[[paste0("demethod",cnt)]]),
        isolate(input[[paste0("limma_normfact",cnt)]]),
        isolate(input[[paste0("limma_fitType",cnt)]]),
        isolate(input[[paste0("normBetween",cnt)]]), sep=",")
    }
  }

  for (i in seq(1:counter)){

    conds <- c(rep(paste(inputconds$cond_name[2*i-1]),
                   length(inputconds$conds_1[[2*i-1]])),
               rep(paste(inputconds$cond_name[2*i]), length(inputconds$conds_1[[2*i]])))
    cols <- c(paste(inputconds$conds_1[[2*i-1]]),
              paste(inputconds$conds_1[[2*i]]))
    params <- unlist(strsplit(inputconds$demethod_params[i], ","))
    withProgress(message = 'Running DE Algorithms', detail = inputconds$demethod_params[i], value = 0, {
      initd <- callModule(debrowserdeanalysis, paste0("DEResults",i), data = data,
                          columns = cols, conds = conds, params = params)
      if (nrow(initd$dat()) > 1){
        inputconds$dclist[[i]] <- list(conds = conds, cols = cols, init_data=initd$dat(),
                                       demethod_params = inputconds$demethod_params[i])
      }
      incProgress(1/counter)
    })
  }

  if(length(inputconds$dclist) <1) return(NULL)

  inputconds$dclist
}




#' getSampleDetails
#'
#' get sample details
#'
#' @param output, output
#' @param summary, summary output name
#' @param details, details ouput name
#' @param data, data
#' @return panel
#' @examples
#'     x <- getSampleDetails()
#'
#' @export
#'
getSampleDetails<- function (output = NULL, summary = NULL, details = NULL, data = NULL) {
  if (is.null(data)) return(NULL)

  output[[summary]]<- renderTable({
    countdata <-  data$count
    samplenums <- length(colnames(countdata))
    rownums <- dim(countdata)[1]
    result <- rbind(samplenums, rownums)
    rownames(result) <- c("# of samples", "# of rows (genes/regions)")
    colnames(result) <- "Value"
    result
  },digits=0, rownames = TRUE, align="lc")

  output[[details]] <- DT::renderDataTable({
    dat <- colSums(data$count)
    dat <- cbind(names(dat), dat)
    dat[, c("dat")] <-  format(
      round( as.numeric( dat[,  c("dat")], digits = 2)),
      big.mark=",",scientific=FALSE)

    if (!is.null(data$meta)){
      met <- data$meta
      dat <- cbind(met, dat[,"dat"])
      rownames(dat) <- NULL
      colnames(dat)[ncol(dat)] <- "read counts"
    }else{
      rownames(dat) <- NULL
      colnames(dat) <- c("samples", "read counts")
    }
    dat
  })
}



#' selectGroupInfo
#'
#' Group info column selection. This can be used in batch effect
#' or coloring the groups in the plots.
#' @param metadata, metadata
#' @param input, input values
#' @param selectname, name of the select box
#' @param label, label of the select box
#' @note \code{selectGroupInfo}
#' @examples
#'     x <- selectGroupInfo()
#' @export
#'
selectGroupInfo <- function(metadata = NULL, input = NULL,
                            selectname = "groupselect",
                            label = "Group info") {
  if (is.null(metadata)) return (NULL)
  lst.choices <- as.list(c("None", colnames(metadata)))
  selectInput(selectname, label = label,
              choices = lst.choices,
              selected = 1)
}


#' addID
#'
#' Adds an id to the data frame being used.
#'
#' @param data, loaded dataset
#' @return data
#' @export
#'
#' @examples
#'     x <- addID()
#'
addID <- function(data = NULL) {
  if (is.null(data)) return (NULL)
  dat1 <- data.frame(data)
  dat1 <- cbind(rownames(data), data)
  colnames(dat1) <- c("ID", colnames(data))
  dat1
}

#' getVariationData
#'
#' Adds an id to the data frame being used.
#'
#' @param inputdata, dataset
#' @param cols, columns
#' @param conds, conditions
#' @param key, gene or region name
#' @return plotdata
#' @export
#'
#' @examples
#'     x <- getVariationData()
#'
getVariationData <- function(inputdata = NULL,
                             cols = NULL, conds = NULL, key = NULL) {
  if (is.null(inputdata)) return (NULL)
  # Pick out the gene with this ID
  vardata <- inputdata[key, ]
  bardata <- as.data.frame(cbind(key, cols,
                                 t(vardata[, cols]), conds) )
  colnames(bardata) <- c("genename", "libs", "count", "conds")
  bardata$count <- as.numeric(as.character(bardata$count))
  data <- rbind(bardata[bardata$conds == levels(bardata$conds)[1], ],
                bardata[bardata$conds == levels(bardata$conds)[2], ])
  data$conds  <- factor(data$conds, levels = unique(data$conds))
  data
}

#' getBSTableUI
#' prepares a Modal to put a table
#' @param name, name
#' @param label, label
#' @param trigger, trigger button for the modal
#' @param size, size of the modal
#' @param modal, modal yes/no
#' @return the modal
#'
#' @examples
#'     x<- getBSTableUI()
#'
#' @export
getBSTableUI<-function(name = NULL,  label = NULL, trigger = NULL, size="large", modal = NULL){
  if (is.null(name)) return (NULL)
  ret <- div(style = "display:block;overflow-y:auto; overflow-x:auto;",
             wellPanel( DT::dataTableOutput(name)))
  if (!is.null(modal) && modal)
    ret <- shinyBS::bsModal(name, label, trigger, size = size, ret)
  ret
}

#' getTableDetails
#'
#' get table details
#' To be able to put a table into two lines are necessary;
#' into the server part;
#' getTableDetails(output, session, "dataname", data, modal=TRUE)
#' into the ui part;
#' uiOutput(ns("dataname"))
#'
#' @param output, output
#' @param session, session
#' @param tablename, table name
#' @param data, matrix data
#' @param modal, if it is true, the matrix is going to be in a modal
#' @return panel
#' @examples
#'     x <- getTableDetails()
#'
#' @export
#'
getTableDetails <- function(output  = NULL, session  = NULL, tablename  = NULL, data = NULL, modal = NULL){
  if (is.null(data)) return(NULL)
  tablenameUI <-  paste0(tablename,"Table")
  output[[paste(tablename, "Download")]] <- downloadHandler(
    filename = function() {
      paste0(tablename,".tsv")
    },
    content = function(file) {
      if(!("ID" %in% names(data)))
        data <- addID(data)
      write.table(data, file, sep = "\t", row.names = FALSE)
    }
  )

  output[[tablename]] <- renderUI({
    ret <- getBSTableUI( session$ns(tablenameUI), "Show Data", paste0("show",tablename), modal = modal)
    if (!is.null(modal) && modal)
      ret <- list( downloadButton(session$ns(paste(tablename, "Download")), "Download"),
                   actionButtonDE(paste0("show",tablename), "Show Data", styleclass = "primary", icon="show"),
                   ret)
    ret
  })

  output[[tablenameUI]] <- DT::renderDataTable({
    if (!is.null(data)){
      DT::datatable(data, extensions = 'Buttons',
                    options = list( server = TRUE,
                                    dom = "Blfrtip",
                                    buttons =
                                      list("copy", list(
                                        extend = "collection"
                                        , buttons = c("csv", "excel", "pdf")
                                        , text = "Download"
                                      ) ), # end of buttons customization

                                    # customize the length menu
                                    lengthMenu = list( c(10, 20,  50, -1) # declare values
                                                       , c(10, 20, 50, "All") # declare titles
                                    ), # end of lengthMenu customization
                                    pageLength = 10))
    }
  })
}

#' push
#'
#' Push an object to the list.
#'
#' @param l, that are going to push to the list
#' @param ..., list object
#' @return combined list
#'
#' @export
#'
#' @examples
#'     mylist <- list()
#'     newlist <- push ( 1, mylist )
push <- function(l, ...) c(l, list(...))

#' round_vals
#'
#' Plot PCA results.
#'
#' @param l, the value
#' @return round value
#' @export
#'
#' @examples
#'     x<-round_vals(5.1323223)
round_vals <- function(l) {
  l <- round(as.numeric(l), digits = 2)
  parse(text = l)
}

#' Buttons including Action Buttons and Event Buttons
#'
#' Creates an action button whose value is initially zero, and increments by one
#' each time it is pressed.
#'
#' @param inputId Specifies the input slot that will be used to access the
#'   value.
#' @param label The contents of the button--usually a text label, but you could
#'   also use any other HTML, like an image.
#' @param styleclass The Bootstrap styling class of the button--options are
#'   primary, info, success, warning, danger, inverse, link or blank
#' @param size The size of the button--options are large, small, mini
#' @param block Whehter the button should fill the block
#' @param icon Display an icon for the button
#' @param css.class Any additional CSS class one wishes to add to the action
#'   button
#' @param ... Other argument to feed into shiny::actionButton
#'
#' @export
#'
#' @examples
#'     actionButtonDE("goDE", "Go to DE Analysis")
#'
actionButtonDE <- function(inputId, label, styleclass = "", size = "",
                           block = FALSE, icon = NULL, css.class = "", ...) {
  if (styleclass %in% c("primary", "info", "success", "warning",
                        "danger", "inverse", "link")) {
    btn.css.class <- paste("btn", styleclass, sep = "-")
  } else btn.css.class = ""

  if (size %in% c("large", "small", "mini")) {
    btn.size.class <- paste("btn", size, sep = "-")
  } else btn.size.class = ""

  if (block) {
    btn.block = "btn-block"
  } else btn.block = ""

  if (!is.null(icon)) {
    icon.code <- HTML(paste0("<i class='fa fa-", icon, "'></i>"))
  } else icon.code = ""
  tags$button(id = inputId, type = "button", class = paste("btn action-button",
                                                           btn.css.class, btn.size.class, btn.block, css.class, collapse = " "),
              icon.code, label, ...)
}


#' getNormalizedMatrix
#'
#' Normalizes the matrix passed to be used within various methods
#' within DEBrowser.  Requires edgeR package
#'
#' @note \code{getNormalizedMatrix}
#' @param M, numeric matrix
#' @param method, normalization method for edgeR. default is TMM
#' @return normalized matrix
#'
#' @examples
#'     x <- getNormalizedMatrix(mtcars)
#'
#' @export
#'
getNormalizedMatrix <- function(M = NULL, method = "TMM") {
  if (is.null(M) ) return (NULL)
  M[is.na(M)] <- 0
  norm <- M
  if (!(method == "none" || method == "MRN")){
    norm.factors <- edgeR::calcNormFactors(M, method = method)
    norm <- edgeR::equalizeLibSizes(edgeR::DGEList(M,
                                                   norm.factors = norm.factors))$pseudo.counts
  }else if(method == "MRN"){
    columns <- colnames(M)
    conds <- columns
    coldata <- prepGroup(conds, columns)
    M[, columns] <- apply(M[, columns], 2,
                          function(x) as.integer(x))
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(M),
                                  colData = coldata, design = ~group)
    dds <- estimateSizeFactors(dds)
    norm <- counts(dds, normalized=TRUE)
  }
  return(norm)
}

#' getCompSelection
#'
#' Gathers the user selected comparison set to be used within the
#' DEBrowser.
#' @param name, the name of the selectInput
#' @param count, comparison count
#' @note \code{getCompSelection}
#' @examples
#'     x <- getCompSelection(name="comp", count = 2)
#' @export
#'
getCompSelection <- function(name = NULL, count = NULL) {
  a <- NULL
  if (count>1){
    a <- list(selectInput(name,
                          label = "Choose a comparison:",
                          choices = c(1:count) ))
  }
  a
}


#' getHelpButton
#' prepares a helpbutton for to go to a specific site in the documentation
#'
#' @param name, name that are going to come after info
#' @param link, link of the help
#' @return the info button
#'
#' @examples
#'     x<- getHelpButton()
#'
#' @export
getHelpButton<-function(name = NULL, link = NULL){
  if (is.null(name)) return(NULL)
  btn <- actionButtonDE(paste0("info_",name),"",icon="info",
                        styleclass="info", size="small")

  HTML(paste0("<a id=\"info_",name,"\" href=\"",link,"\" target=\"_blank\">",
              btn,"</a>"))

}

#' getDomains
#'
#' Get domains for the main plots.
#'
#' @param filt_data, data to get the domains
#' @return domains
#' @export
#'
#' @examples
#'     x<-getDomains()
getDomains <- function(filt_data = NULL){
  if (is.null(filt_data)) return (NULL)
  a <- unique(filt_data$Legend)
  a <- a[a != ""]
  if (length(a) == 1)
    a <- c(a, "NA")
  a
}

#' getColors
#'
#' get colors for the domains
#'
#' @param domains, domains to be colored
#' @return colors
#' @export
#'
#' @examples
#'     x<-getColors()
#'
getColors <- function(domains = NULL){
  if (is.null(domains)) return (NULL)
  colors <- c()
  for ( dn in seq(1:length(domains)) ){
    if (domains[dn] == "NS" || domains[dn] == "NA")
      colors <- c(colors, "#aaa")
    else if (domains[dn] == "Up")
      colors <- c(colors, "green")
    else if (domains[dn] == "Down")
      colors <- c(colors, "red")
    else if (domains[dn] == "MV")
      colors <- c(colors, "orange")
    else if (domains[dn] == "GS")
      colors <- c(colors, "blue")
  }
  colors
}

#' getKEGGModal
#' prepares a modal for KEGG plots
#'
#' @return the info button
#'
#' @examples
#'     x<- getKEGGModal()
#'
#' @export
getKEGGModal<-function(){
  bsModal("modalExample", "KEGG Pathway", "KeggPathway", size = "large",
          div(style = "display:block;overflow-y:auto; overflow-x:auto;",imageOutput("KEGGPlot")))
}

#' getTableModal
#' prepares table modal for KEGG
#'
#' @return the info button
#'
#' @examples
#'     x<- getTableModal()
#'
#' @export
getTableModal<-function(){
  bsModal("modalTable", "Genes in the category", "GeneTableButton", size = "large",
          div(style = "display:block;overflow-y:auto; overflow-x:auto;",
              wellPanel( DT::dataTableOutput("GOGeneTable"))))
}




#' setBatch
#' to skip batch effect correction batch variable set with the filter results
#'
#' @param fd, filtered data
#' @return fd data
#'
#' @examples
#'
#'     x <- setBatch()
#'
#' @export
#'
setBatch <- function(fd = NULL){
  if(!is.null(fd)){
    batchdata <- reactiveValues(count=NULL, meta = NULL)
    batchdata$count <-  fd$filter()$count
    batchdata$meta <-  fd$filter()$meta
    batcheffectdata <- reactive({
      ret <- NULL
      if(!is.null(batchdata$count)){
        ret <- batchdata
      }
      return(ret)
    })
    list(BatchEffect=batcheffectdata)
  }
}

#' setFilter
#' to skip Filter filtd variable set with the upload results
#'
#' @param fd, filtered data
#' @return fd data
#'
#' @examples
#'
#'     x <- setfilter()
#'
#' @export
#'

setFilter <- function(fd = NULL){
  if(!is.null(fd)){
    batchdata <- reactiveValues(count=NULL, meta = NULL)
    batchdata$count <-  fd$load()$count
    batchdata$meta <-  fd$load()$meta
    batcheffectdata <- reactive({
      ret <- NULL
      if(!is.null(batchdata$count)){
        ret <- batchdata
      }
      return(ret)
    })
    list(BatchEffect=batcheffectdata)
  }
}

#' getTabUpdateJS
#' prepmenu tab and discovery menu tab updates
#'
#' @return the JS for tab updates
#'
#' @examples
#'     x<- getTabUpdateJS()
#'
#' @export
getTabUpdateJS<-function(){
  tags$script(HTML( "
                      $(function() {
                      $('#methodtabs').attr('selectedtab', '2')
                      $($('#methodtabs >')[0]).attr('id', 'dataprepMethod')
                      $($('#menutabs >')[0]).attr('id', 'dataprepMenu')
                      for(var i=1;i<=5;i++){
                      $($('#methodtabs >')[i]).attr('id', 'discoveryMethod')
                      }
                      $($('#menutabs >')[1]).attr('id', 'discoveryMenu')
                      $(document).on('click', '#dataprepMethod', function () {
                      if($('#dataprepMenu').attr('class')!='active'){
                      $('#dataprepMenu').find('a').click()
                      }
                      });
                      $(document).on('click', '#dataprepMenu', function () {
                      if($('#dataprepMethod').attr('class')!='active'){
                      $('#dataprepMethod').find('a').click()
                      }
                      });
                      $(document).on('click', '#discoveryMethod', function () {
                      $('#methodtabs').attr('selectedtab', $(this).index())
                      if($('#discoveryMenu').attr('class')!='active'){
                      $('#discoveryMenu').find('a').click()
                      }
                      });
                      $('#discoveryMenu > ').css('display', 'none');
                      $(document).on('click', '#goMain', function () {
                      $('#discoveryMenu > ').css('display', 'block');
                      });
                      $(document).on('click', '#discoveryMenu', function () {
                      $($('#methodtabs >')[ $('#methodtabs').attr('selectedtab')]).find('a').click()
                      });
                      //hide buttons on entrance
                      $('.sidebar-menu > ').css('display', 'none');
                      $('.sidebar-menu > :nth-child(1)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(2)').css('display', 'inline');
                      $(document).on('click', '#Next', function () {
                      $('.sidebar-menu > :nth-child(2)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(3)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(4)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#Batch', function () {
                      $('.sidebar-menu > :nth-child(4)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#goDEFromFilter', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#goDEwithoutFilter', function () {
                      $('.sidebar-menu > :nth-child(2)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(3)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#goDE', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#startDE', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(2)').css('display', 'none');
                      });
                      })
                      "))
}

#' getPCAcontolUpdatesJS
#' in the prep menu we have two PCA plots to show how batch effect correction worked.
#' One set of PCA input controls updates two PCA plots with this JS.
#' @return the JS for tab updates
#'
#' @examples
#'     x<- getTabUpdateJS()
#'
#' @export
getPCAcontolUpdatesJS<-function(){
  tags$script(HTML("
                        var nameInputs = ['pcselx', 'pcsely'];
                        $.each(nameInputs, function (el) {
                                              $(function () {
                                              $(document).on('change keyup', '#batcheffect-beforeCorrectionPCA-'+nameInputs[el], function () {
                                              var value = $(this).val();
                                              $('#batcheffect-afterCorrectionPCA-'+nameInputs[el]).val(value)
                                              $('#batcheffect-afterCorrectionPCA-'+nameInputs[el]).trigger('change');
                                              });
                                              });
                        });

                        var nameDropdowns = ['textonoff', 'legendonoff', 'legendSelect', 'color_pca','shape_pca'];
                        $.each(nameDropdowns, function (el) {
                            $(function () {
                                $(document).on('change', '#batcheffect-beforeCorrectionPCA-'+nameDropdowns[el], function () {
                                    var value = $(this).val();
                                    $('#batcheffect-afterCorrectionPCA-'+nameDropdowns[el])[0].selectize.setValue(value)
                                    $('#batcheffect-afterCorrectionPCA-'+nameDropdowns[el]).trigger('change');
                                });
                            });
                        });
                     $($('#batcheffect-pcacontrols')[0]).css('display', 'none');
                     $(document).on('click','#batcheffect-submitBatchEffect',function () {
                     setTimeout(function () { $($($('#batcheffect-pcacontrols')[0]).children()[1]).children().trigger('click')
                     setTimeout(function () { $($($('#batcheffect-pcacontrols')[0]).children()[0]).children().trigger('click')}, 1000);
                     }, 1000); })
                     "))
}

.initial <- function() {
  req <- function(...){
    reqFun <- function(pack) {
      if(!suppressWarnings(suppressMessages(require(pack, character.only = TRUE)))) {
        message(paste0("unable to load package ", pack))
        require(pack, character.only = TRUE)
      }
    }
    lapply(..., reqFun)
  }
  packs <- c("debrowser", "plotly", "shiny", "jsonlite", "shinyjs", "shinydashboard", "shinyBS")
  req(packs)
}

.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription("debrowser", fields="Version")
  msg <- paste0("DEBrowser v", pkgVersion, "  ",
                "For help: https://debrowser.readthedocs.org/", "\n\n")

  citation <- paste0("If you use DEBrowser in published research, please cite:\n\n",
                     "Alper Kucukural, Onur Yukselen, Deniz M. Ozata, Melissa J. Moore, Manuel Garber\n",
                     "DEBrowser: Interactive Differential Expression Analysis and Visualization Tool for Count Data\n",
                     "BMC Genomics 2019 20:6\n\ndoi:0.1186/s12864-018-5362-x\n")

  packageStartupMessage(paste0(msg, citation))
  .initial()

}





#' debrowsermainplot
#'
#' Module for a scatter, volcano and ma plots that are going to be used
#' as a mainplot in debrowser
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowsermainplot()
#'
debrowsermainplot <- function(input = NULL, output = NULL, session = NULL, conds=NULL, data = NULL) {
  if (is.null(data)) return(NULL)

  plotdata <-  reactive({
    plotData(data, input)
  })
  output$mainplot <- renderUI({
    list(fluidRow(
      column(12,
             shinydashboard::box(
               collapsible = TRUE, title = "Main Plots", status = "primary",
               solidHeader = TRUE,width = NULL,
               draggable = TRUE, plotlyOutput(session$ns("main"),
                                              height=input$plotheight, width=input$plotwidth)
             ))))
  })



  output$mainPlotControlsUI <- renderUI({
    if (input$mainplot == "scatter"){
      a<-levels(factor(conds))[1]
      b<-levels(factor(conds))[2]
      x <- paste0('log10 Norm. Mean(Read Counts) in ',a)
      y <- paste0('log10 Norm. Mean(Read Counts) in ',b)
    }else if  (input$mainplot == "volcano"){
      x <- "log2FC"
      y <- "log10padj"
    }else {
      x <- "A"
      y <- "M"
    }
    list(
      textInput(session$ns('xlab'),'x label', x),
      textInput(session$ns('ylab'),'y label', y)
    )
  })
  selectedPoint <- reactive({
    eventdata <- event_data("plotly_click", source = session$ns("source"))
    key <- ""
    if (!is.null(eventdata$key))
      key <- as.vector(unlist(eventdata$key))

    return(key)
  })

  getSelected  <- reactive({
    keys <- NULL
    selGeneList <- event_data("plotly_selected", source = session$ns("source"))
    if (is.null(selGeneList$key)) return (NULL)
    keys <- as.vector(unlist(selGeneList$key))
    return(keys)
  })

  output$main <- renderPlotly({
    data <- plotdata()$data
    mainScatterNew(input, data, conds, session$ns("source"))
  })

  list( shg = (selectedPoint), shgClicked=(selectedPoint), selGenes=(getSelected))
}

#' getMainPlotUI
#'
#' main plot for volcano, scatter and maplot.
#' @param id, namespace id
#' @note \code{getMainPlotUI}
#' @return the panel for main plots;
#'
#' @examples
#'     x <- getMainPlotUI("main")
#'
#' @export
#'
getMainPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("mainplot"))
}

#' mainScatterNew
#'
#' Creates the main scatter, volcano or MA plot to be displayed within the main
#' panel.
#' @param input, input params
#' @param data, dataframe that has log2FoldChange and log10padj values
#' @param source, for event triggering to select genes
#' @return scatter, volcano or MA plot
#'
#' @examples
#'
#'     x <- mainScatterNew()
#'
#' @export
#'
mainScatterNew <- function(input = NULL, data = NULL, conds=NULL, source = NULL) {
  if ( is.null(data) ) return(NULL)



  p <- plot_ly(source = source, data=data, x=~x, y=~y, key=~key, alpha = 0.8,
               color=~Legend, colors=getLegendColors(unique(data$Legend)),
               type="scatter", mode = "markers",
               width=input$width - 100, height=input$height,
               text=~paste("<b>", ID, "</b><br>",
                           "<br>", "padj=", format.pval(padj, digits = 2), " ",
                           "-log10padj=", round(log10padj, digits = 2),
                           "<br>", "log2FC=", round(log2FoldChange, digits = 2), " ",
                           "foldChange=", round(foldChange, digits = 2),
                           "<br>", sep = " ")) %>%
    plotly::layout(xaxis = list(title = input$xlab),
                   yaxis = list(title = input$ylab)) %>%
    plotly::layout(
      margin = list(l = input$left,
                    b = input$bottom,
                    t = input$top,
                    r = input$right
      ))
  p$elementId <- NULL

  return(p)
}

#' plotData
#'
#' prepare plot data for mainplots
#'
#' @note \code{plotData}
#' @param pdata, data
#' @param input, input
#' @return prepdata
#' @examples
#'     x <- plotData()
#' @export
#'
plotData <- function(pdata = NULL, input = NULL){
  if (is.null(pdata)) return(NULL)
  pdata$key <- pdata$ID
  data_rest <- pdata[ pdata$Legend!="NS",]
  data_NS <- pdata[ pdata$Legend=="NS",]
  backperc <- 10
  if (!is.null(input$backperc))  backperc <- input$backperc
  mainplot <- "scatter"
  if (!is.null(input$mainplot))  mainplot <- input$mainplot

  datapoints <- as.integer(nrow(data_NS) * backperc/ 100)
  if (nrow(data_NS) > datapoints){
    data_rand <- data_NS[sample(1:nrow(data_NS), datapoints,
                                replace=FALSE),]
  }else{
    data_rand  <- data_NS
  }
  plot_init_data <- rbind(data_rand, data_rest)
  plot_init_data$Legend  <- factor(plot_init_data$Legend,
                                   levels = unique(plot_init_data$Legend))

  plot_data <- plot_init_data
  if (mainplot == "volcano") {
    plot_data <- plot_init_data[which(!is.na(plot_init_data$log2FoldChange)
                                      & !is.na(plot_init_data$log10padj)
                                      & !is.na(plot_init_data$Legend)),]
    plot_data$x <- plot_data$log2FoldChange
    plot_data$log10padj[plot_data$log10padj>input$log10padjCutoff] <- input$log10padjCutoff
    plot_data$y <- plot_data$log10padj
  } else if (mainplot == "maplot") {
    plot_data$x <- (plot_init_data$x + plot_init_data$y) / 2
    plot_data$y <- plot_init_data$y - plot_init_data$x
  }
  list( data = (plot_data))
}

#' mainPlotControlsUI
#'
#' Generates the left menu to be used for main plots
#'
#' @note \code{mainPlotControlsUI}
#' @param id, module ID
#' @return mainPlotControls
#' @examples
#'     x <- mainPlotControlsUI("main")
#' @export
#'
mainPlotControlsUI <- function(id) {
  ns <- NS(id)
  list(shinydashboard::menuItem(" Plot Type",
                                startExpanded=TRUE,
                                radioButtons(ns("mainplot"), "Main Plots:",
                                             c(Scatter = "scatter", VolcanoPlot = "volcano",
                                               MAPlot = "maplot"))
  ),
  shinydashboard::menuItem("Main Options",
                           startExpanded=TRUE,
                           sliderInput(ns("backperc"), "Background Data(%):",
                                       min=10, max=100, value=10, sep = "",
                                       animate = FALSE),
                           conditionalPanel(condition <- paste0("input['", ns("mainplot"),"'] == 'volcano'"),
                                            sliderInput(ns("log10padjCutoff"), "Log10 padj value cutoff:",
                                                        min=2, max=100, value=60, sep = "",
                                                        animate = FALSE)
                           ),
                           uiOutput(ns("mainPlotControlsUI"))
  ))

}

#' getLegendColors
#'
#' Generates colors according to the data
#'
#' @note \code{getLegendColors}
#' @param Legend, unique Legends
#' @return mainPlotControls
#' @examples
#'     x <- getLegendColors(c("up", "down", "GS", "NS"))
#' @export
#'

getLegendColors<-function(Legend=c("up", "down", "NS"))
{
  colors <- c()
  for(i in seq(1:length(Legend))){
    if (Legend[i]=="Up"){
      colors <- c(colors, "red")
    }
    else if (Legend[i]=="Down"){
      colors <- c(colors, "blue")
    }
    else if (Legend[i]=="NS"){
      colors <- c(colors, "grey")
    }
    else if (Legend[i]=="GS"){
      colors <- c(colors, "green")
    }
  }
  colors
}

#' getLevelOrder
#'
#' Generates the order of the overlapping points
#'
#' @note \code{getLevelOrder}
#' @param Level, factor levels shown in the legend
#' @return order
#' @examples
#'     x <- getLevelOrder(c("up", "down", "GS", "NS"))
#' @export
#'
getLevelOrder<-function(Level=c("up", "down", "NS"))
{
  levels <- c( "NS", "Up", "Down", "GS")
  for(i in seq(1:length(levels)))
  {
    if(!levels[i] %in% Level){
      levels <- levels[-(i)]
    }
  }
  levels
}

#' generateTestData
#'
#' This generates a test data that is suitable to main plots in debrowser
#' @param dat, DESeq results will be generated for loaded data
#' @return testData
#'
#' @examples
#'     x <- generateTestData()
#'
#' @export
#'
generateTestData <- function(dat = NULL) {
  if (is.null(dat)) return (NULL)
  ##################################################
  columns <- dat$columns
  conds <- dat$conds
  data <- dat$data
  params <-
    #Run DESeq2 with the following parameters
    c("DESeq2", "parametric", F, "Wald")
  non_expressed_cutoff <- 10
  data <- subset(data, rowSums(data) > 10)
  deseqrun <- runDE(data, columns, conds, params)

  met <- as.data.frame(cbind(as.vector(conds), columns))
  colnames(met) <- c("conds", "columns")
  cols1 <- as.vector(met[met$conds==as.character(unique(met$conds)[1]), "columns"])
  cols2 <- as.vector(met[met$conds==as.character(unique(met$conds)[2]), "columns"])

  de_res <- data.frame(deseqrun)
  norm_data <- getNormalizedMatrix(data[, columns])
  rdata <- cbind(rownames(de_res), norm_data[rownames(de_res), columns],
                 log10(rowMeans(norm_data[rownames(de_res),cols1])
                       + 0.1), log10( rowMeans( norm_data[ rownames( de_res ), cols2])
                                      + 0.1), de_res[rownames(de_res),
                                                     c("padj", "log2FoldChange")], 2 ^ de_res[rownames(de_res),                                                                                            "log2FoldChange"], -1 *
                   log10(de_res[rownames(de_res), "padj"]))
  colnames(rdata) <- c("ID", columns, "x", "y", "padj",
                       "log2FoldChange", "foldChange", "log10padj")
  rdata <- as.data.frame(rdata)
  rdata$padj[is.na(rdata$padj)] <- 1

  padj_cutoff <- 0.01
  foldChange_cutoff <- 2


  rdata$Legend <- "NS"
  rdata$Legend[rdata$log2FoldChange > log2(foldChange_cutoff) &
                 rdata$padj < padj_cutoff] <- "Up"
  rdata$Legend[rdata$log2FoldChange <= log2(1 / foldChange_cutoff) &
                 rdata$padj < padj_cutoff] <- "Down"

  dat <- rdata
  dat$M <- rdata$x - rdata$y
  dat$A <- (rdata$x + rdata$y) / 2
  dat
}




#' debrowserIQRplot
#'
#' Module for an IQR plot that can be used in data prep and
#' low count removal modules
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @return IQR
#' @export
#'
#' @examples
#'     x <- debrowserIQRplot()
#'
debrowserIQRplot <- function(input = NULL, output = NULL, session = NULL, data = NULL) {
  if (is.null(data)) return(NULL)
  output$IQR <- renderPlotly({
    getIQRPlot(data, input)
  })
  output$IQRUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("plot"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotlyOutput(session$ns("IQR"),
                                      width = input$width, height=input$height))
  })
}

#' getIQRPlotUI
#'
#' IQR plot UI.
#' @param id, namespace id
#' @note \code{getIQRPlotUI}
#' @return the panel for IQR plots;
#'
#' @examples
#'     x <- getIQRPlotUI("IQR")
#'
#' @export
#'
getIQRPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("IQRUI"))
}

#' IQRPlotControlsUI
#'
#' Generates the controls in the left menu for an IQR plot#'
#' @param id, namespace id
#' @note \code{IQRPlotControlsUI}
#' @return returns the left menu
#' @examples
#'     x <- IQRPlotControlsUI("IQR")
#' @export
#'
IQRPlotControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(id, " - Options"),
                           textInput(ns("breaks"), "Breaks", value = "100" )
  )
}


#' getIQRPlot
#'
#' Makes IQR boxplot plot
#'
#' @param data, count or normalized data
#' @param input, input
#' @param title, title
#'
#' @export
#'
#' @examples
#'     getIQRPlot()
#'
getIQRPlot <- function(data=NULL, input=NULL, title = ""){
  if (is.null(data)) return(NULL)
  data <- as.data.frame(data)
  cols <- colnames(data)
  data[, cols] <- apply(data[, cols], 2,
                        function(x) log10(as.integer(x) + 1))

  data <- addID(data)
  mdata <- melt(as.data.frame(data[,c("ID", cols)]),"ID")
  colnames(mdata)<-c("ID", "samples", "logcount")

  p <- plot_ly(mdata, x = ~samples, y = ~logcount, type = "box",
               width = input$width, height=input$height,
               marker = list(color = 'rgb(8,81,156)',
                             outliercolor = 'rgba(219, 64, 82, 0.6)',
                             line = list(outliercolor = 'rgba(219, 64, 82, 1.0)',
                                         outlierwidth = 2))) %>%
    plotly::layout(title = title,
                   xaxis = list(title = "samples"),
                   yaxis = list(title = "logcount"))
  if (!is.null(input$left))
    p <- p %>% plotly::layout(margin = list(l = input$left,
                                            b = input$bottom,
                                            t = input$top,
                                            r = input$right
    ))
  p$elementId <- NULL
  p
}


#' debrowserbarmainplot
#'
#' Module for a bar plot that can be used in data prep, main plots
#' low count removal modules or any desired module
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @param cols, columns
#' @param conds, conditions
#' @param key, the gene or region name
#' @return density plot
#' @export
#'
#' @examples
#'     x <- debrowserbarmainplot()
#'
debrowserbarmainplot <- function(input, output, session, data = NULL,
                                 cols = NULL, conds=NULL, key=NULL) {
  if(is.null(data)) return(NULL)
  output$BarMainUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("plot"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotlyOutput(session$ns("BarMain"),
                                      height=input$height, width=input$width))
  })
  output$BarMain <- renderPlotly({
    getBarMainPlot(data, cols, conds, key, title = "", input =input)
  })
}

#' getBarMainPlotUI
#'
#' main bar plots UI.
#'
#' @note \code{getBarMainPlotUI}
#' @param id, namespace id
#' @return the panel for Density plots;
#'
#' @examples
#'     x <- getBarMainPlotUI("bar")
#'
#' @export
#'
getBarMainPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("BarMainUI"))
}


#' barMainPlotControlsUI
#'
#' Generates the controls in the left menu for a bar main plot
#'
#' @note \code{barMainPlotControlsUI}
#' @param id, namespace id
#' @return returns the controls for left menu
#' @examples
#'     x <- barMainPlotControlsUI("bar")
#' @export
#'
barMainPlotControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(id, " - Options"),
                           textInput(ns("genename"), "Gene/Region Name", value = "Foxa3" )
  )
}

#' getBarMainPlot
#'
#' Makes Density plots
#'
#' @param data, count or normalized data
#' @param cols, cols
#' @param conds, conds
#' @param key, key
#' @param title, title
#' @param input, input
#' @export
#'
#' @examples
#'     getBarMainPlot()
#'
getBarMainPlot <- function(data=NULL, cols = NULL, conds=NULL, key=NULL, title = "", input = NULL){
  if (is.null(data)) return(NULL)
  vardata <- getVariationData(data, cols, conds, key)
  title <- paste(key, "variation")

  p <- plot_ly(vardata, x = ~libs, y = ~count,
               color=~conds, colors=c("Blue", "Red"),
               type = "bar")
  p <- p %>%
    plotly::layout(title = title,
                   xaxis = list(categoryorder = "array",
                                categoryarray = "conds",
                                title = "Conditions"),
                   yaxis = list(title = "Read Count"),
                   height=input$height, width=input$width,
                   margin = list(l = input$left,
                                 b = input$bottom,
                                 t = input$top,
                                 r = input$right
                   ))
  p$elementId <- NULL
  p
}



#' getLeftMenu
#'
#' Generates the left menu for for plots within the DEBrowser.
#'
#' @param input, input values
#' @note \code{getLeftMenu}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- getLeftMenu()
#' @export
#'
getLeftMenu <- function(input = NULL) {
  if (is.null(input)) return(NULL)
  leftMenu <- list(
    conditionalPanel( (condition <- "input.methodtabs=='panel1'"),
                      getMainPlotsLeftMenu()),
    conditionalPanel( (condition <- "input.methodtabs=='panel2'"),
                      shinydashboard::menuItem(" Plot Type", startExpanded = TRUE,
                                               wellPanel(radioButtons("qcplot",
                                                                      paste("QC Plots:", sep = ""),
                                                                      c(PCA = "pca", All2All = "all2all", Heatmap = "heatmap", IQR = "IQR",
                                                                        Density = "Density")))),
                      getQCLeftMenu(input)),
    conditionalPanel( (condition <- "input.methodtabs=='panel3'"),

                      shinydashboard::menuItem(" Plot Type", startExpanded = TRUE,
                                               wellPanel(radioButtons("goplot", paste("Go Plots:", sep = ""),
                                                                      c(enrichGO = "enrichGO", enrichKEGG = "enrichKEGG",
                                                                        Disease = "disease", compareClusters = "compare")))),
                      actionButton("startGO", "Submit"),
                      getGOLeftMenu()
    ),
    conditionalPanel( (condition <- "input.methodtabs=='panel4'"),
                      shinydashboard::menuItem(" Select Columns", startExpanded=TRUE,
                                               uiOutput("getColumnsForTables")
                      )),
    conditionalPanel( (condition <- "input.methodtabs=='panel5'"),
                      shinydashboard::menuItem("Concordance Analysis - Gene selection", startExpanded=TRUE,
                                               wellPanel(
                                                 uiOutput("getselectionforcon" ),
                                                 uiOutput("getselectionforbiomarker"),
                                                 uiOutput("getselectionforcompound"))),

                      actionButton("conplotselect", "Apply"))
  )
  return(leftMenu)
}

#' getMainPlotsLeftMenu
#'
#' Generates the Main PLots Left menu to be displayed within the DEBrowser.
#'
#' @note \code{getMainPlotsLeftMenu}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- getMainPlotsLeftMenu()
#' @export
#'
getMainPlotsLeftMenu <- function() {
  mainPlotsLeftMenu <- list(
    plotSizeMarginsUI("main",  w=600, h=400),
    shinydashboard::menuItem("Heatmap Options", startExpanded=FALSE,
                             heatmapControlsUI("heatmap"),
                             plotSizeMarginsUI("heatmap", w=550, h=400)),
    plotSizeMarginsUI("barmain", w=550,h=400, t=90),
    plotSizeMarginsUI("boxmain", w=550, h=400, t=90)
  )
  return(mainPlotsLeftMenu)
}

#' getGOLeftMenu
#'
#' Generates the GO Left menu to be displayed within the DEBrowser.
#'
#' @note \code{getGOLeftMenu}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- getGOLeftMenu()
#' @export
#'
getGOLeftMenu <- function() {
  list(
    shinydashboard::menuItem(" Go Term Options", startExpanded=TRUE,
                             textInput("gopvalue", "p.adjust", value = "0.01" ),
                             getOrganismBox(),
                             actionButton("GeneTableButton", "DE Genes"),
                             conditionalPanel( (condition <- "input.goplot=='enrichKEGG'"),
                                               actionButton("KeggPathway", "KeggPathway")),
                             conditionalPanel( ( condition <- "(input.goplot=='enrichGO' ||
            (input.goplot=='compare' && input.gofunc!='enrichDO' &&
            input.gofunc!='enrichKEGG'))" ),
                                               selectInput("ontology", "Choose an ontology:",
                                                           choices =  c( "CC", "MF", "BP"))
                             ),
                             conditionalPanel( ( condition <- "input.goplot!='compare'"),
                                               selectInput("goextplot", "Plot Type:",
                                                           choices =  c("Summary", "Dotplot"))
                             ),
                             conditionalPanel( ( condition <- "input.goplot=='compare'"),
                                               selectInput("gofunc", "Plot Function:",
                                                           choices =  c( "enrichGO", "enrichDO", "enrichKEGG"))
                             )),
    actionButton("goapply","Apply"),
    downloadButton("downloadGOPlot", "Download Plots")
  )

}

#' getQCLeftMenu
#'
#' Generates the left menu to be used for QC plots within the
#' DEBrowser.
#'
#' @param input, input values
#' @note \code{getQCLeftMenu}
#' @return QC left menu
#' @examples
#'     x <- getQCLeftMenu()
#' @export
#'
getQCLeftMenu <- function( input = NULL) {
  if (is.null(input)) return(NULL)
  list(
    shinydashboard::menuItem(" Select Columns", startExpanded=TRUE,
                             uiOutput("columnSelForQC")),
    shinydashboard::menuItem(" QC Options", startExpanded=FALSE,
                             conditionalPanel( (condition <- "input.qcplot=='heatmap'"),
                                               plotSizeMarginsUI("heatmapQC"),
                                               heatmapControlsUI("heatmapQC")),
                             conditionalPanel( condition <- "(input.qcplot=='all2all')",
                                               plotSizeMarginsUI("all2all"),
                                               all2allControlsUI("all2all")
                             ),
                             conditionalPanel( condition <- "(input.qcplot=='Density')",
                                               plotSizeMarginsUI("density"),
                                               plotSizeMarginsUI("normdensity")
                             ),
                             conditionalPanel( condition <- "(input.qcplot=='IQR')",
                                               plotSizeMarginsUI("IQR"),
                                               plotSizeMarginsUI("normIQR")
                             ),
                             getHelpButton("method",
                                           "http://debrowser.readthedocs.io/en/master/heatmap/heatmap.html"),
                             conditionalPanel( (condition <- "input.qcplot=='pca'"),
                                               shinydashboard::menuItem("PCA Options",
                                                                        pcaPlotControlsUI("qcpca")),
                                               plotSizeMarginsUI("qcpca", w=600, h=400, t=0, b=0, l=0, r=0)
                             ))
  )
}

#' getCutOffSelection
#'
#' Gathers the cut off selection for DE analysis
#'
#' @param nc, total number of comparisons
#' @note \code{getCutOffSelection}
#' @return returns the left menu according to the selected tab;
#' @examples
#'     x <- getCutOffSelection()
#' @export
#'
getCutOffSelection <- function(nc = 1){
  compselect <- getCompSelection("compselect", nc)
  compselected<-
    list( conditionalPanel( (condition = "input.dataset!='most-varied' &&
        input.methodtabs!='panel0'"),
                            shinydashboard::menuItem(" Filter",
                                                     #h4("Filter"),
                                                     compselect,
                                                     textInput("padj", "padj", value = "0.01" ),
                                                     textInput("foldChange", "foldChange", value = "2" ),
                                                     actionButton("apply", "Apply")
                            )
    ) )
}

#' getMainPanel
#'
#' main panel for volcano, scatter and maplot.
#' Barplot and box plots are in this page as well.
#'
#' @note \code{getMainPanel}
#' @return the panel for main plots;
#'
#' @examples
#'     x <- getMainPanel()
#'
#' @export
#'
getMainPanel <- function() {
  list(
    fluidRow(column(6,
                    getMainPlotUI("main")
    ),
    column(6,
           getHeatmapUI("heatmap")
    )),
    fluidRow(column(6,
                    getBarMainPlotUI("barmain")),
             column(6,
                    getBoxMainPlotUI("boxmain"))))
}

#' getProgramTitle
#'
#' Generates the title of the program to be displayed within DEBrowser.
#' If it is called in a program, the program title will be hidden
#'
#' @param session, session var
#' @note \code{getProgramTitle}
#' @return program title
#' @examples
#'     title<-getProgramTitle()
#' @export
#'
getProgramTitle <- function(session = NULL) {
  if (is.null(session)) return (NULL)
  DEBrowser <- NULL
  title<-parseQueryString(session$clientData$url_search)$title
  if (is.null(title) || title != "no" )
    DEBrowser <- list(titlePanel("DEBrowser"))
  else
    DEBrowser <- list(titlePanel(" "))
  return(DEBrowser)
}

#' getLoadingMsg
#'
#' Creates and displays the loading message/gif to be displayed
#' within the DEBrowser.
#'
#' @param output, output message
#' @note \code{getLoadingMsg}
#' @return loading msg
#' @examples
#'     x <- getLoadingMsg()
#' @export
#'
getLoadingMsg <- function(output = NULL) {
  addResourcePath(prefix = "www", directoryPath =
                    system.file("extdata", "www",
                                package = "debrowser"))
  imgsrc_full <- "www/images/loading_start.gif"
  imgsrc_small <- "www/images/loading.gif"
  a <- list(
    tags$head(tags$style(type = "text/css", "
            #loadmessage {
            position: fixed;
            top: 0px;
            left: 0px;
            width: 100%;
            height: 100%;
            padding: 5px 0px 5px 0px;
            text-align: center;
            font-weight: bold;
            font-size: 100%;
            color: #000000;
            opacity: 0.8;
            z-index: 100;
            }
            #loadmessage_small {
            position: fixed;
            left: 50%;
            transform: translateX(-50%);
            top: 50px;
            text-align: center;
            opacity: 0.8;
            z-index: 999999;
            }
                             ")),
    conditionalPanel(condition = paste0("$('html').hasClass('shiny-busy')",
                                        "& input.startDE & input.methodtabs=='panel0'"),
                     tags$div(id = "loadmessage",
                              tags$img(src = imgsrc_full
                              ))),
    conditionalPanel(condition =  paste0("$('html').hasClass('shiny-busy')",
                                         "& !(input.startDE & input.methodtabs=='panel0')"),
                     tags$div(id = "loadmessage_small",
                              tags$img(src = imgsrc_small
                              )))
  )
}

#' getLogo
#'
#' Generates and displays the logo to be shown within DEBrowser.
#'
#' @note \code{getLogo}
#' @return return logo
#' @examples
#'     x <- getLogo()
#' @export
#'
getLogo <- function(){
  addResourcePath(prefix = "www", directoryPath =
                    system.file("extdata", "www",
                                package = "debrowser"))
  imgsrc <- "www/images/logo.png"
  a<-list(img(src=imgsrc, align = "right"))
}

#' getStartupMsg
#'
#' Generates and displays the starting message within DEBrowser.
#'
#' @note \code{getStartupMsg}
#' @return return startup msg
#' @examples
#'     x <- getStartupMsg()
#' @export
#'
getStartupMsg <- function() {
  a <- list( column( 12,
                     helpText("Please select a file or load the demo data."),
                     helpText( "For more information;" ),
                     helpText(   a("Quick Start Guide",
                                   href = "http://debrowser.readthedocs.org",
                                   target = "_blank"),
                                 getHelpButton("method", "http://debrowser.readthedocs.org")) ))
}

#' getAfterLoadMsg
#'
#' Generates and displays the message to be shown after loading data
#' within the DEBrowser.
#'
#' @note \code{getAfterLoadMsg}
#' @return return After Load Msg
#' @examples
#'     x <- getAfterLoadMsg()
#' @export
#'
getAfterLoadMsg <- function() {
  a <- list( column( 12, wellPanel(
    helpText( "Please choose the appropriate conditions for DESeq analysis
            and press 'Run DESeq' button in the left menu" ),
    helpText( "To be able to select conditions please click
            'Condition1' or 'Condition2' boxes.
            You can also use delete button to remove the
            samples from the list."))))
}

#' getStartPlotsMsg
#'
#' Generates and displays the starting messgae to be shown once
#' the user has first seen the main plots page within DEBrowser.
#'
#' @note \code{getStartPlotsMsg}
#' @return return start plot msg
#' @examples
#'     x <- getStartPlotsMsg()
#' @export
#'
getStartPlotsMsg <- function() {
  a <- list( conditionalPanel(condition <- "!input.goMain",
                              column( 12,
                                      helpText( "Please choose the appropriate parameters to discover
               more in DE Results" ),
                                      getHelpButton("method", "http://debrowser.readthedocs.io/en/master/quickstart/quickstart.html"))))
}

#' getCondMsg
#'
#' Generates and displays the current conditions and their samples
#' within the DEBrowser.
#'
#' @param dc, columns
#' @param input, selected comparison
#' @param cols, columns
#' @param conds, selected conditions
#' @note \code{getCondMsg}
#' @return return conditions
#' @examples
#'     x <- getCondMsg()
#' @export
#'
getCondMsg <- function(dc = NULL, input = NULL, cols = NULL, conds = NULL) {
  if (is.null(cols) || is.null(conds)) return (NULL)
  num <- input$compselect
  if (is.null(num)) num <- 1
  cnd <- data.frame(cbind(conds, cols))
  params_str <- paste(dc[[as.numeric(num)]]$demethod_params, collapse = ',')
  heatmap_str <-  paste0( "<b>Heatmap Params: Scaled:</b> ", input[['heatmap-scale']],
                          " <b>Centered:</b> ", input[['heatmap-center']],
                          " <b>Log:</b> ", input[['heatmap-log']],
                          " <b>Pseudo-count:</b> ", input[['heatmap-pseudo']])
  a <-list( conditionalPanel(condition <- "input.goMain",
                             shinydashboard::box(
                               collapsible = TRUE, title = "Plot Information", status = "primary",
                               solidHeader = TRUE, width = NULL,
                               draggable = TRUE,
                               style = "overflow-x:scroll",
                               HTML( paste0( "<b>DE Params:</b> ", params_str,
                                             " - <b>Dataset:</b> ", input$dataset," <b>Normalization:</b> ",input$norm_method,
                                             " - ", heatmap_str,
                                             "</br><b>",unique(conds)[1], ":</b> "),
                                     paste(cnd[cnd$conds == unique(conds)[1], "cols"],
                                           collapse =","),
                                     paste0(" vs. ","<b>",unique(conds)[2], ":", "</b> "),
                                     paste(cnd[cnd$conds == unique(conds)[2], "cols"],
                                           collapse =",")),
                               getHelpButton("method",
                                             "http://debrowser.readthedocs.io/en/master/quickstart/quickstart.html#the-main-plots-of-de-analysis")
                             )))
}

#' togglePanels
#'
#' User defined toggle to display which panels are to be shown within
#' DEBrowser.
#'
#' @param num, selected panel
#' @param nums, all panels
#' @param session, session info
#' @note \code{togglePanels}
#' @examples
#'     x <- togglePanels()
#' @export
#'
togglePanels <- function(num = NULL, nums = NULL, session = NULL){
  if (is.null(num)) return (NULL)
  for(i in 0:4){
    if (i %in% nums)
      shinyjs::show(selector =
                      paste0("#methodtabs li a[data-value=panel",i,"]"))
    else
      shinyjs::hide(selector =
                      paste0("#methodtabs li a[data-value=panel",i,"]"))
  }
  if(num)
    updateTabsetPanel(session, "methodtabs",
                      selected = paste0("panel", num))
}


#' getTableStyle
#'
#' User defined selection that selects the style of table to display
#' within the DEBrowser.
#'
#' @param dat, dataset
#' @param input, input params
#' @param padj, the name of the padj value column in the dataset
#' @param foldChange, the name of the foldChange column in the dataset
#' @param DEsection, if it is in DESection or not
#' @note \code{getTableStyle}
#' @examples
#'     x <- getTableStyle()
#' @export
#'
getTableStyle <- function(dat = NULL, input = NULL,
                          padj = c("padj"), foldChange=c("foldChange"), DEsection = TRUE){
  if (is.null(dat)) return (NULL)
  a <- dat
  if(!is.null(padj) && padj != "" && DEsection)
    a <- a %>% formatStyle(
      padj,
      color = styleInterval(c(0, input$padj),
                            c('black', "white", "black")),
      backgroundColor = styleInterval(
        input$padj, c('green', 'white'))
    )
  if(!is.null(foldChange) && foldChange != "" && DEsection)
    a <- a %>%
    formatStyle(
      foldChange,
      color = styleInterval(c(1/as.numeric(input$foldChange),
                              as.numeric(input$foldChange)), c('white', 'black', 'white')),
      backgroundColor = styleInterval(
        c(1/as.numeric(input$foldChange),
          as.numeric(input$foldChange)),
        c('blue', 'white', 'red'))
    )
  a
}

#' textareaInput
#'
#' Generates a text area input to be used for gene selection within
#' the DEBrowser.
#'
#' @param id, id of the control
#' @param label, label of the control
#' @param value, initial value
#' @param rows, the # of rows
#' @param cols, the # of  cols
#' @param class, css class
#' @examples
#'     x <- textareaInput("genesetarea", "Gene Set",
#'         "Fgf21", rows = 5, cols = 35)
#' @export
#'
textareaInput <- function(id, label, value, rows=20, cols=35,
                          class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

#' showObj
#'
#' Displays a shiny object.
#'
#' @param btns, show group of objects with shinyjs
#' @examples
#'     x <- showObj()
#' @export
#'
showObj <- function(btns = NULL) {
  if (is.null(btns)) return (NULL)
  for (btn in seq(1:length(btns)))
    shinyjs::show(btns[btn])
}

#' hideObj
#'
#' Hides a shiny object.
#'
#' @param btns, hide group of objects with shinyjs
#' @examples
#'     x <- hideObj()
#' @export
#'
hideObj <- function(btns = NULL) {
  if (is.null(btns)) return (NULL)
  for (btn in seq(1:length(btns)))
    shinyjs::hide(btns[btn])
}

#' getKEGGModal
#' prepares a helpbutton for to go to a specific site in the documentation
#'
#' @return the info button
#'
#' @examples
#'     x<- getKEGGModal()
#'
#' @export
getKEGGModal<-function(){
  bsModal("modalExample", "KEGG Pathway", "KeggPathway", size = "large",
          div(style = "display:block;overflow-y:auto; overflow-x:auto;",imageOutput("KEGGPlot")))
}

#' getDownloadSection
#'
#' download section button and dataset selection box in the
#' menu for user to download selected data.
#'
#' @param choices, main vs. QC section
#'
#' @note \code{getDownloadSection}
#' @return the panel for download section in the menu;
#'
#' @examples
#'     x<- getDownloadSection()
#'
#' @export
#'
getDownloadSection <- function(choices=NULL) {
  list(conditionalPanel( (condition = "input.methodtabs!='panel0'"),
                         shinydashboard::menuItem(" Data Options",
                                                  selectInput("dataset", "Choose a dataset:",
                                                              choices = choices),
                                                  conditionalPanel( (condition = "input.dataset=='selected'"),
                                                                    selectInput("selectedplot", "The plot used in selection:",
                                                                                choices = c("Main Plot", "Main Heatmap", "QC Heatmap"))),
                                                  selectInput("norm_method", "Normalization Method:",
                                                              c("none", "MRN", "TMM", "RLE", "upperquartile"), selected = "MRN"),
                                                  downloadButton("downloadData", "Download Data"),
                                                  conditionalPanel(condition = "input.dataset=='most-varied'",
                                                                   textInput("topn", "top-n", value = "500" ),
                                                                   textInput("mincount", "total min count", value = "10" )),
                                                  textareaInput("genesetarea","Search",
                                                                "", rows = 5, cols = 35),
                                                  helpText("Regular expressions can be used\n
        Ex: ^Al => Al.., Al$ => ...al"),
                                                  actionButton("getdownload","Submit")
                         )))
}

#' getQCPanel
#'
#' Gathers the conditional panel for QC plots
#'
#' @param input, user input
#' @note \code{getQCSection}
#' @return the panel for QC plots
#'
#' @examples
#'     x <- getQCPanel()
#'
#' @export
#'
getQCPanel <- function(input = NULL) {
  height = "700"
  width = "500"
  if (!is.null(input)) {
    height = input$height
    width = input$width
  }
  qcPanel <- list(
    wellPanel(helpText( "Please select the parameters and press the
                            submit button in the left menu for the plots" ),
              getHelpButton("method",
                            "http://debrowser.readthedocs.io/en/master/quickstart/quickstart.html#quality-control-plots")),
    conditionalPanel(condition = "input.qcplot == 'pca'",
                     getPCAPlotUI("qcpca")),
    conditionalPanel(condition = "(input.qcplot == 'heatmap')",
                     getHeatmapUI("heatmapQC")),
    conditionalPanel(condition = "(input.qcplot == 'IQR')",
                     getIQRPlotUI("IQR"),
                     getIQRPlotUI("normIQR")),
    conditionalPanel(condition = "(input.qcplot == 'Density')",
                     getDensityPlotUI("density"),
                     getDensityPlotUI("normdensity")),
    conditionalPanel(condition = "(input.qcplot == 'all2all')",
                     getAll2AllPlotUI("all2all"))
  )
  return(qcPanel)
}

#' getSelectedCols
#'
#' gets selected columns
#'
#' @param data, all loaded data
#' @param datasetInput, selected dataset
#' @param input, user input params
#'
#' @export
#'
#' @examples
#'     getSelectedCols()
#'
#'
getSelectedCols <- function(data = NULL, datasetInput = NULL, input=NULL){
  if(is.null(data) || is.null(datasetInput)) return(NULL)
  selCols <- NULL
  if (!is.null(input$dataset)){
    selection <- colnames(data)
    if (!is.null(input$col_list))
      selection <- input$col_list

    selection <- selection[selection %in% colnames(data)]

    if (!is.null(selection))
      selCols <- data[rownames(datasetInput), selection]
  }
  return(selCols)
}


#' removeExtraCols
#'
#' remove extra columns for QC plots
#'
#' @param dat, selected data
#'
#' @export
#'
#' @examples
#'     removeExtraCols()
#'
#'
removeExtraCols <- function(dat = NULL){
  rcols <- c(names(dat)[grep("^padj", names(dat))],
             names(dat)[grep("^foldChange", names(dat))],
             names(dat)[grep("^log2FoldChange$", names(dat))],
             names(dat)[grep("^pvalue$", names(dat))],
             names(dat)[grep("^Legend$", names(dat))],
             names(dat)[grep("^Size$", names(dat))],
             names(dat)[grep("^log10padj$", names(dat))],
             names(dat)[grep("^x$", names(dat))],
             names(dat)[grep("^y$", names(dat))],
             names(dat)[grep("^M$", names(dat))],
             names(dat)[grep("^A$", names(dat))],
             names(dat)[grep("^ID$", names(dat))]
  )
  dat <- dat[, !(names(dat) %in% rcols)]
}




#' getSamples
#'
#' Gathers the sample names to be used within DEBrowser.
#'
#' @param cnames, names of the  samples
#' @param index, starting column in a tab separated file
#' @return choices
#' @export
#'
#' @examples
#'     x <- getSamples()
#'
getSamples <- function (cnames = NULL, index = 1) {
  m <- NULL
  if (!is.null(cnames)) {
    cn <- cnames[index:length(cnames)]
    m <- as.list(NULL)
    for (i in seq(cn)) {
      m[i] <- cn[i]
    }
  }
  m
}

#' prepDEOutput
#'
#' Prepares the output data from DE analysis to be used within
#' DEBrowser
#'
#' @param data, loaded dataset
#' @param cols, columns
#' @param conds, conds
#' @param inputconds, inputconds
#' @param i, selected comparison number
#' @param input, input
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDEOutput()
#'
prepDEOutput <- function(data = NULL, cols = NULL,
                         conds = NULL, inputconds=NULL, i=NULL, input = NULL) {
  if (is.null(data)) return (NULL)
  if (length(cols) != length(conds)) return(NULL)
  params <- inputconds$demethod_params[i]
  de_res <- runDE(data, cols, conds, params)
  de_res <- data.frame(de_res)
}


#' applyFilters
#'
#' Applies filters based on user selected parameters to be
#' displayed within the DEBrowser.
#'
#' @param filt_data, loaded dataset
#' @param cols, selected samples
#' @param conds, seleced conditions
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- applyFilters()
#'
applyFilters <- function(filt_data = NULL, cols = NULL, conds=NULL,
                         input = NULL){
  if (is.null(input$padj) || is.null(input$foldChange)
      || is.null(filt_data)) return(NULL)
  compselect <- 1
  if (!is.null(input$compselect) )
    compselect <- as.integer(input$compselect)
  x<- isolate(input[[paste0("conditionname",2*compselect - 1)]])
  y<- isolate(input[[paste0("conditionname",2*compselect)]])

  norm_data <- getNormalizedMatrix(filt_data[, cols],
                                   input$norm_method)
  g <- data.frame(cbind(cols, conds))
  if (length(as.vector(g[g$conds == x, "cols"])) > 1 )
    filt_data$x <- log10(rowMeans(norm_data[,
                                            as.vector(g[g$conds == x, "cols"])]) + 0.1)
  else
    filt_data$x <- log10(norm_data[,
                                   as.vector(g[g$conds == x, "cols"])] + 0.1)
  if (length(as.vector(g[g$conds == y, "cols"])) > 1 )
    filt_data$y <- log10(rowMeans(norm_data[,
                                            as.vector(g[g$conds == y, "cols"])]) + 0.1)
  else
    filt_data$y <- log10(norm_data[,
                                   as.vector(g[g$conds == y, "cols"])] + 0.1)
  filt_data[,cols] <- norm_data

  padj_cutoff <- as.numeric(input$padj)
  foldChange_cutoff <- as.numeric(input$foldChange)
  m <- filt_data
  # Add column which says whether a gene significant or not
  m$Legend <- character(nrow(m))
  m$Size <- character(nrow(m))
  m[, "Size"] <- "40"
  m$Legend <- "NS"
  if (input$dataset == "up" || input$dataset == "up+down" || input$dataset == "selected")
    m$Legend[m$foldChange >= foldChange_cutoff &
               m$padj <= padj_cutoff] <- "Up"
  if (input$dataset == "down" || input$dataset == "up+down" || input$dataset == "selected")
    m$Legend[m$foldChange <= (1 / foldChange_cutoff) &
               m$padj <= padj_cutoff] <- "Down"
  if (input$dataset == "most-varied" && !is.null(cols)) {
    most_varied <- getMostVariedList(m, cols, input)
    m[rownames(most_varied), c("Legend")] <- "MV"
  }
  if (!is.null(input$genesetarea) && input$genesetarea != ""
      && input$methodtabs == "panel1") {
    genelist <- getGeneSetData(m, c(input$genesetarea))
    m[rownames(genelist), "Legend"] <- "GS"
    m[rownames(genelist), "Size"] <- "100"
    tmp <- m["Legend"=="GS", ]
    tmp1 <- m["Legend"!="GS", ]
    m <- rbind(tmp1, tmp)
  }
  m
}

#' getSelectedDatasetInput
#'
#' Gathers the user selected dataset output to be displayed.
#'
#' @param rdata, filtered dataset
#' @param getSelected, selected data
#' @param getMostVaried, most varied data
#' @param mergedComparison, merged comparison data
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- getSelectedDatasetInput()
#'
getSelectedDatasetInput<-function(rdata = NULL, getSelected = NULL,
                                  getMostVaried = NULL, mergedComparison = NULL,
                                  input = NULL) {
  if (is.null(rdata)) return (NULL)
  m <- rdata
  if (input$dataset == "up") {
    m <- getUp(rdata)
  } else if (input$dataset == "down") {
    m <- getDown(rdata)
  } else if (input$dataset == "up+down") {
    m <- getUpDown(rdata)
  } else if (input$dataset == "alldetected") {
    m <- rdata
  } else if (input$dataset == "selected" && !is.null(input$selectedplot)) {
    m <- getSelected
  } else if (input$dataset == "most-varied") {
    m <- rdata[rownames(getMostVaried), ]
  } else if (input$dataset == "comparisons") {
    m <- mergedComparison
  } else if (input$dataset == "searched") {
    m <- getSearchData(rdata, input)
  }
  m
}


#' getMostVariedList
#'
#' Calculates the most varied genes to be used for specific plots
#' within the DEBrowser.
#'
#' @param datavar, loaded dataset
#' @param cols, selected columns
#' @param input, input
#' @return data
#' @export
#'
#' @examples
#'     x <- getMostVariedList()
#'
getMostVariedList <- function(datavar = NULL, cols = NULL, input = NULL){
  if (is.null(datavar)) return (NULL)
  topn <- as.integer(as.numeric(input$topn))
  filtvar <- datavar[rowSums(datavar[,cols]) >
                       as.integer(as.numeric(input$mincount)),]
  cv<-cbind(apply(filtvar, 1, function(x)
    (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))), 1)
  colnames(cv)<-c("coeff", "a")
  cvsort<-cv[order(cv[,1],decreasing=TRUE),]
  topindex<-nrow(cvsort)
  if (topindex > topn) topindex <- topn
  cvsort_top <- head(cvsort, topindex)
  selected_var <- data.frame(datavar[rownames(cvsort_top),])
}


#' getSearchData
#'
#' search the geneset in the tables and return it
#'
#' @param dat, table data
#' @param input, input params
#' @return data
#' @export
#'
#' @examples
#'     x <- getSearchData()
#'
getSearchData <- function(dat = NULL, input = NULL)
{
  if (is.null(dat)) return(NULL)
  if (input$genesetarea != ""){
    dat <- getGeneSetData(dat, c(input$genesetarea))
  }
  dat
}

#' getGeneSetData
#'
#' Gathers the specified gene set list to be used within the
#' DEBrowser.
#'
#' @param data, loaded dataset
#' @param geneset, given gene set
#' @return data
#' @export
#'
#' @examples
#'     x <- getGeneSetData()
#'
getGeneSetData <- function(data = NULL, geneset = NULL) {
  if (is.null(data)) return (NULL)

  geneset1 <- unique(unlist(strsplit(geneset, split="[:;, \t\n\t]")))
  geneset2 <- geneset1[geneset1 != ""]
  if(length(geneset2) > 3)
    geneset2 <- paste0("^", geneset2, "$")

  dat1 <- as.data.frame(data)
  if(!("ID" %in% names(dat1)))
    dat2 <- addID(dat1)
  else
    dat2 <- dat1

  dat2$ID<-factor(as.character(dat2$ID))

  geneset4 <- unique(as.vector(unlist(lapply(toupper(geneset2),
                                             function(x){ sapply(dat2[(grepl(x, toupper(dat2[,"ID"]))), "ID"],
                                                                 as.character) }))))
  retset <- data.frame(dat2[geneset4, ])
  retset
}

#' getUp
#' get up regulated data
#'
#' @param filt_data, filt_data
#' @return data
#' @export
#'
#' @examples
#'     x <- getUp()
#'
getUp <- function(filt_data = NULL){
  if(is.null(filt_data)) return(NULL)
  filt_data[
    filt_data[, "Legend"] == "Up" |
      filt_data[, "Legend"] == "GS", ]
}
#' getDown
#' get down regulated data
#'
#' @param filt_data, filt_data
#' @return data
#' @export
#'
#' @examples
#'     x <- getDown()
#'
getDown <- function(filt_data = NULL){
  if(is.null(filt_data)) return(NULL)
  filt_data[
    filt_data[, "Legend"] == "Down"|
      filt_data[, "Legend"] == "GS", ]
}

#' getUpDown
#' get up+down regulated data
#'
#' @param filt_data, filt_data
#' @return data
#' @export
#'
#' @examples
#'     x <- getUpDown()
#'
getUpDown <- function(filt_data = NULL){
  if(is.null(filt_data)) return(NULL)
  filt_data[
    filt_data[, "Legend"] == "Up" |
      filt_data[, "Legend"] == "Down"|
      filt_data[, "Legend"] == "GS", ]
}

#' getDataForTables
#' get data to fill up tables tab
#'

#' @param input, input parameters
#' @param init_data, initial dataset
#' @param filt_data, filt_data
#' @param selected, selected genes
#' @param getMostVaried, most varied genes
#' @param mergedComp, merged comparison set
#' @param explainedData, pca gene set
#' @return data
#' @export
#'
#' @examples
#'     x <- getDataForTables()
#'
getDataForTables <- function(input = NULL, init_data = NULL,
                             filt_data = NULL, selected = NULL,
                             getMostVaried = NULL,  mergedComp = NULL,
                             explainedData = NULL){
  if (is.null(init_data )) return(NULL)
  if (is.null(filt_data)) filt_data <- init_data
  pastr <- "padj"
  fcstr <- "foldChange"
  dat <- NULL
  if (input$dataset == "alldetected"){
    dat <- getSearchData(filt_data, input)
  }
  else if (input$dataset == "up+down"){
    if (!is.null(filt_data))
      dat <- getSearchData(getUpDown(filt_data), input)
  }
  else if (input$dataset == "up"){
    if (!is.null(filt_data))
      dat <- getSearchData(getUp(filt_data), input)
  }
  else if (input$dataset == "down"){
    if (!is.null(filt_data))
      dat <- getSearchData(getDown(filt_data), input)
  }
  else if (input$dataset == "selected"){
    dat <- getSearchData(selected, input)
  }
  else if (input$dataset == "most-varied"){
    if (!is.null(filt_data)){
      d <- filt_data[rownames(getMostVaried),]
    }else{
      d <- init_data[rownames(getMostVaried),]
    }
    dat <- getSearchData(d, input)
  }
  else if (input$dataset == "comparisons"){
    if (is.null(mergedComp)) return(NULL)
    fcstr<-colnames(mergedComp)[grepl("foldChange", colnames(mergedComp))]
    pastr<-colnames(mergedComp)[grepl("padj", colnames(mergedComp))]
    dat <- getSearchData(mergedComp, input)
  }
  else if (input$dataset == "searched"){
    dat <- getSearchData(init_data, input)
  }
  list(dat, pastr, fcstr)
}


#' getMergedComparison
#'
#' Gathers the merged comparison data to be used within the
#' DEBrowser.
#' @param dc, data container
#' @param nc, the number of comparisons
#' @param input, input params
#' @return data
#' @export
#'
#' @examples
#'     x <- getMergedComparison()
#'
getMergedComparison <- function(dc = NULL, nc = NULL, input = NULL){
  if (is.null(dc)) return (NULL)
  mergeresults <- c()
  mergedata <- c()
  allsamples <- c()
  for ( ni in seq(1:nc)) {
    tmp <- dc[[ni]]$init_data[,c("foldChange", "padj")]

    samples <- dc[[ni]]$cols
    tt <- paste0(isolate(input[[paste0("conditionname", 2*ni-1)]]),".vs.",
                 isolate(input[[paste0("conditionname",2*ni)]]))
    fctt <- paste0("foldChange.", tt)
    patt <-  paste0("padj.", tt)
    colnames(tmp) <- c(fctt,  patt)
    if(ni == 1){
      allsamples <- samples
      mergeresults <- tmp
      mergedata <- dc[[ni]]$init_data[,samples]
    }
    else{
      mergeresults[,fctt] <- character(nrow(tmp))
      mergeresults[,patt] <- character(nrow(tmp))
      mergeresults[rownames(tmp),c(fctt, patt)] <- tmp[,c(fctt, patt)]
      mergeresults[rownames(tmp),patt] <- tmp[,patt]
      mergeresults[is.na(mergeresults[,fctt]),fctt] <- 1
      mergeresults[is.na(mergeresults[,patt]),patt] <- 1
      remaining_samples <- dc[[ni]]$cols[!(samples %in% colnames(mergedata))]
      allsamples <- unique(c(allsamples, remaining_samples))
      mergedata <- cbind(mergedata,  dc[[ni]]$init_data[,remaining_samples])
      colnames(mergedata) <- allsamples
    }
  }
  mergedata[,allsamples] <- getNormalizedMatrix(mergedata[,allsamples], input$norm_method)
  cbind(mergedata, mergeresults)
}

#' applyFiltersToMergedComparison
#'
#' Gathers the merged comparison data to be used within the
#' DEBrowser.
#'
#' @param merged, merged data
#' @param nc, the number of comparisons
#' @param input, input params
#' @return data
#' @export
#'
#' @examples
#'     x <- applyFiltersToMergedComparison()
#'
applyFiltersToMergedComparison <- function (merged = NULL,
                                            nc = NULL, input = NULL)
{
  if (is.null(merged)) return (NULL)
  padj_cutoff <- as.numeric(input$padj)
  foldChange_cutoff <- as.numeric(input$foldChange)
  if (is.null(merged$Legend)){
    merged$Legend <- character(nrow(merged))
    merged$Legend <- "NS"
  }
  for ( ni in seq(1:nc)) {
    tt <- paste0(isolate(input[[paste0("conditionname", 2*ni-1)]]),".vs.",
                 isolate(input[[paste0("conditionname",2*ni)]]))
    merged[which(as.numeric(merged[,c(paste0("foldChange.", tt))]) >=
                   foldChange_cutoff & as.numeric(merged[,c(paste0("padj.", tt))]) <=
                   padj_cutoff), "Legend"] <- "Sig"
    merged[which(as.numeric(merged[,c(paste0("foldChange.", tt))]) <=
                   1/foldChange_cutoff & as.numeric(merged[,c(paste0("padj.", tt))]) <=
                   padj_cutoff), "Legend"] <- "Sig"
  }
  print(head(merged))
  merged
}

#' removeCols
#'
#' remove unnecessary columns
#'
#' @param cols, columns that are going to be removed from data frame
#' @param dat, data
#' @return data
#' @export
#'
#' @examples
#'     x <- removeCols()
#'
removeCols <- function( cols = NULL, dat = NULL) {
  if (is.null(dat)) return (NULL)
  for (colnum in seq(1:length(cols))){
    if (cols[colnum] %in% colnames(dat) )
      dat[, cols[colnum]]<- NULL
  }
  dat
}







#' plotSizeMarginsUI
#'
#' Size and margins module for plotly plots
#'
#' @note \code{plotSizeMarginsUI}
#' @param id, id
#' @param h, height
#' @param w, width
#' @param t, top margin
#' @param b, bottom margin
#' @param l, left margin
#' @param r, right margin
#' @return size and margins controls
#' @examples
#'     x <- plotSizeMarginsUI("heatmap")
#' @export
#'
plotSizeMarginsUI <- function(id, w=800, h=640, t=20, b=100, l=100, r=20) {
  shinydashboard::menuItem(paste0(id, " - Size & Margins"),
                           plotSizeUI(id, w, h),
                           plotMarginsUI(id, t, b, l, r)
  )
}

#' plotSizeUI
#'
#' Size module for plotly plots
#'
#' @note \code{plotSizeUI}
#' @param id, id
#' @param h, height
#' @param w, width
#' @return size and margins controls
#' @examples
#'     x <- plotSizeUI("heatmap")
#' @export
#'
plotSizeUI <- function(id, w=800, h=600){
  ns <- NS(id)
  list(
    checkboxInput(r(ns('size'), "-"), paste0('Plot Size'), value = FALSE),
    conditionalPanel(paste0('input.', r(ns('size'), "-")),
                     sliderInput(ns("width"), "width",
                                 min = 100, max = 2000, step = 10, value = w),
                     sliderInput(ns("height"), "height",
                                 min = 100, max = 2000, step = 10, value = h)
    )
  )
}

r <- function(str, chr)
{
  gsub(chr, '', str)
}

#' plotMarginsUI
#'
#' Margins module for plotly plots
#'
#' @note \code{plotMarginsUI}
#' @param id, id
#' @param t, top margin
#' @param b, bottom margin
#' @param l, left margin
#' @param r, right margin
#' @return size and margins controls
#' @examples
#'     x <- plotMarginsUI("heatmap")
#' @export
#'
plotMarginsUI <- function(id, t=20, b=100, l=100, r=20){
  ns <- NS(id)
  list(
    checkboxInput(r(ns('margins'), "-"), 'Margins', value = FALSE),
    conditionalPanel(paste0('input.', r(ns('margins'), "-")),
                     sliderInput(ns("top"), "Margin Top", min = 0, max = 200, value = t),
                     sliderInput(ns("bottom"), "Margin Bottom", min = 0, max = 200, value = b),
                     sliderInput(ns("left"), "Margin Left", min = 0, max = 200, value = l),
                     sliderInput(ns("right"), "Margin Right", min = 0, max = 200, value = r))
  )
}




#' getBoxMainPlotUI
#'
#' main Box plots UI.
#'
#' @note \code{getBoxMainPlotUI}
#' @param id, namespace id
#' @return the panel for Density plots;
#'
#' @examples
#'     x <- getBoxMainPlotUI("box")
#'
#' @export
#'
getBoxMainPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("BoxMainUI"))
}

#' debrowserboxmainplot
#'
#' Module for a box plot that can be used in DEanalysis main part and
#' used heatmaps
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @param cols, columns
#' @param conds, conditions
#' @param key, the gene or region name
#' @return density plot
#' @export
#'
#' @examples
#'     x <- debrowserboxmainplot()
#'
debrowserboxmainplot <- function(input = NULL, output = NULL, session = NULL, data = NULL,
                                 cols = NULL, conds = NULL, key=NULL) {
  if(is.null(data)) return(NULL)
  output$BoxMain <- renderPlotly({
    getBoxMainPlot(data, cols, conds, key, title="", input)
  })

  output$BoxMainUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("plot"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotlyOutput(session$ns("BoxMain"),
                                      height=input$height, width=input$width))
  })
}

#' BoxMainPlotControlsUI
#'
#' Generates the controls in the left menu for a Box main plot
#'
#' @note \code{BoxMainPlotControlsUI}
#' @param id, namespace id
#' @return returns the controls for left menu
#' @examples
#'     x <- BoxMainPlotControlsUI("box")
#' @export
#'
BoxMainPlotControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(id, " - Options"),
                           textInput(ns("breaks"), "Breaks", value = "100" )
  )
}

#' getBoxMainPlot
#'
#' Makes Density plots
#'
#' @param data, count or normalized data
#' @param cols, cols
#' @param conds, conds
#' @param key, key
#' @param title, title
#' @param input, input
#' @export
#'
#' @examples
#'     getBoxMainPlot()
#'
getBoxMainPlot <- function(data=NULL, cols = NULL, conds=NULL, key=NULL, title = "", input = NULL){
  if (is.null(data)) return(NULL)
  vardata <- getVariationData(data, cols, conds, key)
  title <- paste(key, "variation")
  p <- plot_ly(vardata, x = ~conds, y = ~count,
               color=~conds, colors=c("Blue", "Red"),
               boxpoints = "all", type = "box") %>%
    plotly::layout(title = title,
                   xaxis = list(title = "Conditions"),
                   yaxis = list(title = "Read Count"),
                   height=input$height, width=input$width,
                   margin = list(l = input$left,
                                 b = input$bottom,
                                 t = input$top,
                                 r = input$right
                   ))
  p$elementId <- NULL
  p
}








#' getDensityPlotUI
#'
#' Density plot UI.
#'
#' @param id, namespace id
#' @note \code{getDensityPlotUI}
#' @return the panel for Density plots;
#'
#' @examples
#'     x <- getDensityPlotUI("density")
#'
#' @export
#'
getDensityPlotUI <- function(id) {
  ns <- NS(id)
  uiOutput(ns("DensityUI"))
}

#' debrowserdensityplot
#'
#' Module for a density plot that can be used in data prep and
#' low count removal modules
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param data, a matrix that includes expression values
#' @return density plot
#' @export
#'
#' @examples
#'     x <- debrowserdensityplot()
#'
debrowserdensityplot <- function(input = NULL, output = NULL, session = NULL, data = NULL) {
  if(is.null(data)) return(NULL)
  output$Density <- renderPlotly({
    getDensityPlot(data, input)
  })
  output$DensityUI <- renderUI({
    shinydashboard::box(
      collapsible = TRUE, title = session$ns("plot"), status = "primary",
      solidHeader = TRUE, width = NULL,
      draggable = TRUE,  plotlyOutput(session$ns("Density"),
                                      width = input$width, height=input$height))
  })
}

#' densityPlotControlsUI
#'
#' Generates the controls in the left menu for a densityPlot
#'
#' @note \code{densityPlotControlsUI}
#' @param id, namespace id
#' @return returns the left menu
#' @examples
#'     x <- densityPlotControlsUI("density")
#' @export
#'
densityPlotControlsUI <- function(id) {
  ns <- NS(id)
  shinydashboard::menuItem(paste0(id, " - Options"),
                           textInput(ns("breaks"), "Breaks", value = "100" )
  )
}

#' getDensityPlot
#'
#' Makes Density plots
#'
#' @param data, count or normalized data
#' @param input, input
#' @param title, title
#'
#' @export
#'
#' @examples
#'     getDensityPlot()
#'
getDensityPlot <- function(data=NULL, input = NULL, title = ""){
  if (is.null(data)) return(NULL)
  data <- as.data.frame(data)
  cols <- colnames(data)
  data[, cols] <- apply(data[, cols], 2,
                        function(x) log10(as.integer(x) + 1))

  data <- addID(data)
  mdata <- melt(as.data.frame(data[,c("ID", cols)]),"ID")
  colnames(mdata)<-c("ID", "samples", "density")

  p <- ggplot(data=mdata, aes(x=density)) +
    geom_density(aes(fill = samples), alpha = 0.5) +
    labs(x = "logcount", y = "Density") +
    theme_minimal()
  if (!is.null(input$top))
    p <- p + theme( plot.margin = margin(t = input$top, r =input$right, b =input$bottom, l = input$left, "pt"))
  p <- ggplotly(p, width = input$width, height = input$height)
  p$elementId <- NULL
  p
}



#' debrowserdataload
#'
#' Module to load count data and metadata
#'
#' @param input, input variables
#' @param output, output objects
#' @param session, session
#' @param nextpagebutton, the name of the next page button after loading the data
#' @return main plot
#'
#' @return panel
#' @export
#'
#' @examples
#'     x <- debrowserdataload()
#'
debrowserdataload <- function(input = NULL, output = NULL, session = NULL, nextpagebutton = NULL) {
  if (is.null(input)) return(NULL)
  ldata <- reactiveValues(count=NULL, meta=NULL)
  loadeddata <- reactive({
    ret <- NULL
    if(!is.null(ldata$count)){
      ldata$count <- ldata$count[,sapply(ldata$count, is.numeric)]
      ret <- list(count = ldata$count, meta = ldata$meta)
    }
    return(ret)
  })
  output$dataloaded <- reactive({
    return(!is.null(loadeddata()))
  })
  outputOptions(output, "dataloaded",
                suspendWhenHidden = FALSE)
  observe({
    query <- parseQueryString(session$clientData$url_search)
    jsonobj<-query$jsonobject

    # To test json load;
    # It accepts three parameters:
    # 1. jsonobject=https%3A%2F%2Fdolphin.umassmed.edu%2Fpublic%2Fapi%2F%3Fsource%3Dhttps%3A%2F%2Fbioinfo.umassmed.edu%2Fpub%2Fdebrowser%2Fadvanced_demo.tsv%26format%3DJSON
    # 2. meta=meta=https%3A%2F%2Fdolphin.umassmed.edu%2Fpublic%2Fapi%2F%3Fsource%3Dhttps%3A%2F%2Fbioinfo.umassmed.edu%2Fpub%2Fdebrowser%2Fsimple_meta.tsv%26format%3DJSON
    # 3. title=no
    # The finished product of the link will look like this without metadata:
    #
    # https://127.0.0.1:3427/debrowser/R/?jsonobject=https%3A%2F%2Fdolphin.umassmed.edu%2Fpublic%2Fapi%2F%3Fsource%3Dhttps%3A%2F%2Fbioinfo.umassmed.edu%2Fpub%2Fdebrowser%2Fadvanced_demo.tsv%26format%3DJSON&title=no
    #
    #  With metadata
    #
    #http://127.0.0.1:3427/?jsonobject=https%3A%2F%2Fdolphin.umassmed.edu%2Fpublic%2Fapi%2F%3Fsource%3Dhttps%3A%2F%2Fbioinfo.umassmed.edu%2Fpub%2Fdebrowser%2Fsimple_demo.tsv%26format%3DJSON&meta=https%3A%2F%2Fdolphin.umassmed.edu%2Fpublic%2Fapi%2F%3Fsource%3Dhttps%3A%2F%2Fbioinfo.umassmed.edu%2Fpub%2Fdebrowser%2Fsimple_meta.tsv%26format%3DJSON
    #
    #

    if (!is.null(jsonobj))
    {
      raw <- RCurl::getURL(jsonobj, .opts = list(ssl.verifypeer = FALSE),
                           crlf = TRUE)
      jsondata<-data.frame(fromJSON(raw, simplifyDataFrame = TRUE),
                           stringsAsFactors = TRUE)
      rownames(jsondata)<-jsondata[, 1]
      jsondata<-jsondata[,c(3:ncol(jsondata))]
      jsondata[,c(1:ncol(jsondata))] <- sapply(
        jsondata[,c(1:ncol(jsondata))], as.numeric)
      jsondata <- jsondata[,sapply(jsondata, is.numeric)]
      ldata$count <- jsondata

      metadatatable <- NULL
      jsonmet <-query$meta
      if(!is.null(jsonmet)){
        raw <- RCurl::getURL(jsonmet, .opts = list(ssl.verifypeer = FALSE),
                             crlf = TRUE)
        metadatatable<-data.frame(fromJSON(raw, simplifyDataFrame = TRUE),
                                  stringsAsFactors = TRUE)

      }else{
        metadatatable <- cbind(colnames(ldata$count), 1)
        colnames(metadatatable) <- c("Sample", "Batch")
      }
      ldata$meta <- metadatatable
      input$Next
    }
  })
  observeEvent(input$demo, {
    load(system.file("extdata", "demo", "demodata.Rda",
                     package = "debrowser"))

    ldata$count <- demodata
    ldata$meta <- metadatatable
  })
  observeEvent(input$demo2, {
    load(system.file("extdata", "demo", "demodata2.Rda",
                     package = "debrowser"))
    ldata$count <- demodata
    ldata$meta <- metadatatable
  })

  observeEvent(input$uploadFile, {
    if (is.null(input$countdata)) return (NULL)
    checkRes <- checkCountData(input)

    if (checkRes != "success"){
      showNotification(checkRes, type = "error")
      return(NULL)
    }
    counttable <-as.data.frame(
      try(
        read.delim(input$countdata$datapath,
                   header=T, sep=input$countdataSep,
                   row.names=1, strip.white=TRUE ), TRUE))
    counttable <- counttable[,sapply(counttable, is.numeric)]
    metadatatable <- c()
    if (!is.null(input$metadata$datapath)){
      metadatatable <- as.data.frame(
        try(
          read.delim(input$metadata$datapath,
                     header=TRUE, sep=input$metadataSep, strip.white=TRUE), TRUE))

      checkRes <- checkMetaData(input, counttable)
      if (checkRes != "success"){
        showNotification(checkRes, type = "error")
        return(NULL)
      }
    }
    else{
      metadatatable <- cbind(colnames(counttable), 1)
      colnames(metadatatable) <- c("Sample", "Batch")
    }
    if (is.null(counttable))
    {stop("Please upload the count file")}
    ldata$count <- counttable
    ldata$meta <- metadatatable
  })
  output$nextButton <- renderUI({
    actionButtonDE(nextpagebutton, label = nextpagebutton, styleclass = "primary")
  })
  observe({
    getSampleDetails(output, "uploadSummary", "sampleDetails", loadeddata())
  })
  list(load=loadeddata)
}




#' dataLoadUI
#'
#' Creates a panel to upload the data
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- dataLoadUI("load")
#'
#' @export
#'
dataLoadUI<- function (id) {
  ns <- NS(id)
  list(conditionalPanel(condition =  paste0("!output['", ns("dataloaded"),"']"),
                        fluidRow(
                          fileUploadBox(id, "countdata", "Count Data"),
                          fileUploadBox(id, "metadata", "Metadata")
                        ),
                        fluidRow(column(12,
                                        actionButtonDE(ns("uploadFile"), label = "Upload", styleclass = "primary"),
                                        actionButtonDE(ns("demo"),  label = "Load Demo (Vernia et. al)", styleclass = "primary"),
                                        actionButtonDE(ns("demo2"),  label = "Load Demo (Donnard et. al)", styleclass = "primary")))),
       fluidRow(column(12,
                       conditionalPanel(condition = paste0("output['", ns("dataloaded"),"']"),
                                        uiOutput(ns("nextButton"))
                       ))
       ), br(),
       fluidRow(
         shinydashboard::box(title = "Upload Summary",
                             solidHeader = T, status = "info",
                             width = 12,
                             fluidRow(
                               column(12,
                                      tableOutput(ns("uploadSummary"))
                               )),
                             fluidRow(
                               column(12,div(style = 'overflow: scroll',
                                             DT::dataTableOutput(ns("sampleDetails")))
                               )
                             )
         )
       ))
}


#' fileUploadBox
#'
#' File upload module
#' @param id, namespace id
#' @param inputId, input file ID
#' @param label, label
#' @note \code{fileUploadBox}
#' @return radio control
#'
#' @examples
#'
#'     x <- fileUploadBox("meta", "metadata", "Metadata")
#'
#' @export
#'
fileUploadBox <- function(id = NULL, inputId = NULL, label = NULL) {
  ns <- NS(id)
  shinydashboard::box(title = paste0(label, " File"),
                      solidHeader = TRUE, status = "info",
                      width = 6,
                      helpText(paste0("Upload your '", label," File'")),
                      fileInput(inputId=ns(inputId),
                                label=NULL,
                                accept=fileTypes()
                      ),
                      sepRadio(id, paste0(inputId, "Sep")))
}

#' sepRadio
#'
#' Radio button for separators
#'
#' @param id, module id
#' @param name, name
#' @note \code{sepRadio}
#' @return radio control
#'
#' @examples
#'
#'     x <- sepRadio("meta", "metadata")
#'
#' @export
#'
sepRadio <- function(id, name) {
  ns <- NS(id)
  radioButtons(inputId=ns(name),
               label="Separator",
               choices=c(Comma=',',
                         Semicolon=';',
                         Tab='\t'
               ),
               selected='\t'
  )
}

#' fileTypes
#'
#' Returns fileTypes that are going to be used in creating fileUpload UI
#'
#' @note \code{fileTypes}
#' @return file types
#'
#' @examples
#'     x <- fileTypes()
#'
#' @export
#'
fileTypes <- function() {
  c('text/tab-separated-values',
    'text/csv',
    'text/comma-separated-values',
    'text/tab-separated-values',
    '.txt',
    '.csv',
    '.tsv')
}

#' checkCountData
#'
#' Returns if there is a problem in the count data.
#'
#' @note \code{checkCountData}
#' @param input, inputs
#' @return error if there is a problem about the loaded data
#
#' @examples
#'     x <- checkCountData()
#'
#' @export
#'
checkCountData <- function(input = NULL){
  if (is.null(input$countdata$datapath)) return(NULL)
  tryCatch({
    data <- read.table(input$countdata$datapath, sep=input$countdataSep)
    if (ncol(data) < 3) return ("Error: Please check if you chose the right separator!")
    dups <- data[duplicated(data[,1], fromLast = TRUE),1]
    if (length(dups)>1) return (paste0("Error: There are duplicate entried in  the rownames. (",
                                       paste0(dups, collapse=","),")"))

    return("success")
  }, error = function(err) {
    return (paste0("Error(Count file):",toString(err)))
  }, warning = function(war) {
    return(paste0("Warning(Count file):",toString(err)))
  })
}


#' checkMetaData
#'
#' Returns if there is a problem in the count data.
#'
#' @note \code{checkMetaData}
#' @param input, input
#' @param counttable, counttable
#' @return error if there is a problem about the loaded data
#
#' @examples
#'     x <- checkMetaData()
#'
#' @export
#'
checkMetaData <- function(input = NULL, counttable = NULL){
  if (is.null(counttable) || is.null(input$metadata$datapath)) return(NULL)
  tryCatch({
    metadatatable <- read.table(input$metadata$datapath, sep=input$metadataSep, header=T)
    if (ncol(metadatatable) < 2) return ("Error: Please check if you chose the right separator!")
    met <- as.vector(metadatatable[order(as.vector(metadatatable[,1])), 1])
    count <- as.vector(colnames(counttable)[order(as.vector(colnames(counttable)))])
    difference <- base::setdiff(met, count)
    if (length(difference)>0){
      return(paste0("Colnames doesn't match with the metada table(", paste0(difference,sep=","), ")"))
    }
    return("success")
  }, error = function(err) {
    return (paste0("Error(Matadata file):",toString(err)))
  }, warning = function(war) {
    return (paste0("Warning(Matadata file):",toString(war)))
  })
}




#' deUI
#'
#' Creates a shinyUI to be able to run DEBrowser interactively.
#'
#' @note \code{deUI}
#' @return the panel for main plots;
#'
#' @examples
#'     x<-deUI()
#'
#' @export
#'

deUI <- function() {
  dbHeader <- shinydashboard::dashboardHeader(titleWidth = 250)
  dbHeader$children[[2]]$children <- tags$a(style='color: white;',
                                            id="top_logo" , paste0("DEBrowser v",getNamespaceVersion("debrowser")))
  addResourcePath(prefix = "www", directoryPath = system.file("extdata",
                                                              "www", package = "debrowser"))
  debrowser <- (fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::inlineCSS("
        #loading-debrowser {
        position: absolute;
        background: #000000;
        opacity: 0.9;
        z-index: 100;
        left: 0;
        right: 0;
        height: 100%;
        text-align: center;
        color: #EFEFEF;
    }"),
    # Loading message
    tags$div(h4(paste0("Loading DEBrowser v",getNamespaceVersion("debrowser"))), id = "loading-debrowser",
             tags$img(src = "www/images/initial_loading.gif")),
    tags$head(tags$title(paste0("DEBrowser v",getNamespaceVersion("debrowser"))),
              tags$link(rel = "stylesheet", type = "text/css",
                        href = "www/shinydashboard_additional.css")
    ),
    dashboardPage(
      dbHeader,
      dashboardSidebar(
        width = 250,
        debrowser::getJSLine(),
        uiOutput("loading"),
        tabsetPanel(id = "menutabs", type = "tabs",
                    tabPanel(title = "Data Prep", value = "dataprep", id="dataprep",
                             sidebarMenu(id="DataPrep",
                                         menuItem("Quick Start Guide", icon = icon("user"),
                                                  menuSubItem("Introduction", tabName = "Intro"),
                                                  menuSubItem("Data Assesment", tabName = "assesment"),
                                                  menuSubItem("Data Preparation", tabName = "preparation"),
                                                  menuSubItem("DE Anaylsis", tabName = "deanalysis"),
                                                  menuSubItem("FAQ", tabName ="FAQ")
                                         ),
                                         menuItem("Upload", icon = icon("upload"), tabName = "Upload"),
                                         menuItem("Filter", icon = icon("filter"), tabName = "Filter"),
                                         menuItem("BatchEffect",  icon = icon("align-left"), tabName = "BatchEffect"),
                                         menuItem("CondSelect",  icon = icon("bars"), tabName = "CondSelect"),
                                         menuItem("DEAnalysis", icon = icon("adjust"), tabName = "DEAnalysis"),
                                         menuItem("DEFilter",  icon = icon("code"), tabName = "DEAnalysis",  startExpanded = TRUE,
                                                  uiOutput("cutOffUI"),
                                                  uiOutput("compselectUI"))
                             ),helpText("Developed by ", a("UMMS Biocore.",
                                                           href="https://www.umassmed.edu/biocore/", target = "_blank"))),
                    tabPanel(title = "Discover", value = "discover", id="discover",
                             conditionalPanel(condition = "(output.dataready)",
                                              conditionalPanel( (condition <- "input.methodtabs=='panel1'"),
                                                                debrowser::mainPlotControlsUI("main")),

                                              uiOutput("downloadSection"),
                                              uiOutput('cutoffSelection'),
                                              uiOutput("leftMenu")
                             )
                    ))
      ),
      dashboardBody(
        mainPanel(
          width = 12,
          tags$head(
            tags$style(type = "text/css",
                       "#methodtabs.nav-tabs {font-size: 14px} ")),
          tabsetPanel(id = "methodtabs", type = "tabs",
                      tabPanel(title = "Data Prep", value = "panel0", id="panel0",
                               tabItems(
                                 tabItem(tabName="Intro", debrowser::getIntroText()),
                                 tabItem(tabName="assesment", debrowser::getDataAssesmentText()),
                                 tabItem(tabName="preparation", debrowser::getDataPreparationText()),
                                 tabItem(tabName="deanalysis", debrowser::getDEAnalysisText()),
                                 tabItem(tabName="FAQ",  debrowser::getQAText()),
                                 tabItem(tabName="Upload", debrowser::dataLoadUI("load")),
                                 tabItem(tabName="Filter" ,
                                         conditionalPanel(
                                           (condition <- "input.Next"),
                                           actionButtonDE("goDEwithoutfilter",label="Skip to DE Analysis", styleclass = "primary"),
                                           debrowser::dataLCFUI("lcf"))),
                                 tabItem(tabName="BatchEffect",
                                         conditionalPanel(
                                           (condition <- "input.Batch"),
                                           debrowser::batchEffectUI("batcheffect"))),
                                 tabItem(tabName="CondSelect",
                                         conditionalPanel(
                                           (condition <- "input.goDE || input.goDEFromFilter || input.goDEwithoutfilter"),
                                           debrowser::condSelectUI())),
                                 tabItem(tabName="DEAnalysis",
                                         conditionalPanel(
                                           (condition <- "input.goDE || input.goDEFromFilter || input.goDEwithoutfilter"),
                                           uiOutput("deresUI")))
                               )),
                      tabPanel(title = "Main Plots", value = "panel1", id="panel1",
                               uiOutput("mainmsgs"),
                               uiOutput("mainpanel")),
                      tabPanel(title = "QC Plots", value = "panel2", id="panel2",
                               uiOutput("qcpanel")),
                      tabPanel(title = "GO Term", value = "panel3", id="panel3",
                               uiOutput("gopanel")),
                      tabPanel(title = "Tables", value = "panel4", id="panel4",
                               DT::dataTableOutput("tables"))
                      )


        ),
        debrowser::getTabUpdateJS()
      ))
  )
  )
  debrowser
}

#' PCbrowser: A package to generate the diffential genes
#'
#' This packge in R shiny, provides the function to calculate the DE genes.
#'
#' @docType package
#' @name PCbrowser
"_PACKAGE"

#' Sets up shinyServer to be able to run DEBrowser interactively.
#'
#' @note \code{deServer}
#' @param input, input params from UI
#' @param output, output params to UI
#' @param session, session variable
#' @return the panel for main plots;
#'
#' @examples
#'     deServer
#'
#' @export
#' @importFrom shiny actionButton actionLink addResourcePath column
#'             conditionalPanel downloadButton downloadHandler
#'             eventReactive fileInput fluidPage helpText isolate
#'             mainPanel need numericInput observe observeEvent
#'             outputOptions parseQueryString plotOutput radioButtons
#'             reactive reactiveValues renderPlot renderUI runApp
#'             selectInput shinyApp  shinyServer  shinyUI sidebarLayout
#'             sidebarPanel sliderInput  stopApp  tabPanel tabsetPanel
#'             textInput textOutput titlePanel uiOutput tags HTML
#'             h4 img icon updateTabsetPanel updateTextInput  validate
#'             wellPanel checkboxInput br p checkboxGroupInput onRestore
#'             reactiveValuesToList renderText onBookmark onBookmarked
#'             updateQueryString callModule enableBookmarking htmlOutput
#'             onRestored NS reactiveVal withProgress tableOutput
#'             selectizeInput fluidRow div renderPrint renderImage
#'             verbatimTextOutput imageOutput renderTable incProgress
#'             a h3 strong h2 withMathJax updateCheckboxInput
#'             showNotification updateSelectInput
#' @importFrom shinyjs show hide enable disable useShinyjs extendShinyjs
#'             js inlineCSS onclick
#' @importFrom d3heatmap d3heatmap renderD3heatmap d3heatmapOutput
#' @importFrom DT datatable dataTableOutput renderDataTable formatStyle
#'             styleInterval formatRound
#' @importFrom ggplot2 aes aes_string geom_bar geom_point ggplot
#'             labs scale_x_discrete scale_y_discrete ylab
#'             autoplot theme_minimal theme geom_density
#'             geom_text element_blank margin
#' @importFrom plotly renderPlotly plotlyOutput plot_ly add_bars event_data
#'             hide_legend %>% group_by ggplotly
#' @importFrom gplots heatmap.2 redblue bluered
#' @importFrom igraph layout.kamada.kawai
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom graphics barplot hist pairs par rect text plot
#' @importFrom stats aggregate as.dist cor cor.test dist
#'             hclust kmeans na.omit prcomp var sd model.matrix
#'             p.adjust runif cov mahalanobis quantile as.dendrogram
#'             density
#' @importFrom utils read.csv read.table write.table update.packages
#'             download.file read.delim data install.packages
#'             packageDescription installed.packages
#' @importFrom DOSE enrichDO
#' @importFrom enrichplot gseaplot dotplot
#' @importMethodsFrom DOSE summary
#' @importMethodsFrom AnnotationDbi as.data.frame as.list colnames
#'             exists sample subset head mappedkeys ncol nrow subset
#'             keys mapIds select
#' @importMethodsFrom GenomicRanges as.factor setdiff
#' @importMethodsFrom IRanges as.matrix "colnames<-" mean
#'             nchar paste rownames toupper unique which
#'             as.matrix lapply "rownames<-" gsub
#' @importMethodsFrom S4Vectors eval grep grepl levels sapply t
#' @importMethodsFrom SummarizedExperiment cbind order rbind
#' @importFrom jsonlite fromJSON
#' @importFrom methods new
#' @importFrom stringi stri_rand_strings
#' @importFrom annotate geneSymbols
#' @importFrom reshape2 melt
#' @importFrom Harman harman reconstructData
#' @importFrom clusterProfiler compareCluster enrichKEGG enrichGO
#' @importFrom DESeq2 DESeq DESeqDataSetFromMatrix results estimateSizeFactors
#'             counts
#' @importFrom edgeR calcNormFactors equalizeLibSizes DGEList glmLRT
#'             exactTest estimateCommonDisp glmFit
#' @importFrom shinydashboard dashboardHeader dropdownMenu messageItem
#'             dashboardPage dashboardSidebar sidebarMenu dashboardBody
#'             updateTabItems menuItem tabItems tabItem menuSubItem
#' @importFrom limma lmFit voom eBayes topTable
#' @importFrom sva ComBat
#' @importFrom RCurl getURL
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import shinyBS
#' @import colourpicker
#' @import RColorBrewer
#' @import heatmaply

deServer <- function(input, output, session) {
  options(warn = -1)
  tryCatch(
    {
      if (!interactive()) {
        options( shiny.maxRequestSize = 30 * 1024 ^ 2,
                 shiny.fullstacktrace = FALSE, shiny.trace=FALSE,
                 shiny.autoreload=TRUE, warn =-1)
      }
      # To hide the panels from 1 to 4 and only show Data Prep
      togglePanels(0, c(0), session)

      choicecounter <- reactiveValues(nc = 0)

      output$programtitle <- renderUI({
        togglePanels(0, c(0), session)
        getProgramTitle(session)
      })

      updata <- reactiveVal()
      filtd <- reactiveVal()
      batch <- reactiveVal()
      sel <- reactiveVal()
      dc <- reactiveVal()
      compsel <- reactive({
        cp <- 1
        if (!is.null(input$compselect_dataprep))
          cp <- input$compselect_dataprep
        cp
      })

      observe({
        updata(callModule(debrowserdataload, "load", "Next"))
        updateTabItems(session, "DataPrep", "Upload")
        observeEvent (input$Next, {
          if(!is.null(updata()$load())){
            updateTabItems(session, "DataPrep", "Filter")
            filtd(callModule(debrowserlowcountfilter, "lcf", updata()$load()))
          }
        })
        observeEvent (input$Batch, {
          if(!is.null(filtd()$filter())){
            updateTabItems(session, "DataPrep", "BatchEffect")
            batch(callModule(debrowserbatcheffect, "batcheffect", filtd()$filter()))
          }
        })

        observeEvent (input$goDEFromFilter, {
          if(is.null(batch())) batch(setBatch(filtd()))
          updateTabItems(session, "DataPrep", "CondSelect")
          sel(debrowsercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          choicecounter$nc <- sel()$cc()
        })
        observeEvent (input$goDE, {
          updateTabItems(session, "DataPrep", "CondSelect")
          sel(debrowsercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          choicecounter$nc <- sel()$cc()
        })
        observeEvent(input$goDEwithoutfilter,{
          if(is.null(batch())||is.null(filtd()$filter())) batch(setFilter(updata()))
          updateTabItems(session, "DataPrep", "CondSelect")
          sel(debrowsercondselect(input, output, session,
                                  batch()$BatchEffect()$count, batch()$BatchEffect()$meta))
          choicecounter$nc <- sel()$cc()
        })
        observeEvent (input$startDE, {
          if(!is.null(batch()$BatchEffect()$count)){
            togglePanels(0, c(0), session)
            dc(prepDataContainer(batch()$BatchEffect()$count, sel()$cc(), input))
            updateTabItems(session, "DataPrep", "DEAnalysis")
            buttonValues$startDE <- TRUE
            buttonValues$goQCplots <- FALSE
            hideObj(c("load-uploadFile","load-demo",
                      "load-demo2", "goQCplots", "goQCplotsFromFilter"))
          }
        })

        observeEvent (input$goMain, {
          updateTabItems(session, "methodtabs", "panel1")
          updateTabItems(session, "menutabs", "discover")
          togglePanels(0, c( 0, 1, 2, 3, 4, 5), session)
        })

        output$compselectUI <- renderUI({
          if (!is.null(sel()) && !is.null(sel()$cc()))
            getCompSelection("compselect_dataprep",sel()$cc())
        })

        output$cutOffUI <- renderUI({
          cutOffSelectionUI(paste0("DEResults", compsel()))
        })
        output$deresUI <- renderUI({
          column(12, getDEResultsUI(paste0("DEResults",compsel())))
        })
      })
      output$mainpanel <- renderUI({
        getMainPanel()
      })
      output$qcpanel <- renderUI({
        getQCPanel(input)
      })
      output$gopanel <- renderUI({
        getGoPanel()
      })
      output$cutoffSelection <- renderUI({
        nc <- 1
        if (!is.null(choicecounter$nc)) nc <- choicecounter$nc
        getCutOffSelection(nc)
      })
      output$downloadSection <- renderUI({
        choices <- c("most-varied", "alldetected")
        if (buttonValues$startDE)
          choices <- c("up+down", "up", "down",
                       "comparisons", "alldetected",
                       "most-varied", "selected")
        choices <- c(choices, "searched")
        getDownloadSection(choices)
      })

      output$leftMenu  <- renderUI({
        getLeftMenu(input)
      })
      output$loading <- renderUI({
        getLoadingMsg()
      })
      output$logo <- renderUI({
        getLogo()
      })
      output$startup <- renderUI({
        getStartupMsg()
      })
      output$afterload <- renderUI({
        getAfterLoadMsg()
      })
      output$mainmsgs <- renderUI({
        if (is.null(condmsg()))
          getStartPlotsMsg()
        else
          condmsg()
      })
      buttonValues <- reactiveValues(goQCplots = FALSE, goDE = FALSE,
                                     startDE = FALSE)
      output$dataready <- reactive({
        hide(id = "loading-debrowser", anim = TRUE, animType = "fade")
        return(!is.null(init_data()))
      })
      outputOptions(output, "dataready",
                    suspendWhenHidden = FALSE)

      observeEvent(input$resetsamples, {
        buttonValues$startDE <- FALSE
        showObj(c("goQCplots", "goDE"))
        hideObj(c("add_btn","rm_btn","startDE"))
        choicecounter$nc <- 0
      })

      output$condReady <- reactive({
        if (!is.null(sel()))
          choicecounter$nc <- sel()$cc()
        choicecounter$nc
      })
      outputOptions(output, 'condReady', suspendWhenHidden = FALSE)
      observeEvent(input$goQCplotsFromFilter, {
        if(is.null(batch())) batch(setBatch(filtd()))
        buttonValues$startDE <- FALSE
        buttonValues$goQCplots <- TRUE
        updateTabItems(session, "menutabs", "discover")
        togglePanels(2, c( 0, 2, 4), session)
      })
      observeEvent(input$goQCplots, {
        buttonValues$startDE <- FALSE
        buttonValues$goQCplots <- TRUE
        updateTabItems(session, "menutabs", "discover")
        togglePanels(2, c( 0, 2, 4), session)
      })
      comparison <- reactive({
        compselect <- 1
        if (!is.null(input$compselect))
          compselect <- as.integer(input$compselect)
        dc()[[compselect]]
      })
      conds <- reactive({ comparison()$conds })
      cols <- reactive({ comparison()$cols })
      init_data <- reactive({
        if (buttonValues$startDE && !is.null(comparison()$init_data))
          comparison()$init_data
        else if (!is.null(batch()))
          batch()$BatchEffect()$count
      })
      filt_data <- eventReactive(c(input$apply, input$getdownload), {
        if (!is.null(init_data()) && !is.null(comparison()) && !is.null(input$padj))
          applyFilters(init_data(), cols(), conds(), input)
      })

      selectedQCHeat <- reactiveVal()
      observe({
        if ((!is.null(input$genenames) && input$interactive == TRUE) ||
            (!is.null(input$genesetarea) && input$genesetarea != "")){
          tmpDat <- init_data()
          if (!is.null(filt_data()))
            tmpDat <- filt_data()
          genenames <- ""
          if (!is.null(input$genenames)){
            genenames <- input$genenames
          } else {
            tmpDat <- getSearchData(tmpDat, input)
            genenames <- paste(rownames(tmpDat), collapse = ",")
          }
        }
        if(!is.null(input$qcplot) && !is.null(normdat())){
          if (input$qcplot == "all2all") {
            callModule(debrowserall2all, "all2all", normdat(), input$cex)
          } else if (input$qcplot == "pca") {
            callModule(debrowserpcaplot, "qcpca", normdat(), batch()$BatchEffect()$meta)
          } else if (input$qcplot == "heatmap") {
            selectedQCHeat(callModule(debrowserheatmap, "heatmapQC", normdat()))
          } else if (input$qcplot == "IQR") {
            callModule(debrowserIQRplot, "IQR", df_select())
            callModule(debrowserIQRplot, "normIQR", normdat())
          } else if (input$qcplot == "Density"){
            callModule(debrowserdensityplot, "density", df_select())
            callModule(debrowserdensityplot, "normdensity", normdat())
          }
        }
      })
      condmsg <- reactiveVal()
      selectedMain <- reactiveVal()
      observe({
        if (!is.null(filt_data())) {
          condmsg(getCondMsg(dc(), input,
                             cols(), conds()))
          selectedMain(callModule(debrowsermainplot, "main", conds(), filt_data()))
        }
      })
      selectedHeat <- reactiveVal()
      observe({
        if (!is.null(selectedMain()) && !is.null(selectedMain()$selGenes())) {
          withProgress(message = 'Creating plot', style = "notification", value = 0.1, {
            selectedHeat(callModule(debrowserheatmap, "heatmap", filt_data()[selectedMain()$selGenes(), cols()]))
          })
        }
      })

      selgenename <- reactiveVal()
      observe({
        if (!is.null(selectedMain()) && !is.null(selectedMain()$shgClicked())
            && selectedMain()$shgClicked()!=""){
          selgenename(selectedMain()$shgClicked())
          if (!is.null(selectedHeat()) && !is.null(selectedHeat()$shgClicked()) &&
              selectedHeat()$shgClicked() != ""){
            js$resetInputParam("heatmap-hoveredgenenameclick")
          }
        }
      })
      observe({
        if (!is.null(selectedHeat()) && !is.null(selectedHeat()$shgClicked()) &&
            selectedHeat()$shgClicked() != ""){
          selgenename(selectedHeat()$shgClicked())
        }
      })

      observe({
        if (!is.null(selgenename()) && selgenename()!=""){
          withProgress(message = 'Creating Bar/Box plots', style = "notification", value = 0.1, {
            callModule(debrowserbarmainplot, "barmain", filt_data(),
                       cols(), conds(), selgenename())
            callModule(debrowserboxmainplot, "boxmain", filt_data(),
                       cols(), conds(), selgenename())
          })
        }
      })

      normdat <-  reactive({
        if (!is.null(init_data()) && !is.null(datasetInput())){
          dat <- init_data()
          norm <- c()
          if(!is.null(cols())){
            norm <- removeExtraCols(datasetInput())
          }else{
            norm <- getNormalizedMatrix(dat, input$norm_method)
          }
          getSelectedCols(norm, datasetInput(), input)
        }
      })

      df_select <- reactive({
        if (!is.null(init_data()) && !is.null(datasetInput()) )
          getSelectedCols(init_data(), datasetInput(), input)
      })

      output$columnSelForQC <- renderUI({
        existing_cols <- colnames(removeExtraCols(datasetInput()))
        wellPanel(id = "tPanel",
                  style = "overflow-y:scroll; max-height: 300px",
                  checkboxGroupInput("col_list", "Select col to include:",
                                     existing_cols,
                                     selected=existing_cols)
        )
      })

      selectedData <- reactive({
        dat <- isolate(filt_data())
        ret <- c()
        if (input$selectedplot == "Main Plot"  && !is.null(selectedMain())){
          ret <- dat[selectedMain()$selGenes(), ]
        }
        else if (input$selectedplot == "Main Heatmap" &&  !is.null(selectedHeat())){
          ret <- dat[selectedHeat()$selGenes(), ]
        }
        else if (input$selectedplot == "QC Heatmap" && !is.null(selectedQCHeat())){
          ret <- dat[selectedQCHeat()$selGenes(), ]
        }
        ret
      })

      datForTables <- reactive({
        getDataForTables(input, normdat(),
                         filt_data(), selectedData(),
                         getMostVaried(), mergedComp())
      })

      inputGOstart <- eventReactive(c(input$startGO,input$goapply),{
        if (input$startGO){
          withProgress(message = 'GO Started', detail = "interactive", value = 0, {
            dat <- datForTables()
            getGOPlots(dat[[1]][, isolate(cols())], input)
          })
        }
      })
      observeEvent(input$startGO, {
        inputGOstart()
      })
      output$GOPlots1 <- renderPlot({
        if (!is.null(inputGOstart()$p) && input$startGO){
          return(inputGOstart()$p)
        }
      })
      output$KEGGPlot <- renderImage({
        shiny::validate(need(!is.null(input$gotable_rows_selected),
                             "Please select a category in the GO/KEGG table tab to be able
                to see the pathway diagram")
        )

        withProgress(message = 'KEGG Started', detail = "interactive", value = 0, {

          i <- input$gotable_rows_selected

          pid <- inputGOstart()$table$ID[i]

          drawKEGG(input, datForTables(), pid)
          list(src = paste0(pid,".b.2layer.png"),
               contentType = 'image/png')
        })

      }, deleteFile = TRUE)

      getGOCatGenes <- reactive({
        print(input$gotable_rows_selected)
        if(is.null(input$gotable_rows_selected)) return (NULL)
        org <- input$organism
        dat <- tabledat()
        i <- input$gotable_rows_selected
        genedata <- getEntrezTable(inputGOstart()$enrich_p$geneID[i],
                                   dat[[1]], org)
        dat[[1]] <- genedata
        dat
      })
      output$GOGeneTable <- DT::renderDataTable({
        shiny::validate(need(!is.null(input$gotable_rows_selected),
                             "Please select a category in the GO/KEGG table to be able
                to see the gene list"))
        dat <- getGOCatGenes()
        if (!is.null(dat)){
          DT::datatable(dat[[1]],
                        list(lengthMenu = list(c(10, 25, 50, 100),
                                               c("10", "25", "50", "100")),
                             pageLength = 25, paging = TRUE, searching = TRUE)) %>%
            getTableStyle(input, dat[[2]], dat[[3]], buttonValues$startDE)
        }
      })

      output$getColumnsForTables <-  renderUI({
        if (is.null(table_col_names())) return (NULL)
        selected_list <- table_col_names()
        if (!is.null(input$table_col_list)
            && all(input$table_col_list %in% colnames(tabledat()[[1]])))
          selected_list <- input$table_col_list
        colsForTable <- list(
          wellPanel(id = "tPanel",
                    style = "overflow-y:scroll; max-height: 200px",
                    checkboxGroupInput("table_col_list", "Select col to include:",
                                       table_col_names(),
                                       selected=selected_list)
          )
        )
        return(colsForTable)
      })
      table_col_names <- reactive({
        if (is.null(tabledat())) return (NULL)
        colnames(tabledat()[[1]])
      })
      tabledat <- reactive({
        dat <- datForTables()
        if (is.null(dat)) return (NULL)
        if (nrow(dat[[1]])<1) return(NULL)
        dat2 <- removeCols(c("ID", "x", "y","Legend", "Size"), dat[[1]])

        pcols <- c(names(dat2)[grep("^padj", names(dat2))],
                   names(dat2)[grep("pvalue", names(dat2))])
        if (!is.null(pcols) && length(pcols) > 1)
          dat2[,  pcols] <- apply(dat2[,  pcols], 2,
                                  function(x) format( as.numeric(x), scientific = TRUE, digits = 3 ))
        else
          dat2[,  pcols] <- format( as.numeric( dat2[,  pcols] ),
                                    scientific = TRUE, digits = 3 )
        rcols <- names(dat2)[!(names(dat2) %in% pcols)]
        if (!is.null(rcols) && length(rcols) > 1)
          dat2[,  rcols] <- apply(dat2[,  rcols], 2,
                                  function(x) round( as.numeric(x), digits = 2))
        else
          dat2[,  rcols] <-  round( as.numeric(dat2[,  rcols]), digits = 2)

        dat[[1]] <- dat2
        return(dat)
      })
      output$tables <- DT::renderDataTable({
        dat <- tabledat()
        if (is.null(dat) || is.null(table_col_names())
            || is.null(input$table_col_list) || length(input$table_col_list)<1)
          return (NULL)
        if (!all(input$table_col_list %in% colnames(dat[[1]]), na.rm = FALSE))
          return(NULL)
        if (!dat[[2]] %in% input$table_col_list)
          dat[[2]]= ""
        if (!dat[[3]] %in% input$table_col_list)
          dat[[3]]= ""

        datDT <- DT::datatable(dat[[1]][, input$table_col_list],
                               options = list(lengthMenu = list(c(10, 25, 50, 100),
                                                                c("10", "25", "50", "100")),
                                              pageLength = 25, paging = TRUE, searching = TRUE)) %>%
          getTableStyle(input, dat[[2]], dat[[3]], buttonValues$startDE)
        return(datDT)
      })
      getMostVaried <- reactive({
        dat <- init_data()
        if (!is.null(cols()))
          dat <- init_data()[,cols()]
        getMostVariedList(dat, colnames(dat), input)
      })
      output$gotable <- DT::renderDataTable({
        if (!is.null(inputGOstart()$table)){
          DT::datatable(inputGOstart()$table,
                        list(lengthMenu = list(c(10, 25, 50, 100),
                                               c("10", "25", "50", "100")),
                             pageLength = 25, paging = TRUE, searching = TRUE))
        }
      })
      mergedComp <- reactive({
        dat <- applyFiltersToMergedComparison(
          getMergedComparison(isolate(dc()), choicecounter$nc, input),
          choicecounter$nc, input)
        dat[dat$Legend == "Sig", ]
      })

      datasetInput <- function(addIdFlag = FALSE){
        tmpDat <- NULL
        sdata <- NULL
        if (input$selectedplot != "QC Heatmap"){
          sdata <- selectedData()
        }else{
          sdata <- isolate(selectedData())
        }
        if (buttonValues$startDE) {
          mergedCompDat <- NULL
          if (input$dataset == "comparisons"){
            mergedCompDat <- mergedComp()
          }
          tmpDat <- getSelectedDatasetInput(rdata = filt_data(),
                                            getSelected = sdata, getMostVaried = getMostVaried(),
                                            mergedCompDat, input = input)
        }
        else{
          tmpDat <- getSelectedDatasetInput(rdata = init_data(),
                                            getSelected = sdata,
                                            getMostVaried = getMostVaried(),
                                            input = input)
        }
        if(addIdFlag)
          tmpDat <- addID(tmpDat)
        return(tmpDat)
      }
      output$metaFile <-  renderTable({
        read.delim(system.file("extdata", "www", "metaFile.txt",
                               package = "debrowser"), header=TRUE, skipNul = TRUE)
      })
      output$countFile <-  renderTable({
        read.delim(system.file("extdata", "www", "countFile.txt",
                               package = "debrowser"), header=TRUE, skipNul = TRUE)
      })

      output$downloadData <- downloadHandler(filename = function() {
        paste(input$dataset, "csv", sep = ".")
      }, content = function(file) {
        dat <- datForTables()
        dat2 <- removeCols(c("x", "y","Legend", "Size"), dat[[1]])
        if(!("ID" %in% names(dat2)))
          dat2 <- addID(dat2)
        write.table(dat2, file, sep = ",", row.names = FALSE)
      })

      output$downloadGOPlot <- downloadHandler(filename = function() {
        paste(input$goplot, ".pdf", sep = "")
      }, content = function(file) {
        pdf(file)
        print( inputGOstart()$p )
        dev.off()
      })
    },
    err=function(errorCondition) {
      cat("in err handler")
      message(errorCondition)
    },
    warn=function(warningCondition) {
      cat("in warn handler")
      message(warningCondition)
    })


}


#' startDEBrowser
#'
#' Starts the DEBrowser to be able to run interactively.
#'
#' @note \code{startDEBrowser}
#' @return the app
#'
#' @examples
#'     startDEBrowser()
#'
#' @export
#'
startPCBrowser <- function(){
  if (interactive()) {
    #the upload file size limit is 30MB
    options( shiny.maxRequestSize = 30 * 1024 ^ 2, warn = -1,
             shiny.sanitize.errors = TRUE)
    addResourcePath(prefix = "demo", directoryPath =
                      system.file("extdata", "demo",
                                  package = "debrowser"))
    addResourcePath(prefix = "www", directoryPath =
                      system.file("extdata", "www",
                                  package = "debrowser"))
    environment(deServer) <- environment()

    app <- shinyApp( ui = shinyUI(deUI),
                     server = shinyServer(deServer))
    runApp(app)
  }
}


startPCBrowser()

dir()


