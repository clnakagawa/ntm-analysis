library(shiny)
library(ggplot2)
library(reticulate)
library(dplyr)
library(DT)

# snp distance cutoff for displaying results
cutoff <- 5

use_python("myenv/bin/python")
source_python("cleanFTP.py")
server <- function(input, output, session){
  contents <- eventReactive(input$start, {
    # make files
    dir.create(paste("query/", input$newSample, sep=""))
    write(readLines(input$allSeqs$datapath, warn=F), paste("query/", input$newSample, "/sample.fastq", sep=""))
  
    # build index
    system(paste("sh buildIndex.sh ", input$sp , " ", sub(".*_","",input$sp), sep = ""))
    # process fastq
    system(paste("sh processFQ.sh", 
                 paste("query/", input$newSample, "/sample.fastq", sep=""),
                 input$sp,
                 sub(".*_","",input$sp),
                 input$newSample,
                 sep = " "))
    # move sequence file
    system(paste("cp",
                 paste("fqProcessed/",input$newSample,"/allSeqs.txt",sep=""),
                 paste("query/",input$newSample,sep=""),
                 sep=" "))

    # run test function
    source_python("getTest.py")
    test(input$newSample)
    
    # load data
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    qres <- read.csv(paste("query/",sname,"/totalSNPrefFiles.csv",sep=""))
    qres <- qres[qres$snp < cutoff,]
    qres
  })
  output$resTable <- DT::renderDataTable(datatable(contents(), 
                                                   extensions="Buttons", 
                                                   options=list(dom = 'Bfrtip', 
                                                                buttons=c('copy','csv','excel')))
  )
  
  qcontents <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    qres <- read.csv(paste("query/",sname,"/totalSNPrefFiles.csv",sep=""))
    qres <- qres[qres$snp < cutoff,]
    qres
  })
  output$table <- DT::renderDataTable(datatable(qcontents(), 
                                                extensions="Buttons", 
                                                options=list(dom = 'Bfrtip', 
                                                             buttons=c('copy','csv','excel')))
  )  
  snpplot <- reactive({
    resTab <- contents()
    ggplot(resTab, aes(reorder(sp, snp), snp)) + xlab("species") + 
      geom_bar(stat="identity",  fill='red', width=0.5) +
      ggtitle("SNP Plot") +
      theme(text = element_text(size=20), 
            plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$snpplot <- renderPlot(snpplot())
  
  covplot <- reactive({
    resTab <- contents()
    ggplot(resTab, aes(reorder(sp, snp), cov)) + xlab("species") + 
      geom_bar(stat="identity", fill='blue', width=0.5) +
      ggtitle("Alignment Coverage Plot") +
      theme(text = element_text(size=20), 
            plot.title = element_text(face="bold"), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$covplot <- renderPlot(covplot())
  
  qsnpplot <- reactive({
    resTab <- qcontents()
    ggplot(resTab, aes(reorder(sp, snp), snp)) + xlab("species") + 
      geom_bar(stat="identity",  fill='red', width=0.5) +
      ggtitle("SNP Plot") + 
      theme(text = element_text(size=20), 
            plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$qsnpplot <- renderPlot(qsnpplot())
  
  qcovplot <- reactive({
    resTab <- qcontents()
    ggplot(resTab, aes(reorder(sp, snp), cov)) + xlab("species") + 
      geom_bar(stat="identity", fill='blue', width=0.5) +
      ggtitle("Alignment Coverage Plot") +
      theme(text = element_text(size=20), 
            plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$qcovplot <- renderPlot(qcovplot())
  
  # coverage plots for targets
  p16S <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/16Scoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplot16S <- reactive({
    tab <- p16S()
    ggplot(p16S()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("16S coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scov16S <- renderPlot(scovplot16S())
  
  q16S <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/16Scoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplot16S <- reactive({
    tab <- q16S()
    ggplot(q16S()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("16S coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$cov16S <- renderPlot(covplot16S())
  
  p23S <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/23Scoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplot23S <- reactive({
    tab <- p23S()
    ggplot(p23S()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("23S coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scov23S <- renderPlot(scovplot23S())
  
  q23S <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/23Scoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplot23S <- reactive({
    tab <- q23S()
    ggplot(q23S()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("23S coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$cov23S <- renderPlot(covplot23S())
  
  patpD <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/atpDcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplotatpD <- reactive({
    tab <- patpD()
    ggplot(patpD()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("atpD coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scovatpD <- renderPlot(scovplotatpD())
  
  qatpD <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/atpDcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplotatpD <- reactive({
    tab <- qatpD()
    ggplot(qatpD()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("atpD coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$covatpD <- renderPlot(covplotatpD())
  
  pgroL <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/groLcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplotgroL <- reactive({
    tab <- pgroL()
    ggplot(pgroL()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("groL coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scovgroL <- renderPlot(scovplotgroL())
  
  qgroL <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/groLcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplotgroL <- reactive({
    tab <- qgroL()
    ggplot(qgroL()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("groL coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$covgroL <- renderPlot(covplotgroL())
  
  prpoB <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/rpoBcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplotrpoB <- reactive({
    tab <- prpoB()
    ggplot(prpoB()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("rpoB coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scovrpoB <- renderPlot(scovplotrpoB())
  
  qrpoB <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/rpoBcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplotrpoB <- reactive({
    tab <- qrpoB()
    ggplot(qrpoB()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("rpoB coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$covrpoB <- renderPlot(covplotrpoB())
  
  ptuf <- eventReactive(input$browse, {
    sname <- input$oldSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/tufcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  scovplottuf <- reactive({
    tab <- ptuf()
    ggplot(ptuf()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("tuf coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$scovtuf <- renderPlot(scovplottuf())
  
  qtuf <- eventReactive(input$start, {
    sname <- input$newSample
    if (is.null(sname)) return(NULL)
    tab <- read.table(paste("fqProcessed/",sname,"/tufcoverage.txt",sep=""), header=F)
    colnames(tab) <- c("sp", "position", "depth")
    tab
  })
  covplottuf <- reactive({
    tab <- qtuf()
    ggplot(qtuf()) + 
      geom_bar(aes(x=position,y=depth), stat="identity", fill="skyblue") + 
      ggtitle("tuf coverage") + 
      theme(plot.title=element_text(size=20,face="bold"))
  })
  output$covtuf <- renderPlot(covplottuf())
  
  # kraken code
  options(shiny.maxRequestSize=300*1024^2)
  # so that paths work
  kcontent <- eventReactive(input$check, {
    dir.create(paste("kreports/", input$kSample, sep=""))
    system(paste("kraken2 --db ntmDB --use-names --report kreports/",
                 input$kSample,
                 "/report.kreport --output kreports/",
                 input$kSample,
                 "/output.txt ", 
                 input$fqs$datapath, sep=""))
    kres <- read.delim(paste("kreports/",input$kSample,"/output.txt",sep=""),sep="\t", header = F)
    kres <- kres['V3']
    kres <- kres %>% group_by(V3) %>% summarise(total_count=n()) %>% as.data.frame()
    colnames(kres) <- c("ID", "Count")
    kres[order(kres$Count, decreasing=T),]
  })
  
  output$kTable <- DT::renderDataTable(datatable(kcontent(), 
                                                 extensions="Buttons", 
                                                 options=list(dom = 'Bfrtip', 
                                                              buttons=c('copy','csv','excel'))))
}
