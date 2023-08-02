library(shiny)
library(ggplot2)
library(reticulate)
library(dplyr)
library(DT)
library(reshape2)

# snp distance cutoff for displaying results


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
    qres <- qres[qres$diff < input$cutoff,]
    qres <- qres[,c(1,6,2,5,3)]
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
    print(head(qres))
    qres <- qres[qres$diff < input$qcutoff,]
    qres <- qres[,c(1,6,2,5,3)]
    qres
  })
  output$table <- DT::renderDataTable(datatable(qcontents(), 
                                                extensions="Buttons", 
                                                options=list(dom = 'Bfrtip', 
                                                             buttons=c('copy','csv','excel')))
  )  
  snpplot <- reactive({
    resTab <- contents()
    df = data.frame(species=resTab$sp, total=resTab$diff, snp=resTab$snp, gap=resTab$gdiff)
    df$species = factor(df$species, level=df$species)
    df = melt(df, id.vars=c("species"))
    ggplot(df, aes(species, value, fill=variable)) + geom_bar(stat='Identity',position=position_dodge())+
      ggtitle("SNP Plot") +
      theme(text = element_text(size=20), 
            plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$snpplot <- renderPlot(snpplot())
    
  qsnpplot <- reactive({
    resTab <- qcontents() 
    df = data.frame(species=resTab$sp, total=resTab$diff, snp=resTab$snp, gap=resTab$gdiff)
    df$species = factor(df$species, level=df$species)
    df = melt(df, id.vars=c("species"))
    ggplot(df, aes(species, value, fill=variable)) + geom_bar(stat='Identity',position=position_dodge())+
      ggtitle("SNP Plot") +   
      theme(text = element_text(size=20),
            plot.title = element_text(face="bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5, size=rel(0.8)))
  })
  output$qsnpplot <- renderPlot(qsnpplot())
  
  # coverage plots for targets
  # change this list to set new targets
  targets = list("16S", "23S", "atpD", "groL", "rpoB", "tuf")
  
  # for viewing past runs
  output$plots <- renderUI({
    plot_output_list <- lapply(targets, function(t) {
      plotname <- paste(t, "plot", sep="")
      column(width=3,
        plotOutput(plotname, height=280, width=250)
      )
    })
    do.call(tagList, plot_output_list)
  })
  
  observeEvent(input$browse, {
    for (t in targets) { 
      local({
        target <- t
        plotname <- paste(target, "plot", sep="")
        tab <- read.table(paste("fqProcessed/",input$oldSample,"/",target,"coverage.txt", sep=""), header=F)
        colnames(tab) <- c("sp", "position", "depth")
        output[[plotname]] <- renderPlot(
          ggplot(tab) + geom_bar(aes(x=position, y=depth), stat="identity", fill="skyblue") +
            ggtitle(paste(target, "coverage")) + theme(plot.title=element_text(size=20, face="bold"))
        )
      })
    }    
  })
  
  # for running samples
  output$qplots <- renderUI({
    plot_output_list <- lapply(targets, function(t) {
      plotname <- paste(t, "qplot", sep="")
      column(width=3,
        plotOutput(plotname, height=280, width=250)
      )
    })
    do.call(tagList, plot_output_list)
  })

  observeEvent(input$start, {
    for (t in targets) {
      local({
        target <- t
        plotname <- paste(target, "qplot", sep="")
        tab <- read.table(paste("fqProcessed/",input$newSample,"/",target,"coverage.txt", sep=""), header=F)
        colnames(tab) <- c("sp", "position", "depth")
        output[[plotname]] <- renderPlot(
          ggplot(tab) + geom_bar(aes(x=position, y=depth), stat="identity", fill="skyblue") +
            ggtitle(paste(target, "coverage")) + theme(plot.title=element_text(size=20, face="bold"))
        )
      })
    }
  }, priority=-1)
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
