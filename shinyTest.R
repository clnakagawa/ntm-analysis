library(shiny)
library(ggplot2)
library(reticulate)
library(DT)

setwd("/Users/cln/PycharmProjects/ntm-analysis")
use_python("venv/bin/python")
use_virtualenv("/Users/cln/PycharmProjects/ntm-analysis/venv")
getwd()

cutoff <- 5

ui <- fluidPage(
  titlePanel("NTM ID Dashboard"),
  sidebarLayout(
    sidebarPanel(
      fileInput("allSeqs", "All sequences",
                accept = c("text/fasta",".fasta",".txt",".fna")),
      textInput("newSample", "Name"),
      actionButton("start", "Run"),
      textInput("oldSample", "Browse Results by Sample"),
      actionButton("browse", "Search")
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Run Sample", 
                           dataTableOutput('resTable'), 
                           fluidRow(plotOutput('snpplot', width=500, height=500)),
                           fluidRow(plotOutput('covplot', width=500, height=500))),
                  tabPanel("Search Sample", 
                           dataTableOutput('table'), 
                           fluidRow(plotOutput('qsnpplot', width=500, height=500)),
                           fluidRow(plotOutput('qcovplot', width=500, height=500))))
    )
  )
)

source_python("cleanFTP.py")
server <- function(input, output, session){
  
  # so that paths work
  setwd("/Users/cln/PycharmProjects/ntm-analysis")
  
  contents <- eventReactive(input$start, {
    # run analysis
    dir.create(paste("query/", input$newSample, sep=""))
    write(readLines(input$allSeqs$datapath, warn=F), paste("query/", input$newSample, "/allSeqs.txt", sep=""))
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
      theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size=rel(0.8)))
  })
  output$snpplot <- renderPlot(snpplot())
  
  covplot <- reactive({
    resTab <- contents()
    ggplot(resTab, aes(reorder(sp, snp), cov)) + xlab("species") + 
      geom_bar(stat="identity", fill='blue', width=0.5) +
      theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size=rel(0.8)))
  })
  output$covplot <- renderPlot(covplot())
  
  qsnpplot <- reactive({
    resTab <- qcontents()
    ggplot(resTab, aes(reorder(sp, snp), snp)) + xlab("species") + 
      geom_bar(stat="identity",  fill='red', width=0.5) +
      theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size=rel(0.8)))
  })
  output$qsnpplot <- renderPlot(qsnpplot())
  
  qcovplot <- reactive({
    resTab <- qcontents()
    ggplot(resTab, aes(reorder(sp, snp), cov)) + xlab("species") + 
      geom_bar(stat="identity", fill='blue', width=0.5) +
      theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1, size=rel(0.8)))
  })
  output$qcovplot <- renderPlot(qcovplot())
}

shinyApp(ui = ui, server=server)