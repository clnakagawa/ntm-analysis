library(shiny)
library(ggplot2)
library(reticulate)
library(dplyr)
library(DT)

use_python("myenv/bin/python")


# snp distance cutoff for displaying results
cutoff <- 5

ui <- fluidPage(
  titlePanel("NTM ID Dashboard"),
      navbarPage("Options",
                # using kraken to make tentative ID
                 tabPanel("Kraken Check",
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("fqs", "FASTQ Files",
                                        accept = c(".fastq")),
                              textInput("kSample", "Name"),
                              actionButton("check", "Run"),
                              width=6
                            ),
                            mainPanel(
                              dataTableOutput('kTable') 
                            )
                          )),
                 
                 # Running a sample from fastq for each sequence
                 # input fastq, process fastq to fasta, run fasta
                 tabPanel("Run Sample", 
                          sidebarLayout(
                            sidebarPanel(
                              fileInput("allSeqs", "All sequences",
                                        accept = c(".fastq")),
                              textInput("newSample", "Name"),
                              actionButton("start", "Run"),
                              selectInput(
                                "sp",
                                "Choose Species to Reference",
                                list.dirs(path="mmRef", full.names = F),
                                selected = NULL,
                                multiple = FALSE,
                                selectize = TRUE,
                                width = NULL,
                                size = NULL
                              ),
                              width=6
                            ),
                            mainPanel(
                              dataTableOutput('resTable'), 
                              fluidRow(
                                column(6, plotOutput('snpplot', width=500, height=500)),
                                column(6, plotOutput('covplot', width=500, height=500))
                              ),
                              fluidRow(
                                column(4, plotOutput('cov16S', width=500, height=300)),
                                column(4, plotOutput('cov23S', width=500, height=300)),
                                column(4, plotOutput('covatpD', width=500, height=300))
                              ),
                              fluidRow(
                                column(4, plotOutput('covgroL', width=500, height=300)),
                                column(4, plotOutput('covrpoB', width=500, height=300)),
                                column(4, plotOutput('covtuf', width=500, height=300))
                              )
                            ))),
                 
                 # searching past sample results
                 tabPanel("Search Sample", 
                          sidebarLayout(
                            sidebarPanel(
                              selectInput(
                                "oldSample",
                                "Choose run to view",
                                dir(path="query"),
                                selected = NULL,
                                multiple = FALSE,
                                selectize = TRUE,
                                width = NULL,
                                size = NULL
                              ),
                              actionButton("browse", "Search"),
                              width=6
                            ),
                            mainPanel(
                              dataTableOutput('table'), 
                              fluidRow(
                                column(6, plotOutput('qsnpplot', width=500, height=500)),
                                column(6, plotOutput('qcovplot', width=500, height=500))
                              ),
                              fluidRow(
                                column(4, plotOutput('scov16S')),
                                column(4, plotOutput('scov23S')),
                                column(4, plotOutput('scovatpD'))
                              ),
                              fluidRow(
                                column(4, plotOutput('scovgroL')),
                                column(4, plotOutput('scovrpoB')),
                                column(4, plotOutput('scovtuf'))
                              )
                            )))
  )
)
