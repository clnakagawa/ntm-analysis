library(shiny)
library(ggplot2)
library(reticulate)
library(dplyr)
library(DT)

use_python("/home/carter_nakagawa_health_ny_gov/ntm-analysis/myenv/bin/python")


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
                              dataTableOutput('kTable'), 
                            )
                          )),
                 
                 # Running a sample from fasta's for each sequence
                 # change to run from fastqs?
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
                              fluidRow(plotOutput('snpplot', width=500, height=500)),
                              fluidRow(plotOutput('covplot', width=500, height=500))
                            ))),
                 
                 # searching past sample results
                 tabPanel("Search Sample", 
                          sidebarLayout(
                            sidebarPanel(
                              textInput("oldSample", "Browse Results by Sample"),
                              actionButton("browse", "Search"),
                              width=6
                            ),
                            mainPanel(
                              dataTableOutput('table'), 
                              fluidRow(plotOutput('qsnpplot', width=500, height=500)),
                              fluidRow(plotOutput('qcovplot', width=500, height=500))
                            )))
  )
)
