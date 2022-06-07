#Shiny_Template.R
library(shiny)
library(tidyverse)
source("~/GitHub/Flavi_Sequence_Collections/ShinyServer/CodonAlignment.R") 
CustomHeader <- dashboardHeader(title='datatest')

ui<-fluidPage(
  titlePanel(title = "Zika Research Consortium (ZiRCon) CprME Expression Library"),
  # *Input() functions
  # dateInput("date","Date")  # date
  # dateRangeInput("dateR","Date Range")
  # selectInput("virus","Virus") #virus 
  # radioButtons("virus","Select Virus(es)",choices = unique(flavis$virus)), #virus
  # sliderInput("N","Number of Sequences") #number
  # textInput("seq","Sequence Search") #sequence
  
  # *Output() functions
  fluidRow(
  column(2,
           shiny::plotOutput(outputId = "vpie")
    ),
  column(2,
           dataTableOutput(outputId = "flavis")
    )
  ),
  
fluidRow(
  column(2,
         shiny::plotOutput(outputId = "Tree")
  )
))


server<-function(input, output){
  setwd("~/GitHub/Flavi_Sequence_Collections/ShinyServer/")
  flavisFile<-"Flavis_List.csv"
  flavis<-fread(flavisFile)
  flavis[,ZiRCon_ID:=str_split(Name,"\\|")[[1]][1],by=Name]
  flavis[,accession:=str_split(Name,"\\|")[[1]][2],by=Name]
  flavis[,virus:=str_split(str_split(Name,"\\|")[[1]][3],"_")[[1]][1],by=Name]
  
  #Filter DT
  virus_rows <- reactive({
    flavis[virus %in% input$filter1,   which = TRUE]
  })
  
  '''print(virus_rows)
  filter2_rows <- reactive({
    iris[Sepal.Length < input$filter2, which = TRUE]
  })
  
  filter3_rows <- reactive({
    iris[Sepal.Width > input$filter3,  which = TRUE]
  })
  '''
  
  output$flavis <- renderDataTable(
    #final_rows <- intersect(final_rows,     filter3_rows())
    final_rows <- virus_rows,
    iris[final_rows]
  )
  
  
  output$flavis<-renderDataTable(flavis)
  
  output$vpie<-renderPlot(
    ggplot(flavis)+theme_void()+theme(axis.text.x = element_blank())+
      geom_bar(stat="identity",aes(x="", y=virus, fill=virus))+
      #ggtitle("Viral Species")+
      xlab(paste("N =",length(flavis$virus)))+
      scale_fill_viridis_d(option = "magma")+
      coord_polar(theta = "y"))
    
  output$Tree<-renderPlot(
      ggplot(flavis)+theme_void()+theme(axis.text.x = element_blank())+
        geom_bar(stat="identity",aes(x="", y=virus, fill=virus))+
        #ggtitle("Viral Species")+
        xlab(paste("N =",length(flavis$virus)))+
        scale_fill_viridis_d(option = "magma")+
        coord_polar(theta = "y")
  )
}

shinyApp(ui=ui,server=server)
