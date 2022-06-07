library(shinydashboard)
library(data.table)
library(stringr)
library(msaR)
library(DECIPHER)
library(Biostrings)
library(ape)
library(ggtree)

setwd("~/GitHub/Flavi_Sequence_Collections/ShinyServer/")
flavisFile<-"Flavis_List.csv"#twist csv file (original samples)

flavis<-fread(flavisFile)
newFlDT<-data.table(ZiRConID=flavis[,str_split(Name,"\\|")[[1]][1],by=Name][,2],
  accession=flavis[,str_split(Name,"\\|")[[1]][2],by=Name][,2],
  virus=flavis[,str_split(str_split(Name,"\\|")[[1]][3],"_")[[1]][1],by=Name][,2],
  flavis)
seqs<- Biostrings::DNAStringSet(newFlDT$`Insert sequence`)
names(seqs)<-paste(newFlDT$Name)

CustomHeader <- dashboardHeader(title = 'ZiRCon Database')

ui <- function(request) {
  dashboardPage(
    CustomHeader,
    ## Sidebar content
    dashboardSidebar(
      sidebarMenu(#Set up tabs.
        menuItem("Select Sequences", tabName = "select", icon = icon("dashboard")),
        menuItem("Alignment", tabName = "aaalignment", icon = icon("list-alt")),
        menuItem("Tree", tabName = "tree", icon = icon("th"))
      )
    ),
    ## Body content
    dashboardBody(
      tabItems(
        # First tab content
        tabItem(tabName = "select",#tab for selecting sequences. 
                fluidRow(box(width = 4,
                             checkboxGroupInput(inputId = 'virus',label = 'Virus',choices = sort(unique(newFlDT$virus)))
                             ),
                         box(width = 4,
                             sliderInput(inputId = 'year', value = c(1950,2022), label = 'Year', min=1900,max = 2022,dragRange = T)
                             ),
                         box(width = 4,
                             sliderInput(inputId = 'filter2','filter2:Max.Sepal.Length',min = 0,max = 10,value = 10)
                         )
                        ),
                fluidRow(box(width=9,
                             dataTableOutput('FlaviTable'),collapsible = T
                             )
                )
        ),
                
        
        tabItem(tabName = "aaalignment",
                fluidRow(box(width = 10,
                             msaROutput(outputId = "protmsa")
                ))
                
        ),
        tabItem(tabName = "tree",
                fluidRow(box(width=5,
                             plotOutput(outputId = "tree")
                             ))
                )
        )#end tabitems
      
      
      )#end Dashboard Body
    )
  
}

server<-function(input, output){
  #Filter DT
  virus_rows <- reactive({
    #print(newFlDT[virus.V1 %in% input$virus])
    newFlDT[virus.V1 %in% input$virus]
  })
  
  output$FlaviTable <- renderDataTable(
    #final_rows <- intersect(final_rows,     filter3_rows())
    virus_rows()
  )
  
  MSA=reactive({
    #print(newFlDT[virus.V1 %in% input$virus])
    seqSelect=seqs[newFlDT$virus.V1 %in% input$virus]
    DECIPHER::AlignSeqs(translate(seqSelect))
  })
  
  output$protmsa <- renderMsaR(
    msaR(MSA(),labels = MSA(),labelNameLength = 150,colorscheme = "mae")
  )
  
  
  output$Tree <- renderPlot(
    ggplot(virus_rows())+theme_void()+theme(axis.text.x = element_blank())+
      geom_bar(stat="identity",aes(x="", y=virus.V1, fill=virus.V1))+
      #ggtitle("Viral Species")+
      xlab(paste("N =",length(flavis$virus)))+
      scale_fill_viridis_d(option = "magma")+
      coord_polar(theta = "y")
  )
  
  output$vpie<-renderPlot(
    ggplot(select_rows)+theme_void()+theme(axis.text.x = element_blank())+
      geom_bar(stat="identity",aes(x="", y=virus.V1, fill=virus.V1))+
      #ggtitle("Viral Species")+
      xlab(paste("N =",length(flavis$virus)))+
      scale_fill_viridis_d(option = "magma")+
      coord_polar(theta = "y")
  )
}

shinyApp(ui, server)