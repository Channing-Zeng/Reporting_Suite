

library(shiny)

shinyUI(
  navbarPage("NextSeq500 Run Tools",
             
             
             tabPanel(title="Run Review",
                      selectInput(inputId="txtInRunName",label="Run folder:", choices=list.files("/projects/NGS/projects/DS/PreProcessed/") ),                                       
                      
                      fluidPage(
                        fluidRow(
                          column(4,plotOutput(width=400,height=800,"readPlot")),                                  
                          column(4,plotOutput(width=400,height=800,"BCPlot")),
                          column(4,plotOutput(width=400,height=800,"FilteredPlot"))                           
                        ),
                   
                        fluidRow(
                          titlePanel(title="Mapped reads per organism"),
                          plotOutput(width=1600, height=800,"MetaGenomeMappingPlot")
                        ),
                        fluidRow(
                          titlePanel(title="Mapped reads per chromosome"),
                          plotOutput(width=1600, height=800,"MetaMappingPlot")
                        ),
                        fluidRow(
                          dataTableOutput("DataTableMapping")
                        )
                        
                      )                      
             ),
             tabPanel("Run Submission",
                      selectInput(inputId="txtProjDep",label="Department", choices=c("DS","DSM","RIA","CVMD","ONC","Other")),
                      textInput(inputId="txtProjID",label="Project ID", value="<Input short project name>"),
                      textInput(inputId="txtJIRAID",label="JIRA ID", value="<Input the JIRA Ticket ID>"),              
                      selectInput(inputId="txtRunName",label="Run folder:", choices=list.files("/projects/NGS/projects/DS/fromNextSeq500/") )      
             )
             
  )
)




