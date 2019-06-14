#install.packages("radmixture")
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
# download.file(url = 'https://github.com/wegene-llc/radmixture/raw/master/data/globe4.alleles.RData', destfile = 'globe4.alleles.RData')
# download.file(url = 'http://public-html.biostat.washington.edu/~bsweir/HIBAG/param/Asian-HLA4-hg19.RData', destfile = 'Asian-HLA4-hg19.RData')


library(shiny)
library(HIBAG)
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("HLA分型小工具"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
       fileInput("file1", "Choose txt File",
                multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".txt")),
      tags$hr(),
      selectInput("variable", "请选择一种参考数据模型:",
                  list("亚洲人" = "Asian", 
                       "European" = "European", 
                       "Hispanic" = "Hispanic",
                       "African" = " African"
                  ))
      #checkboxInput("outliers", "Show outliers", FALSE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h4("你的HLA分型结果和准确概率为:"),
      tableOutput(outputId = 'test')
  )
))

# Define server logic required to calculate ancestry
server <- function(input, output) {
  options(shiny.maxRequestSize=30*1024^2)
  output$test <- renderTable({ 
    "You have selected this"
    inFile <- input$file1
    if(is.null(inFile))     
      return(NULL) 
    trained.Rdata <- paste(input$variable, sep = "", "-HLA4-hg19.RData")
    load(trained.Rdata)
    A.model <- hlaModelFromObj(HLA4[['A']])
    B.model <- hlaModelFromObj(HLA4[['B']])
    C.model <- hlaModelFromObj(HLA4[['C']])
    DRB1.model <- hlaModelFromObj(HLA4[['DRB1']])
    DQA1.model <- hlaModelFromObj(HLA4[['DQA1']])
    DQB1.model <- hlaModelFromObj(HLA4[['DQB1']])
    DPB1.model <- hlaModelFromObj(HLA4[['DPB1']])
    cmd = paste('./plink/plink --23file', inFile$datapath, 'make-bed --out',sep = " ", inFile)
    system(cmd, wait = TRUE)
    bed <- paste(inFile, '.bed',sep = "")
    fam <- paste(inFile, '.fam',sep = "")
    bim <- paste(inFile, '.bim',sep = "")
    yourgeno <- hlaBED2Geno(bed.fn=bed, fam.fn=fam, bim.fn=bim)
    # HLA imputation at HLA-A
    
    A.guess <- predict(A.model, yourgeno, type="response+prob",match.type='Position')
    
    
    B.guess <- predict(B.model, yourgeno, type="response+prob",match.type='Position')
    
    
    C.guess <- predict(C.model, yourgeno, type="response+prob",match.type='Position')
    
    DRB1.guess <- predict(DRB1.model, yourgeno, type="response+prob",match.type='Position')
    
    DQA1.guess <- predict(DQA1.model, yourgeno, type="response+prob",match.type='Position')
    
    DQB1.guess <- predict(DQB1.model, yourgeno, type="response+prob",match.type='Position')
    
    DPB1.guess <- predict(DPB1.model, yourgeno, type="response+prob",match.type='Position')
    if(is.null(A.guess))     
      return(NULL) 
    data.frame(A_allele1=A.guess$value$allele1, A_allele2=A.guess$value$allele2, A_prob=A.guess$value$prob,
               B_allele1=B.guess$value$allele1, B_allele2=B.guess$value$allele2, B_prob=B.guess$value$prob,
               C_allele1=C.guess$value$allele1, C_allele2=C.guess$value$allele2, C_prob=C.guess$value$prob,
               DRB1_allele1=DRB1.guess$value$allele1, DRB1_allele2=DRB1.guess$value$allele2, DRB1_prob=DRB1.guess$value$prob,
               DQA1_allele1=DQA1.guess$value$allele1, DQA1_allele2=DQA1.guess$value$allele2, DQA1_prob=DQA1.guess$value$prob,
               DQB1_allele1=DQB1.guess$value$allele1, DQB1_allele2=DQB1.guess$value$allele2, DQB1_prob=DQB1.guess$value$prob,
               DPB1_allele1=DPB1.guess$value$allele1, DPB1_allele2=DPB1.guess$value$allele2, DPB1_prob=DPB1.guess$value$prob)
  })
}
# Run the application 
shinyApp(ui = ui, server = server)