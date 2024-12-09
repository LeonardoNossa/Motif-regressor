library(Biostrings)
library(shiny)
library(shinydashboard)
library(DT)
library(future)
library(future.apply)
library(shinyjs)
library(seqinr)
library(httr)  
library(jsonlite) 

ui <- dashboardPage(
  dashboardHeader(
    title = "Motif Regressor",
    tags$li(
      class = "dropdown",
      actionButton("info", label = NULL, icon = icon("info-circle"),
                   style = "font-size: 20px; color: #2C3E50; border: none; background: none;")
    )
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "Data", icon = icon("upload")),
      menuItem("PFMs to PWMs", tabName = "PFMs_to_PWMs", icon = icon("table")),
      menuItem("Scores", tabName = "Scores", icon = icon("star")),
      menuItem("Regression", tabName = "Regression", icon = icon("chart-line")),
      tags$div(
        style = "position: absolute; bottom: 20px; left: 70px;",
        actionButton("Github", label = NULL, icon = icon("github"),
                     style = "font-size: 20px; color: #2C3E50; border: none; background: none;")
      )
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML(".skin-blue .main-header .logo {background-color: #2C3E50; color: white;}
                       .skin-blue .main-header .logo:hover {background-color: #2C3E50;}
                       .skin-blue .main-header .navbar  {background-color: #C8C8E8;}
                       .skin-blue .main-sidebar {background-color: #C8C8E8;}
                       .main-sidebar .sidebar .sidebar-menu a:hover {background-color: #2C3E50; border-left: #2C3E50;}
                       .main-sidebar .sidebar .sidebar-menu a {color: #4B0082;}
                       .main-sidebar .sidebar .sidebar-menu .active a {background-color: #2C3E50; border-left: #2C3E50;}
                       .skin-blue .sidebar-toggle:hover {background-color: #2C3E50 !important; color: white;}
                       .skin-blue .box.box-primary .box-header {background-color: #2C3E50 !important; color: white;}
                       .box.box-primary {border: 1px solid #2C3E50 !important;}
                       .box.box-success .box-header {background-color: #FF9900 !important;}
                       .box.box-success {border: 1px solid #FF9900 !important;}"))
    ),
    tabItems(
      tabItem(tabName = "Data",
              fluidRow(
                box(
                  title = "Upload Files", status = "primary", solidHeader = TRUE, width = 6,
                  box(status = "primary", width = 12,
                      div(
                        fluidRow(
                          column(6, 
                                 radioButtons("database_choice", "Choose Database:",
                                              choices = c("JASPAR" = "jaspar", "HOCOMOCO" = "hocomoco", "MyData" = "mydata"),
                                              selected = "jaspar")
                          ),
                          column(6, 
                                 conditionalPanel(
                                   condition = "input.database_choice == 'jaspar' | input.database_choice == 'hocomoco'",
                                   textInput("Tax_ID", "Tax_ID")
                                 )
                          )
                        )
                      ),
                      conditionalPanel(
                        condition = "input.database_choice == 'mydata'",
                        fileInput("meme_file", "Add .meme file", multiple = FALSE, accept = c(".meme"))
                      )
                  ),
                  fileInput("gff_file", "Add .gff file", multiple = FALSE, accept = c(".gff")),
                  fileInput("fna_file", "Add .fna file", multiple = FALSE, accept = c(".fna")),
                  actionButton("Conversion_PFM", "Conversion PFM", icon = icon("bolt"))
                ),
                box(
                  title = "Pre-processing", status = "primary", width = 6,
                  p("Number of motifs: ", strong(textOutput("num_motifs", inline = TRUE))),
                  p("Number of sequences: ", strong(textOutput("num_sequences", inline = TRUE))),
                  p("\n"),
                  p("*********************"),
                  p("\n"),
                  p("Total scores: ", strong(textOutput("operations", inline = TRUE))),
                  actionButton("Start_scoring", "Start scoring", icon = icon("star"))
                )
              )
      ),
      tabItem(tabName = "PFMs_to_PWMs",
              fluidRow(
                box(
                  title = "List of PFMs", status = "primary", solidHeader = TRUE, width = 6,
                  DTOutput("List_of_PFMs")
                ),
                box(
                  title = "List of PWMs", status = "primary", solidHeader = TRUE, width = 6,
                  DTOutput("List_of_PWMs")
                )
              ),
              fluidRow(
                box(
                  title = "Selected PFM", status = "success", solidHeader = TRUE, width = 6,
                  tableOutput("Selected_PFM")
                ),
                box(
                  title = "Selected PWM", status = "success", solidHeader = TRUE, width = 6,
                  tableOutput("Selected_PWM")
                )
              )
      ),
      tabItem(tabName = "Scores",
              fluidRow(
                box(
                  title = "Score for all motifs", status = "primary", solidHeader = TRUE, width = 6,
                  DTOutput("Scores")
                ),
                box(
                  title = "Score distribution", status = "primary", solidHeader = TRUE, width = 6,
                  uiOutput("Score_plot")
                )
              )
      ),
      tabItem(tabName = "Regression",
              p("Regression functionality to be implemented")
      )
    )
  )
)




get_jaspar_motif <- function(tax_id) {
  
  motif_id = c()
  url <- glue::glue("https://jaspar.elixir.no/api/v1/species/{tax_id}/")
  query <- list(identifiers = tax_id)
  
  response <- GET(url, query = query)
  content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
  
  while (TRUE){
    
    content <- fromJSON(content(response, as = "text", encoding = "UTF-8"))
    results <- content$results
    urls <- results$url
    url <- content$`next`
    if (length(url) != 0){
      response <- GET(url, query = query)
      motif_id <- c(motif_id, urls)
    } else {
      motif_id <- c(motif_id, urls)
      break
    }
  }
  motif_id <- vapply(X = motif_id, FUN = function(x, split) {return(unlist(strsplit(x, split))[7])} , 
                     FUN.VALUE = character(1), split = "/")
  names(motif_id) <- c()
  
  matrices <- list()
  genes <- c()
  
  for (motif in motif_id) {
    url <- paste0("https://jaspar.elixir.no/api/v1/matrix/",motif,"/")
    response <- GET(url)
    if (status_code(response) == 200) {
      data <- content(response, "text", encoding = 'UTF-8')
      
      json_data <- fromJSON(data)
      
      genes <- c(genes, json_data$name)
      
      df <- data.frame("A"=json_data$pfm$A, "C"=json_data$pfm$C,
                       "G"=json_data$pfm$G, "T"=json_data$pfm$`T`)
      matrices <- append(matrices, list(df))
      
    } else {
      print(paste("Errore nella richiesta:", status_code(response)))
    }
    
  }
  
  names_df <- paste0(genes, "_", motif_id)
  names(matrices) <- names_df
  return(matrices)
  
}

PFM_loader <- function(path) {
  Tab <- read.table(path, skip = 9, fill = TRUE, stringsAsFactors = FALSE)
  
  motif_indices <- grep("MOTIF", Tab$V1)
  motifs <- list()
  
  for (i in seq_along(motif_indices)) {
    start <- motif_indices[i]
    end <- if (i < length(motif_indices)) motif_indices[i + 1] - 1 else nrow(Tab)
    
    motif_name <- Tab[start, 2]
    motif_matrix <- Tab[(start + 2):(end - 1), 1:4]
    colnames(motif_matrix) <- c("A", "C", "G", "T")
    motifs[[motif_name]] <- as.data.frame(apply(motif_matrix, 2, as.numeric), stringsAsFactors = FALSE)
  }
  return(motifs)
}

PFM2PWM <- function(PFMs, background = c(0.25, 0.25, 0.25, 0.25)) {
  PWMs <- list()
  for (name in names(PFMs)) {
    p <- as.matrix(PFMs[[name]]) + 0.01
    p <- p / rowSums(p)
    p <- log2(p / background)
    PWMs[[name]] <- as.data.frame(p)
  }
  return(PWMs)
}

get_sequences <- function(gff,fasta){
  data <- read.table(gff, 
                     sep = "\t", 
                     header = FALSE, 
                     quote = "",
                     comment.char = "#",
                     stringsAsFactors = FALSE)
  
  colnames(data) <- c("seqname","source","feature","start","end","score","strand", "frame","group")
  genome_size <- data[data$feature == "region",]$end
  
  genes <- data[data$feature == "gene",]
  
  genes_fw <- genes[genes$strand == "+",]
  genes_rw <- genes[genes$strand == "-",]
  
  tmp = strsplit(x=genes_fw$group, split=";")
  GeneName_fw = unlist(lapply(X = tmp, function(x)x[1])) 
  GeneName_fw= sub(pattern="ID=gene-",
                   x=GeneName_fw,
                   replacement = "")
  
  tmp2 = strsplit(x=genes_rw$group, split=";")
  GeneName_rw = unlist(lapply(X = tmp2, function(x)x[1]))
  GeneName_rw= sub(pattern="ID=gene-",
                   x=GeneName_rw,
                   replacement = "")
  
  regions_fw <- cbind(start_sequence = ifelse(genes_fw$start-300 >= 1,genes_fw$start-300,(genome_size-300)+genes_fw$start),
                      end_sequence = genes_fw$start-1) 
  regions_fw <- regions_fw[-1,]
  
  regions_rw <- cbind(start_sequence = genes_rw$end+1,
                      end_sequence = ifelse(genes_rw$end+300 <= genome_size,genes_rw$end+300,(genes_rw$end+300)-genome_size))
  
  regions_complete <- rbind(regions_fw,regions_rw)
  rownames(regions_complete) <- c(GeneName_fw[-1], GeneName_rw)
  
  FASTA <- Biostrings::readDNAStringSet(fasta)[[1]]
  FASTA <- as.character(FASTA)
  
  extract_seqs <- function(FASTA,row){
    return(substr(FASTA,row[1],row[2]))
  }
  
  all_seqs <- apply(X = regions_complete, MARGIN = 1, FUN = extract_seqs, FASTA = FASTA)
  return(all_seqs)}

scorer <- function(PWM,dataSeq){
  
  n_seq = length(dataSeq)
  Result <-matrix(0,ncol = 2,nrow= n_seq)
  PWM = matrix(as.numeric(as.matrix(PWM)),ncol = 4)
  
  for (s in (1: n_seq)) {
    UPSTREAM = s2c(dataSeq[s])
    MATRIX<-matrix(0,nrow=4,ncol=length(UPSTREAM))
    mA<-which(UPSTREAM=="A")
    mC<-which(UPSTREAM=="C")
    mG<-which(UPSTREAM=="G")
    mT<-which(UPSTREAM=="T")
    
    MATRIX[1,mA]<-1
    MATRIX[2,mC]<-1
    MATRIX[3,mG]<-1
    MATRIX[4,mT]<-1
    
    n_col = dim(PWM)[2] 
    ScoreMatrix<-matrix(0,ncol = 1, nrow=(length(UPSTREAM)-n_col))
    
    for (i in 1:(dim(MATRIX)[2] - n_col)-1){
      seq = MATRIX[,i:(nrow(PWM))]
      score = sum(diag(PWM %*% seq))
      ScoreMatrix[i,1] <- score
    }
    MaxIndex = which(ScoreMatrix == max(ScoreMatrix))
    Result[s,1] = MaxIndex[1]
    Result[s,2] = max(ScoreMatrix)[1]
  }
  return(Result)
}

AllScores <- function(PWMs, dataSeq){
  plan(multisession, workers = 4) 
  Scores <- future_lapply(PWMs, function(pwm) scorer(pwm, dataSeq))
  Total <- do.call(cbind, lapply(Scores, function(score) score[, 2]))
  plan(sequential)
  rownames(Total) = names(dataSeq)
  return(Total)
}




server <- function(input, output, session) {
  
  observeEvent(input$info, {
    showModal(modalDialog(
      title = "How to use the tool?",
      "Write instructions on how to use the app.",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  PFMs <- reactiveVal(NULL)
  PWMs <- reactiveVal(NULL)
  dataseq <- reactiveVal(NULL)
  Scores <- reactiveVal(NULL)
  

  observeEvent(input$Conversion_PFM, {
    req(input$gff_file, input$fna_file)
    
    if (input$database_choice == "jaspar" && input$Tax_ID != "") {
      PFMs(get_jaspar_motif(input$Tax_ID))
    } else if (input$database_choice == "mydata") {
      req(input$meme_file) 
      PFMs(PFM_loader(input$meme_file$datapath))
    }
    
    req(PFMs())
    PWMs(PFM2PWM(PFMs()))
    dataseq(get_sequences(input$gff_file$datapath, input$fna_file$datapath))
    
    output$num_motifs <- renderText({ length(PFMs()) })
    output$num_sequences <- renderText({ length(dataseq()) })
    output$operations <- renderText({ length(PWMs()) * length(dataseq()) })
    
    output$List_of_PFMs <- renderDT({
      req(PFMs())
      data.frame(Motifs = names(PFMs()))
    }, selection = "single")
    
    output$List_of_PWMs <- renderDT({
      req(PWMs())
      data.frame(Motifs = names(PWMs()))
    }, selection = "single")
  })
  
  observeEvent(input$List_of_PFMs_rows_selected, {
    selected <- input$List_of_PFMs_rows_selected
    req(selected)
    motif_name <- names(PFMs())[selected]
    output$Selected_PFM <- renderTable(PFMs()[[motif_name]], rownames = TRUE)
  })
  
  observeEvent(input$List_of_PWMs_rows_selected, {
    selected <- input$List_of_PWMs_rows_selected
    req(selected)
    motif_name <- names(PWMs())[selected]
    output$Selected_PWM <- renderTable(PWMs()[[motif_name]], rownames = TRUE)
  })
  

  observeEvent(input$Start_scoring, {
    req(PWMs(), dataseq())
    
    Scores(AllScores(PWMs(), dataseq()))
    
    output$Scores <- renderDT({
      req(Scores())
      datatable(as.data.frame(Scores()), options = list(scrollX = TRUE))
    })
  })
  
  observeEvent(input$Github, {
    browseURL("https://github.com/LeonardoNossa/Motif-regressor.git")
  })
}


shinyApp(ui, server)
