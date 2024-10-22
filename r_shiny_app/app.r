library(shiny)
library(ggplot2)
ui <- fluidPage(
  titlePanel("Z-Hunt Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload FASTA or TXT File", accept = c(".fasta", ".fa",".fna",".txt")),
      sliderInput("windowSize", "Window Size:", min = 1, max = 12, value = 6),
      sliderInput("minWindow", "Minimum Window:", min = 1, max = 12, value = 3),
      sliderInput("maxWindow", "Maximum Window:", min = 1, max = 12, value = 9),
      textInput("chromosomeName", "Chromosome Name:"),
      textInput("startSite", "Start Site:"),
      actionButton("runButton", "Run Z-Hunt"),
      downloadButton("downloadData", "Download Results")
    ),
    
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$runButton, {
    req(input$file)
    
    # Validate inputs
    validate(
      need(input$windowSize, "Please set the Window Size."),
      need(input$minWindow, "Please set the Minimum Window."),
      need(input$maxWindow, "Please set the Maximum Window."),
      need(input$chromosomeName != "", "Please enter a Chromosome Name."),
      need(input$startSite != "", "Please enter a Start Site.")
    )
    
    # Extract just the filename from the full path
    filename <- basename(input$file$name)
    
    # Copy the uploaded file to the app's working directory with its original name
    file.copy(input$file$datapath, filename, overwrite = TRUE)
    
    # Define the command to run z-hunt using only the filename
    command <- paste(
      "./zhr",
      input$windowSize,
      input$maxWindow,
      input$minWindow,
      filename,
      input$chromosomeName,
      input$startSite
    )
    
    # Run z-hunt command
    system(command)
    
    # Assuming z-hunt outputs filename.z-score.bedgraph
    resultFile <- paste0(input$file$name, ".Z-SCORE.bedgraph")
    
    # Check if the result file exists before proceeding
    if (file.exists(resultFile)) {
      # Download handler for the bedgraph file
      output$downloadData <- downloadHandler(
        filename = function() {
          resultFile
        },
        content = function(file) {
          file.copy(resultFile, file)
        }
      )
      
      # Plotting (optional)
   
      output$plot <- renderPlot({
        # Read and process the sequence file
        sequence <- readLines(input$file$datapath)
        sequence <- gsub("\\s+", "", sequence)
        sequence <- unlist(strsplit(sequence, ""))
        
        # Create a data frame for nucleotides with positions
        dna_df <- data.frame(Position = seq_along(sequence), Nucleotide = sequence)
        
        # Read z-scores from the result file
        zscores <- read.table(resultFile, header = FALSE)
        
        # Combine nucleotide data with z-scores
        combined_df <- cbind(dna_df, Score = zscores$V4)
        
        # Plotting
        ggplot(combined_df, aes(x = Position, y = Score)) +
          geom_line(color = "blue") +
          geom_point(size = 0.5, color = "blue") +
          geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
          labs(title = "Scores Across Sequence",
               x = "Position",
               y = "Score") +
          theme_minimal()
      })
      
      
    } else {
      showModal(modalDialog(
        title = "Error",
        "The z-hunt did not generate a .bedgraph file. Please check your inputs and try again."
      ))
    }
  })
}

shinyApp(ui = ui, server = server)
