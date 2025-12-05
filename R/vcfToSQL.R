#' VCFtoSQL shiny app
#'
#' @description
#' Launches a Shiny application to convert a VCF file into an SQLite database.
#' The database will contain tables for variants, genotypes, and rejected lines.
#'
#' @import shiny
#' @import shinyFiles
#' @importFrom reticulate source_python
#'
#' @return A Shiny app object
#'
#' @export
vcfToSQL <- function(){

  # Define UI for the application
  ui <- fluidPage(
    titlePanel("Convert a VCF file to SQLite database"),
    verticalLayout(
      shinyFilesButton("vcf_file", "Select VCF input file", "Please select a VCF file", multiple = FALSE),
      textOutput("vcf_path"),
      shinySaveButton("sqlite_file", "Select SQLite output file", "Please select output database file"),
      textOutput("sqlite_path"),
      actionButton("run_convert", "Convert")
    )
  )

  # Define server logic
  server <- function(input, output, session) {

    script_path <- system.file("convert_vcf_to_sql.py", package = "variantsinvestigator")
    if (script_path == "") {
      # Fallback for local testing if package is not installed but loaded via load_all or source
      script_path <- "inst/convert_vcf_to_sql.py"
      if (!file.exists(script_path)) {
         # Try looking one level up just in case (e.g. if working directory is R/)
         script_path <- "../inst/convert_vcf_to_sql.py"
      }
    }
    
    if (file.exists(script_path)) {
        source_python(script_path)
    } else {
        stop("Python script convert_vcf_to_sql.py not found.")
    }

    volumes <- getVolumes()()

    shinyFileChoose(input, "vcf_file",
                                roots = volumes,
                                session = session,
                                filetypes = c("vcf", "gz"))

    shinyFileSave(input, "sqlite_file",
                              roots = volumes,
                              session = session)

    output$vcf_path <- renderText({
      req(input$vcf_file)
      vcf_info <- parseFilePaths(volumes, input$vcf_file)
      if (nrow(vcf_info) == 0) return("")
      as.character(vcf_info$datapath)
    })

    output$sqlite_path <- renderText({
      req(input$sqlite_file)
      db_info <- parseSavePath(volumes, input$sqlite_file)
      if (nrow(db_info) == 0) return("")
      as.character(db_info$datapath)
    })

    observeEvent(input$run_convert, {
      req(input$vcf_file, input$sqlite_file)

      vcf_info <- parseFilePaths(volumes, input$vcf_file)
      db_info  <- parseSavePath(volumes, input$sqlite_file)

      req(nrow(vcf_info) > 0, nrow(db_info) > 0)

      vcf_path <- as.character(vcf_info$datapath)
      db_path  <- as.character(db_info$datapath)

      convertVCFtoSQL(vcf_path, db_path)
      showNotification("Conversion completed!", type = "message")
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)

}
