server <- function(input, output, session) {


  #theme handler
  observe(session$setCurrentTheme(
    if (isTRUE(input$dark_mode)) dark else light
  ))

  #app logo
  output$myImage <- renderImage({
    list(
      src = "/www/VariantInvestigator_logo_blue.png",
      contentType = "image/png",
      width = 75,
      height = 50
    )
  }, deleteFile = FALSE)

  #allowing big files in shiny
  options(shiny.maxRequestSize=3000*1024^2)

  #file input server side processing
  shinyFileChoose(input, 'files', root=c(root='.'), filetypes=c('sqlite', 'txt'))

  #reactive function for database creation
  sqlfunc <- reactive({

    #error message for no file input
    s_file <- parseFilePaths(
      root=c(root='.'), input$files)
    validate(
      need(s_file$datapath != "", message = "Please select a data set"
      )
    )

    #sql database connection
    con <- dbConnect(RSQLite::SQLite(),
                     s_file$datapath
    )

    #first filter and data loading in memory
    sqloutput <-  dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE ?"),
                              params = list(c(paste(input$impact)))
    ) %>%
      dbFetch()

    #genotype table fetching
    Genotype <- dbGetQuery(con, "SELECT AD, DP, GT, GQ, start
                     FROM genotypes")

    #genotype and variants merging in single table
    sqloutput <- merge(sqloutput, Genotype, by.x="start", by.y="start", all.x = TRUE) %>%
      dplyr::mutate(gnomad_af = as.numeric(gnomad_af)) %>%
      relocate("ad", .before = "allele") %>%
      relocate("dp", .before = "ad") %>%
      relocate("gt", .before = "allele") %>%
      relocate("start", .before = "end") %>%
      dplyr::select(-"ac")

    #conversion of certain columns to facilitate browsing
    sqloutput$dp <- as.numeric(sqloutput$dp)

    #labfindings filter
    sqloutput <- sqloutput %>%  mutate(pubmed_count = case_when(
      sqloutput$pubmed == "" ~ 0,
      sqloutput$pubmed == "." ~ 0,
      is.na(sqloutput$pubmed) ~ 0,
      is.null(sqloutput$pubmed) ~ 0,
      TRUE ~ str_count(sqloutput$pubmed, pattern = "&") + 1
    ))

    #applying filters
    sqloutput <- sqloutput %>%
      dplyr::filter(gnomad_af <= input$gnomad |
                      is.na(gnomad_af) == TRUE) %>%
      dplyr::filter(grepl(paste(input$HC, collapse = "|"),
                          lof)) %>%
      dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
      dplyr::filter(grepl(paste0("^", input$searchgene),
                          symbol,
                          ignore.case = TRUE)) %>%
      dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                          consequence)) %>%
      dplyr::filter(grepl(paste(c(input$clinsign), collapse = "|"),
                          clin_sig)) %>%
      dplyr::filter(grepl(paste(c(input$genotype), collapse = "|"),
                          gt)) %>%
      dplyr::filter(case_when(
        input$labf != 100 ~ pubmed_count <= input$labf,
        input$labf == 100 ~ is.numeric(pubmed_count)
      ))


    sqloutput <- sqloutput %>%
      dplyr::filter(gq >= input$gq | is.na(gq)==TRUE)

    sqloutput <- sqloutput %>%
      dplyr::filter(dp >= input$dp[1] | is.na(dp) == TRUE)


    sqloutput <- sqloutput %>%
      separate_wider_delim("sift",
                           delim = "(",
                           names = c("sift_prediction","sift_values"),
                           too_few = "align_start") %>%
      separate_wider_delim("polyphen",
                           delim = "(",
                           names = c("polyphen_prediction","polyphen_values"),
                           too_few = "align_start")

    #conversion of sift column to numeric
    sqloutput$sift_values <- gsub(")", "", sqloutput$sift_values) %>%
      as.numeric(sqloutput$sift_values)

    #same for polyphen
    sqloutput$polyphen_values <- gsub(")", "", sqloutput$polyphen_values) %>%
      as.numeric(sqloutput$polyphen_values)


    sqloutput <- sqloutput %>%
      dplyr::filter(sift_values <= input$SIFT |
                      is.na(sift_values) == TRUE) %>%
      dplyr::filter(polyphen_values >= input$PolyPhen |
                      is.na(polyphen_values) == TRUE) %>%
      dplyr::select(-"pubmed")


    sqloutput <- sqloutput %>%  mutate(Varsome = paste0(sqloutput$chromosome, ":", sqloutput$start, ":", sqloutput$ref, ":", sqloutput$alt))

    sqloutput <- sqloutput %>%
      relocate(symbol) %>%
      relocate(consequence, .before = id) %>%
      relocate(clin_sig, .before = chromosome) %>%
      relocate(gq, .after = gt) %>%
      rowwise() %>%
      mutate(ref = case_when(nchar(ref) > 10 ~
                               paste(str_sub(ref, 1, 10), "..."),
                             nchar(ref) <= 10 ~ ref)) %>%
      mutate(alt = case_when(nchar(alt) > 10 ~
                               paste(str_sub(alt, 1, 10), "..."),
                             nchar(alt) <= 10 ~ alt))

    sqloutput$ad <- gsub(",", ", ", sqloutput$ad)

    sqloutput <- sqloutput %>%
      mutate(chr = gsub("chr", "", chromosome))

    sqloutput$chr <- sqloutput$chr %>%
      as.numeric(sqloutput$chr)

    sqloutput

  })

  #variants table output
  output$table <- DT::renderDataTable({

    #correct display of variants table
    data_vars <- sqlfunc() %>%
      dplyr::mutate(Varsome = paste0("<a href = 'https://varsome.com/variant/hg38/", Varsome, "' target='blank' >", id,  "</a>")) %>%
      dplyr::mutate(Varsome = map(Varsome, gt::html)) %>%
      relocate(Varsome, .before = clin_sig) %>%
      dplyr::select(-"gene") %>%
      dplyr::select(-"id") %>%
      dplyr::group_by(chr, start) %>%
      dplyr::arrange(chr) %>%
      dplyr::select(-"chr")



    #retrieve filtered dataset for download
    filtered <- reactive({
      get("data_vars")
    })

    #download handler server
    output$down <- downloadHandler(
      filename = function() {
        paste0("filt_vars", ".xlsx")
      },
      content = function(file) {
        # Write the dataset to the `file` that will be downloaded
        write.table(filtered(), file, sep = " ", row.names = FALSE)
      }
    )

    DT::datatable(
      data = data_vars,
      class = 'cell-border stripe',
      plugin = "ellipsis",
      filter = 'top', selection = 'none',
      height = "100%", width = "100px",
      extensions = 'Buttons',
      options = (
        list(
          buttons = c('copy', 'csv', 'excel')
        ))
    )
  }, server = TRUE)

  #tissue informations tab
  output$dropdownbutt2 <- renderUI({

    #ordering genes in alphabetical order for convenience
    data_vars <- sqlfunc() %>%
      dplyr::select(symbol) %>%
      arrange(symbol)

    pickerInput(
      "gene2",
      "Gene expression by tissue",
      choices = c(unique(data_vars$symbol)),
      options = c(`actions-box` = TRUE),
      multiple = TRUE
    )


  })

  output$notefilter <- renderUI({
    sliderInput(
      inputId = "note",
      "NoTE (Normalised Tissue Expression) filter",
      min = 0,
      max = 1,
      step = 0.005,
      value = 0
    )
  })

  output$tissue_info <- DT::renderDataTable({




    tvars <- sqlfunc() %>%
      mutate(var_position = paste0(chromosome, ":", start, "_", ref, "/", alt)) %>%
      dplyr::select("symbol", "var_position") %>%
      rename("symbol" = "gene")

    tissues_results <- dbConnect(RSQLite::SQLite(), "data/express_data.sqlite")
    vars_expression <- dbSendQuery(tissues_results, paste("SELECT *
                          FROM expression_data
                          WHERE gene LIKE ?"),
                                   params = list(c(paste(input$gene2)))
    ) %>%
      dbFetch()




    #creation of the NoTE score
    vars_expression <- merge(tvars, vars_expression) %>%
      arrange(Tissue) %>%
      dplyr::group_by(Tissue, .add = TRUE)  %>%
      mutate(maxTissue = max(nTPM)) %>%
      mutate(nTissue = nTPM/maxTissue) %>%
      ungroup()

    vars_expression <- vars_expression %>%
      group_by(gene, .add = TRUE) %>%
      mutate(ppsum = sum(nTissue, na.rm = TRUE)) %>%
      ungroup()



    # creazione dello score BiTE
    vars_expression <- vars_expression %>%
      mutate(
        cutoff_met = case_when(nTPM < 5 ~ 0,
                               nTPM >= 5 ~ 1)
      )

    vars_expression <- vars_expression %>%
      group_by(gene) %>%
      mutate(BiTE = sum(cutoff_met, na.rm = TRUE)) %>%
      ungroup()

    vars_expression <- vars_expression %>%
      mutate(NoTE = ppsum/max(ppsum)) %>%
      dplyr::arrange(desc(NoTE)) %>%
      dplyr::select("gene", "var_position", "NoTE", "BiTE") %>%
      dplyr::distinct() %>%
      dplyr::filter(NoTE >= input$note |
                      is.na(NoTE) == TRUE)





    DT::datatable(
      data = vars_expression,
      class = 'cell-border stripe',
      plugin = "ellipsis",
      filter = 'top', selection = 'none',
      height = "100%", width = "100px"
    )

  })

  #gene selector
  output$dropdownbutt <- renderUI({

    #ordering genes in alphabetical order for convenience
    data_vars <- sqlfunc() %>%
      dplyr::select(symbol) %>%
      arrange(symbol)

    #selector of gene for graphical representation
    pickerInput(inputId = 'gene',
                label = 'Genes chosen',
                choices = c("",
                            unique(data_vars$symbol)
                )
    )
  })

  output$symbol_num <- renderDT({

    #number of variants per gene table
    data_vars <- sqlfunc() %>%
      dplyr::select(symbol) %>%
      dplyr::group_by(symbol) %>%
      dplyr::summarise(total_vars = dplyr::n(),
                       .groups = 'drop') %>%
      dplyr::arrange(desc(total_vars))


    DT::datatable(
      data = data_vars,
      class = 'cell-border stripe',
      plugin = "ellipsis",
      filter = 'top', selection = 'none',
      height = "100%", width = "100px"
    )

  }, server=TRUE)

  #gene rendering functions
  output$genestbl <- renderPlot({
    req(input$gene)

    sqlfunc() %>%
      dplyr::filter(symbol == input$gene) %>%
      gene_information_integration()
  })
}
