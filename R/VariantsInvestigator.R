#' VariantsInvestigator app core
#'
#' @description
#' Launches the main Shiny application for investigating variants from the generated SQLite database.
#'
#' @import bslib
#' @import DBI
#' @importFrom DT datatable
#' @import GenomicFeatures
#' @import GenomicRanges
#' @import gt
#' @importFrom Gviz AnnotationTrack IdeogramTrack GenomeAxisTrack GeneRegionTrack plotTracks
#' @import IRanges
#' @import reticulate
#' @import RSQLite
#' @import shiny
#' @importFrom dplyr mutate filter select arrange relocate group_by summarise rename ungroup distinct case_when n rowwise
#' @importFrom tidyr separate_wider_delim
#' @importFrom purrr map
#' @importFrom utils write.table
#' @importFrom stringr str_count str_sub str_detect
#' @import tibble
#' @import shinyFiles
#' @import shinyWidgets
#'
#' @return A Shiny app object
#'
#' @export
VariantsInvestigator <- function(){

  ui <- page_sidebar(
    sidebar = sidebar(
      div(img(
        src = "VariantInvestigator_logo_blue.png",
        height = 130,
        width = 200,
        style = "margin:1px 1px"
      ), ""),


      #standard theme
      theme = get_light_theme(),

      #night mode switch
      materialSwitch
      (
        "dark_mode",
        "Dark mode",
        status = "primary"
      ),

      #input file
      shinyFilesButton
      (
        'files',
        label='Upload file',
        title='Please select an SQL file',
        multiple=FALSE,
        buttonType = "primary"
      ),

      #lof filter
      div
      (
        pickerInput
        (
          "HC",
          "LoF variants",
          choices = c("HC"),
          selected = c(""),
          multiple = TRUE
        ),
        shiny::tags$small(style = "color: gray; font-size: 12px;", "High-confidence loss-of-function variants")
      ),

      #impact filter
      pickerInput
      (
        inputId = "impact",
        label = "Variant impact",
        choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
        selected = c("HIGH"),
        options = c(`actions-box` = TRUE),
        multiple = TRUE
      ),

      #clinical significance filter
      pickerInput
      (
        inputId = "clinsign",
        label = "ClinVar classifications",
        choices = c("benign",
                    "uncertain_significance",
                    "pathogenic",
                    "drug_response",
                    "association",
                    "affects",
                    "other",
                    "conflicting",
                    "risk_factor",
                    "not_provided",
                    "protective",
                    "association_not_found",
                    "confers_sensitivity",
                    ".",
                    "NA"),
        selected = c("pathogenic", "conflicting", "NA"),
        options = c(`actions-box` = TRUE),
        multiple = TRUE
      ),

      #genotype filter
      pickerInput
      (
        inputId = "genotype",
        label = "Genotype",
        choices = c("0/1",
                    "1/1",
                    "1/0",
                    "0|1",
                    "1|1",
                    "1|0"),
        selected = c("0/1",
                     "1/1",
                     "1/0",
                     "0|1",
                     "1|1",
                     "1|0"),
        options = c(`actions-box` = TRUE),
        multiple = TRUE
      ),

      #labfindings
      div(
        sliderInput
        (
          inputId = "labf",
          "Maximum number of publications",
          min = 0,
          max = 100,
          step = 1,
          value = 100
        ),
        shiny::tags$small(style = "color: gray; font-size: 12px;", "When choosing 100, all variants are included (even if they have more than 100 publications)")
      ),


      #genotype quality filter
      sliderInput
      (
        inputId = "gq",
        "Min genotype quality",
        min = 0,
        max = 99,
        step = 1,
        value = 10
      ),

      #consequence filter
      pickerInput(
        inputId = "conseq",
        label = "Variant consequences",
        choices =
          c(
            "transcript_ablation",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
            "feature_elongation",
            "feature_truncation",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "splice_donor_5th_base_variant",
            "splice_region_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
            "incomplete_terminal_codon_variant",
            "start_retained_variant",
            "stop_retained_variant",
            "synonymous_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "coding_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "regulatory_region_variant",
            "intergenic_variant",
            "sequence_variant",
            "NA"
          ),
        selected =
          c(
            "transcript_ablation",
            "splice_acceptor_variant",
            "splice_donor_variant",
            "stop_gained",
            "frameshift_variant",
            "stop_lost",
            "start_lost",
            "transcript_amplification",
            "feature_elongation",
            "feature_truncation",
            "inframe_insertion",
            "inframe_deletion",
            "missense_variant",
            "protein_altering_variant",
            "splice_donor_5th_base_variant",
            "splice_region_variant",
            "splice_donor_region_variant",
            "splice_polypyrimidine_tract_variant",
            "incomplete_terminal_codon_variant",
            "start_retained_variant",
            "stop_retained_variant",
            "synonymous_variant",
            "coding_sequence_variant",
            "mature_miRNA_variant",
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant",
            "intron_variant",
            "NMD_transcript_variant",
            "non_coding_transcript_variant",
            "coding_transcript_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "TFBS_ablation",
            "TFBS_amplification",
            "TF_binding_site_variant",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "regulatory_region_variant",
            "intergenic_variant",
            "sequence_variant",
            "NA"
          ),
        options = c(`actions-box` = TRUE),
        multiple = TRUE
      ),

      #chromosome choice filter
      pickerInput
      (
        inputId = "chr",
        label = "Chromosomes in analysis",
        choices = c(paste0("chr", as.character(c(1:22))),
                    "chrX",
                    "chrY",
                    "chrM"),
        selected =
          c(
            paste0(
              "chr", as.character(c(1:22))),
            "chrX",
            "chrY",
            "chrM"
          ),
        options = list(`actions-box` = TRUE),
        multiple = TRUE
      ),

      #gnomad filter
      sliderInput
      (
        inputId = "gnomad",
        "GnomAD allele frequency",
        min = 0,
        max = 1,
        step = 0.005,
        value = 0.5
      ),

      #allele depth filter
      sliderInput
      (
        inputId = "dp",
        "Read depth", #Minimum depth read
        min = 0, max = 200,
        value = 0
      ),

      #sift filter
      sliderInput
      (
        inputId = "SIFT",
        "SIFT value", #Maximum SIFT value
        min = 0,
        max = 1,
        step = 0.01,
        value = 1
      ),

      #polyphen filter
      sliderInput
      (
        inputId = "PolyPhen",
        "PolyPhen value", #Minimum
        min = 0,
        max = 1,
        step = 0.01,
        value = 0
      ),

      #gene search filter
      div(
        textInput
        (
          "searchgene",
          "Gene name(s)",
          value = ""
        ),
        shiny::tags$small(style = "color: gray; font-size: 12px;", "Use `|` separator for multiple genes")
      ),

      #download button widget
      downloadButton("down", "Download Dataset", class = "btn-primary")
    ),



    #tabs
    navset_tab(

      #variant table tab
      nav_panel(
        "Variants",
        DT::DTOutput(outputId = "table"),
        fillable = FALSE,

      ),

      #tissue informations tab
      nav_panel(
        "Tissues expression scores",
        uiOutput("dropdownbutt2"),
        uiOutput("notefilter"),
        DT::DTOutput(outputId = "tissue_info")

      ),

      #genome view tab
      nav_panel(
        "Affected transcripts",
        uiOutput("dropdownbutt"),


        plotOutput("genestbl")
      ),

      #number of mutations per gene tab
      nav_panel(
        "Mutations on each gene",
        DT::DTOutput(outputId = "symbol_num")

      )
    )
  )

  server <- function(input, output, session) {


    #theme handler
    observe(session$setCurrentTheme(
      if (isTRUE(input$dark_mode)) get_dark_theme() else get_light_theme()
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
        roots=c(root='.'), input$files)
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

    output$symbol_num <- DT::renderDT({

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

  shinyApp(ui, server)
}
