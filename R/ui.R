ui <- page_sidebar(
  sidebar = sidebar(
    div(img(
      src = "VariantInvestigator_logo_blue.png",
      height = 130,
      width = 200,
      style = "margin:1px 1px"
    ), ""),


    #standard theme
    theme = light,

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
      DTOutput(outputId = "table"),
      fillable = FALSE,

    ),

    #tissue informations tab
    nav_panel(
      "Tissues expression scores",
      uiOutput("dropdownbutt2"),
      uiOutput("notefilter"),
      DTOutput(outputId = "tissue_info")

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
      DTOutput(outputId = "symbol_num")

    )
  )
)
