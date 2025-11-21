# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

gene_information_integration <- function(res){

  res$strand[res$strand == -1] <- "-"
  res$strand[res$strand == 1] <- "+"

  refgenome <- readRDS("data/refgenome_model.rds")

  refgenome <- refgenome %>%
    dplyr::filter(grepl(paste0("^", res$gene),
                        gene,
                        ignore.case = TRUE)) %>%
    dplyr::filter(is.na(exon)==FALSE & is.na(transcript)==FALSE)

  chr <- as.character(unique(res$chromosome))

  atrack <- AnnotationTrack(res)

  itrack <- IdeogramTrack(chromosome = chr, genome = "hg38")

  GR_res <- GRanges(
    seqnames = unique(res$chromosome),
    IRanges(start = res$start, end = res$end),
    strand = res$strand,
    gene_name = res$symbol,
    gene_id = res$gene,
    ref = res$ref,
    alt = res$alt,
    conseq = res$consequence,
    pos = res$start,
    id = res$id
  )

  gtrack <- GenomeAxisTrack()

  variants <- AnnotationTrack(GR_res, chromosome = chr, name = "Variants",
                              strand = GR_res@strand,
                              group = GR_res$id,
                              genome = "hg38",
                              fill = "red", color = "red")

  # refgenome_subset <- subsetByOverlaps(refgenome, GR_res)
  genetrack <- GeneRegionTrack(refgenome, chromosome = chr, name = as.character(res$symbol[1]),
                               trancriptAnnotation = "transcript")

  getOption("Gviz.scheme")
  scheme <- getScheme()
  scheme$GeneRegionTrack$fill <- "blue"
  scheme$GeneRegionTrack$col <- NULL
  scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
  addScheme(scheme, "myScheme")
  options(Gviz.scheme = "myScheme")

  plotTracks(list(genetrack, variants),
             background.panel = "#FFFEDB",
             background.title = "black", groupAnnotation = "group", just.group = "left")
}

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
