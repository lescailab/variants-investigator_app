
#' Object defining the app's light mode ui layout
#'
#' @name light
#'
#'
#' @format ## `R object`
#' An object describing the ui light mode's colour patterns
#' \describe{
#'  \item{light}{a bstheme object containing information on theme and colours of the application's light mode}
#'  }
#'
#'

"light"


#' Object defining the app's dark mode ui layout
#'
#' @name dark
#'
#'
#' @format ## `R object`
#' An object describing the ui dark mode's colour patterns
#' \describe{
#'  \item{dark}{a bstheme object containing information on theme and colours of the application's dark mode}
#'  }
#'
#'

"dark"



#' Dataset containing exome references for gene matching
#'
#' @name refgenome_model
#'
#'
#' @format ## `tibble`
#' A tibble with 10 columns
#' \describe{
#'  \item{chromosome}{the chromosome containing the feature}
#'  \item{start}{start position of the feature}
#'  \item{end}{end position of the feature}
#'  \item{width}{number of nucleotides of the feature}
#'  \item{strand}{strand on which the feature is positioned}
#'  \item{feature}{type of data regarding the gene, either transcript, exon or gene}
#'  \item{gene}{Ensembl's standard gene reference}
#'  \item{exon}{Ensembl's standard exon reference}
#'  \item{transcript}{Ensembl's standard transcript reference}
#'  \item{symbol}{gene name}
#'  }
#' @source created with TxDb.Hsapiens.UCSC.hg19.knownGene R package
#'

"refgenome_model"

