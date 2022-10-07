#' Fetching SNP positions based on their names
#' 
#' @param snp_names - character vector with names of the SNPs
#' @param attribs - which attributes do you want to fetch? Check `listAttributes`
#'   in {biomaRt} package
#' @param asGRanges - should the result be returned in GRanges format? (default:
#'   FALSE)
#' @param genome_ver - by default, function will fetch information from the
#'   latest genome release; if any specific version is required, please provide
#'   it here (e.g., "37" will give GRCh37)
#'   
#' @return data.frame or GRanges object with the fetched information; columns 
#'   are named as the given attributes
grabSNPpositions <- function(
    snp_names,
    attribs,
    asGRanges = FALSE,
    genome_ver = NULL
){
  if(!requireNamespace('biomaRt')){
    stop("Couldn't find 'biomaRt' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('biomaRt')
  if(!requireNamespace('regioneR')){
    stop("Couldn't find 'regioneR' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('regioneR')

  if(is.null(genome_ver)){
    ensembl <- useEnsembl(
      biomart = "snps",
      dataset = "hsapiens_snp"
    )
  } else {
    ensembl <- useEnsembl(
      biomart = "snps",
      dataset = "hsapiens_snp",
      GRCh = genome_ver
    )
  }
  
  filter <- "snp_filter"
  
  snps <- getBM(
      attributes = attribs,
      filters = filter,
      values = snp_names,
      mart = ensembl
    )

  if(asGRanges){
    if(nrow(snps) == 0){
      return(NULL)
    }
    snps <- toGRanges(
      as.data.frame(
        snps %>%
          mutate(seqnames = paste0("chr", chr_name)) %>%
          dplyr::select(
            seqnames,
            start = chrom_start,
            end = chrom_end,
            everything()
          ) %>%
          distinct()
      )
    )
  }
  return(snps)
}

#' Fetching genes positions based on their names or IDs
#' 
#' @param gene_symbols - character vector with names or IDs of the genes
#' @param filter_name - which filter is used? Check `listFilters` function in
#'   the {biomaRt} package for help (default: "hgnc_symbol")
#' @param attribs - which attributes do you want to fetch? Check `listAttributes`
#'   in {biomaRt} package; default: `c('chromosome_name', 'start_position',
#'   'end_position', 'hgnc_symbol')`
#' @param asGRanges - should the result be returned in GRanges format? (default:
#'   FALSE)
#' @param genome_ver - by default, function will fetch information from the
#'   latest genome release; if any specific version is required, please provide
#'   it here (e.g., "37" will give GRCh37)
#'   
#' @return data.frame or GRanges object with the fetched information; columns 
#'   are named as the given attributes
grabGenesPositions <- function(
    gene_symbols,
    filter_name = "hgnc_symbol",
    attribs = c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
    asGRanges = FALSE,
    genome_ver = NULL
){
  if(!requireNamespace('biomaRt')){
    stop("Couldn't find 'biomaRt' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('biomaRt')
  if(!requireNamespace('regioneR')){
    stop("Couldn't find 'regioneR' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('regioneR')

  if(is.null(genome_ver)){
    ensembl <- useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl"
    )
  } else {
    ensembl <- useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      GRCh = genome_ver
    )
  }
  genes <- getBM(
      attributes = attribs,
      filters = filter_name,
      values = gene_symbols,
      mart = ensembl
    )

  if(asGRanges){
    if(nrow(genes) == 0){
      return(NULL)
    }
    genes <- toGRanges(
      as.data.frame(
        genes %>%
          mutate(seqnames = paste0("chr", chromosome_name)) %>%
          dplyr::select(
            seqnames,
            start = start_position,
            end = end_position,
            everything()
          ) %>%
          distinct()
      )
    )
  }
  return(genes)
}
