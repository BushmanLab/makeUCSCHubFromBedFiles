#' grab reads from sites within a region
#' Rscript ~/makeUCSCHubFromBedFiles/grabSites.R -w . -r chr1:1-100 -o tmp.csv
#' 

options(stringsAsFactors=FALSE)

#' set all argumentgs for the script
#' @return list of argumentgs
#' @example set_args()
#'          set_args(c("~", "-g=test"))
#' Rscript ~/intSiteUploader/intSiteUploader.R 
set_args <- function(...) {
    ## arguments from command line
    suppressMessages(library(argparse))
    
    parser <- ArgumentParser(description="generate bed file from a sample folder")
    parser$add_argument("-w", "--workDir", nargs=1,
                        default='.',
                        help="Result sample folder from intSiteCaller")
    parser$add_argument("-r", "--region", nargs=1,
                        default='chrX:48,425,222-48,436,669',
                        help="Result sample folder from intSiteCaller")
    parser$add_argument("-o", "--output", nargs=1,
                        default='region.csv',
                        help="Result sample folder from intSiteCaller")
    args <- parser$parse_args(...)
    
    args$workDir <- normalizePath(args$workDir, mustWork=TRUE)
    
    args$region <- gsub(',', '', args$region)
    
    return(args)
}
##set_args()
args <- set_args()
print(args)

libs <- c("dplyr", "stats", "stringr", "GenomicRanges",
          "BiocGenerics", "parallel", "IRanges", "GenomeInfoDb")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))


region2gr <- function(regs) {
    df <- str_match(regs, "(chr.*):(\\d+)-(\\d+)")
    df <- data.frame(chr=df[,2],
                     start=as.integer(df[,3]),
                     end=as.integer(df[,4]))
    gr <- makeGRangesFromDataFrame(df)
    return(gr)
}

target <- region2gr(args$region)


rawSitesFile <- list.files(args$workDir, pattern="rawSites.RData",
                           recursive=TRUE)

sitesInTg <- lapply(rawSitesFile, function(f) {
                        load(f)
                        hits <-  findOverlaps(allSites, target)
                        idx <- as.data.frame(hits)$queryHits
                        return(as.data.frame(allSites[idx]))
                    } )
sitesInTg <- plyr::rbind.fill(sitesInTg)

sitesInTg$width <- NULL
sitesInTg$revmap <- NULL
sitesInTg$pairingID <- NULL

write.csv(sitesInTg, file=args$output, row.names=FALSE, quote=FALSE)

q()

