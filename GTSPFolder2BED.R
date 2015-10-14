#' Generate bed file from a GTSP sample folder
#' Rscript ~/makeUCSCHubFromBedFiles/GTSPFolder2BED.R GTSP0515-10
#' To generate all bed files, run
#' for i in $(ls -d GTSP*/); do echo ${i%%//}; Rscript ~/makeUCSCHubFromBedFiles/GTSPFolder2BED.R ${i%%//}; done

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
    parser$add_argument("workDir", nargs='?',
                        default='.',
                        help="Result sample folder from intSiteCaller")
    args <- parser$parse_args(...)
    
    args$workDir <- normalizePath(args$workDir, mustWork=TRUE)
    
    return(args)
}
##set_args()
args <- set_args()
print(args)

libs <- c("dplyr", "RMySQL", "stats", "GenomicRanges",
          "BiocGenerics", "parallel", "IRanges", "GenomeInfoDb")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))


## load uniq hits sites
load(file.path(args$workDir, "sites.final.RData"))
load(file.path(args$workDir, "allSites.RData"))

sites <- data.frame(
    "sampleID"=sites.final$sampleName,
    "chr"=as.character(seqnames(sites.final)),
    "strand"=as.character(strand(sites.final)),
    "position"=start(flank(sites.final, -1, start=T)),
    "idx"=seq(sites.final))

gtspid <- unique(sites.final$sampleName)
stopifnot(length(gtspid)==1)
stopifnot(grepl(gtspid, args$wordDir))

revmap <- sites.final$revmap
allBreakpoint <- start(flank(allSites, -1, start=F))

bp <- plyr::ldply(seq(revmap), function(i)
    data.frame(
        "idx"=i,
        "breakpoint"=allBreakpoint[revmap[[i]]] ) )

bp <- bp %>% dplyr::group_by(idx, breakpoint) %>%
    dplyr::mutate(count=n()) %>% 
        dplyr::distinct()

sites.uniq <- merge(sites, bp, by="idx")


## load multi hits sites
load(file.path(args$workDir, "multihitData.RData"))

sites.multi <- data.frame(
    "chr"=as.character(seqnames(multihitData[[1]])),
    "strand"=as.character(strand(multihitData[[1]])),
    "position"=start(flank(multihitData[[1]], -1, start=T)),
    "breakpoint"=start(flank(multihitData[[1]], -1, start=FALSE)))
sites.multi <- dplyr::distinct(sites.multi)

allsites <- plyr::rbind.fill(sites.uniq, sites.multi)
allsites <- dplyr::select(allsites, chr, strand, position, breakpoint)

## write to file
bed <- plyr::arrange(allsites, chr, position, strand, breakpoint)

bed <- data.frame(chrom=as.character(bed$chr),
                  chromStart=as.integer(pmin(bed$position, bed$breakpoint)),
                  chromEnd=as.integer(pmax(bed$position, bed$breakpoint)),
                  name=args$workDir,
                  score=500,
                  strand=as.character(bed$strand))

fileName <- paste(gtspid, "bed", sep=".")

write.table(bed, file=fileName,
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
message("File written to ", fileName)


