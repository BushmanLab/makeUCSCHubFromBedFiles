bedfile <- list.files(".", "*.bed")
Notes <- ifelse( grepl("old", bedfile), "old", "new")
Sample <- sub(".bed", "", bedfile)
freeze <- "hg18"
hub <- "vectorTrimBetaThal"

csv <- data.frame(Sample, bedfile, Notes, freeze, hub)
write.csv(csv, "bed.csv", row.names=FALSE, quote=FALSE)
system("cat bed.csv")


options(stringsAsFactors=F)

get_args <- function(...) {
    suppressMessages(library(argparse))
    
    codeDir <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))
    if( length(codeDir)!=1 ) codeDir <- list.files(path="~", pattern="makeUCSCHubFromBedFiles$", recursive=TRUE, include.dirs=TRUE, full.names=TRUE)
    codeDir <- head(codeDir, 1)
    stopifnot(file.exists(file.path(codeDir, "makeUCSChubFromBedFiles.R")))
    
    p <- ArgumentParser(description="Generate bed file from GTSP")
    p$add_argument("gtspid", nargs='+', default='GTSP0000')
    p$add_argument("-f", "--freeze", type="character", nargs=1,
                   default="hg18",
                   help="hg18, hg19, etc")
    p$add_argument("-c", "--codeDir", type="character", nargs=1,
                   default=codeDir,
                   help="Directory of code")
    
    ##args <- p$parse_args(commandArgs(trailingOnly=TRUE))
    args <- p$parse_args(...)
    
    return(args)
}
## get_args(c("GTSP001", "GTSP002"))
stopifnot(get_args(c("GTSP001", "GTSP002"))$gtspid == c("GTSP001", "GTSP002"))
args <- get_args(c("GTSP0877", "GTSP0878"))
args <- get_args()

print(args)

#' required packages
libs <- c("stats", "methods", "RMySQL")
null <- suppressMessages(sapply(libs, library, character.only=TRUE))

#' increase output width to console width
wideScreen <- function(howWide=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]])[2]) {
    options(width=as.integer(howWide))
}
wideScreen()

