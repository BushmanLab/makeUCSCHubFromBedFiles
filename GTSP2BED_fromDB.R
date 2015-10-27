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
    
    args <- p$parse_args(...)
    
    return(args)
}
## get_args(c("GTSP001", "GTSP002"))
stopifnot(get_args(c("A", "B"))$gtspid == c("A", "B"))
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


processGTSP <- function(gtspid) {
    stopifnot(length(gtspid)==1)
    message("\nProcessing\t", gtspid, "\n")
    
    ## check if file exist and permission .my.cnf 
    stopifnot(file.exists("~/.my.cnf"))
    stopifnot(file.info("~/.my.cnf")$mode == as.octmode("600"))
    
    ## initialize connection to database
    ## ~/.my.cnf must be present
    junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
    dbConn <- dbConnect(MySQL(), group="intsites_miseq.read") 
    stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)
    
    allSampleName <- suppressWarnings( dbGetQuery(dbConn, "SELECT * FROM samples") )
    
    ## get all replicates
    replicates <- subset(allSampleName,
                         grepl(paste0("^",gtspid), sampleName) &
                         refGenome==args$freeze)
    if( nrow(replicates)==0 ) return() 
    
    ## get sampleInfo
    sql <- sprintf("SELECT * FROM specimen_management.gtsp WHERE SpecimenAccNum = '%s'", gtspid)
    sampleInfo <- suppressWarnings( dbGetQuery(dbConn, sql) )
    stopifnot( nrow(sampleInfo)==1)
    colnames(sampleInfo) <- tolower(colnames(sampleInfo))
    
    sampleIDin <- sprintf("(%s)", paste(unique(replicates$sampleID), collapse=","))
    
    ##get unique sites
    sql <- sprintf("SELECT DISTINCT *
                FROM samples JOIN sites
                ON samples.sampleID = sites.sampleID
                JOIN pcrbreakpoints
                ON pcrbreakpoints.siteID = sites.siteID 
                WHERE samples.refGenome = '%s' AND
                samples.sampleID in %s", args$freeze, sampleIDin)
    message(sql)
    sites.uniq <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
    sites.uniq <- sites.uniq[, !duplicated(colnames(sites.uniq))]
    
    ##get multihit sites
    sql <- sprintf("SELECT *
                FROM samples JOIN multihitpositions
                ON samples.sampleID = multihitpositions.sampleID
                JOIN multihitlengths
                ON multihitpositions.multihitID = multihitlengths.multihitID
                WHERE samples.refGenome = '%s' AND
                samples.sampleID in %s",  args$freeze, sampleIDin)
    message(sql)
    sites.multi <- suppressWarnings( dbGetQuery(dbConn, sql) ) 
    sites.multi <- sites.multi[, !duplicated(colnames(sites.multi))]
    sites.multi$breakpoint <- ifelse(sites.multi$strand=="+",
                                     sites.multi$position+sites.multi$length-1,
                                     sites.multi$position-sites.multi$length+1)
    
    ##output bed file
    needed <- c("chr", "position", "breakpoint", 
                "sampleName", "refGenome", "strand")
    
    bed <- rbind(subset(sites.uniq, select=needed),
                 subset(sites.multi, select=needed))
    
    bed <- plyr::arrange(bed, chr, position, strand)
    
    ##bed$name <- gtspid
    bed$name <- rep(paste(sampleInfo$patient, sampleInfo$timepoint, sampleInfo$celltype, sep=":"),
                    nrow(bed))
    bed$score <- rep(400, nrow(bed))
    
    bed <- plyr::arrange(bed, chr, position, strand)
    
    bed <- data.frame(chrom=as.character(bed$chr),
                      chromStart=as.integer(pmin(bed$position, bed$breakpoint)),
                      chromEnd=as.integer(pmax(bed$position, bed$breakpoint)),
                      name=as.character(bed$name),
                      score=rep(500, nrow(bed)),
                      strand=as.character(bed$strand))
    
    fileName <- paste(gtspid, sampleInfo$patient, sampleInfo$timepoint, sampleInfo$celltype, "bed", sep=".")
    fileName <- paste(gtspid, "bed", sep=".")
    
    write.table(bed, file=fileName,
                quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    message(fileName,"\t", nrow(bed), "\trows")
    
    return(list(fileName=fileName, reads=nrow(bed)))
}

res <- sapply(args$gtspid, processGTSP)
message()
print(as.data.frame(t(res)))

