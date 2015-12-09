## UCSC track doc 
## http://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html
##
## tools needed
## UCSC tools 
## wget -e robots=off -r -nH --cut-dirs=2 --no-parent --reject="index.html*" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
## bedtools https://github.com/arq5x/bedtools2/releases
## 
message("checking required external programs, mysql, UCSC tools, bedtools...")
app <- list(mysql="mysql",
            bedToBigBed="~/opt/UCSC/bedToBigBed",
            genomeCoverageBed="~/opt/bedtools/genomeCoverageBed",
            bedGraphToBigWig="~/opt/UCSC/bedGraphToBigWig",
            bedSort="~/opt/UCSC/bedSort")
null <- sapply(app, function(x) system2("which", x, stdout=FALSE))
if( any(null!=0) ) stop(paste(app[null!=0], collapse=" "), " not availabe")

library(BiocParallel)
ncore <- multicoreWorkers()
BP <- MulticoreParam(ncore)

options(stringsAsFactors=FALSE)

#### input ####
#' @param csvFile, csv file, see below
#' @param hubName, a name for the hub
#' @param freeze, reference genome, hg18, hg19 etc
csvFile <- "bed.csv"
if( length(commandArgs(trailingOnly=TRUE))==1 ) csvFile <- commandArgs(trailingOnly=TRUE)[1]

if( !file.exists(csvFile) ) { message("
Usage: Rscript makeUCSChubFromBedFiles.R file.csv

A sample csv file looks like the following
sample,bedfile,notes,freeze,hub
GTSP0001,GTSP0001.bed,sample0001,hg18,GENE 
GTSP0002,GTSP0002.bed,sample0002,hg18,GENE
GTSP0003,GTSP0003.bed,sample0003,hg18,GENE
")
                              message("\n", csvFile, " not found")
                              q(status=1) }

#### processing files ####
sampleInfo <- read.csv(csvFile)
colnames(sampleInfo) <- tolower(colnames(sampleInfo))
hubName <- unique(sampleInfo$hub)
freeze <- unique(sampleInfo$freeze)
stopifnot(all(file.exists(sampleInfo$bedfile)))
stopifnot(length(freeze)==1)
stopifnot(length(hubName)==1)

## set color for track by notes
colors <- c("235,127,127", "127,127,235")
sampleInfo$color <- colors[as.integer((as.factor(sampleInfo$notes)))%%length(colors)+1]

message("\nProcessing bed files...")
write.table(sampleInfo, "", quote=FALSE, sep="\t", row.names=FALSE)

## get $freeze.genome.sizes from UCSC
cmd <- sprintf("%s --user=genome --host=genome-mysql.cse.ucsc.edu -sNA -e 'select chrom, size from %s.chromInfo' > %s.genome.sizes", app$mysql, freeze, freeze) 
message("\n",cmd)
stopifnot(system(cmd)==0)
stopifnot(file.exists(paste0(freeze, ".genome.sizes")))

## sort bed files
message("\nSort bed files...")
cmd <- sprintf("%s %s %s.sort.bed", 
               app$bedSort,
               sampleInfo$bedfile, 
               sampleInfo$bedfile)
##null <- sapply(cmd, function(x) {message(x); system(x)} )
null <- bplapply(cmd, function(x) {message(x); system(x)}, BPPARAM=BP)
stopifnot(all(null==0))

## convert sorted bed file to bigbed
message("\nconvert bed files to bigBed files...")
cmd <- sprintf("%s %s.sort.bed %s.genome.sizes %s.bb",
               app$bedToBigBed,
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
##null <- sapply(cmd, function(x) {message(x); system(x)} )
null <- bplapply(cmd, function(x) {message(x); system(x)}, BPPARAM=BP)
stopifnot(all(null==0))

## convert sorted bed file to bedGraph file
message("\nconvert bed files to bigGraph files...")
cmd <- sprintf("%s -bg -i %s.sort.bed -g %s.genome.sizes > %s.bedGraph",
               app$genomeCoverageBed,
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
##null <- sapply(cmd, function(x) {message(x); system(x)} )
null <- bplapply(cmd, function(x) {message(x); system(x)}, BPPARAM=BP)
stopifnot(all(null==0))

## convert bedGraph file to bigwig
message("\nconvert bigGraph files to bigWig files...")
cmd <- sprintf("%s %s.bedGraph %s.genome.sizes %s.bw",
               app$bedGraphToBigWig,
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
##null <- sapply(cmd, function(x) {message(x); system(x)} )
null <- bplapply(cmd, function(x) {message(x); system(x)}, BPPARAM=BP)
stopifnot(all(null==0))

#### put files in hub folder, generate trackDB_bw.txt and trackDB_bb.txt ####
message("\nCreating hub ",hubName, "...")
dir.create(file.path(hubName, freeze), recursive=TRUE)
stopifnot(all(file.copy(sprintf("%s.bb", sampleInfo$bedfile), file.path(hubName, freeze))))
stopifnot(all(file.copy(sprintf("%s.bw", sampleInfo$bedfile), file.path(hubName, freeze))))

## bigWig, for displaying coverage
trackDbFile <- file.path(hubName, freeze, "trackDb_bw.txt")
trackDb <- data.frame(
    track = sampleInfo$sample,
    bigDataUrl = basename(paste0(sampleInfo$bedfile, ".bw")),
    shortLabel = sampleInfo$sample,
    longLabel = paste(sampleInfo$sample, sampleInfo$notes),
    type = "bigWig",
    color = sampleInfo$color,
    visibility = "full",
    windowingFunction = "maximum",
    configurable = "on" )
trackDb.list <- lapply(1:nrow(trackDb), function(i){ 
    c(paste(colnames(trackDb), trackDb[i,]), "")
} )
message("TrackDb for bw files: ", trackDbFile)
write(unlist(trackDb.list), file=trackDbFile)

## bigBed, for displaying individual reads
trackDbFile <- file.path(hubName, freeze, "trackDb_bb.txt")
trackDb <- data.frame(
    track = sampleInfo$sample,
    bigDataUrl = basename(paste0(sampleInfo$bedfile, ".bb")),
    shortLabel = sampleInfo$sample,
    longLabel = paste(sampleInfo$sample, sampleInfo$notes),
    type = "bigBed",
    color = sampleInfo$color,
    visibility = "pack",
    windowingFunction = "maximum",
    configurable = "on" ,
    maxHeightPixels = "100:24:8" )
trackDb.list <- lapply(1:nrow(trackDb), function(i){ 
    c(paste(colnames(trackDb), trackDb[i,]), "")
} )
message("TrackDb for bb files: ", trackDbFile)
write(unlist(trackDb.list), file=trackDbFile)


#### generate genomes.txt ####
## Use bb for reads
## use bw for coverage
message("\ngenomes files that point to track files")
message(file.path(hubName, "genomes_bb.txt"))
write(sprintf("genome %s \ntrackDb %s/trackDb_bb.txt", freeze, freeze),
      file=file.path(hubName, "genomes_bb.txt") )
message(file.path(hubName, "genomes_bw.txt"))
write(sprintf("genome %s \ntrackDb %s/trackDb_bw.txt", freeze, freeze),
      file=file.path(hubName, "genomes_bw.txt") )

#### generate hub.txt ####
message("\nUCSC hub entry files")
message("bigBed display: ", file.path(hubName, "hub_bb.txt"))
write(sprintf("hub %s_hub \nshortLabel %s \nlongLabel %s \ngenomesFile genomes_bb.txt \nemail yinghua@mail.med.upenn.edu",
        hubName, hubName, hubName),
      file=file.path(hubName, "hub_bb.txt") )
message("bigWig display: ", file.path(hubName, "hub_bw.txt"))
write(sprintf("hub %s_hub \nshortLabel %s \nlongLabel %s \ngenomesFile genomes_bw.txt \nemail yinghua@mail.med.upenn.edu",
              hubName, hubName, hubName),
      file=file.path(hubName, "hub_bw.txt") )
message("\nHub ", hubName, " generated, host the folder somewhere and point to the hub file with\n", 
        sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&hubClear=[URLOFhub.txt]", freeze),
        "\nNote the difference between &hubUrl= &hubClear=")

#### host hub ###
hostOnMicrob32 <- function() {
    message("\nHosting on microb32...")
    hublink <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=chr17%%3A1-78774742&hubClear=http://bushmanlab.org/ucsc/ywu/%s", freeze, file.path(hubName, c("hub_bb.txt", "hub_bw.txt")))
    write(hublink, file=file.path(hubName, "hub_links.txt"))
    
    hubHost <- "microb32.med.upenn.edu:~/pub"
    cmd <- sprintf("scp -r %s %s", hubName, hubHost)
    message(cmd)
    system(cmd)
    message("\nUCSC links are saved in ", file.path(hubName, "hub_links.txt"))
    message(paste(hublink, collapse="\n"))
}

hostOnMicrob215 <- function() {
    message("\nHosting on microb215...")
    hublink <- sprintf("http://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=chr17%%3A1-78774742&hubClear=https://microb215.med.upenn.edu/Download/data/share/UCSCTracksGTKJKOOKJHGTF/%s", freeze, file.path(hubName, c("hub_bb.txt", "hub_bw.txt")))
    write(hublink, file=file.path(hubName, "hub_links.txt"))
    
    hubHost <- "microb215.med.upenn.edu:~/Sites/share/UCSCTracksGTKJKOOKJHGTF"
    cmd <- sprintf("scp -r %s %s", hubName, hubHost)
    message(cmd)
    system(cmd)
    message("\nUCSC links are saved in ", file.path(hubName, "hub_links.txt"))
    message(paste(hublink, collapse="\n"))
}

##hostOnMicrob32()
hostOnMicrob215()

#### clean up ####
message("\nClean up...")
null <- file.remove(paste0(sampleInfo$bedfile, ".sort.bed"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bedGraph"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bw"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bb"))
null <- file.remove(paste0(freeze, ".genome.sizes"))

message("\nUse private/incognito mode to test the links.")


