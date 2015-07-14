## UCSC track doc 
## http://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html

## tools needed
## UCSC tools 
## wget -e robots=off -r -nH --cut-dirs=2 --no-parent --reject="index.html*" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
## bedtools http://bedtools.readthedocs.org/en/latest/content/installation.html
## https://github.com/arq5x/bedtools2/releases
message("checking required external programs, mysql, UCSC tools, bedtools...")
stopifnot( system("which mysql")==0 )
stopifnot( system("which ~/opt/UCSC/bedToBigBed")==0 )
stopifnot( system("which ~/opt/bedtools/genomeCoverageBed")==0 )
stopifnot( system("which ~/opt/UCSC/bedGraphToBigWig")==0 )

options(stringsAsFactors=FALSE)

#### input ####
#' @param csvFile, csv file, see below
#' @param hubName, a name for the hub
#' @param freeze, reference genome, hg18, hg19 etc
csvFile <- "bed.csv"
hubName <- "GENE"
freeze <- "hg18"

#### processing files ####
## the csv should look something like this, bedfile should point to the BED files.
##sample    bedfile    notes
##GTSP0001 GTSP0001.bed sample0001
##GTSP0002 GTSP0002.bed sample0002
##GTSP0003 GTSP0003.bed sample0003


sampleInfo <- read.csv(csvFile)
colnames(sampleInfo) <- tolower(colnames(sampleInfo))
stopifnot(file.exists(sampleInfo$bedfile))
sampleInfo <- subset(sampleInfo, file.exists(bedfile))

message("\nProcessing bed files...")
write.table(sampleInfo, "", quote=FALSE, sep="\t", row.names=FALSE)

## get $freeze.genome.sizes from UCSC
cmd <- sprintf("mysql --user=genome --host=genome-mysql.cse.ucsc.edu -sNA -e 'select chrom, size from %s.chromInfo' > %s.genome.sizes", freeze, freeze) 
message("\n",cmd)
stopifnot(system(cmd)==0)
stopifnot(file.exists(paste0(freeze, ".genome.sizes")))

## sort bed files
message("\nSort bed files...")
cmd <- sprintf("sort -S1G -k1,1 -k2,2n %s > %s.sort.bed", 
               sampleInfo$bedfile, 
               sampleInfo$bedfile)
null <- sapply(cmd, function(x) {message(x); system(x)} )
stopifnot(all(null==0))

## convert sorted bed file to bigbed
message("\nconvert bed files to bigBed files...")
cmd <- sprintf("~/opt/UCSC/bedToBigBed %s.sort.bed %s.genome.sizes %s.bb", 
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
null <- sapply(cmd, function(x) {message(x); system(x)} )
stopifnot(all(null==0))

## convert sorted bed file to bedGraph file
message("\nconvert bed files to bigGraph files...")
cmd <- sprintf("~/opt/bedtools/genomeCoverageBed -bg -i %s.sort.bed -g %s.genome.sizes > %s.bedGraph",
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
null <- sapply(cmd, function(x) {message(x); system(x)} )
stopifnot(all(null==0))

## convert bedGraph file to bigwig
message("\nconvert bigGraph files to bigWig files")
cmd <- sprintf("~/opt/UCSC/bedGraphToBigWig %s.bedGraph %s.genome.sizes %s.bw",
               sampleInfo$bedfile, 
               freeze,  
               sampleInfo$bedfile)
null <- sapply(cmd, function(x) {message(x); system(x)} )
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
    bigDataUrl = paste0(sampleInfo$bedfile, ".bw"),
    shortLabel = sampleInfo$sample,
    longLabel = paste(sampleInfo$sample, sampleInfo$notes),
    type = "bigWig",
    color = c("235,127,127", "127,127,235")[1:nrow(sampleInfo) %% 2 + 1],
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
    bigDataUrl = paste0(sampleInfo$bedfile, ".bb"),
    shortLabel = sampleInfo$sample,
    longLabel = paste(sampleInfo$sample, sampleInfo$notes),
    type = "bigBed",
    color = c("235,127,127", "127,127,235")[1:nrow(sampleInfo) %% 2 + 1],
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
hostOnMicrob98 <- function() {
    message("\nHosting on microb98...")
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

hostOnMicrob215()

#### clean up ####
message("\nClean up...")
null <- file.remove(paste0(sampleInfo$bedfile, ".sort.bed"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bedGraph"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bw"))
null <- file.remove(paste0(sampleInfo$bedfile, ".bb"))
null <- file.remove(paste0(freeze, ".genome.sizes"))

message("\nUse private/incognito mode to test the links.")


#### sample track configuration file ####
## [yinghua@microb237 makeUCSChubFromBed]$ tree CYS
## CYS
## ├── genomes_bb.txt
## ├── genomes_bw.txt
## ├── hg18
## │── ├── GTSP0689.bed.bb
## │── ├── GTSP0689.bed.bw
## │── ├── GTSP0691.bed.bb
## │── ├── GTSP0691.bed.bw
## │── ├── GTSP0692.bed.bb
## │── ├── GTSP0692.bed.bw
## │── ├── trackDb_bb.txt
## │── └── trackDb_bw.txt
## ├── hub_bb.txt
## ├── hub_bw.txt
## └── hub_links.txt


