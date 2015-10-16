bedfile <- list.files(".", "*.bed")
Notes <- ifelse( grepl("old", bedfile), "old", "new")
Sample <- sub(".bed", "", bedfile)
freeze <- "hg18"
hub <- "vectorTrimBetaThal"

csv <- data.frame(Sample, bedfile, Notes, freeze, hub)
write.csv(csv, "bed.csv", row.names=FALSE, quote=FALSE)
system("cat bed.csv")

