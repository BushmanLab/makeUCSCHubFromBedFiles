library("dplyr")
library("R.utils")

options(dplyr.width = Inf)

cmd <- "Rscript ~/geneTherapyPatientReportMaker/check_patient_gtsp.R | grep bthal | grep pPa"

res <- system(cmd, intern=TRUE)
res <- plyr::ldply(strsplit(res, "\t"))
colnames(res) <- c("trial", "patient", "cell", "time", "gtsp", "replicate", "freeze", "sex")

cmd <- paste("Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromDB.R",
             paste(unique(res$gtsp), collapse=" "))
system(cmd)


process_sample <- function(trial, patient, cell, time, gtsp) {
    fileName <- unique(paste(trial, patient, cell, time, "bed", sep="."))
    cmd <- paste("cat", paste(paste0(unique(gtsp), ".bed"), collapse=" "),
                 ">",  fileName)
    system(cmd)
    return(fileName)
}

bedfiles <- (res  %>% group_by(trial, patient, cell, time) %>%
             summarize(fn=process_sample(trial, patient, cell, time, gtsp)))



csv <- (mutate(bedfiles,
               Sample=paste(trial, patient, cell, time, sep="."),
               bedfile=fn,
               Notes=Sample,
               freeze="hg18",
               hub="bthal") %>%
        ungroup %>%
        select(Sample,bedfile,Notes,freeze,hub) )

csv$nrow <-  sapply(csv$bedfile, countLines)

csv <- subset(csv, nrow != 0 )

write.csv(csv, file="bed.csv", row.names=FALSE, quote=FALSE)       

cmd <- "Rscript ~/makeUCSCHubFromBedFiles/makeUCSChubFromBedFiles.R bed.csv"
system(cmd)

