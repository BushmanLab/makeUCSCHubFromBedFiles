for gtsp in `cut -d, -f1 CAR.csv | grep -iv sample`; do
    echo $gtsp;
    Rscript GTSP2BED.R $gtsp
done

Rscript ~/makeUCSCHubFromBedFiles/makeUCSChubFromBedFiles.R CAR.csv

