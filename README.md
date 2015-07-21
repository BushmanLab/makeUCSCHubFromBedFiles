### makeUCSCHubFromBedFiles
This script does the following things to make the bed files into a customer hub for the UCSC genome browser.

1. sort bed files,
2. convert sorted bed files to bigbed format,
3. calculate coverage from the sorted bed files to get bedGraph files,
4. convert bedGraph files to bigWig files.

Both ```bigbed``` and ```bigwig``` can be displayed on the genome browser with the distinction:  
- ```bigbed``` is more for showing individual reads,  
- ```bigWig``` is more towards the coverage.

In order to make switch of display easy, two hubs are generated:
```
1. hub_bb.txt ---> genomes_bb.txt ---> trackDb_bb.txt ---> the bigbed files,
2. hub_bw.txt ---> genomes_bw.txt ---> trackDb_bw.txt ---> the bigwig files,
```
so that one only needs to change ```bb``` to ```bw``` to switch the display.

### Code example
Change the following variables to adapt to your data and then run the script.
```
Rscript makeUCSChubFromBedFiles.R bed.csv

A csv file should look like this
Sample,bedfile,Notes,freeze,hub
GTSP0001,GTSP0001.bed,mock,hg18,GENE
GTSP0002,GTSP0002.bed,real,hg18,GENE
GTSP0003,GTSP0003.bed,mock,hg18,GENE
GTSP0004,GTSP0004.bed,real,hg18,GENE
GTSP0005,GTSP0005.bed,real,hg18,GENE

# The tracks will be displayed by the order given in the csv file.
# Short labels are given by sample column.
# Long labels are given by sample notes columns.
# The colors alternate acording to the notes column.
# Freeze column should only have one value.
# Hub column should only have one value.
```

In order to generate `GTSP####.bed`, use `Rscript GTSP2BED.R GTSP####`. A sript to make the hub without the bed files could be:
```
for gtsp in `cut -d, -f1 bed.csv | grep -iv sample`; do
    echo $gtsp;
    Rscript GTSP2BED.R $gtsp
done
Rscript ~/makeUCSCHubFromBedFiles/makeUCSChubFromBedFiles.R CAR.csv
```

### Output example
```
GENE
├── genomes_bb.txt
├── genomes_bw.txt
├── hg18
│   ├── GTSP0001.bed.bb
│   ├── GTSP0001.bed.bw
│   ├── GTSP0002.bed.bb
│   ├── GTSP0002.bed.bw
│   ├── GTSP0003.bed.bb
│   ├── GTSP0003.bed.bw
│   ├── GTSP0004.bed.bb
│   ├── GTSP0004.bed.bw
│   ├── GTSP0005.bed.bb
│   ├── GTSP0005.bed.bw
│   ├── trackDb_bb.txt
│   └── trackDb_bw.txt
├── hub_bb.txt
├── hub_bw.txt
└── hub_links.txt
```

### Link example
Please open in private/incognito mode.
```
[yhwu GENE]$ cat hub_links.txt
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr17%3A1-78774742&hubClear=http://bushmanlab.org/ucsc/ywu/GENE/hub_bb.txt
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=chr17%3A1-78774742&hubClear=http://bushmanlab.org/ucsc/ywu/GENE/hub_bw.txt
```

### Host
One needs a public host that the genome browser can pull data from without password. Either http or https is fine. If your apache server allows override the following ```.htaccess``` file enables directory listing and public access.
```
Options +Indexes
IndexOptions IgnoreCase FancyIndexing FoldersFirst NameWidth=* DescriptionWidth=* SuppressHTMLPreamble
IndexIgnore header.html footer.html favicon.ico .htaccess *.php ..
HeaderName header.html

Order allow,deny
Allow from env=allow
Satisfy any
Allow from all
Require all granted
```

### Link to hub url
The links can be given by 
```
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hubClear=[URLTOTHEHUB.txt]
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hubUrl=[URLTOTHEHUB.txt]
```
where ```hg18``` should be changed to the reference genome for the bed files. With
- ```hubClear```, only this hub is displayed,
- ```hubUrl```, this hub is added to existing costomer tracks.  

It is always a good idea to check the link from private/incognito mode.  
If you like a specific display, save it to a session and link to it by the session's link. 

### Requirements
1. UCSC genome browser binary tools available from http://hgdownload.cse.ucsc.edu/admin/exe/
2. bedtools from https://github.com/arq5x/bedtools2/releases
3. mysql  

To donwload all the UCSC bonaries, use the follow command:
```
wget -e robots=off -r -nH --cut-dirs=2 --no-parent --reject="index.html*" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
```

