### makeUCSCHubFromBedFiles
This script does the following things to make the bed files into a costomer hub for the UCSC genome browser.

1. sort bed files,
2. convert sorted bed files to bigbed format,
3. calculate coverage from the sorted bed files to get bedGraph files,
4. convert bedGraph files to bigWig files.

Both ```bigbed``` and ```bigwig``` can be displayed on the genome browser with the distinction:
1. ```bigbed``` is more for showing individual reads,
2. ```bigWig``` is more towards the coverage.

In order to make switch of display easy, two hubs are generated:
1. hub_bb.txt ---> genomes_bb.txt ---> trackDb_bb.txt ---> the bigbed files,
2. hub_bw.txt ---> genomes_bw.txt ---> trackDb_bw.txt ---> the bigbed files,
so that one only need to change ```bb``` to ```bw``` to switch the display.

### code example
The script is easy to follow. Just adapt it to your data.

### host
One needs a public host that the genome browser can pull data from, either http and https is fine. But, the files can not be password protected. If your apache server allows override the following ```.htaccess``` file enables directory listing and public access.
```
Options +Indexes
IndexOptions IgnoreCase FancyIndexing FoldersFirst NameWidth=* DescriptionWidth=* SuppressHTMLPreamble
IndexIgnore header.html footer.html favicon.ico .htaccess *.php *.css .. .ftpquota .DS_Store icons *.log *,v *,t .??* *~ *#
HeaderName header.html

Order allow,deny
Allow from env=allow
Satisfy any
Allow from all
Require all granted
```

### link
The links can be given by 
```
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hubClear=[URLTOTHEHUB.txt]
http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&hubUrl=[URLTOTHEHUB.txt]
```
where ```hg18``` should be changed to the reference genome of the bed files. With
- ```hubClear``` only that hub is displayed,
- ```hubUrl``` that hub is added to existing costomer tracks.  

It is always a good idea to check the link from private/incognito mode. If you like a specific display, save it to a session and then show it to collaberators. Don't just save the link. 

### requirements
1. UCSC genome browser binary tools available from http://hgdownload.cse.ucsc.edu/admin/exe/
2. bedtools from https://github.com/arq5x/bedtools2/releases
3. mysql
To donwload all the UCSC bonaries, use the follow command:
```
wget -e robots=off -r -nH --cut-dirs=2 --no-parent --reject="index.html*" http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
```


