# treehugger.pl
A pipeline to grab, label, align and trim sets of homologous proteins using a query and calculate a phylogeny <br />
it runs along with fungison.pl using the same database

# Installation
Just download the script and place it in the same folder fungison.pl and the /bin directory  <br /><br />

# Dependencies:
blast+<br />
      `sudo apt install blast+`<br />
muscle<br />
      `sudo apt install muscle`<br />
trimal<br />
```
wget https://github.com/inab/trimal/archive/refs/heads/trimAl.zip
unzip trimAl.zip
cd trimAl
cd source
make
iqtree
sudo apt install muscle
```

The /bin/GENOMES folder is requiered <br />
The /bin/GENOMES.IDs file should be available <br />

# Usage
```
USAGE: perl treehugger.pl <OPTIONS>
OPTIONS:\n
-q FILE.query   	|QUERY FILE, [a file with .query extension}
-e 0.0000001		|E-VALUE CUTOFF, [a number]
-s 200	        	|BIT-SCORE CUTOFF [a number]
-x n or -F FORMATDB	|FORMAT THE DATABASE  ['no' is the recommeded option or 'FORMATDB']
-a yes/no		|ALIGNS THE HIT SEQUENCES WITH MUSCLE, [yes/no]
-t yes/no		|TRIMS THE HIT SEQUENCES WITH MUSCLE, [yes/no]
-p yes/no		|CALCULATES A PHYLOGENETIC TREE [yes/no]

Examples:
perl treehugger.pl  -q gluc.query  -e 0.00001 -s 100 -x no -a yes -t yes -p yes
perl treehugger.pl  -q gluc.query  -e 0.00001 -s 100 -x FORMATDB -a yes -t yes -p yes
```
