#Created by Pablo Cruz-Morales, sept 2022
#A pipeline to grab, label, align and trim sets of homologous proteins using a query and calculate a phylogeny 
#it runs along with fungison.pl using the same database
#the genomes database must be in /fungison/bin/genomes
#the GENOMES.Ids file must be in /fungison/bin/

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

print "USAGE: perl treehugger.pl <OPTIONS>\n\n";
print "OPTIONS:\n\n";
print "-q FILE.query   	|QUERY FILE, [a file with .query extension}\n";
print "-e 0.0000001		|E-VALUE CUTOFF, [a number]\n";
print "-s 200	        	|BIT-SCORE CUTOFF [a number]\n";
print "-x n or -F FORMATDB	|FORMAT THE DATABASE  ['no' is the recommeded option or 'FORMATDB']\n";
print "-a yes/no		|ALIGNS THE HIT SEQUENCES WITH MUSCLE, [yes/no]\n";
print "-t yes/no		|TRIMS THE HIT SEQUENCES WITH MUSCLE, [yes/no]\n";
print "-p yes/no		|CALCULATES A PHYLOGENETIC TREE [yes/no]\n\n";

my $query; my $evalue; my $score; my $database;
my $formatting; my @fastas; my $cont; my $dbname; my @blastouts;
my @columns; my $entry; my $line; my $original; my $genomename;
my $names; my $peg; my $id; my $hits; my $alignment; my $trimming; my $phylogeny;

GetOptions(
'q=s' => \$query,
'e=s' => \$evalue,
's=s' => \$score,
'x=s' => \$formatting,
'a=s' => \$alignment,
't=s' => \$trimming,
'p=s' => \$phylogeny,

) or die "Missing parameters\n";

	if ($query=~/.+query/) {
	print "QUERY= $query\n";
	}
	else {
	die "Missing argument -q FILE.query [QUERY FILE, a file with .query extension]\n\n"; 
	}
	if ($evalue) {
	print "EVALUE= $evalue\n";
	}
	else {
	die "Missing argument -e 0.0000001 [E-VALUE CUTOFF, a number]\n\n"; 
	}
	if ($score) {
	print "BITSCORE=$score\n";
	}
	else {
	die "Missing argument -s 200 [BIT-SCORE CUTOFF a number]\n\n";
	}

	#formatting database
	if ($formatting eq "FORMATDB"){
		@fastas=qx'ls ./bin/GENOMES/*.faa';
		print "The full database will be formatted\n";
		foreach (@fastas){
		    chomp $_;
		    system "makeblastdb -in $_ -dbtype prot -parse_seqids -logfile makeblastdb.log";
		    $cont++;
		    print "$cont $_ has been added to the blast database\n";
		}
	}
	elsif ($formatting=~/no/) {
	print "\nYou selected the current database, no DB formatting is requiered\n";
	}
	else {
	die "Missing argument -x n or -x FORMATDB [FORMAT THE DATABASE SELECTED WITH THE -d OPTION, 'no' or 'FORMATDB]\n\n";
	}
	
	#blasting on the database
	print "The blast search will be performed with evalue $evalue and score cutoff of $score\n";
	$query=~s/\.query//;
	open OUT, ">$query.hitlist";
	@fastas=qx'ls ./bin/GENOMES/*.faa';
	foreach (@fastas){
		chomp $_;
		$dbname= "$_";
		@blastouts= `blastp -query $query.query -db $_ -evalue $evalue -outfmt 6`;
		foreach (@blastouts){

			@columns=split(/\s+/,$_);
			if ($columns[11] >= "$score"){
				$entry="$columns[1]";
				#$entry=~s/fig://;
				print OUT "$entry\n";
				print "$entry\n";
			}
		}
	}

	#grabbing the non-redundant IDs
	system "sort $query.hitlist|uniq > $query.uniq";
	print "Grabbing the sequences with blastdbcmd\n";
	foreach(@fastas){
		chomp $_;
		$dbname= "$_";
		system "blastdbcmd -db $dbname -dbtype prot -entry_batch $query.uniq  -logfile blastdbcmd.log >> $query.txt";
	}

	open FASTA, ">$query.fasta";
	open FILE, "$query.txt" or die "Raw sequences file $query.txt not found\n";
	while ($line=<FILE>){
	if ($line=~/>fig:/){
        	$line=~/(>|fig\:)(\d+\.\d+)(\.peg\.\d+)/;
        	$original="$line";
        	$genomename="$2";
        	$peg="$3";
        	$peg=~s/\.//g;
	        open GENOMESID, "./bin/GENOMES.IDs" or die "I cant find the GENOMES.IDs file\n";
        	while ($id=<GENOMESID>){
              	  if ($id=~/$genomename/){
              		  $id=~/(.+)\t(.+)\t(.+)\t(.+)/;
              		  $names="$3";
              		  $names=~s/\-//g;
              		  $names=~s/\;//g;
              		  $names=~s/\(//g;
              		  $names=~s/\)//g;
              		  $names=~s/\'//g;
              		  $names=~s/\=//g;
              		  $names=~s/\.1//g;
              		  $names=~s/\.2//g;
              		  $names=~s/ATCC /ATCC/g;
              		  $names=~s/NRRL /NRRL/g;
              		  $names=~s/NBRC /NBRC/g;
              		  $names=~s/DSM /DSM/g;
              		  $names=~s/subsp //g;
              		  $names=~s/strain //g;
              		  $names=~s/\s+/\_/g;
              		  print FASTA ">$names\_$peg\n";
              	  }
        }
}
else {print FASTA $line;}
}

#cleaning and moving reporting so far
$hits=`grep ">" -c $query.fasta`;
chomp $hits;
print "\nfound $hits hit sequences, they are in the $query.fasta file\n\n";

#aligning
	if ($alignment) {
		if ($alignment=~/yes/) {
			print "Aligning the sequences in $query.fasta\n";
			system "muscle -in $query.fasta -out $query.aln -quiet";
			print "The alignments are in $query.aln\n";
		}
		if ($alignment=~/no/) {
			print "You selected not aligning the sequences\n";
			print "All done, have a great day\n";	
		}
	}
	else {
		print "The sequences will not be aligned\n";
		print "All done, have a great day\n";
		die "Missing argument -a yes/no [yes or no]\n\n";
	}

#trimming
#trimal is needed, get it here: https://github.com/inab/trimal

	if ($trimming) {
		if ($trimming=~/yes/) {
			print "Trimming the sequences in $query.aln\n";
			system "trimal -in $query.aln -automated1 -out $query.trimmed";
			print "The trimmed aligment is in $query.trimmed\n";
		}
		if ($trimming=~/no/) {
			print "You selected not trimming the sequences\n";
			print "All done, have a great day\n";	
		}
	}
	else {
		print "The sequences will not be trimmed\n";
		print "All done, have a great day\n";		
		die "Missing argument -t yes/no [yes or no]\n\n";
	}

#making a phylogeny
	if ($phylogeny) {
		if ($phylogeny=~/yes/) {
			print "A phylogenetic tree will be constructed using $query.trimed\n";
			system "iqtree -s $query.trimmed -m TEST -alrt 10000 -bb 10000 -quiet";
			print "The phylogenetic tree is in $query.trimmed.contre\n";
		}
		if ($phylogeny=~/no/) {
			print "You selected not constructing a phylogeny\n";
			print "All done, have a great day\n";	
		}
	}
	else {
		print "A phylogeny will not be constructed\n";
		print "All done, have a great day\n";		
		die "Missing argument -p yes/no [yes or no]\n\n";
	}

system "rm *.hitlist blastdbcmd.log";
system "rm *.uniq *.perf  $query.txt";
system "rm *.bionj *.ckp.gz *.iqtree *.mldist *.model.gz *.splits.nex *.treefile";
system "mkdir $query\_Treehugger_results";
system "mv *.contree ./$query\_Treehugger_results";
system "mv *.trimmed ./$query\_Treehugger_results";
system "mv *.aln ./$query\_Treehugger_results";
system "mv *.fasta ./$query\_Treehugger_results";
system "mv *.log ./$query\_Treehugger_results";
print "\nAll outputs are in the $query\_Treehugger_results folder\n";

print "All done, have a great day\n\n";
