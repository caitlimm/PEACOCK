#!/usr/local/bin/perl -w
use strict;

#generatefeatures.pl

#ARGV[0] is out, e.g., v1featuresHCT116
#ARGV[1] is the file of new links to be featured e.g., total_cislinks if starting from scratch
#ARGV[2] is Colon_Transverse_eQTL_Zstats.txt
#ARGV[3] is highH3K27ac_HCT116_links
#ARGV[4] is H3K4me1_HCT116_links
#ARGV[5] is H3K4me3_HCT116_genes
#ARGV[6] is closest__genes
#ARGV[7] is yes/no for p300
#ARGV[8] is intronic_enh
#ARGV[9] is pantherGeneList2020.txt
#ARGV[10] is pre_eQTL_Colon_Transverse_links_coeff_abs
#ARGV[11] is H3K27ac_HCT116_genes
#ARGV[12] is EP300_HepG2 if applicable

#To make the data features file:
#Enhancer	gene	f1	f2	f3	f4	f5	f6	f7	f8	f9	f10

#F1: H3K27ac (enhancer): 0/1
#F2: H3K4me1 (enhancer): score
#F3: H3K4me3 (promoter): 0/1
#F4: eQTL: combined Z score Colon Transverse
#F5: eQTL coeff absolute value CS
#F6: nearest gene? 0/1
#F7: intronic?
#F9: H3K27ac genes
#F10: p300 enhancers

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);

my $start_time = [Time::HiRes::gettimeofday()];

use POSIX qw(ceil floor);

my ($enhancer, $gene, $score, $cell, $link, $mark, %one, %two, %three, %seven, $cownt, $zstat, %five, %all, $count, @array, $panth, $ensg, %panther, $enh, %intronic , %links, $coeff, %coeff, %absCS, %new, $sum, $names, @names, %p300, %fiveCS, $entry);

my $p300=$ARGV[7];

open OUT, "> $ARGV[0]" or die $!; 

open (GEN,'<',$ARGV[11]) or die $!; 
#HUMAN|HGNC=24517|UniProtKB=Q9Y3T9       3.29504 H3K27ac HCT116
while (my $line=<GEN>){
	chomp $line;
	($gene, $score)=((split /\t/, $line)[0,1]);
	$new{$gene}{$score}=1;
}
close GEN;

if ($p300 eq "yes"){
	open (ABS,'<',$ARGV[12]) or die $!; 
	#41      175.40901       p300    HepG2
	while (my $line=<ABS>){
		chomp $line;
		($enh, $coeff)=((split /\t/, $line)[0,1]);
		$p300{$enh}=$coeff;
	}
close ABS;
}

open (NEW,'<',$ARGV[10]) or die $!; 
#EH37E0109680	HUMAN|HGNC=12781|UniProtKB=Q93097	30	3	-0.0324175
while (my $line=<NEW>){
	chomp $line;
	($enh, $gene, $coeff)=((split /\t/, $line)[0,1,4]);
	$link="$enh\t$gene";
	$absCS{$link}=$coeff;
}
close NEW;

open (PAN,'<',$ARGV[9]) or die $!;
#HUMAN|HGNC=33906|UniProtKB=Q075Z2       BSPH1,ENSG00000188334   Binder of sperm protein homolog 1;BSPH1;ortholog        BINDER OF SPERM PROTEIN HOMOLOG 1 (PTHR22918:SF12)    Homo sapiens
while (my $line=<PAN>){
	chomp $line;
	($panth, $names)=((split /\t/, $line)[0,1]);
	@names=(split /\,/, $names);
	foreach my $gene (@names){
		$panther{$gene}=$panth;
	}
}
close PAN;

open (INT,'<',$ARGV[8]) or die $!; 
#chr1    100104444       100104686       2       chr1    100104415       100104761       ENSG00000007952.13
while (my $line=<INT>){
	chomp $line;
	($enh, $ensg)=((split /\t/, $line)[3,7]);
	$ensg=((split /\./, $ensg)[0]);
	$intronic{$enh}=1;
	if ($ensg ~~ %panther){
		$gene=$panther{$ensg};
		$link="$enh\t$gene";
		$links{$link}=1;
	}
}
close INT;

open (SIX,'<',$ARGV[1]) or die $!;
#1       HUMAN|Gene=LPPR4|UniProtKB=Q7Z2D5
while (my $line=<SIX>){
	chomp $line;
	($enhancer, $gene)=((split /\t/, $line)[0,1]);
	$link="$enhancer\t$gene";
	$all{$link}=1;
}
close SIX;

$count=0;
open (FOUR,'<',$ARGV[2]) or die $!;
#enhancer        gene    tissue  k       chistat pchisq  Zstat
#28137   HUMAN|Gene=LPPR2|UniProtKB=Q96GM1       23      3       7.59517258439882        0.269286755388118       1.12510537168738
while (my $line=<FOUR>){
	chomp $line;
	$count++;
	next if ($count==1);
	($enhancer, $gene, $zstat)=((split /\t/, $line)[0,1,6]);
	$link="$enhancer\t$gene";
	$fiveCS{$link}=$zstat;
}
close FOUR;

open (ONE,'<',$ARGV[3]) or die $!;
#3       9_state_high       H3K27ac HepG2
while (my $line=<ONE>){
	chomp $line;
	($enhancer, $score, $mark, $cell)=(split /\t/, $line);
	$one{$enhancer}=1;
}
close ONE;

open (TWO,'<',$ARGV[4]) or die $!;
#2       4.92004 H3K4me1 HepG2
while (my $line=<TWO>){
	chomp $line;
	($enhancer, $score, $mark, $cell)=(split /\t/, $line);
	$two{$enhancer}{$score}=1;
}
close TWO;

$cownt=0;
foreach my $keys (sort keys %two){
	foreach my $keyz (sort keys %{$two{$keys}}){
		$cownt++;
	}
}
#print scalar keys %two;
#print " $cownt\n";

open (THREE,'<',$ARGV[5]) or die $!;
#HUMAN|Ensembl=ENSG00000005189|UniProtKB=H3BM72  9_state_high    H3K4me3 HepG2
while (my $line=<THREE>){
	chomp $line;
	($gene, $score)=((split /\t/, $line)[0,1]);
	$three{$gene}=1;
}
close THREE;

open (SEV,'<',$ARGV[6]) or die $!;
#1       HUMAN|HGNC=15846|UniProtKB=Q9NP74       closest_genes   HepG2
while (my $line=<SEV>){
	chomp $line;
	($enhancer, $gene)=((split /\t/, $line)[0,1]);
	$link="$enhancer\t$gene";
	$seven{$link}=1; 
}
close SEV;

$count=0;

if ($p300 eq "yes"){
	print OUT "enhancer\tgene\tH3K27ac_enhancer\tH3K4me1\tH3K4me3\teQTL_Zscore\tnearest\tintronic\tcoeff_abs\tH3K27ac_promoter\tp300\n";
}
elsif ($p300 eq "no"){
	print OUT "enhancer\tgene\tH3K27ac_enhancer\tH3K4me1\tH3K4me3\teQTL_Zscore\tnearest\tintronic\tcoeff_abs\tH3K27ac_promoter\n";
}
else{
	print "UHOH!";
}

foreach my $links (sort keys %all){
	print OUT "$links\t";
	($enhancer, $gene)=(split /\t/, $links);
	if ($enhancer ~~ %one){ #H3K7ac enhancer yes/no high
		print OUT "1\t";
	}
	else {
		print OUT "0\t";
	}
	if ($enhancer ~~ %two){ #H3K4me1
		$count=0;
		foreach my $entry (keys %{$two{$enhancer}}){
			$count++;
			if ($count==1){
			print OUT "$entry\t";
			}
		}
	}
	else {
		print OUT "0\t";
	}
	if ($gene ~~ %three){ #H3K4me3
		print OUT "1\t";
	}
	else{
		print OUT "0\t";
	}
	if ($links ~~ %fiveCS){ #eqtl z score
		$entry=$fiveCS{$links};
		chop $entry;
		print OUT "$entry\t";
	}
	else {
		print OUT "0\t";
	}
	if ($links ~~ %seven){ #nearest
		print OUT "1\t";
	}
	else {
		print OUT "0\t";
	}
	if ($links ~~ %links){
		print OUT "1\t"; #the enhancer is intronic to its target gene
	}
	else{
		if ($enhancer ~~ %intronic){
			print OUT "2\t"; #the enhancer is intronic to a different gene
		}
		else{
			print OUT "0\t"; #the enhancer is not intronic
		}
	}
	if ($links ~~ %absCS){	#eqtl coeff abs value
		$coeff=$absCS{$links};
		print OUT "$coeff\t"; 
	}
	else{
		print OUT "0\t";
	}
	if ($gene ~~ %new){	#H3K27ac gene value
		$count=0;
		$sum=0;
		foreach my $key (keys %{$new{$gene}}){
			$count++;
			$sum=($sum+$key);
		}
		$coeff=($sum/$count);
		print OUT "$coeff\t"; #the averaged gene value is added to the end of the line
	}
	else{
		print OUT "0\t"; #the gene value given is zero
	}
	if ($p300 eq "yes"){	#p300
		if ($enhancer ~~ %p300){
			$coeff=$p300{$enhancer};
			print OUT "\t$coeff\n";
		}
		else {
			print OUT "\t0\n";
		}
	}
	elsif ($p300 eq "no"){
		print OUT "\n";
	}
	else{
		print "PROBLEM";
	}
}

close OUT;

my $diff = Time::HiRes::tv_interval($start_time);

print "\n$diff\n";