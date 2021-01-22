### This pipeline explains how to generate the dataset ascribing machine learning features of interest (to me) to all cis enhancer-gene links involving PEREGRINE enhancers and PANTHER genes. 

Here are the features:  
F1: H3K27ac (enhancer): 0/1  
F2: H3K4me1 (enhancer): score  
F3: H3K4me3 (promoter): 0/1  
F4: p300 binding (enhancer): score  
F5: eQTL: combined Z score  
F6: eQTL averaged coefficient (absolute values of original coefficients)  
F7: H3K27ac (promoters): score  
F8: nearest gene? 0/1  
F9: intronic? 0=no, 1=yes, target gene 2=yes, different gene  

### Feature 1: H3K27ac for enhancers
Download the "high H3K27ac for [cell line]" bed files in hg19 build, e,g., https://www.encodeproject.org/annotations/ENCSR370SPA/ for HepG2.
Use this Perl script to make a file with only the elements that had high peaks:
```
$ perl highstate.pl ENCFF386PZF.bed highstateout enhancer
```
See which PEREGRINE enhancers meet the 50% overlap threshold for these peaks:
```
$ bedtools intersect -wa -wb -f 0.5 -F 0.5 -e -a PEREGRINEenhancers -b highstateout > int
$ perl intersectenhancers.pl int highH3K27ac_HepG2_links H3K27ac HepG2 yes
```
The output file has the following format:
```
$ head highH3K27ac_HepG2_links
3       .       H3K27ac HepG2
12      .       H3K27ac HepG2
13      .       H3K27ac HepG2
41      .       H3K27ac HepG2
```

### Feature 2: H3K4me1 for enhancers
Download the “replicated peaks” bed file from ENCODE ChIP-seq experiments targeting H3K4me1 in hg19 build.  
Use a 50% overlap threshold to see which of the PEREGRINE enhancers overlap with these peaks:
```
$ bedtools intersect -wa -wb -f 0.5 -F 0.5 -e -a PEREGRINEenhancers -b ENCFF428FAW.bed > H3K4me1_HepG2_int
```
Use this perl script to output the enhancers that meet the 50% threshold for overlap with the H3K27ac peaks.
```
$ perl intersectenhancers.pl H3K4me1_HepG2_int H3K4me1_HepG2_links H3K4me1 HepG2 no
```
The output file format is:
```
$ head H3K4me1_HepG2_links
2       4.92004 H3K4me1 HepG2
3       5.38432 H3K4me1 HepG2
8       6.21267 H3K4me1 HepG2
13      7.09743 H3K4me1 HepG2
```

### Feature 3: H3K4me3 for promoters
Download the "high H3K4me3 for [cell line]" bed files in hg19 build, e,g., https://www.encodeproject.org/annotations/ENCSR368ISU/ for HepG2.
Use this Perl script to make a file with only the elements that had high peaks:
```
$ perl highstate.pl ENCFF579BAM.bed highstateout promoter
```
See which genes meet the 50% overlap threshold for these peaks:
```
$ bedtools intersect -wa -wb -f 0.5 -F 0.5 -e -a first500bpgenes -b highstateout > int
$ perl H3K4me3binding.pl int H3K4me3_HepG2_genes H3K4me3 pantherGeneList2020.txt HepG2
```
This is the output file format:
```
$ head H3K4me3_HepG2_genes
HUMAN|Ensembl=ENSG00000005189|UniProtKB=H3BM72  9_state_high    H3K4me3 HepG2
HUMAN|Ensembl=ENSG00000111780|UniProtKB=H0YIV9  9_state_high    H3K4me3 HepG2
HUMAN|Ensembl=ENSG00000141979|UniProtKB=B8ZZF3  9_state_high    H3K4me3 HepG2
HUMAN|Ensembl=ENSG00000144785|UniProtKB=F8W031  9_state_high    H3K4me3 HepG2
HUMAN|Ensembl=ENSG00000179240|UniProtKB=A0A0U1RQW1      9_state_high    H3K4me3 HepG
```
### Feature 4: p300 binding of enhancers
Download the bed file with “optimal” thresholded peaks in hg19 build from ENCODE's ChIP-seq experiments targeting p300. (e.g., https://www.encodeproject.org/experiments/ENCSR000BLW/)
The 7th column is the level of enrichment for the peak. 8th column is -log10(pval). 9th column is -log10(qval) for FDR. Use the enrichment peak amount.  
Find the level of overlap between these peaks and the PEREGRINE enhancers.
```
$ bedtools intersect -wa -wb -f 0.5 -F 0.5 -e -a PEREGRINEenhancers -b ENCFF547POU.bed > EP300_1
$ perl intersectenhancers.pl EP300_1 EP300_HepG2 p300 HepG2 no
```
The output file has this format:
```
$ head EP300_HepG2
41      175.40901       p300    HepG2
54      66.27693        p300    HepG2
106     33.32324        p300    HepG2
```
### Feature 5: eQTL combined Z score
Download the eQTL all_pairs file from GTEx (https://gtexportal.org/home/datasets) for the tissue of interest (click on GTEx Analysis V7 on the lefthand menu and select the relevant tissue under the Tissue-Specific All SNP Gene Associations section) and use this perl script to separate it by chromosomes (this is a speed/memory tip):
```
$ perl splitv4.pl Artery_Coronary.allpairs.txt
```
Then run the following script in R (but make sure to change the file names to match whatever tissue you're working on and the directory to wherever the files are located).
```
eqtl.R
```
Then execute this file (but make sure to change the file names to match whatever tissue you're working on and the directory to wherever the files are located):
```
$ ./all_eqtl_bash.sh
```
which contains this for each chromosome:
```
perl eqtlprocessv2.pl newdata1.txt Cells_EBV-transformed_lymphocytes_eqtl_all
perl pantherIDweed.pl Cells_EBV-transformed_lymphocytes_eqtl_all pantherGeneList2020.txt Cells_EBV-transformed_lymphocytes_eqtls leftovers
perl eliminateoverlapeqtlv2.pl Cells_EBV-transformed_lymphocytes_eqtls exons_genes.txt Cells_EBV-transformed_lymphocytes_eqtl2 overlaps unmatched
#rm newdata1.txt
bedtools intersect -wa -wb -a PEREGRINEenhancers -b Cells_EBV-transformed_lymphocytes_eqtl2 > Cells_EBV-transformed_lymphocytes_intersect
perl eqtllinks.pl Cells_EBV-transformed_lymphocytes_intersect links_chr1Cells_EBV-transformed_lymphocytes_eqtl tissuetable_03302020 pantherGeneList2020.txt
perl linksDBeqtl.pl links_chr1Cells_EBV-transformed_lymphocytes_eqtl linksDBeqtlchr1Cells_EBV-transformed_lymphocytes
```
At the end, the links_chr()Liver_eqtl files are concatenated into file Liver_links_eqtl.  
The linksDBeqtlchr()Liver files are concatenated into file Liver_links_eqtlnum.  
```
cat links_chr*Breast_Mammary_Tissue_eqtl > Breast_links_eqtl
cat linksDBeqtlchr*Breast > Breast_links_eqtlnum
```
Then run the following perl script.
```
$ perl eqtl.pl Liver_links_eqtl Liver_links_eqtlnum pre_eQTL_Liver_links
```
Then run eQTL_Zstats.R to get the final output file, Liver_eQTL_Zstats.txt.
```
$ head Liver_eQTL_Zstats.txt
enhancer        gene    tissue  k       chistat pchisq  Zstat
EH37E0330081    HUMAN|HGNC=575|UniProtKB=Q9Y587 30      1       2.81633570234595        0.244591        1.18666248410109
EH37E0437000    HUMAN|HGNC=7809|UniProtKB=P08138        30      1       1.77906392514693        0.410848        0.943150021244481
EH37E0144960    HUMAN|HGNC=4693|UniProtKB=Q16774        30      2       1.25198213681103        0.869468172454407       0.55946003807489
EH37E0594856    HUMAN|HGNC=15791|UniProtKB=Q9H490       30      1       2.65776484772833        0.264773        1.15277162693404
```
### Feature 6: eQTL averaged coefficient (absolute values of original coefficients)
Download the eQTL all_pairs file from GTEx (https://gtexportal.org/home/datasets) for the tissue of interest (click on GTEx Analysis V7 on the lefthand menu and select the relevant tissue under the Tissue-Specific All SNP Gene Associations aectioin) and use this perl script to separate it by chromosomes (this is a speed/memory tip):
```
$ perl splitv4.pl Artery_Coronary.allpairs.txt
```
Then run the following script in R (but make sure to change the file names to match whatever tissue you're working on and the directory to wherever the files are located).
```
eqtl.R
```
Then execute this file (but make sure to change the file names to match whatever tissue you're working on and the directory to wherever the files are located):
```
$ ./all_eqtl_bash.sh
```
which contains this for each chromosome:
```
perl eqtlprocessv2.pl Cells_EBV-transformed_lymphocytescoeffnewdata1.txt Cells_EBV-transformed_lymphocytes_eqtl_all 
perl pantherIDweed.pl Cells_EBV-transformed_lymphocytes_eqtl_all pantherGeneList2020.txt Cells_EBV-transformed_lymphocytes_eqtls leftovers 
perl eliminateoverlapeqtlv2.pl Cells_EBV-transformed_lymphocytes_eqtls exons_genes.txt Cells_EBV-transformed_lymphocytes_eqtl2 overlaps unmatched 
bedtools intersect -wa -wb -a PEREGRINEenhancers -b Cells_EBV-transformed_lymphocytes_eqtl2 > Cells_EBV-transformed_lymphocytes_intersect
perl eqtllinks.pl Cells_EBV-transformed_lymphocytes_intersect links_chr1Cells_EBV-transformed_lymphocytes_eqtlcoeff tissuetable_03302020 pantherGeneList2020.txt 
```
At the end, the links_chr()Breast_Mammary_Tissue_eqtlcoeff files are concatenated into file Breast_links_eqtlcoeff.  
```
cat links_chr*Breast_Mammary_Tissue_eqtlcoeff > Breast_links_eqtlcoeff
```
Then run the following perl script.
```
$ perl eqtlcoeff.pl Liver_links_eqtlcoeff Liver_links_eqtlnum pre_eQTL_Liver_links_coeff_abs yes
```
The output file has the format:
```
$ head pre_eQTL_Liver_links_coeff_abs
EH37E0617026    HUMAN|HGNC=1911|UniProtKB=Q13112        30      3       0.0466321523333333
EH37E0489983    HUMAN|HGNC=18753|UniProtKB=Q8NI29       30      2       0.4264625
EH37E0086258    HUMAN|HGNC=11732|UniProtKB=Q96S53       30      1       0.0273263
```
### Feature 7: H3K27ac for promoters
Download the “replicated peaks” bed file in hg19 build from ENCODE ChIP-seq experiments targeting H3K27ac (e.g., https://www.encodeproject.org/experiments/ENCSR000AMO/).  
The file promoters is the promoters of the genes.
```
$ bedtools intersect -wa -wb -f 0.5 -F 0.5 -e -a promoters -b ENCFF825NBZ.bed > H3K27acgenes_HepG2_int
$ perl p300bindinggenes.pl H3K27acgenes_HepG2_int H3K27ac_HepG2_genes H3K27ac pantherGeneList2020.txt HepG2
```
The output file format is:
```
$ head H3K27ac_HepG2_genes
HUMAN|HGNC=24517|UniProtKB=Q9Y3T9       3.29504 H3K27ac HepG2
HUMAN|HGNC=24023|UniProtKB=Q6TDP4       14.20247        H3K27ac HepG2
HUMAN|HGNC=24023|UniProtKB=Q6TDP4       3.29504 H3K27ac HepG2
```
### Feature 8: to find the closest gene to each enhancer
(This has already been generated; you can just take the outfile closest__genes but here is the process by which it was generated.)
```
$ sort -k1,1 -k2,2n PEREGRINEenhancers > sortedCRE
$ sort -k1,1 -k2,2n genes.bed > sortedgenes
$ bedtools closest -a sortedCRE -b sortedgenes > closestgenes
$ perl closestgenes.pl closestgenes closest__genes pantherGeneList2020.txt genes.txt
```
The output file format is:
```
$ head closest__genes
1       HUMAN|HGNC=15846|UniProtKB=Q9NP74       closest_genes   HepG2
10      HUMAN|HGNC=40020|UniProtKB=E7ERA6       closest_genes   HepG2
1000    HUMAN|HGNC=29408|UniProtKB=B1AJZ9       closest_genes   HepG2
10000   HUMAN|HGNC=19622|UniProtKB=Q8NG78       closest_genes   HepG2
10001   HUMAN|HGNC=8477|UniProtKB=Q15620        closest_genes   HepG2
```
### Feature 9: intronic? 0=no, 1=yes, target gene 2=yes, different gene
(This has already been generated; you can just take the outfile intronic_enh but here is the process by which it was generated.)  
First we need to organize the exons file exons_genes.txt so that we have all the exons grouped together by gene and in order:
```
$ perl exons.pl exons_genes.txt sorted_exons
```
Now we can pull the introns:
```
$ perl introns.pl sorted_exons intronz
$ sort intronz| uniq > introns
```
Now we can use bedtools intersect to get the enhancers that are intronic.
```
$ bedtools intersect -wa -wb -f 0.33 -a PEREGRINEenhancers -b introns > intronic_enh
```
The output file has the format:
```
$ head intronic_enh
chr1    100104444       100104686       2       chr1    100104415       100104761       ENSG00000007952.13
chr1    100113534       100114036       3       chr1    100106314       100117160       ENSG00000007952.13
chr1    100169951       100170218       7       chr1    100166438       100172334       ENSG00000036530.4
chr1    100169951       100170218       7       chr1    100166438       100172984       ENSG00000036530.4
chr1    100422351       100422596       8       chr1    100404050       100441812       ENSG00000152034.6
```
### To make a new features dataset:
This script takes all of the above generated files and creates a features dataset that is ready to be loaded intoR for analysis. The top example is for a cell line which did not have p300 ChIP-seq data and the bottom example is of a cell line that did have p300 ChIP-seq data.
```
$ perl generatefeatures.pl v1featuresHCT116 total_cislinks Colon_Transverse_eQTL_Zstats highH3K27ac_HCT116_links H3K4me1_HCT116_links H3K4me3_HCT116_genes closest__genes no intronic_enh pantherGeneList2020.txt pre_eQTL_Colon_Transverse_links_coeff_abs H3K27ac_HCT116_genes
$ perl generatefeatures.pl v1featuresK562 total_cislinks Whole_Blood_eQTL_Zstats.txt highH3K27ac_K562_links H3K4me1_K562_links H3K4me3_K562_genes closest__genes yes intronic_enh pantherGeneList2020.txt pre_eQTL_Whole_Blood_links_coeff_abs H3K27ac_K562_genes EP300_K562
```
