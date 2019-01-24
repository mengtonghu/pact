10-10-18
Outlines
1.First c
#PCAT
##1.Download dataset from `https://www.med.unc.edu/pgc/results-and-downloads`
Chose ANGST(Anxiety Neuro Genetics Study)
Only need `SNPID`, `CHR`, `Effect`,`P.value`. I can use `awk` function to extract columns before using R
The positive and negative sign of the effect can tell us the sign of z score
negative effect -- negative z score
postive effect -- positive z score
##2.convert Pvalue to z score
```R
setwd("/Users/mengtonghu/Desktop/Fall2018/Honors/angst")
angst<-read.csv("anxiety.meta.full.cc.tbl",sep="\t")
angst$z<-qnorm(1 - (angst$P.value/2)) 
#one tail or two tail?? one tail due to null hypothesis
angst$sign <- 1
angst$sign[angst$Effect < 0] <- -1
angst$z<-angst$z*angst$sign
```
##3. 
3.1 Use Bystro to locate SNPS on genes using 1000 genome data
Bystro can import data from AWS S3, having trouble creating basket.
10-16-18
Jingjing gave me access to 1000 genome data in her lab on hgcc. I tried to download it and upload it on to Bystro using pc in the library.
I did not understand the results.
3.2 
I followed the instruction on this website `http://www.gettinggeneticsdone.com/2011/06/mapping-snps-to-genes-for-gwas.html`
Table Browser:
`http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=694099781_ujY6T7YeDrJi9dhRfMn4uPmLhBNK&clade=mammal&org=Human&db=hg38&hgta_group=varRep&hgta_track=hg38Patch11&hgta_table=0&hgta_regionType=genome&position=1%3A11102837-11267747&hgta_outputType=primaryTable&hgta_outFileName=`

10-23-18
Installing bedtools
First install homebrew:
```shell
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" < /dev/null 2> /dev/null
```
install bget:
```shell
brew install wget
```
follow instruction on bedtools website to install bedtools
downloaded the latest version from github
```shell
make
```
compare genes overlaping rate using BEDtools
```
intersectBed -a SNPs.bed -b entrezgenes.bed -wa -wb | awk '{printf("%s\t%s\n", $4, $8);}'
```
We decided to use HUGO.txt 
To get the correlation matrix, use plink r. First filter out snps in ped for gene 1, use Plink to get the correlation matrix for gene1. Repeat this for all genes. Try write a bash script to do this.
When start to play around, try chromosome 22 which has the least number of genes.
```plink
--r
```
10-31-18
1.5kb up stream
shell/ bash language 
I've ran the methods 
https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/
###step 1: mapping snps to hugo.txt
The code written in python only compares position information. For different chromosomes, I need to first extract both angst.table and hugo.txt for each chromosome
##R-codes for extracting chromosomes in 'anxiety.tables':
``` R
angst<-read.csv("anxiety.table",sep=" ")
angst$z<-qnorm(1 - (angst$P.value/2)) 
#one tail or two tail??
angst$sign <- 1
angst$sign[angst$Effect < 0] <- -1
angst$z<-angst$z*angst$sign
angst<-angst[,1:8]
df_list<- split(angst, as.factor(angst$CHR))
write.table(df_list$`1`,'angst1.table',row.names=FALSE,quote=FALSE)
write.table(df_list$`2`,'angst2.table',row.names=FALSE,quote=FALSE)
#repeat the above steps for all 22 chromosomes
```

##R-codes for extracting chromosomes in 'hugo.txt':
```R
hugo<-read.table("hugo.txt",sep="\t")
df_list22<- split(hugo, as.factor(hugo$V1))
write.table(df_list2$`1`,'hugo1.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`2`,'hugo2.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`3`,'hugo3.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`4`,'hugo4.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`5`,'hugo5.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`6`,'hugo6.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`7`,'hugo7.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`8`,'hugo8.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`9`,'hugo9.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`10`,'hugo10.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`11`,'hugo11.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`12`,'hugo12.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`13`,'hugo13.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`14`,'hugo14.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`15`,'hugo15.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`16`,'hugo16.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`17`,'hugo17.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`18`,'hugo18.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`19`,'hugo19.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`20`,'hugo20.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`21`,'hugo21.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
write.table(df_list2$`22`,'hugo22.txt',row.names=FALSE,col.names =FALSE,quote=FALSE,sep="\t")
```

##actual mapping
```python
import pandas as pd
data2 = pd.read_csv('angst22.table', sep=" ", header=0)
hugo=pd.read_csv('hugo22.txt',sep='\t',header=None)
hugo.columns = ["CHR", "LB", "UB", "GENE"]
newDF = pd.DataFrame()
data2[['BP']] = data2[['BP']].apply(pd.to_numeric)
import numpy as np


for row2 in hugo.iterrows():

    newDF = data2.loc[data2['BP'].apply(lambda x: x >= row2[1][1]-5000)]
    newDF = newDF.loc[newDF['BP'].apply(lambda x:x <=row2[1][2])].copy()
    newDF['id'] = row2[1][3]
    success = newDF.to_csv('/Users/mengtonghu/Desktop/Fall2018/Honors/matching/success22.csv', mode='a',header=False)            
```

###using 1000 genomes browser for LD matrix

`success22.csv` has an additional columns of DNA names
go to this website for the starting and ending position for each chromosome
https://www.ncbi.nlm.nih.gov/variation/tools/1000genomes/
on the search place, type in Chr22, and the results show that 
Chr22: 16256332-46639653
22:42,089,417..42,089,518
Use this website to finish the conversion from VCF to PED
http://grch37.ensembl.org/Homo_sapiens/Tools/VcftoPed
Name for this job:Chr22
Region Lookup:22:16256332-46639653
4:122868000-122946000
22:16256332-46639653

phase3? or phase 1
 population: CEU TSI FIN GBR IBS
 http://www.internationalgenome.org/category/population/

http://www.internationalgenome.org/category/tabix/

#Testing for gene POTEH
```R
succcess22[1:5,]
```
```shell
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:16260678-16286465 >chr22.vcf

plink --vcf chr22.vcf --recode --out chr22

```
`chr22.map`, `chr22.ped` contains genotypes fro 1000genome for snps that are in this region as gene POTEH
Because we only care about spns that are identified by `success22.txt`, I manually wrote a script `snps.txt` contains the snps name of the snps that we want to extract from 1000 genome

```shell
plink --file chr22 --recode --extract snps.txt --out chr22filter
plink --file chr22filter --r --out correlation
```

##Using P-ACT.R
http://csg.sph.umich.edu/boehnke/p_act/README.txt
http://csg.sph.umich.edu/boehnke/p_act/p_ACT_meta_0.1/p_act_meta.php

https://github.com/gkichaev/PAINTOR_V3.0/wiki/2a.-Computing-1000-genomes-LD
https://data.broadinstitute.org/mpg/snpsnap/database_download.html

According to the input file format, we want to rewrite the snps into 0,1,2 as allele counts. To do that in plink, we will need a recode.txt
To make the recode. txt, I first created a frequency file
```shell
plink --file chr22filter --freq --out chr22 

```
awk '{print $2,$4}' chr22.frq >chr22filter_recode.txt
```
plink --file chr22filter --recodeA --recode-allele chr22filter_recode.txt --out chr22filter_recode
awk '{print $7,$8,$9,$10,$11}' chr22filter_recode.raw >chr22filter_recode_transpose.txt
```
```R
library(tidyverse)
table<-read.table("chr22filter_recode_transpose.txt",sep=" ")
a<-str_split_fixed(t(table[1,]),"_",2)
colnames(table)<-a[,1]
table=table[-1,]
table_out<-t(table)
write.table(table_out,"geno1.csv",quote=FALSE,sep=",",col.names=FALSE)

```

use A2 as reference allele

##Download vcf tools
`https://github.com/vcftools/vcftools`
```
git clone https://github.com/vcftools/vcftools.git
cd vcftools
brew install autoconf
brew install automake
./autogen.sh
brew install pkg-config
```
do brew install pkg-config
and then run ./autogen.sh again
and then ./configure
```
make
sudo make install

mv vcftools /Users/mengtonghu
export PERL5LIB=/Users/mengtonghu/vcftools/src/perl/
```
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr1.phase1.projectConsensus.genotypes.vcf.gz 17:1471000-1472000 . | head -20 | vcf-subset -c HG00110 > HG00110.vcf
```
How to include only the individuals from european population (CEU,GBR,FIN,TSI,IBS)?
`http://www.internationalgenome.org/data-portal/sample`
Choose the population and download the list of data
```
 awk '{print $1}' /Users/mengtonghu/Desktop/Fall2018/Honors/matching/1000-genome-indiviudal-information.txt >europeanids.txt
```
`europeanids.txt` contains the European population
```
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr1.phase1.projectConsensus.genotypes.vcf.gz  22:16260678-16286465 . | head -20 | vcf-subset -c europeanids.txt> gene1.vcf
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:16260678-16286465 . | head -20 | vcf-subset -c europeanids.txt| bgzip -c> gene1.vcf
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:16260678-16286465 >chr22.vcf
```
Manually deleted all individuals listed below
No such column in the VCF file: "HG00107"
No such column in the VCF file: "HG00181"
No such column in the VCF file: "HG00132"
No such column in the VCF file: "HG00105"
No such column in the VCF file: "HG00115"
No such column in the VCF file: "HG00145"
No such column in the VCF file: "HG00157"
No such column in the VCF file: "HG00379"
No such column in the VCF file: "HG00371"
No such column in the VCF file: "HG00365"
No such column in the VCF file: "HG00380"
No such column in the VCF file: "HG00368"
NA06985
NA11840
NA12234
NA11881
NA12762
NA12005
NA11832
NA12813
NA12156
NA12828
NA12414
NA12878
NA12892
NA12891
NA12760
NA12776
NA20526
NA20514

```R
awk '{print $1,$1}' /Users/mengtonghu/Desktop/Fall2018/Honors/matching/1000-genome-indiviudal-information.txt >europeanids.txt
plink --vcf chr22.vcf --recode --keep europeanids.txt --extract snps_POTEH.txt  --out POTEH_euro
```
obtained the ped file `POTEH_euro.ped` and fam file `POTEH_euro.fam`

To run PACT
Script P_ACT.R
Input: space-delimited
1. `pvalues.txt:`
1 rs148271204 0.926
1 rs146675715 0.5352
1 rs139144787 0.5357
1 rs141944226 0.08782
1 rs150703810 0.1121
2. `genotype.txt`
    One column for each genotype score. 
    One row for each individual plus a header row of genotype labels.  
    Missing values must be coded as NA.
    Not necessary if only one marker is considered.

    plink --file POTEH_euro --recode A-transpose --recode-allele reference_alleles_POTECH.txt --out POTEH_recode

    paste <(cut -f1 -d " " reference_alleles_POTEH.txt) <(cut -f2 -d " "  reference_alleles_POTEH.txt | tr '[:lower:]' '[:upper:]')
    12-3
    I met with Karen today and I realize there is no need to consider about reference_alleles

1-8
We need a julia code that can get vcf file for each gene
```Julia
using DataFrames, CSV

```
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:16260678-16286465 >chr22.vcf

tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz 22:43083155-43117037 >A4GALT.vcf

22:43083155-43117037

##4.Creating both `pvalues.txt` and `genotyp.txt`
###4.1 Use julia to get pvals.txt
```julia
using DataFrames, CSV
location = CSV.read("gene_22.txt",delim=' ')
for i in 1:size(location)[1]
    bound = string("22:",location[:start][i],"-", location[:end][i])
    @show bound
    file= string(location[:GENE][i],".vcf")
     @show file
    #script= "tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr22.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz $bound > $file"
    #run(`sh -c $script`)
    file_euro=string(location[:GENE][i],"_euro")
   run(`plink --vcf $file --recode --keep europeanids.txt --maf 0.05 --out $file_euro`)
    file_transpose=string(location[:GENE][i],"_euro_trans")
   run(`plink --file $file_euro --recode A-transpose --out $file_transpose`)
    file_transpose_r=string(location[:GENE][i],"_euro_trans.traw")
    transpose = CSV.read(file_transpose_r,delim='\t')
    dfl = stack(transpose,7:size(transpose)[2])
    dfnew = unstack(dfl, :variable, :SNP, :value)
    out=dfnew[:,2:size(dfnew)[2]]
    output= string(location[:GENE][i],"_pact.txt")
    CSV.write(output,out, delim=' ')
end
```
obtain `_pact.txt` for each gene on chormosome 22

Error in mvt(lower = lower, upper = upper, df = 0, corr = corr, delta = mean,  : 
  only dimensions 1 <= n <= 1000 allowed
60:TOP3B
167:MTMR3
213:SYN3
215:LARGE
411:TBC1D22A
412:FAM19A5

After getting rid of the snps with minor allele frequency<0.05, the below error goes away
diag(.) had 0 or NA entries; non-finite result is doubtful[1] 95
 Show Traceback
Error in probval.GenzBretz(algorithm, n, df, lower, upper, infin, corr, : NA/NaN/Inf in foreign function call (arg 6)
159
217:HMGXB4
218:
230:MYH9
348:WBP2NL
350:NAGA
351:FAM109B
352:C22orf32
355:CYP2D7P1
356:TCF20

More than 1000 snps for some genes:
Check the gene:
Check with Claudia: using Julia
Check 0.05 minor allele frequency. If they don't help, don't need to do this step.
eleminate monorphic: 
R code- dubugging. -a browser-

error catching
error flag- not completely break

Complete `python hmt1.py` `python hmt2.py` on hgcc
First time using Python3.6
```shell
module load Anaconda3
conda create --clone base --name mhu
conda install -c conda-forge -python-levenshtein
```

at node00, to check the status of work, type in `qstat`
sumbit a job
```
qsub -q b.q -cwd -j y -M mengtong.hu@emory.edu /home/mhu/honors/hmt_job3.sh
```

one snp might match up with 2 genes

scp /Users/mengtonghu/Desktop/Honors/angst/*.py mhu@hgcc.genetics.emory.edu:/home/mhu/honors/

scp /Users/mengtonghu/Desktop/Honors/angst/*.sh mhu@hgcc.genetics.emory.edu:/home/mhu/honors/

scp mhu@hgcc.genetics.emory.edu:/home/mhu/honors/success_3.txt .
Hmt0818@$2017
```r
snps<-read.table("EG302_01_in_1MDuov3.table",header=FALSE)
snps1<-t(snps)
 write.table(snps1, file = "snps_header_EG302_01.txt", row.names = FALSE,quote = FALSE,  append = FALSE, col.names = FALSE, sep = ",")
```