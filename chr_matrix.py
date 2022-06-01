import numpy as np
import time
import csv
import numpy as np

gwas_chr_snps=[]

with open("p.1gwasSNPs_chr1.txt","r") as in_handle:                                       #CHANGE CHROMOSOME
    for line in in_handle:
        sline = line.strip().split("\t")
        gwas_chr_snps.append(sline)

num_snps = len(gwas_chr_snps) #number comes from quantity of common snps retrieved from new vcf
chr_snps = np.asarray(gwas_chr_snps, dtype=str)
num_samples = 9600 #total number of people in study

print("number of gwas snps p<.1: ", num_snps)

snp=0
matrix=[]
with open("/labs/mpsnyder/kainj/gvcfgenotyper.9600.2019-02-16.chr1.filt.norm.QC.vcf","r") as in_handle:                          #CHANGE CHROMOSOME
    for line in in_handle:
        if line.startswith("#"):
            continue #SKIP HEADER LINES
        sline = line.strip().split("\t")
        allele=sline[2]
        if allele not in chr_snps:
            continue #SKIP SNPS WITH P>0.1 IN ALS GWAS
        snp = snp + 1
        sample=1
        #print("SNP: ", snp, " allele: ", allele)
        x=allele+" "
        for info in sline[9:]: #sline[9:] is correct for starting at the first sample "0/0:....."
            genotype = info[0:3] #select first 3 characters from input value per sample for genotype
            if genotype in ["0/0", "./0"]:
                x=x+'0 ' #SET ALLELE COUNT TO 0 (SNP absent)
            elif genotype in ["0/1", "./1"]:
                x=x+'1 ' #SET ALLELE COUNT TO 1 (heterozygous)
            else:
                x=x+'2 ' #SET ALLELE COUNT TO 2 (homozygous)
        matrix.append(x)
        sample = sample + 1
print("done")
np.savetxt('gwas.1_matrix_chr1_filtnorm.txt', matrix, delimiter='', fmt="%s")                               #CHANGE CHROMOSOME
print("saved")
