# Polygenic Risk Scores
Polygenic Risk Scores, or PRS, consists of taking whole-genome genotyping information and their association to a phenotype in a population to predict the individual-level genetic risk of having or developing the phenotype. Instead of looking at specific SNPs that passed genome-wide significance, PRS, like LDSC, make use of all SNPs in the genome for its calculation. There is a very nice review paper from [Choi, Mak, and O'Reilly, 2020, Nat. Protoc](https://doi.org/10.1038/s41596-020-0353-1) in which every single detail of PRS calculation is very well explained. Actually, they are the developers of [PRSice-2](10.1093/gigascience/giz082), one of the most widely used tools for PRS calculation. Briefly, PRS are calculated using base data, that is a GWAS of large sample size for a phenotype of interest, and a target data, that is genotype-level data from an independent, tipycally much smaller, cohort. The effect sizes derived from the base data (which come from the association testing against the phenotype of interest) are used to calculate the PRS in the target data. The outcome PRS model will be able to predict the phenotype in the target sample with varying degree of accuracy and precision. During the PRS calculation, a clumping + thesholding approach is followed. First, all common SNPs between base and target data are LD clumped to select only independent signals. Then, a different *P*-value threshold is applied to select the one that integrates the number of SNPs that better predicts the phenotype. Therefore, the best model fit will have ideally a significant *P*-value with a high R<sup>2</sup> statistic. Besides the two papers linked above, there are two very nice tutorials that one can follow to perform PRS in a controlled and reproducible way and ensuring hgih quality standards:

-  [Basic Tutorial for Polygenic Risk Score Analysis](https://choishingwan.github.io/PRS-Tutorial/)
-  [PRSice-2 Tutorial](https://choishingwan.github.io/PRSice/)

Before calculating the PRS, we need to make sure that our base and target data fullfil the quality standards, I will just name the QC steps:
## QC of base data
- **Heritability**. Your GWAS should have an h<sup>2</sup><sub>SNP</sub> > 0.05 (5% heritability). Otherwise we could reach misleading conclusions from the application of PRS. Check [LDSC](https://github.com/cdiazmun/MDAgwas/tree/main/ldsc) section for more info on how to calculate this. 
- **Effect alleles**. One should be sure of which of the alleles is the effect allele and the other allele. If you have produced the GWAS you're using as base data, then you should be fine. If you're downloading GWAS from GWASCatalog you should also be fine (99% of times I would say). For the latter it may be a good idea to homogenize the alleles anyway to the 1000G reference panel of the corresponding population/ancestry. Check [postGWAS-QC](https://github.com/cdiazmun/MDAgwas/tree/main/postGWAS-QC) section on how to perform this using GWASinspector. 
-  **Rare variants**. For PRS calculation, only common variants are tipycally employed. If you don't have the minimum allele frequency (MAF) column in your GWAS sumstats, calculate it from the effect allele frequency (EAF) and use it to filter out MAFs < 0.01.
``` r
library(dplyr)
library(data.table)

setwd("PAth/to/your/working/directory/")

# Read the GWAS meta-analysis, our base data
df <- fread("gwas-meta-analysis.txt.gz")

# Ensure we are only using common SNPs, so that the PRS are reliable
df1 <- df %>% 
  mutate(MAF = if_else(EFF_ALL_FREQ <= 0.5, EFF_ALL_FREQ, 1-EFF_ALL_FREQ)) %>% 
  filter(MAF >= 0.01)  %>%
  # I created a column named SNP, because we will use this column to find common markers between base and target data
  mutate(SNP = paste(CHR, POSITION, OTHER_ALL, EFFECT_ALL, sep = ":"))
```
-  **Duplicated SNPs**. Remove them although PRSice-2 will detect them and remove them anyway.
```r
# Take only SNPs with rsIDs --> In my case, I had variants without an rsID assigned to them, therefore these will be taken as duplicated, while they're not
df2 <- df1 %>% 
  filter(grepl("rs", ID))
# Check if there are duplicated values, by rsID and by SNP
df_all_dups <- df2[duplicated(df2$ID) | duplicated(df2$ID, fromLast = TRUE), ]
df_all_dups_SNP <- df1[duplicated(df1$SNP) | duplicated(df1$SNP, fromLast = TRUE), ]

# Remove duplicated SNPs 
df3 <- df1 %>% 
  filter(!(SNP %in% df_all_dups$SNP)) %>% 
  filter(!(SNP %in% df_all_dups_SNP$SNP))
```
-  **Ambiguous SNPs**. Remove them although PRSice-2 will detect them and remove them anyway. We should remove them because they may be strand-ambiguous and difficult to flip if there are inconsistencies between the base and the target (especially when the allele frequencies are close to 0.5).
``` r
df4 <- df3 %>% 
  filter(!(EFFECT_ALL == "A" & OTHER_ALL == "T")) %>% 
  filter(!(EFFECT_ALL == "T" & OTHER_ALL == "A")) %>% 
  filter(!(EFFECT_ALL == "G" & OTHER_ALL == "C")) %>% 
  filter(!(EFFECT_ALL == "C" & OTHER_ALL == "G"))
```
-  **Transform BETA to OR**. Do this if your phenotype is binary (cases/controls). We also rename and select only the necessary columns. See the tutorials to check these columns and their names.
``` r
base <- df4 %>% 
  mutate(OR = exp(EFFECT)) %>% 
  rename(BP = POSITION, A1 = EFFECT_ALL, A2 = OTHER_ALL, SE = STDERR, N = N_TOTAL, P = PVALUE) %>% 
  select(-EFF_ALL_FREQ, -EFFECT, -ID)

fwrite(base, "base_files/base.data.nondup.nonamb.txt.gz", 
       quote = F, row.names = F, sep = "\t", compress = "gzip")
```

## QC of target data
-  **Genotyping and sample QC**. If not done already, perform the QC that is normally executed on PLINK binary files (MAF, HWE, missing rates, *etc*). I will write more details about this in another section.
-  If necessary, transform PLINK2 files (pgen, psam, pvar) into **PLINK binary files** (bim, bed, fam).
``` bash
plink2  \
	--pfile genotype_data \ # prefix of pgen, psam, pvar files
	--make-bed \
  --out \ # genotype_data # prefix of bim, bed, fam files
```
-  **LiftOver**. If your target data is not in the same genome build as your base data, you should perform liftover. The chain files can be downloaded [here](https://hgdownload.soe.ucsc.edu/downloads.html#liftover). The R packages `GenomicRanges` and `rtracklayer` can be used to liftover your data. Check [postGWAS-QC](https://github.com/cdiazmun/MDAgwas/tree/main/postGWAS-QC) section on how to perform this.
``` r
# Read in bim file 
bim <- fread("genotype_data.bim") %>%
 rename(CHR = V1,
        SNP = V2,
        CM = V3,
        POSITION = V4,
        EFFECT_ALL = V5,
        OTHER_ALL = V6)

# Perform liftOver from 38 to 37
library(rtracklayer)
library(GenomicRanges)

chain <- import.chain("hg38ToHg19.over.chain")

bim38 <- bim %>% 
  mutate(CHR = paste0("chr", CHR))

bim38_ranges <- makeGRangesFromDataFrame(df = bim38,
                                         keep.extra.columns = T,
                                         ignore.strand = T,
                                         seqnames.field = "CHR",
                                         start.field = "POSITION",
                                         end.field = "POSITION")

bim37 <- liftOver(x = bim38_ranges, chain = chain) %>% 
  unlist() %>% 
  as.data.table()

efficiency <- nrow(bim37)/nrow(bim38)*100
print(paste0("Percentage of markers lift-overed: ", efficiency))

# Check if there are duplicated values --> 0, no duplicated variants
bim_dups <- bim37[duplicated(bim37$SNP), ]

# Generate the target data
target <- bim37[, -c(3:5)] %>% 
  dplyr::rename(CHR = seqnames, POSITION = start) %>% 
  mutate(CHR = as.numeric(sapply(strsplit(basename(as.character(CHR)), "chr"), `[`, 2))) %>% 
  filter(!(is.na(CHR))) %>% # These NA are because some variants are in the chrX or chrY when liftovered
  mutate(newSNP = paste(CHR, POSITION, OTHER_ALL, EFFECT_ALL, sep = ":"))
```
-  Check that most of the SNPs between base and target data are in **common** and that the alleles are **homogenized** (so no need to flip them).
``` r
# Merge summary statistic with target 
info <- df4 %>% 
  inner_join(target, by = c("SNP" = "newSNP"))

# Get SNPs that have the same alleles across base and target --> ALL, we don't need to flip alleles
info.match <- subset(info, EFFECT_ALL.x == EFFECT_ALL.y & OTHER_ALL.x == OTHER_ALL.y)
```
-  Save the chain file to update the PLINK binary files.
```r
write.table(target, "target_files/target_chain_file.txt", row.names = F, col.names = F, quote = F, sep = "\t")
```
-  Update PLINK binary files for GRCh7. Create three files with the first column as the old SNP names and the second column as (1) CHR, (2) POSITION, and (3) new SNP names. And a final file containing only the variants to keep (old SNP names to be able to select them from the original bim file). 
``` bash
awk '{print $3, $1}' target_chain_file.txt > updated_chr.txt
awk '{print $3, $2}' target_chain_file.txt > updated_pos.txt
awk '{print $3, $7}' target_chain_file.txt > updated_names.txt
awk '{print $1}' updated_names.txt > keep_vars.txt
plink \
	--bfile target_data \
	--extract keep_vars.txt \
	--make-bed \
	--out tmp
plink \
	--bfile tmp \
	--update-chr updated_chr.txt 2 1 \
	--update-map updated_pos.txt 2 1 \
	--make-bed \
	--out tmp_chr_pos
plink \
	--bfile tmp_chr_pos \
	--update-name  updated_names.txt 2 1 \
	--make-bed \
	--out target_data_37
```

## PRSice-2
Download the executable from [here](https://choishingwan.github.io/PRSice/).

``` bash
unzip PRSice_linux.zip

Rscript ~/opt/prsice-2/PRSice.R \ # Soft links do not work
	--prsice ~/opt/prsice-2/PRSice_linux \
	--base base_files/base_data.txt.gz \
	--target target_files/target_data_37 \
	--binary-target T \ # Binary target set to TRUE
	--pheno target_files/pheno_file.txt \
	--cov target_files/covariates_file.txt \ # It must include covariates like age and sex and the PCs all together in a single file
	--stat OR \ # Tell the name of the column for the genetic effects
	--or \ # Tell that we're dealing with Odds Ratio
	--out phenotype \ # The name of the output file
	--thread 12 # Don't know if necessary, it's really fast
```
For some reason, it found two (2) duplicated SNPs in the target data, so it selected them so we can re-run without them:
``` bash
Rscript ~/opt/prsice-2/PRSice.R \ 
	--prsice ~/opt/prsice-2/PRSice_linux \
	--base base_files/base_data.txt.gz \
	--target target_files/target_data_37 \
	--binary-target T \ 
	--pheno target_files/pheno_file.txt \
	--cov target_files/covariates_file.txt \ 
	--stat OR \ 
	--or \ 
	--out phenotype \ 
	--thread 12 \
	--extract phenotype.valid # This is the file that automatically generates PRSice
```
It was expecting that the --pheno file contained FIDs, but we can tell the script to ignore this
``` bash
Rscript ~/opt/prsice-2/PRSice.R \ 
	--prsice ~/opt/prsice-2/PRSice_linux \
	--base base_files/base_data.txt.gz \
	--target target_files/target_data_37 \
	--binary-target T \ 
	--pheno target_files/pheno_file.txt \
	--cov target_files/covariates_file.txt \ 
	--stat OR \ 
	--or \ 
	--out phenotype \ 
	--thread 12 \
	--extract phenotype.valid \
	--ignore-fid
```

## Per-sample PRS
The best PRS model will be produced with a subset of the GWAS SNPs, which are the ones that help the most to predict the risk of the disease. We can then use that set of SNPs to calculate the polygenic score of an individual for which we have the genotype available. If we have the genotypes of hundreds of thousando of individuals we can calculate the score for all of them and stratify them according to their genetic risk to suffer the disease under study.

We are going to work with a case in which the PRS for a Phenotype has already been calculated in an independent cohort (not included in our GWAS meta-analysis). That PRS model is based on ~50K variants. We need to assign a score for the presence/absence of each of the risk alleles for this set of variants in each individual from the UK Biobank (UKB). This score will be weighted by the beta, so the genetic effects of each risk allele. Besides stratifying the individuals by their associated PRS, we can also use this per-sample PRS to test for association with ICD10 codes to see if individuales with a higher score there is also higher prevalence of other conditions. This kind of analysis tell us how genetically similar are different conditions. It's another layer on top of genetic correlation or Mendelian randomization analyses. 

This pipeline is based on plink2 and some data handling with R. We need the following input files:
-	`weights`: A data frame with the ~50K variants obtained from the PRS and their betas.
-	`genotypes`: We will use imputed data from the UKB in `.pgen/.psam/.pvar` format, split per chromosome. 

### Check input files
The first thing we need to do is to inspect our `weights` file, which will look something like this:

|CHR	|BP	|A1	|A2	|SNP	|ID	|BETA	
|-------|-------|-------|-------|-------|-------|-------
|1	|100066515	|T	|C	|1:100066515:C:T	|rs1234567	|0.005516434
|1	|100097419	|C	|T	|1:100097419:T:C	|rs1234567	|0.004740695
|1	|100116253	|T	|G	|1:100116253:G:T	|rs1234567	|0.003353760
|1	|100177526	|A	|G	|1:100177526:G:A	|rs1234567	|0.011832565
|1	|100196224	|T	|C	|1:100196224:C:T	|rs1234567	|0.015933962

In this case, `A1` will be the risk allele. We will need the `ID` column to extract the variants from the `genotype` files and the `SNP`column for the last step when calculating the scores. We need to check if the genome build of the `weights` file and the `genotype` files is the same and how are the variants reported in the `genotype` files. As I'm saying, in this case it's reported with the rsID, that's why we use the `ID` column. 

### Extract the PRS variants
Use the `weights` file to extract only the rsIDs, without headers, and save it to `rsids.prs.txt`. Then run:
``` bash
for i in {1..22};
	do plink2 --bgen /path/to/Genotypes/Imputation_field-22828/ukb22828_c${i}_b0_v3.bgen \
	ref-first \
	--sample /path/to/Genotypes/Imputation_field-22828/ukb22828_c${i}_b0_v3_s487203.sample \
	--extract /path/to/Weights/rsids.prs.txt \
	--make-pgen \
	--out ukb_prs_genotypes_per_chromosome/ukb.prs.genotypes.ref1.c${i};
done
```
### Merge the files
If there are some rsIDs that correspond to multiallelic variants, the merging will fail because it will encounter more than one variant for the same rsID. To avoid this happening we have different options but the easiest is to simply reset all variants to `CHR:POSITION:REF:ALT`: 
``` bash
for i in {1..22};
	do plink2 --pfile ukb_prs_genotypes_per_chromosome/ukb.prs.genotypes.ref1.c$i \
	--set-all-var-ids @:#:\$r:\$a \
	--make-pgen \
	--out ukb_prs_genotypes_per_chromosome_setvars/ukb.prs.genotypes.ref1.setvars.c$i;
done
```
Then we create a list with all files we want to merge:
``` bash
ll ukb_prs_genotypes_per_chromosome_setvars/*.pgen | cut -d " " -f 10 | cut -d "." -f 1,2,3,4,5,6 > merge_list_setvars.txt
```
And we finally merge the files:
``` bash
plink2 --pmerge-list merge_list_setvars.txt --make-pgen --out ukb.prs.genotypes.ref1.setvars.allchr
```
### Calculate the polygenic scores
We are going to calculate the score with `plink2` default parameters. This means we are going to use the **_additive model_** to calculate the allelic risk effects (as done with most GWAS). For what I could read, `plink2` would check the risk allele (A1) in the .pvar UKB genotype file andflip the alelle if necessary to calculate the score. The `ukb.prs.genotypes.ref1.setvars.allchr` generated should look like this:

|#CHROM  |POS     |ID      |REF     |ALT
|--------|--------|--------|--------|---
|1       |835499  |1:835499:A:G    |A       |G
|1       |890104  |1:890104:G:A    |G       |A
|1       |902069  |1:902069:T:C    |T       |C
|1       |1038796 |1:1038796:G:A   |G       |A
|1       |1219941 |1:1219941:G:A   |G       |A

And the `weights` file to be used in this step should only contain three columns (be aware we changed `SNP` by `ID` in the column name to match the `.pvar` file):

|ID      |A1      |BETA
|--------|--------|----
|1:100066515:C:T |T       |0.00551643416434031
|1:100097419:T:C |C       |0.00474069531077363
|1:100116253:G:T |T       |0.00335376003369533
|1:100177526:G:A |A       |0.0118325647205171
|1:100196224:C:T |T       |0.0159339621510669

Once we have both files like this, we can calculate the scores:
``` bash
plink2 --pfile ukb.prs.genotypes.ref1.setvars.allchr --score SNP_A1_BETA_PRS_WEIGHTS.txt 1 2 3 header --out prs.ukb.all.samples
```
The output should look like this:

|#FID    |IID     |ALLELE_CT       |NAMED_ALLELE_DOSAGE_SUM |SCORE1_AVG
|--------|--------|----------------|------------------------|----------
|xxx     |xxx     |110500  |22454.789       |4.84725e-05
|xxx     |xxx     |110500  |23448.194       |5.28227e-05
|xxx     |xxx     |110500  |23588.847       |6.82782e-05
|xxx     |xxx     |110500  |23704.081       |7.14413e-05
