# Post-GWAS QC

Thistutorial drives you through the quality check (QC) of several independent genome-wide association study (GWAS) for the same phenotype. A real-world scenario would be that of working with many cohorts that supply us with the GWAS summary statistics files for the trait under study. There are several bioinformatic tools to perform a GWAS ([PLINK](10.1086/519795), [BOLT-LMM](10.1038/ng.3190), [SAIGE](https://doi.org/10.1038/s41588-018-0184-y), or [REGENIE](https://doi.org/10.1038/s41588-021-00870-7)) and even more specificities (*e.g.,* covariates used, type of summary statistics reported) so that the final file that we receive can vary a lot from cohort to cohort. There are many things to check and even more GWAS files to perform these checks in. If we jump over critical issues (lacking essential columns, wrong information, non-matching numbers, huge deviations in some stats, *etc*, ...) there are a few things that we need to solve *prior* to perform the actual QC. Although one can argue that the following steps are actually part of the quality control...

## Homogenization

When dealing with many cohorts, I consider important to have the less heterogenity as possible, starting with the file names. Having a standard file name for every GWAS file will make our life much easier and our workflow much faster and efficient. In the example below, I just read the raw GWAS summary statistics file and a table containing all the information about every GWAS (path, sex, ancestry, N, trait, *etc*), and I output the same file, renamed (automatically) and with the new column "N" added, if needed:

``` r
# Import the libraries we will need for specific functions
library(readxl) # for read_excel()
library(data.table) # for fread() and fwrite()
library(dplyr) # for data manipulation

trait1 <- "your_phenotype"

# Get the summary file
gwas.sum <- read_excel("path/to/the/Excel_file.xlsx") %>%
  filter(Trait == trait1)

# Get the file with the whole path of the files to rename
code1 <- "pmb"

file.names <- subset(gwas.sum, Code == code1, select = File_complete_win)
file.names <- file.names[[1]]

# Get the ancestries
ancestry <- subset(gwas.sum, Code == code1, select = Ancestry)
ancestry <- ancestry[[1]]

# Get the sex
sex <- subset(gwas.sum, Code == code1, select = Sex)
sex <- sex[[1]]

# Get the extension
compr <- gwas.sum$Compression[gwas.sum$Code == code1]
compr <- compr[1]


# Run if you ONLY need to rename
for (i in 1:length(file.names)){
  # Set the output path
  out.path <- paste0("/path/to/output/folder/", trait1, "/sumstats/", ancestry[i])
  # Set the output name
  if (compr == TRUE) {
    out.name <- paste("sumstat", code1, ancestry[i], sex[i], trait1, "txt", "gz", sep = ".")
  } else {
    out.name <- paste("sumstat", code1, ancestry[i], sex[i], trait1, "txt", sep = ".")
  }
  print(paste0("copying ", file.names[i], " to ", paste(out.path, out.name, sep = "/")))
  file.copy(from = file.names[i],
            to = paste(out.path, out.name, sep = "/"))
}


# Run if you need to add an extra column with the total N
n_total <- subset(gwas.sum, Code == code1, select = N_total)
n_total <- n_total[[1]]

for (i in 1:length(file.names)){
  # Set the output path
  out.path <- paste0("/path/to/output/folder/", trait1, "/sumstats/", ancestry[i])
  # Set the output name
  out.name <- paste("sumstat", code1, ancestry[i], sex[i], trait1, "txt", sep = ".")
  # Read the GWAS summary statistics
  # Doesn't matter if it's compressed or not, delimiter, it just reads it
  print(paste0("Processing ", file.names[i], "..."))
  print("Reading file...")
  df <- fread(file.names[i])
  print("Adding column...")
  df$N <- rep(n_total[i], nrow(df))
  print(paste0("Writing file to ", paste(out.path, out.name, sep = "/")))
  fwrite(df, paste(out.path, out.name, sep = "/"), col.names = TRUE, row.names=FALSE, quote=FALSE, sep = "\t")
  print("---------------------------------------------------")
}
```

Again, I just prompt the Rscript above in case you find the code useful for yourself. But to be able to use this script completely, you need to organize the data accordingly. 

## Allele frequency
We need to check how is the allele frequency reported. Sometimes is the "effect allele frequency (EAF)" and others the "minimum allele frequency (MAF)". This is not much of a problem because if it's the MAF, `GWASinspector` will flip or switch the variants accordingly to match the reference dataset. However, it's recommended that they report to us the EAF. And if they report both, use the EAF column (A1FREQ, EAF, AF_ALLELE1, *etc*). In any case, it is normally reported from 0 to 1. But there are some cases in which they can report this as a percentage. In that case we need to simply divide that column by 100. Just leave the Rscript to do it iteratively over many files, but it's just a simple thing to do. 

``` r
library(data.table)
library(dplyr)

filepath <- "/path/to/input/files/"

filenames <- list.files(full.names = TRUE,
                        path = filepath,
                        pattern = "*.txt.gz")

out.path <- "/path/to/output/folder/"

for (i in filenames) {
  # Set output name. Yes, I know it is redundant, but it doesn't consume time and it works so I don't touch it.
  split.name <- unlist(strsplit(basename(i), "[.]"))
  input.name <- paste(split.name[1:5], collapse = ".")
  print(paste0("Processing ", input.name, "..."))
  output.name <- paste(input.name, "txt.gz", sep = ".")
  # Read the GWAS summary statistics and transform the EAF column.
  # Doesn't matter if it's compressed or not, delimiter, it just reads it.
  df <- fread(i) %>% 
    mutate(EAF = EAF/100)
  # Write the output file
  fwrite(df, paste0(out.path, output.name), col.names = TRUE, row.names=FALSE, quote=FALSE, sep = "\t", compress = "gzip")
  # Print everything.
  print("------------------------------")
}
q()
```

## Column headers of the alleles

As mentioned before, there are a few softwares that can be used to perform a GWAS and each of them outputs a tab- or space-spearated file containing the summary statistics, but under different nomenclature. Later, we will use `GWASinspector` to standardize the common (and most important) columns but we need to avoid missinterpretation by this tool when reading GWAS from different softwares name equally different things. In this case, we decided to rename "Allele2" and "Allele1" from some GWAS produced with SAIGE to "ALT" and "REF", respectively:

``` r
zcat headers/sumstat.cohort.ancestry.sex.trait.txt.gz | sed '1s/Allele2/ALT/g; 1s/Allele1/REF/g' > sumstat.cohort.ancestry.sex.trait.txt
gzip sumstat.cohort.ancestry.sex.trait.txt
```

## Flipped variants
Sometimes it can happen that the alleles are shifted and the EAF is therefore reported for the "other allele" instead of the "effect allele". In this case we need to do two things to manually flip the variants:

Swap the name of the alleles:
``` bash
zcat flipped/sumstat.cohort.ancestry.sex.trait.txt.gz | sed '1s/ALLELE0/ALLELETEMP/g; 1s/ALLELE1/ALLELE0/g; 1s/ALLELETEMP/ALLELE1/g' > sumstat.cohort.ancestry.sex.trait.txt
gzip sumstat.cohort.ancestry.sex.trait.txt
```
Recalculate the effect allele frequency as 1 - the reported allele frequency:
``` r
# Import and export as shown in one of the previous scripts
# Mutate the column
df <- df %>% mutate(A1FREQ = 1-A1FREQ)
```

Actually, `GWASinspector` can do that automatically for every variant, but it only does so with those variants that it finds in the reference dataset (which varies from ~70% to 99.99% in most of the cases).


## LiftOver
What we also need to do is to check the genome build that has been used to perform the GWAS. In this pipeline, we will always use the human genome build GRCh37 because most of our summary statistics are done with that genome build. But we have one GWAS performed with the GRCh38 (the most up-to-date one). In order to classify appropriately every variant, we need to perform the so-called liftover. We can do that using a few `R-packages`, namely `GenomicRanges` and `rtracklayer`. You will need to check and adapt the script according to your own GWAS colnames and content:

``` r
library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicRanges)

# Import the GWAS sum stats file
gwas38 <- fread("sumstat.cohort.ancestry.sex.trait.txt.gz")

# Import the chain file, to be able to compare and swap between variants from one genome build to the other
# It can be found here: wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
chain <- import.chain("/path/to/hg38ToHg19.over.chain")
chain 

# Be sure chromosomes are coded same way in the sumstats, 
# In this case, for example we need to add "chr" at the beginning and chrX is coded as chr23, instead

# Reformat the chromosome column (watch out to use the appropriate column names)
gwas38[, CHR := paste0("chr", CHR)] %>% 
  .[, CHR := ifelse(CHR =="chr23","chrX",CHR)]

# Convert to GRanges object (watch out to use the appropriate column names)
gwas38_ranges <- makeGRangesFromDataFrame(df = gwas38
                                        , keep.extra.columns = T
                                        , ignore.strand = T
                                        , seqnames.field = "CHR"
                                        , start.field = "BP"
                                        , end.field = "BP")

# Perform the liftover
gwas37 <- liftOver(gwas38_ranges, chain) %>% unlist() %>% as.data.table()
# Check the percentage of SNPs that have been converted 
nrow(gwas37)/nrow(gwas38)*100 

# Modify a bit to accommodate to GWASinspector expected input
# Remove unused columns
gwas.final <- gwas37[,-c(3:5)] 
# Change column names
colnames(gwas.final)[1:2] <- c("CHR","BP")
# Go back to a numeric indicator for the chromosomes
gwas.final$CHR <- sapply(strsplit(as.character(gwas.final$CHR), "chr"), `[`, 2)
# Go back to 23 as the nomenclature for the chromosome X
gwas.final$CHR[gwas.final$CHR == "X"] <- 23

# write.file
fwrite(gwas.final, "sumstat.cohort.ancestry.sex.trait.txt.gz", sep = "\t")
```

If you're interested on an iterative way of running this script, you can do as in the examples before. However, for a reason I ignore (after much investiation), it can accumulate a lot of RAM usage when running iteratively, so I wouldn't do it for more than 3-4 files at a time. 

## GWASinspector
The main objective of this step is to make sure that we only include variants that have a high quality, in other words, variants that we can trust for any downstream analysis. As a summary, we will filter out the following variants:

- **Minimum Allele Frequency (MAF) < 0.01**. We will exclude rare variants and keep only common variants. Rare variants are difficult to discover in, especially, cohorts with a small sample size.
- **Imputation score (INFO) < 0.8**. To only include imputed variants that we can highly trust. If the summary statistics file does not have a column indicating the `INFO`, check if this filter was already applied upstream of the analysis.
- **Hardy-Weinberg Equilibrium (HWE) > 1e-6**. To only keep variants that are in HWE, as otherwise we would be including variants that are influenced by population stratification or other evolutionary forces that will only contribute to increase the risk of false positive discovery or add noise to the overall analysis. This is typically checked before performing the association analysis, so in practice this is irrelevant here, but included in the tool we will use.
- **Call rate > 0.95**. To include variants that are present in at least 95% of the individuals in the cohort. This step is typically performed before the association test, so in practice this is irrelevant here, but included in the tool we will use.
- **Difference in Allele Frequency (AF) > 0.2**. We will exclude variants whose AF is very different from the ones reported in the reference dataset. This could entail issues with population stratification that alter the expected AF for a given variant.
- **Multiallelic variants**. We will exclude multiallelic variants because they are usually characterized by being of bad quality and poorly annotated.

We will also check for:
- **Population stratification** - QQ plots
- **Allele frequency differences** - Query*vs*Reference plots
- **Standard Error of SNP effect sizes** - SE*vs*Nmax plots
- **Distribution of AF** - Skewness*vs*Kurtosis plots

There are a few R packages that we can use for this matter:
1.  [GWASinspector](https://doi.org/10.1093/bioinformatics/btaa1084) - [download](http://gwasinspector.com/)
2.  [EasyQC](https://doi.org/10.1089/cmb.2017.0186) - [download](https://bio.tools/easyqc#!)
3.  [QCGWAS](https://doi.org/10.1093/bioinformatics/btt745) - [download](https://cran.r-project.org/web/packages/QCGWAS/index.html)

We will be using the first option. You can access the [vignettes](https://cran.r-project.org/web/packages/GWASinspector/vignettes/GWASinspector.html) and the [manual](http://gwasinspector.com/references/Introduction_to_GWASinspector.pdf). Most of the scripts you'll see below are coming from there. I will detail everything here anyway, so that we are all in the same page. 

### Prepare the files 
The first thing we need to do is to install the R package (`install.packages("GWASinspector")`). Then we will start getting the necessary input files. Create a general folder as well as a different folder for the reference files.
``` bash
mkdir gwas_inspector gwas_inspector/dir_references
```
Download the genotype reference files. There are a few options here to use as reference. Since we are only using populations of European ancestry, we will directly download that one. Change if necessary. The files need to be decompressed to be used, in the `dir_references` folder.
``` bash
wget http://www.gwasinspector.com/references/1000GENOMES-p_3_EUR.sqlite.gz
gunzip 1000GENOMES-p_3_EUR.sqlite.gz
```
Prepare the environment and files to run GWASinspector.
``` r
library(GWASinspector)
system_check()
get_headerTranslation('~/gwas_inspector/dir_references')
get_config('~/gwas_inspector')
```
With the initial R script you prepare the additional input files that serve as parameter configuration file (`config.ini`) and a translator of the thousands of different ways of putting a column name for the same thing (`alt_headers.txt`). Change the `config.ini` file to set the proper directories and input file locations and adapt a few parameters. The file is too large to print it all here, but the most impotant things to be considered are the following:
```
# Folder containing the input file(s) (do not use a trailing slash)
dir_data = /complete/path/to/sum_stats

# Folder where the output will be placed. default value is 'dir_data' folder (do not use a trailing slash)
dir_output = /complete/path/to/gwas_inspector/QC_output

# Folder containing reference and header files. default value is 'dir_data' folder (do not use a trailing slash)
dir_references = /complete/path/to/gwas_inspector/dir_references

allele_ref_std_population = EUR

gzip_final_dataset = TRUE
html_report = TRUE
object_file = TRUE
add_column_multiallelic = TRUE
add_column_AFmismatch = TRUE
add_column_HQ = TRUE
add_column_rsid = TRUE
add_column_AF = TRUE
add_column_hid = TRUE
ordered = TRUE
make_plots = TRUE

[filters]
# Threshold values for the high-quality (HQ) variant selection
# Variants that do not meet or exceed all four of those values will be excluded from several QC tests.
# The filters are for allele-frequency (HQfilter_FRQ), HWE p-value (HQfilter_HWE), callrate (HQfilter_cal) & imputation quality (HQfilter_imp) respectively.
HQfilter_FRQ = 0.01
HQfilter_HWE = 1e-6
HQfilter_cal = 0.95
HQfilter_imp = 0.8

# The threshold for the difference between reported and reference allele-frequency values.
# SNPs for which the difference exceeds the threshold are counted and reported. default value = 0.15
threshold_diffEAF = 0.20
```

Good, now turn to add a few extra lines to the `alt_headers.txt`, especially if we're working with GWAS summary statistics that have been produced using different softwares and by different organizations. A quick way of checking the column names and leave them in a file to check everytime we want is the following:
``` bash
for i in sumstats/*.txt; do echo $i >> colnames.txt && cat $i | head -n 1 >> colnames.txt; done
```
I ended up adding the following lines to the `alt_headers.txt` file:
```
EFFECT_ALL      ALT
OTHER_ALL       REF
EFF_ALL_FREQ    AF_ALT
EFF_ALL_FREQ    AF.alt
N_TOTAL num
PVALUE  P_BOLT_LMM_INF
```

### Run the QC
We can finally run the QC, using the following R script:

``` r
## load the package
library(GWASinspector)

## import the QC-configuration file 
job <- setup_inspector("/absolute/path/to/gwas_inspector/config.ini")

## check the created instance
## input result files that will be inspected are also displayed
job

## run the algorithm 
job <- run_inspector(job)

## check the results
## comprehensive report and result file are already saved in the output folder
result_inspector(job)
```

And executing it in your command line:
``` bash
Rscript gwas.inspector.script.cdm.R
```

To re-run the QC, the only thing you need to change is the `config.ini`, the rest remains unchanged. This is the step that takes the biggest amount of time. The output from GWASinspector includes a bunch of files with individual and collective reports and plots. We can explore them and check a few statistics as part of the quality control procedure:

- **Population stratification** - QQ plots. We should see that the distribution of Expected *vs* Observed *P*-values follows the expected pattern, probably deviating at the end, which mostly correspond to the variants that are trully associated to the trait under study. The colors indicate the AF. 

- **Allele frequency differences** - Observed*vs*Expected plots. We should see that most variants follow the expected distribution, matching the reference allele frequency (AF). It is possible to see that some of them deviate, which will be the ones we need to filter out. We should check that most of the variants are in the expected distribution. Otherwise there's something wrong with the data (it can be that the reported AF is completely opposite to the one in the reference, which would mean that the alleles are shifted in the summary statistics file). 

- **Standard Error of SNP effect sizes** - SE*vs*Nmax plots. We should see that there's no deviation from a straight line starting from the axis (0,0) and following the diagonal. As sample size increases, you would expect precision to increase proportionally. As such, you expect all your results to be on a roughly diagonal line. This is checked to detect potential issues with trait transformations contrasting the study-specific standard errors with the sample size. We can also check the effect size distribution. This serves the same purpose as the precision plot: with increased sample size you expect a reduced effect-size range because precision improves. If a file has notably wider or narrower range than expected, this strongly suggests that its results are not comparable to the others.

- **Distribution of AF** - Skewness*vs*Kurtosis plots. Kurtosis is a measure of how well a distribution matches a Gaussian distribution. A Gaussian distribution has a kurtosis of 0. Negative kurtosis indicates a flatter distribution curve, while positive kurtosis indicates a sharper, thinner curve. Skewness is a measure of distribution asymmetry. A symmetrical distribution has a skewness of 0. A positive skewness indicates a long tail towards higher values, while a negative skewness indicates a long tail towards lower values. Ideally, one expects both the skewness to be around 0 and kurtosis to be below 10.

### Filter based on the QC analysis
Then we need to filter for those markers that have been flagged as high quality (HQ =  1) by `GWASinspector`, markers whose AF do not deviate > 20% from the reference panel (highDiffEAF = 0 or NA) and only biallelic markers (MULTI_ALLELIC = 0). All those SNPs that are not found in the reference dataset, will have this column filled with `NA`, as it's not possible to calculate the difference in AF. We can do that automatically for all files using the following R script:

``` r
library(dplyr)
library(data.table)

# Do not use a trailing slash
path <- "path/to/the/output/folder/of/GWASinspector/"
setwd(path)

filenames <- list.files(path, pattern="txt.gz", full.names = TRUE)

samplenames1 <- sapply(strsplit(basename(filenames), "qc_"), `[`, 2)
samplenames1

check.tab <- NULL

for (i in samplenames1) {
        print(paste("Processing", i, "cohort..."))
        df <- fread(paste0("qc_", i))
        df1 <- df %>% filter(HQ == 1)
        df2 <- df1 %>% filter(highDiffEAF %in% c(0, NA))
        df3 <- df2 %>% filter(MULTI_ALLELIC == 0)
        check.tab <- rbind(check.tab, data.frame(Name=i, Input=nrow(df), HQ=nrow(df1), highDiffEAF=nrow(df2), MultiAllelic=nrow(df3), Pass=(nrow(df3)*100/nrow(df))))
        fwrite(df3, paste0("hq.", i), col.names = TRUE, row.names=FALSE, quote=FALSE, sep = "\t", compress = "gzip")
}

print(check.tab)

q(save = 'no')
```

And executing it in your command line:
``` bash
Rscript filter.HQ.R /path/to/the/output/folder # Do not use a trailing slash 
```

This will take some time, since the files are typically huge and reading and writing takes a lot of time. The script outputs a summary table to follow up how many variants we lost at each step of the filtering process.


**Don't forget** to gzip all the files that you don't need anymore. And actually, for the next step, the tool that we will use can work with gzipped files, so let's gzip the new files we just created. A faster way of compressing files than `gzip` is `pigz` which can parallelize if you supply more threads. As a good practice, always check how many threads are available and if you're going to use a lot for a long time (*is not the case here*) let the other users know you will do so, you're not the only one working on the server :).
``` bash
pigz filtered_sumstats/*.txt -p 8
```
