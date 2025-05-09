# Meta-analysis
This step is performed to merge all the summary statistics files coming from different cohorts and perform a meta-analysis. There are a few methods to do this:

- **Fixed effect model**. Assumes error variances are equal across the cohorts. The main bioinformatic tool using this method is [METAL](https://doi.org/10.1093/bioinformatics/btq340).
- **Random effect model**. Tests for heterogeneity in the cohorts through the Cochran's Q-test. It excludes markers whose effect-size estimates are too different between the cohorts, and therefore too uncertain to consider for downstream analysis. [MANTRA](https://onlinelibrary.wiley.com/doi/10.1002/gepi.20630), [METASOFT](), and [POPCORN](https://www.cell.com/ajhg/fulltext/S0002-9297(16)30135-5).
- **Meta-regression**. It uses genome-wide metrics of diversity between populations. Allelic effects are modelled in a linear regression with axes of genetic variation. This method is applied in [MR-MEGA](https://academic.oup.com/hmg/article/26/18/3639/3976569).
- **Heterogeneity modelling**. Models allele frequency and LD differences while accounting for heterogeneity of marginal effect sizes across populations. This method is applied in [MAMA](https://www.biorxiv.org/content/10.1101/2021.04.23.441003v1)

## Single-ancestry meta-analysis
We will use a software commonly used for this matter:
- [METAL](https://doi.org/10.1093/bioinformatics/btq340) - [Quick Start Guide](https://genome.sph.umich.edu/wiki/METAL_Quick_Start)

As with `GWASinspector`, `METAL` runs using a configuration file that we need to adapt to accomodate our data. We will use an inverse-variance weighted method to perform the meta-analysis (`SCHEME STDERR`), which is suitable for our kind of data, in which the transformation to calculate effect sizes should be very similar. These were the main options, check the file itself to see more: 

```
SCHEME   STDERR

AVERAGEFREQ ON
MINMAXFREQ ON

CUSTOMVARIABLE	N_TOTAL
LABEL N_TOTAL as N_TOTAL

MARKER	varID
WEIGHT	N_TOTAL
ALLELE	EFFECT_ALL OTHER_ALL # Very important the order in which you set here the effect allele and the other allele (!!!)
FREQ	EFF_ALL_FREQ
EFFECT	EFFECT
STDERR	STDERR
PVAL	PVALUE

PROCESS file/to/QC_output/filtered_sumstats/hq.qc.sumstat.ancestry.sex.trait.txt.gz

PROCESS [...] (All the other cohorts)

OUTFILE meta.het.hq.qc.sumstat.ancestry.sex.trait .tbl

ANALYZE HETEROGENEITY
```

``` bash
metal metal.config.txt
```

Subsequently, we will run another Rscript to retrieve the chromosome (`CHR`) and position (`BP`) columns from the `varID` column and filter for variants present in two or more cohorts and with no heterogeneity of effects between cohorts (`HetPVal` > 0.05, calculated using a Cochran's Q test). These conditions may change depending on the number and type of cohorts you count with. Another valid approach is to filter for variants present in a minimum % of the total populations (total N, from all cohorts combined). In case we are processing several meta-analysis files, we can process them iteratively:

``` r
library(dplyr)
library(data.table)

# Ask for prompt arguments
args <- commandArgs(trailingOnly = TRUE) 

# Read the input files
path <- args[1]
setwd(path)

# Set the file name and output name
filenames <- list.files(full.names = TRUE,
                        pattern = "tbl.gz")
print(paste0("Processing ", length(filenames), " files."))

# Create an empty object so that we can add the information of every iteration to the same object
check.tab <- NULL

for (i in filenames) {
  print(paste0("Processing ", i, "..."))
  # Set the input and output name
  split.name <- unlist(strsplit(basename(i), "[.]"))
  input.name <- paste(split.name[1:4], collapse = ".")
  output.name <- paste("filtered", input.name, "txt.gz", sep = ".")
  # Read the table
  df <- fread(i)
  # Filtering of the meta-analysis file to only keep markers that:
  # The marker is present in at least two cohorts
  df1 <- df %>% filter(HetDf > 0)
  # Show no heterogeneity effect (Cochran's Q-test p-value > 0.05)
  df2 <- df1 %>% filter(HetPVal > 0.05)
  # Create the chromosome and position columns from the info at the varID column
  df2$CHR <- sapply(strsplit(basename(as.character(df2$MarkerName)), ":"), `[`,1)
  df2$BP <- sapply(strsplit(basename(as.character(df2$MarkerName)), ":"), `[`,2)
  # Append the new data to check.tab
  check.tab <- rbind(check.tab, data.frame(Name=i, Input=nrow(df), Presence.2.cohorts=nrow(df1), Heterogeneity=nrow(df2), Pass=(nrow(df2)*100/nrow(df))))
  # Calculate and print LambdaGC
  chisq <- qchisq(1-df2$`P-value` , 1)
  gc <- median(chisq)/qchisq(0.5, 1)
  print(paste0("LambdaGC: ", gc))
  # Export the filtered dataframe into a tab-separated file
  fwrite(df2, output.name, col.names = TRUE, row.names=FALSE, quote=FALSE, sep = "\t", compress = "gzip")
}

# Print the check table
print("------------------------------") 
print(check.tab)
print("------------------------------")

```

Now we can run the script in our command line:
``` bash
Rscript filter.meta.iterative.R /absolute/path/to/the/meta-analysis/files/
```

Now we have an "unfiltered" meta-analysis and a filtered one. The first thing we need to do before jumping into the functional annotation is to see if we get rid of any interesting signal by applying these filters. We should at this point compare the genome-wide significant signals between the filtered and unfiltered meta-analysis summary statistics. Once we have FUMA correctly set for reading our summary statistics file, it's very easy to submit a new job (in this case for the unfiltered data) in case we want to check the genes that are mapped to the region which we (hypothetically) drop after the filtering step. 

## Multi-ancestry meta-analysis
In the case we have GWAS from different populations and we are interested in meta-analyzing them to obtain a multi-ethnic picture of the trait, we need to make use of multi-ancestry meta-analysis tools. In this case we will use:

- [MR-MEGA](https://academic.oup.com/hmg/article/26/18/3639/3976569) - [Documentation](https://genomics.ut.ee/en/tools)

Tu run MR-MEGA, we need to edit the summary statistics from GWASinspector to remove empty columns, as this tool fails ro recognize them properly and starts mixing one with another. If you want to run this from your command line, just change the input and output path to arguments you can specify as done in previous scripts:

``` r
# Import libraries
library(data.table)
library(dplyr)
# Set input path
in.path <- "/path/to/input/files/for/meta-analysis/"
setwd(in.path)
#Set output path
out.path <- "/path/to/output/meta-analysis/file/"
# List all files
filenames <- list.files(in.path, pattern="hq.sumstat.", full.names = TRUE)
print(filenames)
# Get the names of the samples
samplenames1 <- basename(filenames)
print(samplenames1)
# Iterate over them to select only the columns we want
for (i in samplenames1) {
  print(paste("Processing", i, "cohort..."))
  df <- fread(i) 
  df1 <- df[, c("varID", "CHR", "POSITION", "EFFECT_ALL", "OTHER_ALL", "EFF_ALL_FREQ", "PVALUE", "EFFECT", "STDERR", "N_TOTAL")]
  fwrite(df1, paste0(out.path, i), col.names = TRUE, row.names=FALSE, quote=FALSE, sep = "\t", compress = "gzip")
}
```

We actually only kept the columns that MR-MEGA needs. We also need to generate the `MR-MEGA.in` file, which contains a row per file (full path). Then we can run MR-MEGA as follows:
``` bash
mr-mega \
-i MR-MEGA.in \
-o output.name \
--no_std_names \
--name_pos POSITION \
--name_chr CHR \
--name_n N_TOTAL \
--name_se STDERR \
--name_beta EFFECT \
--name_eaf EFF_ALL_FREQ \
--name_nea OTHER_ALL \
--name_ea EFFECT_ALL \
--name_marker varID \
--pc xxx \
--qt # Only if it is a quantitative trait\
--debug
```

For the choice of the numbers of PC, we use these general rules: 
1. **From their documentation**: PC count must be < cohort count - 2. If 8 cohorts have been used in the analysis, then the maximum number of PC-s can be 5.
2. **From Reedik Mägi**: PC count should be <= population count - 1. Therefore, if 3 populations (*i.e.,* ancestries, ethnicities) have been used in the analysis, then the number of enough PCs to explain the variance between populations can be two.
3. **After running MR-MEGA for the first time**, use the matrix of PC coordinates that appear at the end of the log file and plot them. Visually inspect how many PCs are enough to separate accordingly the populations.
4. **After running MR-MEGA for the first time**, check the column `Ncohort` to see the total number of cohorts where the marker was present. The higher the number of PCs, the more strict will be MR-MEGA in the number of required cohorts in which a marker needs to be present. We don't want to be too strict here, as we would loose a lot of statistical power (many markers will drop).

## rsID assignment
There are a few ways of (re)assigning the rsID to each of the variants meta-analysed.

1. Using the information from the individual GWAS if those contained the rsID and adding them by matching the markers by CHR:POSITION:OTHER_ALL:EFFECT_ALL. This is an easy solution but one should be carefult to be sure that no mistakes are introduced (*e.g.,* duplicated or incorrectly assigned rsIDs).
2. Assigning *de novo* the rsIDs using the [dbSNP](https://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz) reference database and the [Human Reference Genome 1000G](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz). The pipeline runs using diverse `R` packages, `SAMtools`, `tabix`, and `BCFtools`. The first step will be to download the reference files. We need to further prepare the files, by indexing the `dbSNP` and decompressing the `Human Reference Genome` accordingly. We will then get the information from the summary statistics file (`tbl`) to generate the `vcf` with the exact header and number of columns that is expected by `BCFtools` for the comparison later on. Next, we proceed to compress the file in a compatible format with `SAMtools` using `bgzip`, we index it using `tabix`, and then we can use that file to (1) Left-align and normalize indels, check if REF alleles match the reference, split multiallelic sites into multiple rows; recover multiallelics from multiple rows and (2) annotate the variants with the rsIDs. Finally, we only need to join the newly generated table with the rsIDs and the meta summary statistics from before.

**NOTE:** The best option is to run the option (2) and then if there are some unassigned variants, try option (1) to fill them. Still, it's very much possible that you don't end up assigning the 100% of them. We have two options here then, or continuing the downstream annotation with all variants (assigned or unassigned to an rsID) or filter only those assigned by an rsID, which, to me, seems a bit unnecessary, but it's up to you.

**NOTE 2:** I did not include the scripts for option (2) because the original lines of code were from a previous collaborator and I just automated them and made them easy to run in a single workflow. If interested on them, you can contact me at cristiandiaz93@gmail.com.


## Manhattan plot and QQ plot
There's a general R package for this named [qqman](https://cran.r-project.org/web/packages/qqman/) that is very easy to use and there is also a [vignette](https://r-graph-gallery.com/101_Manhattan_plot.html) from r-graph-gallery with very clear explanations on how to generate these plots. As a note, I will just mention that I produced the QQ plot using all the variants in the final meta-file and the Manhattan plot using only variants with *p*-value < 0.05, as usually perform elsewhere. I included a few R scripts to generate fast and easy manhattan plots. Again, you only need to grant the path to make it work:

``` r
library(qqman) 
library(data.table) 
library(dplyr) 

# Ask for prompt arguments
args <- commandArgs(trailingOnly = TRUE) 

# Set the variables
path <- args[1]
setwd(path)
filenames <- list.files(pattern = "*.txt.gz", full.names = TRUE)

# Set the outpath
outpath <- paste0(path, "manhattan_plots")
outnames <- sub("\\.txt\\.gz$", "", basename(filenames))

for (i in 1:length(filenames)) {
  # Read the GWAS file
  df <- fread(filenames[i])
  # Make sure these columns are numeric
  df$CHR <- as.numeric(df$CHR)
  df$POSITION <- as.numeric(df$POSITION)
  # Make the Q-Qplot.
  jpeg(paste0(outpath, "qq.plot.", outnames[i], ".jpg"), width = 8, height = 6, units = "in", res = 300)
  qq(df$PVALUE)
  dev.off()
  # Filter for SNPs with a p-value <= 0.05 to reduce the size of the df Also, Manhattan plots are normally plotted with variants p<0.05
  df1 <- df %>% filter(PVALUE <= 0.05)
  # Fill the ID column with the marker names whenever there is no rsID assigned
  df1$ID[df1$ID == "."] <- df1$SNP[df1$ID == "."]
  # Create a character vector with the SNPs to highlight
  high.snps <- df1$ID[df1$PVALUE <= 5e-8]
  # If there are significant snps
  if (length(high.snps > 0)) {
    write.table(high.snps, file = paste0(outpath, outnames[i], ".high.snps.txt"), 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    # Make the Manhattan plot. 
    jpeg(paste0(outpath, "manhattan.", outnames[i], ".jpg"), width = 16/1.5, height = 9/1.5, units = "in", res = 300)
    manhattan(df1, chr = "CHR", bp = "POSITION", snp = "ID", p = "PVALUE", 
              col = c("#cccccc", "#969696"), highlight = high.snps, 
              suggestiveline = -log10(5e-6), annotatePval = 5e-8, annotateTop = TRUE)
    dev.off()
  } else {
    # Make the Manhattan plot. 
    jpeg(paste0(outpath, "manhattan.", outnames[i], ".jpg"), width = 16/1.5, height = 9/1.5, units = "in", res = 300)
    manhattan(df1, chr = "CHR", bp = "POSITION", snp = "ID", p = "PVALUE", 
              col = c("#cccccc", "#969696"), highlight = FALSE, 
              suggestiveline = -log10(5e-6), annotatePval = 5e-8, annotateTop = FALSE)
    dev.off()
  }
  
  print(paste0("DONE with ", outnames[i], "!"))
}
```






