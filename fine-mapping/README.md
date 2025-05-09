# Fine-mapping

There is a very nice review on this topic which I used to write this section: [From genome-wide associations to candidate causal variants by statistical fine-mapping](https://www.nature.com/articles/s41576-018-0016-z). 

The general strategy is to use the GWAS list of SNPs associated with a trait to identify regions of interest.Each region is the visually explored for its LD structure and for genes known to be mapped to the region. The selected SNPs are further evaluated for their likely function based on publicly available genomic annotation. The overarching goal is to determine which variants are most likely to be functional and to quantify the strength of evidence. To do that, additional statistical methods are needed in order to discriminate likely functional variants from variants that are merely correlated with the functional variants. 

A limitation of the leadSNP is that there is a reasonable chance that it is not the causal variant. True associations re not likely to result in the smallest *p*-values, in part owing to the small effect sizes of the variants on complex traits. These findings emphasize the impotance of caution when considering the lead SNP as likely causal and the importance of fine-mapping in order to identify the causal variant or variants. A number of approaches have been used to perform fine-mapping. We present three main strategies that have been used in the literature: heuristic methods, penalized regression models and Bayesian methods:

1. **Heuristic fine-mapping approaches**. Examine the correlation (r<sup>2</sup>) among the SNPs surrounding a lead SNP retaining as potentially causal only those SNPs with an r<sup>2</sup> above a threshold. They do not account for the joint effects of the SNPs on the trait, and they do not give an objective measure of the confidence that a SNP is causal but rather rely on somewhat arbitrary thresholds and subjective interpretations of correlations among SNPs.
2. **Penalized regression models**. These models simultaneously perform estimation of SNP effect sizes and SNP selection (using P values to determine whether a SNP should be included in a model) into a model by shrinking small effect estimates towards zero. Penalized models tend to result in sparse models, selecting only one or a few SNPs belonging to a group of correlated SNPs. This can result in a good prediction model that includes non- causal SNPs and excludes a causal SNP when they are highly correlated.
3. **Bayesian methods**. Bayesian methods have been specifically designed for fine- mapping, offering advantages over heuristic and penalized regression approaches. Bayesian models also try to determine which SNPs have non-zero effect sizes (regression $\beta$-values) on a trait. Bayesian inference focuses on the probability of a specific hypothesis or specific model, thus providing probabilistic interpretation of models of interest. A model for fine-mapping can be represented by an indicator variable for each SNP, with values of 1 for causal and 0 for not, and by organizing these indicators for all SNPs of interest in vector *c*. For *m* SNPs, there are 2<sup>*m*</sup> possible *c* vectors (hence, 2<sup>*m*</sup> possible models), ranging from all values of *c* equal to 0 for no SNPs causal to all values equal to 1 for all SNPs causal. Using Bayes' formula, the prior probability of a model can be combined with the likelihood of the data *D* (trait and SNPs) to compute the posteriror probability of a specified model M<sub>*c*</sub>. There are several ways to specify the prior probability for a model, such as assuming that variants are independent and equally likely to be causal or assuming a fixed number of causal variants out of the total variants. The posterior probabilities for different models can be used to determine the posterior probability of including each SNP in any of the models (posterior inclusion probability (PIP)), as well as determining the minimum set of SNPs required to capture the likely causal SNPs (credible sets). Ranking SNPs by their PIP is a convenient way to select putative causal SNPs.


**Notes:**

- **Prior probability**: In Bayesian probability theory, the probability distribution assigned to parameters of interest, specified to represent prior knowledge of their values before observing the data.
- **Posterior probability**: In Bayesian probability theory, the updated probability distribution of parameters of interest, conditional on the observed data.
- **Posterior inclusion probability**: The marginal probability that a SNP is included in any causal model, conditional on the observed data, thereby providing weight of evidence that a SNP should be included as potentially causative.
- **Credible sets**: the minimum set of SNPs that contains all causal SNPs with probability α.

## PolyFun

In this section we are going to focus on [PolyFun](https://doi.org/10.1038/s41588-020-00735-5), which allows you to incorporate prior-probabilities of the causality of the SNPs and it uses two Bayesian fine-mapping tools:
- [FINEMAP](https://doi.org/10.1093/bioinformatics/btw018) - [Website](http://www.christianbenner.com/)
- [SuSiE](https://doi.org/10.1111/rssb.12388) - [Website](https://stephenslab.github.io/susie-paper/)

### LD matrices

This is tpypically a correlation matrix between the SNPs of the region of interest that can be calculated directly from the sample (genotype) data or from an external reference panel (*i.e.,* 1000G). Since we are conducting a meta-analysis we can't obtain this matrix from the genotype data, so we are going to retrieve this matrix from the 1000G reference panel. The LD matrix must contain SNPs included in the region surrounding the lead SNP identified by FUMA. I used the "start" and "end" region also identified by FUMA to flank the regions of interest for each significant locus. With that, we can extract the list of rsIDs from each region of interest:

``` r
# Import the summary statistics file
sumstats <- fread(path/to/your/meta-analysis/file")
# Import a data frame containing every significant locus identified
genriskloc <- read_excel("path/to/your/metadata/excel/file")

# Create function to get the variants from each significant locus
getLocusVariants <- function(df, chr, start, end){
  region <- df %>% 
    filter(Chromosome == chr) %>% 
    filter(Position >= start & Position <= end)
  snplist <- region$ID
  return(snplist)
}

# Create a for loop to generate a file per significant locus containing the list of variants
for (i in 1:nrow(genriskloc)){
  chr <- genriskloc$chr[i]
  start <- genriskloc$start[i]
  end <- genriskloc$end[i]
  a <- getLocusVariants(sumstats, chr, start, end)
  write.table(a, paste0("W:/USERS/cdmunoz/stool_traits/fine_mapping/", genriskloc$Name[i]), row.names=FALSE, quote=FALSE, sep = "\t")
}
```

With that list of variants we can then obtain a correlation matrix using the 1000G reference panel. To do that, I chose for the [LDmatrix](https://ldlink.nih.gov/?tab=ldmatrix) web-based tool, from the [LDlink](https://ldlink.nih.gov/?tab=home) set of tools of the NIH. I selected the EUR populations and downloaded the correlation matrix. With the correlation matrix downloaded I imported it into the R space to work with it:

``` r
r2_matrix <- as.matrix(read.table("path/to/r2_matrix/from/your/region/of/interest"))
r2_matrix[is.na(r2_matrix)] <- 0
```

If you count with UKB data within the individual GWAS you performed the meta-analysis with, it is also possible to get the LD matrix using the UKB genotype data. We can download the UKB precomputed LD matrices from the Amazon Web Services (AWS). To do that, we need to install `awscli`:

``` bash
# With the python installation of conda
pip install awscli 
# Verify Installation
aws --version
# List Files in the S3 Bucket
aws s3 ls --no-sign-request s3://broad-alkesgroup-ukbb-ld/
# Download Files from S3 Bucket
aws s3 cp --no-sign-request s3://broad-alkesgroup-ukbb-ld/ . --recursive
```

There are way too many files. We need to do this only for the regions of interest. In this database, every file spans a genomic region of 3 Mbp. In order to download only the necessary files for us, I created this Rscript in which I use the start and end position identified by FUMA for every significant locus:

``` r
library(tidyverse)
library(readxl)

# Read meta-data
meta <- read_excel("path/to/your/metadata/excel/file", sheet = 3)

# Set options to avoid scientific notation
options(scipen = 999)

# Get the round number around the start of the region
meta$start_round <- meta$start - (meta$start %% 1000000)
meta$start_round <- meta$start_round +1
meta$end_round <- meta$start_round + 3000000

meta$section <- paste0("chr", meta$chr, "_", meta$start_round, "_", meta$end_round)

meta2 <- meta[34] %>% 
  unique() %>% 
  mutate(section_gz = paste(section, "gz", sep = ".")) %>% 
  mutate(section_npz = paste(section, "npz", sep = "."))

write_delim(meta2[2], "path/to/output/the/sections_to_download_gz.txt", col_names = FALSE, quote = "none", delim = "\n")
write_delim(meta2[3], "path/to/output/the/sections_to_download_npz.txt", col_names = FALSE, quote = "none", delim = "\n")
```

Then I used these files to download all necessary LD matrices:

``` bash
for i in `cat path/to/output/the/sections_to_download_gz.txt`;
do aws s3 cp --no-sign-request s3://broad-alkesgroup-ukbb-ld/UKBB_LD/$i ./;
done

for i in `cat path/to/output/the/sections_to_download_npz.txt`;
do aws s3 cp --no-sign-request s3://broad-alkesgroup-ukbb-ld/UKBB_LD/$i ./;
done
```

### Non-functional fine-mapping

To standardize the sumstats file for any analysis done with polyfun, we first need to munge the summary statistics file:

``` bash
python munge_polyfun_sumstats.py 
--sumstats input_sumstats/meta.analysis.ancestry.sex.trait.txt.gz 
--out munged_sumstats/meta.analysis.ancestry.sex.trait.parquet 
--min-info 0 
--min-maf 0
```

With the sumstats prepared and the LD matrices downloaded for the specific region of interest, we can continue the non-functional fine-mapping pipeline. I used the same Excel file than before to get all the options for the polyfun script so I can run every region consecutively in a for loop with a bash script that you can find at the end of this section. Down here you can find a simple example for only one region (CHR, START, and END should be changed by the actual numbers):

``` bash
python finemapper.py 
--ld /path/to/UKB_LD_matrices/CHR_START_END 
--sumstats munged_sumstats/meta.analysis.ancestry.sex.trait.parquet
--chr CHR
--start START
--end END 
--method susie 
--max-num-causal 10 
--out susie.UKB.CHR.START.END.gz 
--n XXX 
--non-funct
--susie-outfile susie.object.UKB.CHR.START.END.RDS
--no-sort-pip
--verbose
```

Afterwards I ran the fine-mapper script using finemap instead of susie to compare the results:

``` bash
python finemapper.py 
--ld /path/to/UKB_LD_matrices/CHR_START_END
--sumstats munged_sumstats/meta.analysis.ancestry.sex.trait.parquet
--chr CHR
--start START
--end END 
--method finemap 
--finemap-exe ~/bin/finemap 
--max-num-causal 10 
--out polyfun_finemap_output/finemap.UKB.CHR.START.END.gz
--n xxx 
--non-funct
--verbose
```

In my case, the results are very much comparable between FINEMAP and SuSiE. FINEMAP outputs, however, the probability of the credible set to consist on 0, 1, 2, ..., 10 variants, which SuSiE does not generate. This is very useful to know which credible set of variants to report. It is generated as standard output, so not in the log file, so you need to check it right away (or save the standard output to a text file). 

### PolyFun for all regions of interest in a go

``` bash
# We first need to define empty variables. One for each flag we want to include
meta_flag=''
verbose='false'

# We print how to use the script by default
print_usage() {
  printf "Usage: [-m] <meta-data>\n" $0
}

# We use the package getopts, which is embedded in bash.
# The column (:) indicates that the flag must be followed by a file.
# What we put before the ')' is the actual flag we need to write in the command line when running the script.
# The actual variable is stored when granting the "$OPTARG".
while getopts 'm:v' flag; do
  case "${flag}" in
    m) meta_flag="true"
       meta_val="$OPTARG";;
    v) verbose='true' ;;
    ?)   printf "Usage: [-i] <input> [-m] <meta-data>\n" $0
          exit 2;;
  esac
done

outpath="output_name";
ldpath="/path/to/UKB_LD_matrices/"

for i in {1..26}; do
  # Obtain the variables for the python script
  the_chromosome=$(cut -f 4 $meta_val | tail -n +2 | awk "NR==$i");
  the_start=$(cut -f 7 $meta_val | tail -n +2 | awk "NR==$i");
  the_end=$(cut -f 8 $meta_val | tail -n +2 | awk "NR==$i");
  the_start_round=$(cut -f 32 $meta_val | tail -n +2 | awk "NR==$i");
  the_end_round=$(cut -f 33 $meta_val | tail -n +2 | awk "NR==$i");
  the_sumstats=$(cut -f 31 $meta_val | tail -n +2 | awk "NR==$i");
  the_section=$(cut -f 34 $meta_val | tail -n +2 | awk "NR==$i");
  the_ld="$ldpath$the_section";
  the_trait=$(cut -f 24 $meta_val | tail -n +2 | awk "NR==$i");
  the_ancestry= # Set one of EUR, EAS, SAS, AMR, AFR
  if [[ $the_trait = "TRAIT_1" ]]; then
    n=XXX
  elif [[ $the_trait = "TRAIT_2" ]]; then
    n=XXX
  elif [[ $the_trait = "TRAIT_3" ]]; then
    n=XXX
  else
    echo "Error: trait not found"
  fi
  echo "--------------------------------------------------------------------------"
  echo "Running fine-mapping on $the_trait.$the_ancestry.$the_chromosome.$the_start.$the_end"
  echo "--------------------------------------------------------------------------"
  python finemapper_modified.py\
  --ld /path/to/UKB_LD_matrices/$the_chromosome\_$the_start_round\_$the_end_round \
  --sumstats $the_sumstats\
  --chr $the_chromosome\
  --start $the_start\
  --end $the_end\
  --method susie\
  --max-num-causal 10\
  --out $outpath/susie.$the_trait.$the_ancestry.$the_chromosome.$the_start.$the_end.txt\
  --n $n\
  --non-funct\
  --susie-outfile $outpath/susie.object.$the_trait.$the_ancestry.$the_chromosome.$the_start.$the_end.RDS\
  --allow-missing\
  --no-sort-pip\
  --verbose;
done
```

## MR-MEGA

## MR-MEGA Bayesian fine-mapping

In case you are dealing with multi-ancestry meta-analysis of GWAS, MR-MEGA counts with an inner algorithm to calculate the probability of each variant to be the causative variant for the phenotype you are testing the association for, accounting for the heterogeneity that is correlated with ancestry. The output is called the Bayes'factor in favour of association. We need to calculate the natural logarithm of the Bayes'factor for each variant, select a window size (1Mbp) and from the set of variants in that window (locus) then derive a 95% (or 99%) credible set for the association signal by: (i) ranking all variants according to their Bayes’ factor; and (ii) including ranked variants until their cumulative posterior probability of driving the association attains or exceeds 0.95 (or 0.99). We accomplished that by using a simple `Rscript`. We just need the summary statistics for the multi-ancestry meta-analysis (the output from MR-MEGA) and a list of loci to finemap.

``` r
library(data.table)
library(dplyr)
library(tidyr)

# MR-MEGA
# Read the sumstats
sumstats <- fread("path/to/your/multi.ancrestry.meta.analysis.txt.gz")
# Read the metadata containing the information of every locus we want to finemap with the respective start and end coordinates
meta <- read.table("significant.loci.meta.data.multi.sf.txt", header = TRUE, sep = "\t")
# Get the lead SNPs
leadSNPs <- meta$rsID
# Define the function to get the credible set
GetCredSet <- function(df, v) {
  # Set an empty dataframe
  credible.sets <- NULL
  # Iterate over the rows, which are sorted by their PIP, in descending order
  for (i in 1:nrow(df)) {
    # Include the row if the cumulative PIP is under the 95 %
    if (v[i] < 0.95) {
      credible.sets <- rbind(credible.sets, df[i]) 
      # Include the row if the cumulative PIP goes over 95 % the first time
    } else if (v[i] > 0.95) {
      credible.sets <- rbind(credible.sets, df[i])
      # Exit the loop when the condition is met
      break 
    }
  }
  return(credible.sets)
}

cs.mega <- NULL
# Iterate over all loci of interest
for (i in leadSNPs){
  pos <- sumstats$POSITION[sumstats$ID == i]
  chr <- sumstats$CHR[sumstats$ID == i]
  locus <- sumstats %>% 
    filter(CHR == chr) %>% 
    filter(POSITION > (pos-1000000), POSITION < (pos+1000000)) %>% 
    mutate(START = pos-1000000) %>% 
    mutate(BF = exp(LNBF)) %>% 
    mutate(CUMBF = sum(BF)) %>% 
    mutate(PIP = BF/CUMBF) %>% 
    arrange(desc(PIP)) %>% 
    mutate(CUMPIP = cumsum(PIP))
  top <- GetCredSet(locus, locus$CUMPIP)
  cs.mega <- rbind(cs.mega, top)
}

causals.mega <- cs.mega %>% 
  filter(PIP >= 0.5) 

nvar.mega <- cs.mega %>% 
  mutate(uniqID = paste(CHR, START, sep = "-")) %>% 
  group_by(uniqID) %>% 
  count(uniqID)

write.table(cs.mega, "output_mr-mega/credible.sets.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
write.table(causals.mega, "output_mr-mega/high.causal.variants.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
write.table(nvar.mega, "output_mr-mega/number.variants.credible.set.tsv", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")
```
