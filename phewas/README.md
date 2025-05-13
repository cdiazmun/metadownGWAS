# Phenome-wide association study (PheWas)

A PheWAS is nothing else than to take a variant and see if it has been previously associated to other phenotypes by inspecting databases of thousands of GWAS already published. There are different web services and tools that rely on different sources that are used as databases to perform the query searches to. There used to be a tool very widely used for this, named PhenoScanner (coming as an R package), but the server is down since March 2024 and there's no foreseen improvement in the situatuion. Therefore, I have used the following three tools:

1. [GWAS ATLAS](https://atlas.ctglab.nl/). It has 600 GWAS from the UKB and 4,756 in total that are community-added and maintained. You need to do the searches variant per variant.
2. [Open Targets Genetics](https://genetics.opentargets.org/). It has +12,000 GWAS from UKB, FinnGen, and GWAS Catalog. You need to do the searches variant per variant. It also gives you eQTL information.
3. [IEU OpenGWAS](https://gwas.mrcieu.ac.uk/). It has +13,000 GWAS from UKB, FinnGen, Biobank of Japan, EBI and own curated GWAS. You can install the R package ieugwasr and perform batch searches using a list of variants. Also, you can use this package to extract the proxies in the nearby region and perform the PheWAS not only on the leadSNPs but also on the proxies in LD.

There is also a tool from the NIH LDlink resource, named [LDtrait](https://ldlink.nih.gov/?tab=ldtrait) that is also worth highlighting, as it automatically recognizes variants in LD and use them to detect associated traits. 

## The ieugwasr package
To use this package, you first need to set it up, as it needs to connect to the server using an ID. The searches that we launch to their server will be limited by ID, but under my experience, the limit is veeeery large.

``` r
install.packages("ieugwasr")

# Add OPENGWAS_JWT=<token> to your .Renviron file.
usethis::edit_r_environ()

# To check that your token is being recognised
ieugwasr::get_opengwas_jwt()

# To check that your token is working
library(ieugwasr)
user()

# get API status
api_status()
```

Then we are going to calculate the proxies for the leadSNPs of interest:

``` r
# Import the libraries needed
library(dplyr)
library(data.table)
library(ieugwasr)

# Set the working directory
setwd("/path/to/your/working/directory/")

# Define a function to detect the SNPs in LD with your lead SNP
ld_proxy <- function(df, leadsnp, ancestry, flank){
  # Define the locus by obtaining the genomic coordinates and flanking region
  locus_chr <- df$CHR[df$ID == leadsnp]
  locus_start <- df$POSITION[df$ID == leadsnp] - flank
  locus_end <- df$POSITION[df$ID == leadsnp] + flank
  locus_region <- df %>% 
    filter(CHR == locus_chr, POSITION >= locus_start, POSITION <= locus_end)
  print(paste0("Number of SNPs in region: ", nrow(locus_region)))
  # Calculate the LD correlation matrix (signed, an not squared)
  r_matrix <- as.data.frame(ld_matrix(variants = locus_region$ID, pop = "EAS", with_alleles = FALSE))
  print(paste0("Number of SNPs in matrix: ", nrow(r_matrix)))
  # Obtain the r2 information for our lead SNP
  r2_leadsnp <- r_matrix[leadsnp] ^ 2
  print(paste0("Length of r2 vector: ", nrow(r2_leadsnp)))
  # Reformat the dataframe to get the IDs and the r2 info
  proxies <- r2_leadsnp %>% 
    mutate(ID = rownames(r2_leadsnp))
  colnames(proxies)[1] <- "R2"
  # print(proxies$R2[proxies$ID == "rs4975013"])
  # Get the proxies, defined by r2 >= 0.8 
  proxies <- proxies %>% 
    filter(R2 >= 0.8)
  # Join the two dataframes and filter in variants with p-value < 1e-5
  final <- locus_region %>% 
    inner_join(proxies, by = "ID") %>% 
    filter(PVALUE <= 1e-5)
  return(final)
}

# Import the GWAS meta-analysis 
sumstats <- fread("/path/to/your/meta.analysis.ancestry.sex.trait.txt.gz")

# An example using one of the significant locus
locus1 <- ld_proxy(df, "rs123456789", "EUR", 100000)
write.table(locus1, "proxies_significant_loci/locus1.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
```

Once we have all the proxy files (*i.e.,* the sumstats for all variants in LD for each locus), we can run the PheWAS for the leadSNPs and their proxies:

``` r
# Import the libraries needed
library(readxl)
library(tidyverse)
library(ieugwasr)

# Set the working directory
setwd("/path/to/your/working/directory/")

# Get the file names of all significant loci
filenames <- list.files("proxies_significant_loci/", pattern = ".txt", full.names = TRUE)

loci <- sapply(strsplit(basename(filenames), split = ".", fixed = TRUE), `[`, 1)

phewas <- NULL

for (i in 1:length(filenames)){
  # Import proxies, per locus
  df <- read.table(filenames[i], header = TRUE)
  # Define the leadSNP and the nearest gene
  LEADSNP <- df$ID[df$PVALUE == min(df$PVALUE)]
  GENE <- loci[i]
  # Get the list of proxy variants to perform the PheWAS
  proxies <- df$ID
  # Perform the PheWAS
  pw <- phewas(variants = proxies, 
               pval = 5e-8)
  if (nrow(pw) > 0){
    # Collapse common traits to avoid redundancy and assign all this data to a new column named leadSNP
    pw1 <- pw %>%
      group_by(trait) %>% 
      slice_min(p) %>% 
      mutate(leadSNP = LEADSNP, Gene = GENE)
    # Concatenate all novel loci in a single data frame
    phewas <- rbind(novel_phewas, pw1)
  }
  else {
    print(paste0("No signals found for ", LEADSNP))
  }
}

write.table(phewas, "phewas/summary_phewas_proxies_opengwas.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
```

## GWAS ATLAS, OpenTargetGenetics, and LDTrait

To use these tools, you need to go locus by locus using their website. Easier to do, slower to execute and analyse.
