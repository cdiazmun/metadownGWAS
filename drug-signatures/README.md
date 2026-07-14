# Drug analysis
There are several ways of performing this kind of analysis. The objective is to find drugs that are approved or in development stage that can target the genes or pathwyas that the genes we have associated with the phenotype under study are part of.  or by performing an unbiased pipeline that uses the summary stats to perform a kind of gene set enrichment analysis (GSEA) for targeting enriched functions or pathways. 

## Manual query in specific databases
We can do that by manually inspecting the genes or pathways of interest in dedicated databases, such as:
- [DrugBank](https://go.drugbank.com/)
- [Therapeutic Target Database](https://db.idrblab.net/ttd/)
- [Open Targets Platform](https://platform.opentargets.org/)
- [Drug-Gene Interaction Database](https://dgidb.org/)

## Genes-tissue signature and connectivity mapping

The translational potential of GWAS findings for drugs discovery can be assessed using a two-step computational pipeline translating genetic variation associated with the phenotype under study into a candidate gene expression signature that can be used to screen public expression databases for compounds able to mimic or oppose such signature.

A gene-level association analysis can be performed using [S-MultiXcan](https://doi.org/10.1371/journal.pgen.1007889), a tool that combines GWAS summary statistics and gene-expression prediction models for multiple cells/tissues. This tool is embedded in [CTG-VL](https://vl.genoma.io/). We should restrict the search to relevant tissues for the disease. For example, tissue and cell-type analysis of IBS showed that there is an increased heritability in the central nervous system and at the same time, there were cell types within the enteric nervous system that were highlighted as some genes were idfferentially expressed there. This can be used to direct our research question toward specific GTEx tissues, such as, brain and gastrointestinal tissues. 

*__Important:__* If you're going to run several tissues where you expect difference results (*i.e.,* brain and colon), run them separately. Here is important to establish the focus of your study (are you interested in drugs targeting genes dysregulated in the brain only? everywhere in the body? in the GI tract?). The research question matters a lot here. 

The output of `S-MultiXcan` should have the following columns (March 2026):

- `gene`: gene identifier used in the transcriptome prediction model (typically an Ensembl gene ID).
- `gene_name`: gene symbol (HGNC gene name) corresponding to `gene`.
- `pvalue`: significance p-value of the S-MultiXcan association.
- `n_models`: number of tissue models available for the gene and included in the multivariate analysis.
- `n_indep`: number of independent components retained after dimensionality reduction (PCA/SVD) of predicted expression across tissues.
- `p_i_best`: smallest (most significant) p-value among all single-tissue S-PrediXcan associations.
- `t_i_best`: tissue corresponding to the best single-tissue S-PrediXcan association.
- `p_i_worst`: largest (least significant) p-value among all single-tissue S-PrediXcan associations.
- `t_i_worst`: tissue corresponding to the worst single-tissue S-PrediXcan association.
- `eigen_max`: In the PCA decomposition of predicted expression, the maximum eigenvalue.
- `eigen_min`: In the PCA decomposition of predicted expression, the minimum eigenvalue.
- `eigen_min_kept`: In the PCA decomposition of predicted expression, the minimum eigenvalue retained after PCA/SVD truncation (i.e., after removing near-zero components).
- `z_min`: minimum single-tissue S-PrediXcan Z-score observed across all tissues for the gene.
- `z_max`: maximum single-tissue S-PrediXcan Z-score observed across all tissues for the gene.
- `z_mean`: mean of the single-tissue S-PrediXcan Z-scores across all available tissues.
- `z_sd`: standard deviation of the single-tissue S-PrediXcan Z-scores across all available tissues.
- `status`: reports whether the computation completed successfully; if not, contains the error or warning message explaining the issue.

You should set a threshold to select the most significant gene-tissue signatures:
- `pvalue`: I applied a false discovery rate (FDR) correction ofr multiple tests and selected only significant ones P<sub>FDR</sub><0.05: The signature has a high level of confidence. 
- `z-mean`: I applied a threshold to select altered gene expression across datasets using the effects observed |z-mean| > 2.0: Te signature has a heavy impact on phenotype risk.

We can finally use that list of genes that are down/over-regulated in gastro-relevant tissues as a query for [iLincs](https://doi.org/10.1038/s41467-022-32205-3), which performs a connectivity map analysis to match this list of genes with pharmacological compounds that also alter their expression:

<img width="611" height="66" alt="image" src="https://github.com/user-attachments/assets/6048395c-4da4-48ea-8b42-6984d070c56e" />

In general, a perturbagen connectivity analysis is performed in order to highlight compounds that could cause perturbations in gene expression in concorance or discordance with the queried signature. More specifically, in iLincs we can perform these two different analysis:
1) __Signature-to-Signature:__ Signatures in iLincs libraries that are connected to the current signature. Micro-level. This is a pairwise comparison in which iLincs computes a signed similarity measure (called concordance) between two ranked gene vectors: (1) Genes that move in the same direction (increase similarity) and (2) genes that move in the opposite direction (decrease similarity).
   - __Comparison__: Your query signature _vs_ One specific reference signature from the iLINCS database. Each reference signature corresponds to one perturbagen, one dose, one time point, one cell line, one experiment.
   - __Positive signature-to-signature connection__: The specific perturbational condition produces a transcriptional state similar to your phenotype signature. Phenocopying conditions, mechanisms, disease states.
   - __Negative signature-to-signature connection__: That specific perturbational condition reverses your phenotype-like signature. Potential therapeutic effects, mechanistic antagonism.
2) __Signature-to-Perturbation:__ The perturbagen connectivity analysis compares the query signature to all signatures for a given perturbagen as a group, thus extending the pair-wise connectivity analysis to account for diversity of responses in different cellular contexts. The question here is if across all experimental contexts in which this compound was tested, does this perturbagen tend to mimic or reverse my query signature. This is a compound-level question, not an experiment-level question.
   - __Comparison__: Your query signature _vs_ ALL signatures associated with a given perturbagen, pooled together.
   - __Positive signature-to-signature connection__: 		Across many contexts, this compound tends to _reproduce_ the phenotype-associated transcriptional pattern. Useful for understanding biology, not therapy.
   - __Negative signature-to-signature connection__: Across many cellular contexts, this compound tends to consistently _reverse_ your phenotype-associated signature. This is the strongest evidence for therapeutic relevance in connectivity mapping.

<img width="644" height="320" alt="image" src="https://github.com/user-attachments/assets/233ea888-4c80-4bbb-98c6-cbcc28402d49" />

The list of unique perturbagens can then be classified according to the Anatomical-Therapeutic-Chemical (ATC) classes from the [World Health Organization Collaborating Centre for Drug Statistics Methodology](atcddd.fhi.no/atc_ddd_index/). This list of classified drugs can then be subjected to an enrichment analysis to identify overrepresented classes at the three first ATC levels. We can do that as follows:

<img width="695" height="682" alt="image" src="https://github.com/user-attachments/assets/82928f5d-2847-4616-b568-79e8b57ddb47" />


Import libraries and set the working directory.
``` r
library(tidyverse)
library(readxl)

setwd("path/to/your/working/directory")
```

To classify the iLincs perturbagen compounds with ATC codes, we need to get the database already formatted as a tsv for easy manipulation. I left [here](https://github.com/cdiazmun/MDAgwas/blob/main/Drug-analysis/WHO%20ATC-DDD%202024-07-31.csv) a version of the database from 2024. But I have seen that the [original repository](https://github.com/fabkury/atcd) I used to get it from has updated versions, check it out.

``` r
# Read the list of ATC drugs and their categories
who <- read.csv("WHO ATC-DDD 2024-07-31.csv")
# Some entries ar repeated because of different doses, just get the uniques
who1 <- who %>% 
  distinct(atc_code, .keep_all = TRUE) %>% 
  select(atc_code, atc_name)
# Change names to UPPERCASE to match properly
who1$atc_name <- toupper(who1$atc_name)
```

Then we need to incorporate the list of perturbagen from iLincs. In this case, I downloaded the results from the __Signature-to-Perturbation__ analysis:

``` r
# Read the list of perturbagens
df <- read_excel("path/to/your_phenotype_ConnectedPerturbations_LIB_5_year_month_day_h_m_s.xlsx") %>% 
  mutate(fdr = p.adjust(p = pValue, method = "fdr"), 
         bonferroni = p.adjust(p = pValue, method = "bonferroni"),
         trait = "PHENO",
         Perturbagen = toupper(Perturbagen))
```

Join both data frames and use the column `atc_code` to generate new columns with the different ATC levels:

``` r
df.atc <- df %>% 
  left_join(who1, by = c("Perturbagen" = "atc_name")) %>% 
  mutate(ATC_1 = sapply(strsplit(atc_code, ""), `[`, 1),
         ATC_2 = if_else(is.na(ATC_1) == TRUE, NA, paste0(sapply(strsplit(atc_code, ""), `[`, 1), 
                                                          sapply(strsplit(atc_code, ""), `[`, 2),
                                                          sapply(strsplit(atc_code, ""), `[`, 3))),
         ATC_3 = if_else(is.na(ATC_1) == TRUE, NA, paste0(ATC_2, 
                                                          sapply(strsplit(atc_code, ""), `[`, 4))))

# Check how many perturbagens are classified
length(df.atc$Perturbagen[is.na(df.atc$atc_code)])

# Select those with ATC code assigned
df.atc1 <- df.atc %>% 
  filter(!(is.na(atc_code)))

# Write down the tables
write.table(df.atc, "path/to/iLincs/phenotype.ilincs.atc.txt",
            row.names = F, col.names = T, quote = F, sep = "\t", na = "")
write.table(df.atc1, "path/to/iLincs/phenotype.ilincs.atc.classified.txt",
            row.names = F, col.names = T, quote = F, sep = "\t", na = "")
```

Now we need to incorporate the list of all iLincs perturbagens, which you can find [here](https://github.com/cdiazmun/MDAgwas/blob/main/Drug-analysis/sigSearch_LIB_5_2026_4_20_15_5_13.xls) for download. This is the database that we will use to compare our list with and check if there is any group of compounds enriched in our dataset.  

``` r
# Import the list of all iLincs perturbagens 
pert <- fread("path/to/iLincs/sigSearch_LIB_5_2026_4_20_15_5_13.xls")
# Check the data frame
head(pert)
# Remove duplicated compounds tested at different concentrations
pert1 <- pert %>% 
  group_by(Perturbagen) %>% 
  slice_min(order_by = Concentration) %>% 
  ungroup() %>% 
  mutate(Perturbagen = toupper(Perturbagen))
# Check the data frame again
head(pert1)
# Join the data frame with that containing the ATC codes
pert2 <- pert1 %>% 
  left_join(who1, by = c("Perturbagen" = "atc_name"))
# Check how many perturbagens are classified
length(pert2$Perturbagen[is.na(pert2$atc_code)])
# Use the column atc_code to generate new columns with the different ATC levels, as before
pert3 <- pert2 %>% 
  filter(!(is.na(atc_code))) %>% 
  mutate(ATC_1 = sapply(strsplit(atc_code, ""), `[`, 1),
         ATC_2 = paste0(sapply(strsplit(atc_code, ""), `[`, 1), 
                        sapply(strsplit(atc_code, ""), `[`, 2),
                        sapply(strsplit(atc_code, ""), `[`, 3)),
         ATC_3 = paste0(ATC_2, 
                        sapply(strsplit(atc_code, ""), `[`, 4)))
# Check 
table(pert3$ATC_1)
```

We are going to finally compare the percentage of ATC codes (compounds) assigned to each category:

``` r
# Reshape the phenotype data for joining later
df.atc_1 <- na.omit(df.atc) %>% 
  group_by(ATC_1) %>% 
  summarize(n = n()) %>% 
  mutate(perc = n / sum(n) * 100,
         trait = "PHENOTYPE")
# Reshape the iLincs data for joining later
lincs.atc_1 <- pert3 %>% 
  group_by(ATC_1) %>% 
  summarize(n = n()) %>% 
  mutate(perc = n / sum(n) * 100,
         trait = "iLincs")
# Join both data frames
atc_1 <- rbind(df.atc_1, lincs.atc_1) %>% 
  group_by(trait) %>% 
  mutate(rank = rank(-n, ties.method = "first")) %>% 
  ungroup() %>% 
  mutate(trait = factor(trait, levels = c("PHENOTYPE", "iLincs"))
# Calculate the total number of compounds per trait and 
atc_1_stats <- atc_1 %>% 
  select(-perc, -rank) %>% 
  group_by(trait) %>% 
  mutate(Total = sum(n)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = trait, values_from = c(Total, n))
# Transform NA values to 0, otherwise the Fisher test fails
atc_1_stats$n_PHENOTYPE[is.na(atc_1_stats$n_PHENOTYPE)] <- 0
atc_1_stats$Total_PHENOTYPE[is.na(atc_1_stats$Total_PHENOTYPE)] <- 0
# Run Fisher's Exact Test for each category and compute Odds Ratio (OR)
atc_1_res <- atc_1_stats %>% 
  rowwise() %>% 
  mutate(fisher_PHENOTYPE = list(fisher.test(matrix(c(n_PHENOTYPE, 
                                                 n_iLincs-n_PHENOTYPE, 
                                                 Total_HS-n_PHENOTYPE, 
                                                 Total_iLincs-n_iLincs-(Total_PHENOTYPE-n_PHENOTYPE)), 
                                                 nrow = 2))),
         p_PHENOTYPE = fisher_PHENOTYPE$p.value,
         or_PHENOTYPE = fisher_PHENOTYPE$estimate[[1]]) %>% 
  ungroup() %>% 
  mutate(fdr_PHENOTYPE = p.adjust(p_PHENOTYPE, method = "BH")) %>% 
  select(-fisher_PHENOTYPE)
# Write down the results table 
write.table(atc_1_res, "path/to/iLincs/atc_1_enrichment_analysis.txt",
            row.names = F, col.names = T, quote = F, sep = "\t")
```

Then you can repeat this analysis but with ATC level 2 and 3.

``` r
# ATC level 2
df.atc_2 <- na.omit(df.atc) %>% 
  group_by(ATC_2) %>% 
  summarize(n = n()) %>% 
  mutate(perc = n / sum(n) * 100,
         trait = "PHENOTYPE")
# ATC level 3
df.atc_3 <- na.omit(df.atc) %>% 
  group_by(ATC_3) %>% 
  summarize(n = n()) %>% 
  mutate(perc = n / sum(n) * 100,
         trait = "PHENOTYPE")
```
