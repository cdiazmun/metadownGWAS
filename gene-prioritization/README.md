# Gene prioritization
Similar to the [functional annotation](https://github.com/cdiazmun/metadownGWAS/tree/main/functional-annotation) analysis, gene prioritization is performed to add biological meaning to the results obtained and try to gather as much evidence as possible of the mechanistic pathways that may be behind the genetic variant-trait association. This is done by prioritizing genes, tissues, and pathways. The ultimate objective is to find specific genes affected by the genetic variants associated to our trait so that its expression is deregulated in a tissue of interest (brain for psychiatric disorders, artieries for cardiovascular disorders, GI tract for digestive disorders, *etc*). There are different tools to perform this, which apply different methodologies and data sources to assess gene-level associations, functional relevance and regulatory mechanisms. 

## MAGMA
[MAGMA: Generalized Gene-Set Analysis of  GWAS Data](https://dx.plos.org/10.1371/journal.pcbi.1004219)

[CTGLAB](https://cncr.nl/research/magma/)

**Type:** Gene-based association test. 

**Input:** GWAS summary statistics. 

**Method:** Multiple linear regression models to account for LD. F-tests for gene *p*-values. 

**Output:** Genes significantly associated with the trait. 

## fastBAT
[Fast set-based association analysis  using summary data from GWAS  identifies novel gene loci for human  complex traits](https://www.nature.com/articles/srep32894)

[YangLab](https://yanglab.westlake.edu.cn/software/gcta/#fastBAT)

**Type:** Gene-based association test.

**Input:** GWAS summary statistics. 

**Method:** Calculates *p*-values from an approximated distribution of the $\chi$<sup>2</sup> statistics.

**Output:** Genes significantly associated with the trait. 

## DEPICT
[Biological interpretation of genome-wide  association studies using predicted gene functions](https://www.nature.com/articles/ncomms6890)

[GitHub](https://github.com/perslab/depict)

**Type:** Gene prioritization, pathway analysis and tissue/cell type enrichment analysis.

**Input:** GWAS summary statistics. 

**Method:** Predicts gene function based on co-regulation and expression data across tissues. It assigns each gene a probability to be part of a (previously reconstituted) gene set which is characterized by a certain functionality. Using these precomputed gene functions and a set of trait-associated loci, DEPICT assesses whether any of the 14,461 reconstituted gene sets are significantly enriched for genes in the associated loci, and prioritizes fenes that share predictd functions with genes from the other associated loci more often than expected by chance. In addition, DEPICT utilizes a set of 37,427 human microarrays to identify tissue/cell types in which genes from associated loci are highly expressed. 

**Output:** Prioritized genes, enriched pathways, and tissues. 

## SMultiXcan
[Integrating predicted transcriptome from  multiple tissues improves association  detection](https://dx.plos.org/10.1371/journal.pgen.1007889)

[GitHub](https://github.com/hakyimlab/MetaXcan)

**Type:** Transcriptome-wide association study (TWAS). Gene level. 

**Input:** GWAS summary statistics + eQTL data (*e.g.,* GTEx)

**Method:** Predicts gene expression across tissues and correlates with phenotype. PrediXcan tests the hypothesis that genetic variants affect the phenotypes through the regulation of gene expression traits. To do that, it correlated genetically predicted gene expression and the phenotype with the idea that causal genes are likely to show a significant association. Linear prediction models of expression using genetic variation in the vicinity of the gene are trained in reference transcriptome datasets such as GTEx. MultiXcan tests the joint effects of gene expression variation from different tissues. 

**Output:** Genes whose predicted expression is associated with the trait. 

## SMR
[Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets](https://www.nature.com/articles/ng.3538)

[YangLab](https://yanglab.westlake.edu.cn/software/smr/)

**Type:** Transcriptome-wide association study (TWAS). Variant level. 

**Input:** GWAS summary statistics + eQTL data (*e.g.,* GTEx)

**Method:** Tests whether gene expression mediates GWAS associations. Under the assumption of either causality (where the effect of a genetic variant on a trait is mediated by gene expression) or pleiotropy (where a genetic variant has direct effects on both a trait and gene expression), SMR is equivalent to MR analysis if genotype, gene expression and phenotype data are available from the same sample and that the power of SMR an be increased by orders of magnitude if *b<sub>zx<sub>* (effect of genetic variant on gene expression) and *b<sub>zy<sub>* (effect of genetic variant on the trait) re estimated separately from two independent samples with very large sample sizes. SMR is unable to distinguish between causality and pleiotropy regardless of whether the effect of *z* (genetic variant) on *x* (gene expression) is direct or mediated by a latent variable. 

**Output:** Genes with evidence of causality via expression.

