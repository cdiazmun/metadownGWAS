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

**Type:** Gene prioritization and pathway enrichment

**Input:** GWAS summary statistics. 

**Method:** Predicts gene function based on co-regulation and expression data across tissues. 

**Output:** Prioritized genes, enriched pathways, and tissues. 

## SMultiXcan
[Integrating predicted transcriptome from  multiple tissues improves association  detection](https://dx.plos.org/10.1371/journal.pgen.1007889)

[Tool website]()

**Type:** Transcriptome-wide association study (TWAS). Gene level. 

**Input:** GWAS summary statistics + eQTL data (*e.g.,* GTEx)

**Method:** Predicts gene expression across tissues and correlates with phenotype. 

**Output:** Genes whose predicted expression is associated with the trait. 

## SMR
[Integration of summary data from GWAS and eQTL studies predicts complex trait gene targets](https://www.nature.com/articles/ng.3538)

[Tool website]()

**Type:** Transcriptome-wide association study (TWAS). Variant level. 

**Input:** GWAS summary statistics + eQTL data (*e.g.,* GTEx)

**Method:** Tests whether gene expression mediates GWAS associations. 

**Output:** Genes with evidence of causality via expression.

