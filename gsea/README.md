
# Gene Set Enrichment Analysis (GSEA)

This section has been elaborated between Francisco Heredia (a former collaborator) and me.

This kind of analysis consists on taking the genes that mapped to the significant variants in the GWAS (*via* eQTL, physically, or chromatin interaction mapping) and using them as a query to look for enriched biological functions and/or pathways. Therefore, the list of genes we start with will depend on the defined locus flanking region, proxies definition, and the databases and parameters defined by FUMA during the [functional annotation](https://github.com/cdiazmun/metadownGWAS/tree/main/functional-annotation) step. Furthermore, we can select to run GSEA using only protein coding genes, all type of genes mapped by FUMA and/or mapped genes _via_ specific mapping conditions. The selection of genes is usually slected case-by-case. With the list of genes we can then use a wide range of tools that rely their annotation on different databases and that consist of different algorithms to define enriched functions and/or pathways. These are the most widely used:

- [enrichR](https://maayanlab.cloud/Enrichr/)
- [Gene Network](https://genenetwork.nl/)
- [PASCAL](https://www2.unil.ch/cbg/index.php?title=Pascal)
- [gProfiler](https://biit.cs.ut.ee/gprofiler/gost)
- [GARFIELD](https://www.ebi.ac.uk/birney-srv/GARFIELD/)
- [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html)

## Enrichr

The following section is based on a very nice article that details the different analyses that can be performed in Enrichr: https://doi.org/10.1002/cpz1.90

Enrichr is a popular gene set search engine that offers a large collection of gene sets and libraries for conducting GSEA. It contains approximately 400,000 annotated gene sets organized into around 300 libraries, covering various biological functions such as pathways, diseases, and transcription factor regulation. To use Enrichr, you can visit its website (https://maayanlab.cloud/Enrichr/), where you can upload lists of human or mouse genes. Enrichr accepts genes from H. sapiens or M. musculus encoded in the Entrez Gene Symbol or the HGNC Gene Symbol format. Input files must be in the BED file format or TXT file format (genes must be on their own row). These genes are compared against numerous gene set libraries of known biological function, such as pathways, diseases, or gene sets regulated by transcription factors. Matched gene sets are ranked using various methods to measure how well your input genes match the gene sets in each library.

The GSEA results are initially displayed as a grid of bar charts, showing a summary of the top enriched terms for each gene set library in each category. You can select different categories of gene set libraries to explore. We usually focus on the pathways (Reactome and KEGG) and ontologies (GO Biological Process, GO Molecular Function, Human Phenotype Ontology and GO Cellular component) categories. By clicking the "Table"

option, you can view the results in a downloadable table format. Enrichr also provides several visualization options, including Network and Clustergram views, though these may not be available for all gene set libraries.


## GeneNetwork

The GeneNetwork tool (https://www.genenetwork.nl/) is based on a comprehensive analysis of large-scale RNA-seq data to predict gene functions and their potential associations with disease. Using data from over 31,000 RNA-seq samples across multiple tissues and cell types, GeneNetwork identifies patterns of co-regulation between genes. This helps to prioritize genes that may be associated with specific phenotypes, even if they lack prior biological annotations.

The method is based on the concept that genes that cause similar diseases often share similar molecular functions or are involved in related biological processes. By analyzing patterns of gene co-expression, GeneNetwork assigns scores to genes, ranking them based on their likelihood of contributing to a given phenotype. GeneNetwork integrates Human Phenotype Ontology (HPO) terms describing patient symptoms to refine its predictions. This allows it to identify subsets (clusters) of genes with similar functions from the user-provided list that can be analyzed independently.


## PASCAL

Pascal (Pathway Scoring Algorithm) is a computational tool designed to analyze GWAS data and extract biological insights by calculating gene and pathway scores from SNP-based summary statistics. It uses statistical methods to integrate SNP-level association data and map them to genes and pathways.

The tool uses two main gene scoring methods: the sum of association signals across SNPs and the maximum signal within a gene region. For pathways, Pascal uses a parameter-free scoring method that avoids arbitrary thresholds, ensuring that all contributing genes are considered. This is complemented by a strategy that accounts for linkage disequilibrium (LD) to avoid inflated scores from correlated genes.


## g:Profiler

g:Profiler is a web-based tool (https://biit.cs.ut.ee/gprofiler/gost) designed for GSEA that integrates up-to-date data from multiple sources, including Gene Ontology (GO), KEGG, Reactome, and Human Phenotype Ontology (HPO).

The core functionality, g:GOSt, maps genes to known functional information sources and detects statistically significantly enriched terms by applying multiple testing corrections.

## GARFIELD

GARFIELD (GWAS Analysis of Regulatory or Functional Information Enrichment with LD correction) is a computational tool designed to integrate GWAS findings with regulatory and functional genomic annotations to find features relevant to a phenotype of interest.

The workflow begins with GWAS summary statistics and genomic annotation data as inputs. It then prunes SNPs with LD r2 > 0.1 and then annotates them based on functional information overlap. It then quantifies enrichment using odds ratios (OR) at different GWAS p-value cutoffs and assesses their significance using generalized linear model testing, taking into account minor allele frequency, distance to the nearest transcription start site, and number of LD proxies (r2 > 0.8). Within this framework, it also accounts for various sources of confounding.

GARFIELD provides a way to assesess the enrichment of association analysis signals in 1005 features extracted from ENCODE, GENCODE and Roadmap Epigenomics projects, including genic annotations, chromatin states, histone modifications, DNaseI hypersensitive sites and transcription factor binding sites, among others, in several publicly available cell lines.


## clusterProfiler

clusterProfiler is an R-based tool designed to discover biological insights from high-throughput omics data using functional enrichment analysis. It contains functional features of coding and non-coding genomic data from thousands of species with updated gene annotations. The tool performs both overrepresentation analysis and GSEA using databases such as Gene Ontology (GO), KEGG pathways, Reactome, Disease Ontology and WikiPathways. It also simplifies the interpretation of GO results by eliminating redundancy.

clusterProfiler allows the creation of different plots using ggplot2. All this and more is well explained in this paper (https://doi.org/10.1016/j.xinn.2021.100141) and in its bioconductor manual (https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf).
