# Colocalization

If you want to understand the genetic relationship between two traits, you probably did before this a [genetic correlation](https://github.com/cdiazmun/metadownGWAS/ldsc/) analysis or maybe a bi-directional [Mendelian randomization](https://github.com/cdiazmun/metadownGWAS/mendelian-randomization/) analysis to check the causality between the two traits. If both analyses gave results that make you think that the two traits under study are somehow related (simply correlated, confounded, or causally-related) there is one more analysis one can perform to know more about the detailed mechanisms behind this relationship and that is the colocalization analysis. When one performs a MR study, one makes use of IVs, which are nothing else than genetic variants associated to the exposure, the outcome or both. This already gives you a hint that these associated variants, and therefore associated loci, may shared a common causal effect on both traits. Another possibility is that the same locus is associated to the two different traits but the causal variant, and therefore the mechanism behind the association, is different. The colocalization analysis will calculate the probabilities of a given locus to be associated to trait1 (H1) to trait2 (H2) or to both of them with different causal variants (H3) or the same causal variant (H4).

To know more about colocalization and the R package we use for this analysis you can check [the original paper](https://doi.org/10.1371/journal.pgen.1004383), the [GitHub](https://github.com/chr1swallace/coloc) or the [vignettes](https://chr1swallace.github.io/coloc/index.html).

## Pipeline

*TBD - coming soon*
