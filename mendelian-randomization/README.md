# Mendelian Randomization

Mendelian randomization (MR) is a method in which genetic variants (SNPs) are used as instrumental variables (IVs) to infer causality between two traits. The IVs need to be strongly associated with trait A (exposure) and associated with trait B (outcome) only through their effects on trait A. Because of that, a key step in MR analysis is the selection of valid or strong IVs (*i.e.,* removal of weak IVs). For this, there are a number of approaches, from LD clumping to avoid including variants in LD with another variant that influences the outcome through a different mechanism, to horizontal pleiotropy (when a gene/genetic variant affects two different phenotypes independently) estimation and removal of outliers. 

Genetic variants are excellent instruments because they are not affected by confounders like behavioral/environmental factors, the direction is always from SNP to the phenotype/trait, there is little measurement error or bias, if the desired SNP is not available a marker in LD is enough to perform the analysis, and because of the wide availability of genetic data (GWAS). In this case we will try to see if genome-wide significan SNPs (IVs) associated to trait1 (exposure) affect a variety of diseases (outcomes; *e.g.,* trait2, trait3, trait4). 

## TwoSampleMR

We will use a widely used R package to perform MR, named [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/index.html). If you click there, you'll find the information on how to install it. This package uses the [IEU GWAS database](https://gwas.mrcieu.ac.uk/) to obtain a myriad of GWAS, eQTLs, metabolites, and protein levels. It contains GWAS data from UKB, FinnGen, Biobank of Japan, EBI, and IEU (~50K GWAS). The list of resources can be found at the end of this section. 

Load the libraries.
``` r
library(TwoSampleMR) # To perform MR
library(ieugwasr) # To access the IEU publicly available GWAS (to be used as exposures or outcomes)
library(data.table) # To read the summary stats file 
library(dplyr) # To manipulate data
```
Set the working directory and read the sumstats.
``` r
setwd("/path/to/your/working/directory/")
sumstats <- fread("/path/to/your/meta.analysis.ancestry.sex.trait.txt.gz")
```
### Exposures
The first thing we need to do is to select the IVs from the sumstats by performing LD clumping. We are first going to filter for genome-wide significant SNPs so that the `ld_clump()` function has less work to do (it improves in orders of magnitude the speed of this step). We are going to set the r<sup>2</sup> to 0.001, the window size to 10 Mbp and the *p*-value threshold to 5x10<sup>-08</sup>. We also need to select a reference population from the 1000G reference panel.
``` r
high.snps <- sumstats %>% 
  filter(PVALUE <= 5e-8)

high.clump <- ld_clump(dplyr::tibble(rsid=high.snps$ID, pval=high.snps$PVALUE, id="trait1"),
                        clump_r2 = 0.001, clump_kb = 10000, clump_p = 5e-8, pop = "EUR")
```
We then need to accomodate the data frame to what the package expects as exposure object.
``` r
exposure <- high.clump %>% 
  left_join(high.snps, by = c("rsid" = "ID")) %>% 
  select(-PVALUE, -SNP, -LNBF) %>% 
  rename(pval.exposure = pval, samplesize.exposure = N_TOTAL, 
         chr.exposure = CHR, se.exposure = STDERR, beta.exposure = EFFECT,
         pos.exposure = POSITION, id.exposure = id, SNP = rsid, 
         effect_allele.exposure = EFFECT_ALL, other_allele.exposure = OTHER_ALL,
         eaf.exposure = EFF_ALL_FREQ) %>% 
  mutate(exposure = "Trait1", mr_keep.exposure = TRUE, 
         pval_origin.exposure = "reported", data_source.exposure = "in-house")
```
Once we have the exposure set with a list of independent IVs we are going to calculate the F-statistic for each variant to test the strength of the association between the IV and the exposure. Hence, we are going to remove weak variants that can bias our analysis on the outcome. As a general rule, an F-statistic > 10 indicates that the level of weak instrument bias is likely to be small.
``` r
exposure_purged <- exposure %>% 
  mutate(F = beta.exposure^2/se.exposure^2) %>% 
  filter(F >= 10)
```
### Outcomes
To select the outcomes, we are going to examine all available outcomes and select the ones we are interested on testing. If you want to do a pre-selection in case you have too many related traits, there is a column named priority that you can use to get a selection of them. 
```r
ao <- available_outcomes()
ao1 <- ao %>% filter(priority == 1)

# Here you should imagine that "trait" is something like "major depresive disorder" so that you look for that using grepl()
trait2 <- subset(ao$id, grepl("trait2", ao$trait, ignore.case = TRUE))
trait3 <- subset(ao$id, grepl("trait3", ao$trait, ignore.case = TRUE))

# And so on and so on with the list of diseases/traits you want to use as outcome
outcome_list <- c(trait2, trait3, ...)
```
Now that we have a list of outcomes IDs, we can extract the IVs from the exposure on the outcomes.
``` r
outcome <- extract_outcome_data(snps=exposure_purged$SNP, outcomes = outcome_list)
```
If we wish to import a GWAS that is not in the list of available outcomes, we can do so as well.
``` r
outside_GWAS <- read_outcome_data(
  snps = exposure_purged$SNP,
  filename = "path/to/the/outside/GWAS/sumstats.trait4.txt.gz",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf_col = "A1FREQ",
  pval_col = "P")

outside_GWAS <- outside_GWAS %>% 
  mutate(outcome = "Trait4", samplesize.outcome = XXX) %>% 
  mutate(id.outcome = "Author et al., XXX") %>% 
  select(-pval_origin.outcome)

# In case that we want to join outcomes from IEU and impoted ones
cols_to_keep <- colnames(outside_GWAS)
outcome_dat <- rbind(outcome[, cols_to_keep], outside_GWAS)
```

### Harmonise and perform MR
Next step is to harmoise the data, to be sure that the alleles are in the same direction. This means that the effect of a SNP on the exposure and the effect of that SNP on the outcome must each correspond to the same allele. Palindromic variants can give problems here, but there are [ways to solve it if needed](https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html)
``` r
dat <- harmonise_data(exposure_purged, outcome_dat)
```
Now we are ready to perform MR on this list of IVs
``` r
res <- mr(dat)
res
```
If you explore the data frame, you will see different methods that have been performed: Inverse variance weighted (IVW), MR Egger, Simple mode, Weighted mode, and Weighted median. 

- **Inverse variance weighted (IVW):**  It assumes the inexistence of horizontal pleiotropy, so it assumes that all IVs are valid. It uses the Wald ratio, so the ratio between the beta describing the association of the variant and the outcome by the beta describing the association of the variant and the exposure. These calculated ratios are weighted by the inverse of the variance of the association between the variant and the outcome. This is done to give more weight to stronger instruments (higher precision).
- **Mode-based estimate:** The mode-based methods assume that the most common causal estimate across the genetic variants is the true causal effect. In other words, if most of your genetic variants are valid instruments (*i.e.,* they satisfy the instrumental variable assumptions), then their estimates will "cluster" around the correct value. It is robust to the presence of invalid instruments (*e.g.,* variants affected by horizontal pleiotropy) as long as the majority of valid instruments cluster around the true effect. Can be weighted (by IVW) or unweighted.
  - **Simple mode:** Assumes all variants contribute equally, regardless of instrument strength.
  - **Weighted mode:** Gives more weight to stronger instruments (variants with higher precision).
- **Median-based estimate:** The weighted median method assumes that at least 50% of the information comes from valid instruments. It’s robust to outliers or invalid instruments, as it focuses on the middle of the distribution rather than the average. Even if some variants are invalid, as long as most (weighted) information comes from valid instruments, the median estimate will be reliable.
- **MR Egger:** A method developed for two-sample MR settings that combines Wald ratio together into a meta-regression (with an intercept and slope parameter) to estimate the causal effect adjusted for any directional pleiotropy. Unlike standard MR methods like IVW, MR-Egger doesn’t assume that all genetic variants are perfectly valid instruments. The intercept of a MR-Egger regression provides an indication of horizontal pleiotropy when it is not null.

Here I make some `ggplots` to check how the overall results look like and decide on the next steps. Since we can select several GWAS for a single trait, we would want to check if there is consistency between the outcomes in the direction of causality and the significance. Ideally, we would have an outcome that concludes that trait A causes trait B with statistical significance in as many methods as possible. 

### Sensitivity analysis
When we have a desired outcome to report, we also need to check a few things that will add robustness to the analysis and confidence in the interpretation of the results. The first thing we will do is to check the heterogeneity through Cochran's Q statistic. If the Q statistic much larger than its degrees of freedom (*i.e.,* the number of IVs minus 1), this provides evidence for heterogeneity and invalid IVs.
``` r
het <- mr_heterogeneity(dat)
het
```
The horizontal pleiotropy can be checked by inspecting the egger intercept as mentioned before in the **MR Egger** method description.
``` r
pleio <- mr_pleiotropy_test(dat)
pleio
```
Another test that we can do to add more evidence in the results is Mendelian Randomization Pleiotropy RESidual Sum and Outlier (MR-PRESSO). MR-PRESSO is an extension of the IVW method, which attempts to perform outlier removal based on their contributions to heterogeneity, which is achieved through a simulation approach. It gives you a *p*-value for the MR test with the "raw" data and after the "removal of outliers". By inspecting the R object, we can see which variants contribute the most to the heterogeneity and decide to remove them from the analysis or not. 
``` r
library(MRPRESSO)
presso <- run_mr_presso(dat, NbDistribution = 1000)
presso[[1]]
```
In MR it is assumed that the instruments influence the exposure first and then the outcome through the exposure. But sometimes this is difficult to evaluate. The causal direction between the hypothesised exposure and outcomes can be tested using the Steiger test. It calculates the variance explained in the exposure and the outcome by the instrumenting SNPs, and tests if the variance in the outcome is less than the exposure.
``` r
out <- directionality_test(dat)
out
```
To obtain the MR estimates using each of the SNPs singly we can do a single SNP analysis. The method used by default is the Wald ratio. We can use this information to detect which variant contributes the most to the overall analysis and detect the presence of outliers if the effect estimates differ a lot from all others, in one way or the other. 
``` r
res_single <- mr_singlesnp(dat)
```
A complementary analysis for the single SNP MR is the leave-one-out analysis. In this case, MR is performed again but leaving out each SNP in turn, to identify if a single SNP is driving the association. By default the method used is the inverse variance weighted method.
``` r
res_loo <- mr_leaveoneout(dat)
```
Finally, we can choose to generate a report in which all MR analyses, sensisitivty analyses, and plots are generated. This function will create a separate report file for every exposure-outcome combination that is present in the `dat` object.
``` r
mr_report(dat, output_path = xxx, output type = "html" or "md")
```

## Resources

Guides: https://mrcieu.github.io/TwoSampleMR/articles/introduction.html

Manual: https://mrcieu.r-universe.dev/TwoSampleMR/doc/manual.html

Tutorial (from user): https://andrewslabucsf.github.io/MR-tutorial/scripts/mr.html

MR Dictionary: https://mr-dictionary.mrcieu.ac.uk/


## References

Davey Smith & Shah, 2003 - [‘Mendelian randomization’: can genetic  epidemiology contribute to understanding  environmental determinants of disease?](https://doi.org/10.1093/ije/dyg070)

Davey Smith & Hemani, 2014 - [Mendelian randomization: genetic anchors for causal inference in epidemiological studies](https://doi.org/10.1093/hmg/ddu328)

Bowden *et al*., 2015 - [Mendelian randomization with invalid  instruments: effect estimation and bias  detection through Egger regression](10.1093/ije/dyv080)

Hemani *et al*., 2017 - [Orienting the causal relationship between  imprecisely measured traits using GWAS  summary data](https://doi.org/10.1371/journal.pgen.1007081)
