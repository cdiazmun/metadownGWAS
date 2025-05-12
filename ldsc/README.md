# Genetic correlation & heritability estimation

Genomic inflation (calculated previously using $\chi$<sup>2</sup> statistics) is highly influenced by polygenicity. Through linkage disequilibrium score (LDSC) regression, one can better estimate the true inflation and polygenic signal.  LDSC estimates can also be used to determine the genetic heritability of the trait under study (given by the slope of the LDSC regression) and, if two summary statistics are provided, the genetic correlation (cross-trait LDSC) between both traits. To know more about LD Score regressions please see:

- [LDSC](https://www.nature.com/articles/ng.3211)
- [Genetic correlations using LDSC](https://www.nature.com/articles/ng.3406)

We will be following the instructions from their [Wiki](https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation). Everything is very nicely explained there, so I won't detail much here. As a note, I installed ldsc as a `conda` environment, which was quite straightforward. So, we are going to analyse all the statistics based on pre-computed LD Scores of European populations (see the wiki to download the files). Because of that, we need the list of rsIDs for each of our variants to compare with the database. 

## Munge sumstats
The first thing to do is to reformat the summary statistics file to accomodate the `ldsc` format and to only keep markers for which we have LD information from the reference LD files. We need to make sure we don't have 2 columns with information about the SNPs ID, as is the case here. We have the `MarkerName` column with info about the `CHR:BP:ALT:REF` on the one hand and the `ID` column with the rsIDs on the other. We need to keep only the column with the rsIDs, otherwise a `ValueError: Found 2 different SNP columns` will occur. Finally, there is also the flag `--merge-alleles` to only keep only markers included in the reference HapMap3 list of SNPs (without the MHC). However, this option is recommended when there's little trust on the set of markers that have been imputed and therefore used for the analysis. In our case we can trust our markers. 

``` bash
munge_sumstats.py
  --sumstats path/to/your/input/meta.analysis.ancestry.sex.trait.txt.gz
  --signed-sumstats EFFECT,0 # Name of the column and indication that is the BETA (0) and not the OR (1)
  --N-col N_TOTAL
  --a1 EFFECT_ALL
  --a2 OTHER_ALL
  --p PVALUE
  --snp ID
  --frq EFF_ALL_FREQ
  --out path/to/your/munged_sumstats/meta.analysis.ancestry.sex.trait
```

The markers removed because of missing values are those that do not have an rsID (*e.g.,* markers not present in the 1000G reference files used for the rsID assignment). Just a few variants are removed by the minimum allele frequency (MAF) filter while many more are removed because of ambiguity/not SNPs (INDELs). Then, other variants, not many, are removed because they are not in a number of samples equal to the 90th percentile of total sample size / 2. If too many variants are removed in this step, you may consider to adapt the number using the corresponding flag. Also, if you already applied a filtering step after the meta-analysis to remove variants not present in xxx % of the total sample size, you also may consider to avoid this step, as it is not necessary anymore. 

In case you want to run the upper script for several files, you can do that with a simple `for` loop in your command line (you can apply this strategy to barely any program you want to run):

``` bash
for i in path/to/your/input/*.txt.gz;
  do bn=$(basename ${i%.txt.gz});
  munge_sumstats.py --sumstats $i --signed-sumstats EFFECT,0 --N-col N_TOTAL --a1 EFFECT_ALL --a2 OTHER_ALL --p PVALUE --snp ID --frq EFF_ALL_FREQ --out path/to/your/munged_sumstats/$bn;
done
```

## Heritability
Once we have the summary statistics formatted, we can initiate the estimation of the heritability (h<sup>2</sup><sub>SNP</sub>) of the trait. First, we need to download the reference files for the markers LD. The code snippet below corresponds to the EUR population. We would need to do the same for the EAS population. For any other population we would need to generate ourselves the LD matrix from the genotype data. There is a guide [here](https://github.com/bulik/ldsc/wiki/LD-Score-Estimation-Tutorial). 

``` bash
# LD scores:
 wget https://zenodo.org/records/7768714/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1

# Weights:
wget https://zenodo.org/records/7768714/files/1000G_Phase3_weights_hm3_no_MHC.tgz?download=1
```

Once we have the reference files correctly downloaded, we can initiate the main ldsc command:

``` bash
ldsc.py --h2 path/to/your/munged_sumstats/meta.analysis.ancestry.sex.trait.sumstats.gz 
	--ref-ld-chr /path/to/your/database/EUR/baselineLD/baselineLD.
	--w-ld-chr /path/to/your/database/EUR/baselineW/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
	--out heritability/meta.analysis.ancestry.sex.trait.h2
```

In the reference panel there are ~1M SNPs so it's normal to only retain that amount. You should check that (1) the estimated h<sup>2</sup><sub>SNP</sub> is as expected for your trait (in case you have an idea of that), the Lambda GC is not too high (ideally close to 1.0 and not going higher than ~1.2), the mean $\chi$<sup>2</sup> and intercepts are also within the renge and that he ratio is also within the "fine" range (below 20%). To know more about these ranges, see the guide I included before. 


## Genetic correlations
We can estimate the genetic relatedness between two or more traits by calculating the so-called genetic correlations (r<sub>g</sub>). The command is very similar to the one used to calculate the h<sup>2</sup><sub>SNP</sub> and the LD reference files are the same. In the example below we are calculating the r<sub>g</sub> between the three traits. It always calculates the r<sub>g</sub> of the first trait against all the others after the first comma. 

``` bash
ldsc.py
  --rg path/to/your/munged_sumstats/meta.analysis.ancestry.sex.trait1.sumstats.gz,path/to/your/munged_sumstats/meta.analysis.ancestry.sex.trait2.sumstats.gz,path/to/your/munged_sumstats/meta.analysis.ancestry.sex.trait3.sumstats.gz
  --ref-ld-chr /path/to/your/database/EUR/baselineLD/baselineLD.
  --w-ld-chr /path/to/your/database/EUR/baselineW/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
  --out genetic_correlations/trait1_vs_trait_2_3
```

We need to differentiate here between the `h2_obs` and the h<sup>2</sup><sub>SNP</sub> calculated using the `--h2` flag. These will not be the same since `h2_obs` it's calculated using only SNPs shared between the two traits you're calculating the genetic correlation for. It will always be better to calculate the Heritability of a given trait by itself. 

## Batch genetic correlations (website based)
Using the website-based service of [CTG-VL](https://vl.genoma.io/) we can calculate the genetic correlations of our trait of interest toward a set of more than 1,500 publicly available GWAS in their database. After running it, we can download the results and have a look at them. For reporting them, we typically do the following:
1. Filter for h<sup>2</sup><sub>SNP</sub> Z-score (`h2_obs`/`h2_se`) higher than 2, which would correspond to a significant heritability estimation of the trait (*p* < 0.05).
2. Calculate the False Discovery Rate (FDR). I normally do so in R with the `p.adjust(df$p-value, method = "fdr")` function.
3. Filter for FDR < 0.05.


