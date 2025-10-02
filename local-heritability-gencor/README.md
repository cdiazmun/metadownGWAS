# Local heritability and genetic correlations

You may have performed a genetic correlation analysis between your trait of interest and some clinically relevant traits or traits with a shared aetiology or symptomatology. Or you maybe did a batch genetic correlation analysis with everything that is out there. These analysis calculate the genetic similarity at whole-genome level, which is interesting to find connections between disorders and investigate a bit more on what is driving that connection. But one can't enfer from these analysis which is the exact mechanism or underlying genetic variant that is driving this correlation between the two traits. One option would be to do a colocalization analysis to identify common signals. But there is another way, which consists of calculating the heritability and genetic correlation between two traits in a "local" fashion by splitting the genome in ~2,000 pieces and testing these 2,000 independently. There are quite some tools to perform that, but we will be using these ones:

- [MiXeR](https://www.nature.com/articles/s41467-019-10310-0). 2019. Very visual interpretation of the results. MiXeR is more focused on trait architecture rather than direct local correlation, but can complement other tools in broader analyses.
- [LAVA](https://www.nature.com/articles/s41588-022-01017-y). 2022. Less restrictive than MiXeR, quite slow to run, but easy to interpret results.
- [HDL-L](https://www.nature.com/articles/s41588-025-02123-3). 2025. Similar to LAVA but stricter to reduce false positives. HDL-L is emerging as a strong alternative to LAVA, especially for large-scale analyses due to its computational efficiency and robustness.

## Comparison of tools for local genetic correlation analysis
LD estimation is a critical factor across all methods. Poor LD reference panels can lead to inflated type I error rates or biased estimates.

| Tool         | Strengths                                                                 | Limitations                                                                 |
|--------------|---------------------------------------------------------------------------|------------------------------------------------------------------------------|
| **LAVA**     | - Models bivariate and multivariate local genetic correlations<br>- Accounts for sample overlap<br>- Offers conditional modeling (partial correlation, regression)<br>- Widely used and documented | - Sensitive to LD reference panel quality<br>- Slower than HDL-L<br>- May yield false positives in some settings |
| **HDL-L**    | - High accuracy and consistency<br>- Operates on semi-independent LD blocks<br>- ~50x faster than LAVA<br>- Robust across diverse traits | - Newer tool with less widespread adoption<br>- Requires careful LD block definition |
| **MiXeR**    | - Estimates polygenicity and discoverability<br>- Models shared and unique genetic architecture<br>- Useful for global and local trait architecture | - Not primarily designed for local correlation<br>- Less granular than LAVA or HDL-L |
| **ρ-HESS**   | - Simple and interpretable<br>- Early method using summary statistics | - Sensitive to LD estimation errors<br>- Limited scalability and flexibility |
| **SUPERGNOVA** | - Improved LD modeling over ρ-HESS<br>- Handles sample overlap better<br>- Strong performance in simulations | - Sensitive to LD matrix estimation<br>- May yield inconsistent results across traits or regions |

## MiXeR
This tool works really as a black box. The "only" thing you need to do is to install it following these [instructions](https://github.com/precimed/mixer?tab=readme-ov-file#install), download the required [references](https://github.com/comorment/mixer/blob/main/README.md), and reformat the summary statistics [like this](https://github.com/precimed/mixer?tab=readme-ov-file#gwas-summary-statistics-format). It is all build as a docker container or singularity. This was my first time dealing with Docker soit took me some time, but it is actually quite straightforward.

Docker was already installed in the system, so the first thing I had to do was to give my user permission to run docker:
``` bash
sudo usermod -aG docker $USER
```
Once I can use it, we need to download the Docker version of MiXeR:
``` bash
docker pull ghcr.io/precimed/gsa-mixer:2.2.1
```
I can the set up the environment by running `source setup_mixer_env.sh`. This script will connect (by mounting) the docker container and my working and reference directory so that everything is accessible. 
``` bash
export DOCKER_RUN="docker run -v $PWD:/path/to/workdir -v /path/to/databases/mixer:/path/to/databases/mixer -w /path/to/workdir"
export MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:latest python /tools/mixer/precimed/mixer.py"
```
The last thing we need to do is to set the arguments, also as environment variables:
``` bash
export MIXER_COMMON_ARGS="--ld-file /path/to/databases/mixer/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld --bim-file /path/to/databases/mixer/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim --threads 16"
```
Once everything is setup and worked good, then we can run the mixer script for cross-trait analysis:
``` bash
${MIXER_PY} fit1 $MIXER_COMMON_ARGS --trait1-file input_sumstats/trait1.mixer.sumstats.gz --out output/trait1.fit
${MIXER_PY} fit1 $MIXER_COMMON_ARGS --trait1-file input_sumstats/trait2.mixer.sumstats.gz --out output/trait2.fit

```

