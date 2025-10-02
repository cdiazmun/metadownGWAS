# Local heritability and genetic correlations

You may have performed a genetic correlation analysis between your trait of interest and some clinically relevant traits or traits with a shared aetiology or symptomatology. Or you maybe did a batch genetic correlation analysis with everything that is out there. These analysis calculate the genetic similarity at whole-genome level, which is interesting to find connections between disorders and investigate a bit more on what is driving that connection. But one can't enfer from these analysis which is the exact mechanism or underlying genetic variant that is driving this correlation between the two traits. One option would be to do a colocalization analysis to identify common signals. But there is another way, which consists of calculating the heritability and genetic correlation between two traits in a "local" fashion by splitting the genome in ~2,000 pieces and testing these 2,000 independently. There are quite some tools to perform that:

- [MiXeR](https://www.nature.com/articles/s41467-019-10310-0). Published in 2019, one of the first ones for this.
- [LAVA](https://www.nature.com/articles/s41588-022-01017-y). Less restrictive than MiXeR, quite slow to run, but easy to interpret results.
- [HDL-L](https://www.nature.com/articles/s41588-025-02123-3). Similar to LAVA but stricter to reduce false positives.


