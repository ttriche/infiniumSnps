# infiniumSnps

SNP genotypes for Infinium methylation chips, mostly derived from GIAB calls for NA12878.    

The calls are identical on hm450 and hmEPIC (I checked, you can too: use the rgSet in minfiDataEPIC to get 3 replicates of GM12878 on EPIC; grab the ENCODE hm450 IDATs from [GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40699&format=file) and compare them.  They're identical) so I concentrated on the 58 unique SNPs common to both platforms. 

"Score" is the call from fitting a three-component Gaussian mixture model to a matrix of SNP beta values:

```r
library(minfiDataEPIC)
GM12878.EPIC.SNPs <- getSnpBeta(RGsetEPIC)

library(mclust)
mfit <- Mclust(as.numeric(GM12878.EPIC.SNPs), G=3)
mcalls <- t(apply(GM12878.EPIC.SNPs, 1, function(x) predict(mfit, x)$classification))
GM12878$score <- rowMeans(mcalls)[names(GM12878)] # they are all identical, of course, but it never hurts to check
```

Coordinates and ref/alt alleles for the SNPs are GRCh37p13 from SNPlocs.Hsapiens.dbSNP144.GRCh37 (some are multiallelic; in those cases I used the GM12878 genotype where available along with its score to determine what fluoresced as U or M). 

There are 18 rs probes common to EPIC and 450 where the base corresponding to the M allele and that of the U allele could not be determined.  Their M and U columns are therefore marked as NA.  When I have more time I'll get back to these. 

The other 40 probes have high-confidence GM12878 genotypes from Genome In A Bottle (NIST) lined up next to their hm450/hmEPIC genotype calls.  Anyone who wants to figure out the remaining 18 is welcome to submit a pull request, although it would be helpful if you indicate how (i.e. "IMR90 genotype", "H1 genotype", whatever).  

At some point this will go into minfi.
