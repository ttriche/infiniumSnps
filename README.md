# infiniumSnps

SNP genotypes for Infinium methylation chips, mostly derived from GIAB calls for NA12878.    

The calls are identical on hm450 and hmEPIC (I checked, you can too: use the rgSet in minfiDataEPIC to get 3 replicates of GM12878 on EPIC; grab the ENCODE hm450 IDATs from [GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE40699&format=file) and compare them.  They're identical) so I concentrated on the 58 unique SNPs common to both platforms. 

```r
library(minfi) 
GM12878 <- makeGRangesFromDataFrame(
    read.csv(file="https://raw.githubusercontent.com/ttriche/infiniumSnps/master/GM12878.csv", row.names=1), 
  keep=TRUE)
message("Please link to https://github.com/ttriche/infiniumSnps in your methods section,")
message("if you use this data in a paper or to develop your own methods. Thank you. --t")
```
*Please link to https://github.com/ttriche/infiniumSnps in your methods section,    
if you use this data in a paper or to develop your own methods. Thank you. --t*

```r
head(GM12878, 2)
# GRanges object with 2 ranges and 6 metadata columns:
```

|ID|seqnames|ranges|strand|score|genotype|ref|alt|U|M|
|--|--------|------|------|-----|--------|---|---|-|-|
|rs3936238|chr1|[ 4031586,  4031586]|\*|3|A:A|A|G|G|A|
|rs877309|chr1|[11412265, 11412265]|\*|2|A:G|A|G|A|G|

```r
seqinfo: 18 sequences from an unspecified genome; no seqlengths
```

"Score" (above, after "strand") is the call from fitting a 3-component mixture model to the SNP beta values:

```r
library(minfiDataEPIC)
GM12878.EPIC.SNPs <- getSnpBeta(RGsetEPIC)

library(mclust)
mfit <- Mclust(as.numeric(GM12878.EPIC.SNPs), G=3)
mcalls <- t(apply(GM12878.EPIC.SNPs, 1, function(x) predict(mfit, x)$classification))
GM12878$score <- rowMeans(mcalls)[names(GM12878)] # they are all identical, of course, but it never hurts to check
```

Coordinates and ref/alt alleles for the SNPs are GRCh37p13 from SNPlocs.Hsapiens.dbSNP144.GRCh37 (some are multiallelic; in those cases I used the GM12878 genotype where available along with its score to determine what fluoresced as U or M).  In cases where GIAB has no well-supported call for a position, I left the $genotype set to NA. 

There are 18 rs probes common to EPIC and 450 where the base corresponding to the M allele and that of the U allele could not be determined, either because GM12878 is heterozygous at that position or because I don't have a genotype from GIAB at the SNP. The M and U columns for these SNPs are therefore marked as NA.  When I have more time I'll get back to these. 

The other 40 probes have high-confidence GM12878 genotypes from Genome In A Bottle (NIST) lined up next to their hm450/hmEPIC genotype calls.  Anyone who wants to figure out the remaining 18 is welcome to submit a pull request, although it would be helpful if you indicate how (i.e. "IMR90 genotype", "H1 genotype", whatever).  

If you have the SNP beta values in your GenomicRatioSet's metadata, for example, you can bolt the two together:

```r
library(minfiDataEPIC) 
source("https://raw.githubusercontent.com/ttriche/minfi/master/R/tcgaPipeline.R")
GM12878_EPIC_SNPs <- metadata(tcgaPipeline(RGsetEPIC))$SNPs[names(GM12878), ] 
mcols(GM12878) <- cbind(mcols(GM12878)[, c("M","U")], GM12878_EPIC_SNPs)
head(GM12878,2)
# GRanges object with 2 ranges and 5 metadata columns:
```

|       ID|seqnames|ranges|strand|U|M|X200144450018_R04C01|X200144450019_R07C01|X200144450021_R05C01|
|---------|--------|------|------|-|-|--------------------|--------------------|--------------------|
|rs3936238|    chr1|[ 4031586,  4031586]|\*|A|G|0.0534766974|0.0504425326|0.0607401706|
| rs877309|    chr1|[11412265, 11412265]|\*|G|A|0.5237327350|0.5314140191|0.5298953202|

```r
seqinfo: 18 sequences from an unspecified genome; no seqlengths
```
If you do things this way, i.e. bolt your actual SNP probe intensities and/or "calls" onto the side of the GRanges, you can perform subsetByOverlaps queries against e.g. GWAS data, or liftOver the sites to a different genome, or the usual stunts.  You don't have to, of course, but it may make your life easier, especially if "barcoding" samples by shared SNPs. 

At some point this will most likely go into minfi.  If you use it prior to then, kindly cite or acknowledge the source.
