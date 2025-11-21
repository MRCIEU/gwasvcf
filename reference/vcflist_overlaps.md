# Reduce list of VCFs to intersecting regions

Reduce list of VCFs to intersecting regions

## Usage

``` r
vcflist_overlaps(vcflist, chrompos)
```

## Arguments

- vcflist:

  List of VCF objects, or list of VCF filenames, or mix of VCF objects
  and filenames

- chrompos:

  Either vector of chromosome and position ranges e.g. "1:1000" or
  "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.

## Value

List of VCFs
