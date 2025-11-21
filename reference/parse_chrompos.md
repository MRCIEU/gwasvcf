# Parse chromosome:position

Takes data frame or vector of chromosome position ranges and parses to
granges object

## Usage

``` r
parse_chrompos(chrompos, radius = NULL)
```

## Arguments

- chrompos:

  Either vector of chromosome and position ranges e.g. "1:1000" or
  "1:1000-2000", or data frame with columns `chrom`, `start`, `end`.

- radius:

  Add radius to the specified positions. Default = NULL

## Value

GRanges object
