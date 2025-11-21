# Set plink binary location

Set plink binary location

## Usage

``` r
set_plink(path = "")
```

## Arguments

- path:

  If "" (default), then will use the MRCIEU/genetics.binaRies to get
  binaries that are appropriate for the detected operating system.
  Otherwise, provide the path to the plink binary. If NULL then will set
  the option to NULL.

## Value

NULL, sets option 'tools_plink'
