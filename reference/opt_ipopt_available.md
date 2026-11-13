# Is the IPOPT backend available?

optcompR compiles an optional binding to the IPOPT C interface when
`pkg-config` finds IPOPT at installation (for example after
`brew install ipopt` or `apt install coinor-libipopt-dev`). Reinstall
optcompR after installing IPOPT to enable `solver = "ipopt"`.

## Usage

``` r
opt_ipopt_available()
```

## Value

`TRUE` when the IPOPT backend was compiled in.

## Examples

``` r
opt_ipopt_available()
#> [1] TRUE
```
