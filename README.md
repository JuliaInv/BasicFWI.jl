[![Build Status](https://travis-ci.org/JuliaInv/FWI.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/FWI.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/FWI.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/FWI.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/regr4ahc7xrtopox?svg=true)](https://ci.appveyor.com/project/lruthotto/fwi-jl)


# FWI.jl - Full Waveform Inversion in Julia

# Requirements

This package is inteded to use with Julia versions 0.4.x.

This package is an add-on for [`jInv`](https://github.com/JuliaInv/jInv.jl), which needs to be installed. 

# Installation

In julia type:

``` 
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/FWI.jl","FWI")
Pkg.test("FWI")
```
