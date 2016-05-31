[![Build Status](https://travis-ci.org/JuliaInv/FWI.jl.svg?branch=master)](https://travis-ci.org/JuliaInv/BasicFWI.jl)
[![Coverage Status](https://coveralls.io/repos/github/JuliaInv/FWI.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaInv/BasicFWI.jl?branch=master)
[![Build status](https://ci.appveyor.com/api/projects/status/regr4ahc7xrtopox?svg=true)](https://ci.appveyor.com/project/lruthotto/basicfwi-jl)


# BasicFWI.jl - Basic Full Waveform Inversion in Julia.

# Requirements

This package is inteded to use with Julia versions 0.4.x.

This package is an add-on for [`jInv`](https://github.com/JuliaInv/jInv.jl), which needs to be installed. This is a basic FWI package, mostly used for teaching.

# Installation

In julia type:

``` 
Pkg.clone("https://github.com/JuliaInv/jInv.jl","jInv")
Pkg.clone("https://github.com/JuliaInv/BasicFWI.jl","BasicFWI")
Pkg.test("BasicFWI")
```
