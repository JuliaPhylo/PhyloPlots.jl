```@setup installation
using PhyloPlots, RCall
```

# Installation

This will assume you have installed Julia and PhyloNetworks. For information on how to 
install them, see
[here](https://crsl4.github.io/PhyloNetworks.jl/dev/man/installation/#Installation)

## Installation of the package PhyloPlots

PhyloPlots depends on PhyloNetworks, and has further dependencies
like [Gadfly](http://gadflyjl.org/stable/) and
[RCall](https://github.com/JuliaInterop/RCall.jl)

To install in the Julia REPL, enter package mode with `]`, and:

```
add PhyloPlots
```
Or in julian mode:

```@repl installation
using Pkg
Pkg.add("PhyloPlots")
```
