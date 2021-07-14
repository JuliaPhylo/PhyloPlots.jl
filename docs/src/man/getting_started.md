```@setup getting_started
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Getting started

## Plotting
To demonstrate, we will plot the simple network: `(A,((B,#H1),(C,(D)#H1)));`

To start plotting, use the packages:

```@repl getting_started
using PhyloNetworks
using PhyloPlots
```
Then read the topology
```@repl getting_started
net = readTopology("(A,((B,#H1),(C,(D)#H1)));")
```
and call `plot`, using :R for full funtionality:

!!! note "About the two versions of plot()" 
    One is made using Gadfly (`plot(net)`), the other using RCall (`plot(net, :R)`).
    Since only the RCall version supports icytree style trees and 
    the `style = :fulltree` parameter, we will only use the RCall.

(To the funtion's full documentation, see here: [`plot`](@ref))

This will draw the following plot.

```@example getting_started
R"svg"(figname("gettingstarted.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/gettingstarted.svg)
