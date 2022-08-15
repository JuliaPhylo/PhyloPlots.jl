```@setup getting_started
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Getting started

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
and call `plot`, as shown below.

```@example getting_started
R"svg"(figname("gettingstarted.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/gettingstarted.svg)

For the function's full documentation, see here: [`plot`](@ref).

!!! note "version history"
    Compared to v0.3, v1 does not support the Gadfly-based plots,
    and uses small-case-only argument names.

    The v0.3 syntax `plot(net, :R; ...)` still works in v1.0 but is
    **deprecated**, and will be removed in a future release. For example,
    you can still use `plot(net, :R; showNodeNumber=true)`,
    but you should instead start using `plot(net; shownodenumber=true)`.
