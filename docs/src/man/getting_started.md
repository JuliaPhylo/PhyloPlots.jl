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
net = readnewick("(A,((B,#H1),(C,(D)#H1)));")
```
and call `plot`, as shown below.

```@example getting_started
R"svg"(figname("gettingstarted.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readnewick("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/gettingstarted.svg)

For the function's full documentation, see here: [`plot`](@ref).

!!! note "version history"
    - By default, v3 uses a `:majortree` style and `:both` hybrid edges curved.
      Instead, v2 used the `:fulltree` style by default,
      and had no option to curve edges.
    - The v0.3 syntax `plot(net, :R; ...)` still worked in v1 but was
      deprecated, and then removed in v2. For example, we could still use
      `plot(net, :R; showNodeNumber=true)` in v1, but in v2 (and later)
      we have to use instead `plot(net; shownodenumber=true)`.
    - Compared to v0.3, v1 does not support the Gadfly-based plots,
      and uses small-case-only argument names.
