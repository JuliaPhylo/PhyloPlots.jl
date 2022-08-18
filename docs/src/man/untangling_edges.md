```@setup untangling
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Untangling the network

This plot may not be the easiest to read, as the hybrid edge crosses over C's
edge:

![example1](../assets/figures/gettingstarted.svg)

To fix this, we can to rotate C and D's edges around their parent node.

First we need to know the number of this parent node. By showing node numbers
with the `shownodenumbers = true` option, we can find the number of the node
whose child edges we should rotate.

```@example untangling
R"svg"(figname("untangling1.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, shownodenumber=true);
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/untangling1.svg)

As we can see, rotating edges around node `-5` will make for a prettier network.

```@example untangling
R"svg"(figname("untangling2.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
rotate!(net, -5)
plot(net)
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/untangling2.svg)


This may seem unnecessary for a small network as shown, but it is a useful tool for plotting
large networks.
