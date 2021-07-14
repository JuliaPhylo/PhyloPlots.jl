```@setup untangling
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Untangling the network

## Rotating edges to untangle them

This plot may not be the easiest to read, as the hybrid edge crosses over C's 
edge:

![example1](../assets/figures/gettingstarted.svg)

To fix this, we can to rotate C and D around their parent node. 

First we need to know the number of this parent node. By showing node numbers
with the `showNodeNumbers = true` parameter, we can find the number of the node 
we should rotate.

```@example untangling
R"svg"(figname("untangling1.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, showNodeNumber=true);
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/untangling1.svg)

As we can see, rotating node `-5` will make for a prettier network.

```@example untangling
R"svg"(figname("untangling2.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
rotate!(net, -5)
plot(net, :R)
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/untangling2.svg)


This may seem unnecesary for a small network as shown, but it is a usefull tool for plotting 
large networks.
