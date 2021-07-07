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
R"svg"(figname("gettingstarted1.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/gettingstarted1.svg)

## Rotating networks for more readable plots

This plot may not be the easiest to read, as the hybrid edge crosses over C's 
edge. To fix this, we can to rotate C and D around their parent node. 

First we need to know the number of this parent node. By showing node numbers
with the `showNodeNumbers = true` parameter, we can find the number of the node 
we should rotate.

```@example getting_started
R"svg"(figname("gettingstarted2.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, showNodeNumber=true);
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/gettingstarted2.svg)

As we can see, rotating node `-5` will make for a prettier network.

```@example getting_started
R"svg"(figname("gettingstarted3.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
rotate!(net, -5)
plot(net, :R)
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/gettingstarted3.svg)


This may seem unnecesary for a small network as shown, but it is a usefull tool for plotting 
large networks.

## Adding labels

!!! note
    For demonstration purposes, I will walk through the process of adding labels to edges, 
    with notes on how to do the same for nodes in parentheses.

To add labels on edges (or nodes), we need to know their numbers. We can use the 
`showNodeNumbers = true` parameter for this. (Use `showEdgeNumbers = true` to see node numbers).

```@example getting_started
R"svg"(figname("gettingstarted4.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, showEdgeNumber=true);
R"dev.off()" # hide
nothing # hide
```
![example4](../assets/figures/gettingstarted4.svg)

We will need to define a DataFrame with two collumns of information: the number of the edge (or 
node), and the label that goes on it, like this:

| Number | Label            |
|--------|------------------|
| 1      | "My first edge"  |
| 2      | "My second edge" |

After including the DataFrames package, we can define it as so:
```@repl
using DataFrames
DataFrame(Number=[1, 2], Label=["My first edge", "My second edge"])
```
Using this dataframe as input to the `edgeLabel` parameter puts the text on the correct nodes:
```@example getting_started
R"svg"(figname("edge_labels_example.svg"), width=5, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, edgeLabel=DataFrame(Number=[1, 2], Label=["My first edge", "My second edge"]));
R"dev.off()" # hide
nothing # hide
```
![example5](../assets/figures/edge_labels_example.svg)

## Hybrid edge styles

We can use the `style` parameter to visualize minor hybrid edges as simple lines, unlike the icytree.org style visualization. `style` is by default `:fulltree`, but by switching it 
to `:majortree`, we can draw minor hybir edges as simple lines.

```@example getting_started
R"svg"(figname("style_example.svg"), width=4, height=4) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, style=:majortree);
R"dev.off()" # hide
nothing # hide
```
![example7](../assets/figures/style_example.svg)


## Using edge lengths

We can use the `useEdgeLength=true` parameter to draw a plot that uses the network's edge lengths to determine the lengths of the lines. For this, we'll use a network that can be found [here](https://github.com/nkarimi/Adansonia_HybSeq/blob/master/trait-evolution/BestH1_372g_calibrated.tre).

```@example getting_started
R"svg"(figname("edge_len_example.svg"), width=6, height=6) # hide
net = readTopology("((Smi165:1.6261423761885154,Pcr070:1.6261423761885154):0.0345640579033647,(((#H18:0.23725347915651637::0.12440503333556951,(Adi001:0.2497235156848997,(Adi003:0.22937089525518586,Adi002:0.22937089525518586):0.02035262042971384):0.36375851550073335):1.3644833483534315e-9,((((Asu001:0.3743115977039037,(Aga001:0.24655560654629088,Aga002:0.24655560654629088):0.1277559911576128):0.0019169543252129813)#H18:0.10531506732808882::0.8755949666644305,(Aza135:0.35341877971983404,((Aza037:0.2766129548639923,(Ama018:0.2198029990484657,Ama006:0.2198029990484657):0.056809955815526614):0.030764173677767633,(Ape009:0.2219108012530084,Ape001:0.2219108012530084):0.08546632728875153):0.04604165117807412):0.12812483963737142):0.016691432891123632,(Aru001:0.2294802905252357,Aru127:0.2294802905252357):0.2687547617230934):0.11524698030178727):0.3865179674498836,Age001:1.0):0.6607064340918802);") # hide

plot(net, :R, useEdgeLength=true);
R"dev.off()" # hide
nothing # hide
```
![example8](../assets/figures/edge_len_example.svg)