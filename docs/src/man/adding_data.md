```@setup adding_data
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Adding data

In this section, we will look over ways of adding extra information or data to a plot.

## Adding labels

!!! note
    For demonstration purposes, I will walk through the process of adding labels to edges, 
    with notes on how to do the same for nodes in parentheses.

To add labels on edges (or nodes), we need to know their numbers. We can use the 
`showNodeNumbers = true` option for this. (Use `showEdgeNumbers = true` to see node numbers).

```@example adding_data
R"svg"(figname("adding_data1.svg"), width=3, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),((C)#H1, D)));") # hide
plot(net, :R, showEdgeNumber=true);
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/adding_data1.svg)

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
Using this dataframe as input to the `edgeLabel` option puts the text on the correct nodes:
```@example adding_data
R"svg"(figname("edge_labels_example.svg"), width=4, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, :R, edgeLabel=DataFrame(Number=[1, 2], Label=["My first edge", "My second edge"]));
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_labels_example.svg)

## Adding other data using R

We can use the return values of [`plot`](@ref) to get some information on the coordinates of different elemnts of the plot. Using this, we can add any other information we want.

The [`plot`](@ref) function returns the following tuple:
```
(xmin, xmax, ymin, ymax, node_x, node_y, node_yB, node_yE, 
edge_xB, edge_xE, edge_yB, edge_yE, ndf, edf)
```
See the documentation for descriptions on every element: [`plot`](@ref)

## Side clade bars example

Here's some example code that adds bars to denote clades in the margin:

```@example adding_data
R"svg"(figname("side_bars.svg"), width=4, height=4) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
plot(net, :R, xlim=(1, 10))
R"segments(c(9, 9, 9), c(0.8, 7.8, 9.8), c(9, 9, 9), c(7.2, 9.2, 10.2))"
R"text(c(9.5, 9.5, 9.5), c(4, 8.5, 10), c(\"C\", \"B\", \"A\"))"
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/side_bars.svg)

First, we read the topology, and plot the graph normally. Below is the 
return value, which will be usefull.
```julia
plot(net, :R)
```

Looking at `xmin` and `xmax` in the return value, we can see that the x 
dimentions are about `(0.3, 9)`. To give us extra space to work with, we can 
set `xlim` to `(1,10)`, forcing the dimentions to be wider.

```julia
plot(net, :R, xlim=(1, 10));
```

Knowing the coordinates, we can now add more information to the plot through 
`RCall`. Here, I use the `segments` and `text` functions to add side bars with 
text on them. Don't forget to use backslashes to excape quotes in the RCall 
text.

```
R"segments(c(9, 9, 9), c(0.8, 7.8, 9.8), c(9, 9, 9), c(7.2, 9.2, 10.2))"
R"text(c(9.5, 9.5, 9.5), c(4, 8.5, 10), c(\"C\", \"B\", \"A\"))"
```