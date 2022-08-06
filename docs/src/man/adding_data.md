```@setup adding_data
using PhyloNetworks, PhyloPlots, RCall, DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Adding data

In this section, we look over ways of adding extra information or data to a plot.

## Adding labels

!!! note
    For demonstration purposes, we walk through the process of adding labels to edges,
    with notes on how to do the same for nodes in parentheses.

To add labels on edges (or nodes), we need to know their numbers. We can use the
`showedgenumbers = true` option for this. (Use `shownodenumbers = true` to see node numbers).

```@example adding_data
R"svg"(figname("adding_data1.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]); # hide
net = readTopology("(A,((B,#H1),((C)#H1, D)));") # hide
plot(net, showedgenumber=true);
plot(net, showedgenumber=true, edgenumbercolor="red4");
R"dev.off()" # hide
nothing # hide
```
![example1](../assets/figures/adding_data1.svg)

Edge numbers are shown in grey by default (to avoid mistaking them
for edge lengths), but their color can be adjusted as shown above.

We then need to define a DataFrame with two columns of information: the number of the edge (or
node), and the label that goes on it, like this:

| number | label            |
|--------|------------------|
| 1      | "my first edge"  |
| 2      | "my second edge" |

After including the DataFrames package, we can define it as so:
```@repl
using DataFrames
DataFrame(number=[1,2], label=["my first edge","my second edge"])
```
Using this data frame as input to the `edgelabel` option (`nodelabel` for nodes)
puts the text on the correct edges:
```@example adding_data
R"svg"(figname("edge_labels_example.svg"), width=4, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, edgelabel=DataFrame(number = [1,2],
                              label = ["my first edge", "my second edge"]),
          edgelabelcolor = "orangered");
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_labels_example.svg)

## Adding other annotations using R

We can use the return values of [`plot`](@ref) to get information on the coordinates of
different elements of the plot. Using this, we can add any other information we want.

The [`plot`](@ref) function returns the following tuple:
```
(xmin, xmax, ymin, ymax,
 node_x, node_y, node_yB, node_yE,
 edge_xB, edge_xE, edge_yB, edge_yE,
 nodedataframe, edgedataframe)
```
See the documentation for descriptions of these elements: [`plot`](@ref)

## Side clade bars example

Here's example code that adds bars to denote clades in the margin:

```@example adding_data
R"svg"(figname("side_bars.svg"), width=4, height=4) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readTopology("(((((((1,2),3),4),5),(6,7)),(8,9)),10);");
plot(net, xlim=(1,10))
using RCall # to send any R command, to make further plot modifications
R"segments"([9, 9, 9], [0.8, 7.8, 9.8], [9, 9, 9], [7.2, 9.2, 10.2])
R"text"([9.5, 9.5, 9.5], [4, 8.5, 10], ["C", "B", "A"])
R"dev.off()" # hide
nothing # hide
```
![example3](../assets/figures/side_bars.svg)

Let's break this down step by step.
First, we read the topology, and plot the graph normally. `plot` actually returns
a value, from which we can get useful information.
Below, we store the plot output in `res`, then check its first two values
because they contain the default range of the x axis; `xmin` and `xmax`.

```@example adding_data
res = plot(net);
res[1:2]
```

Looking at `xmin` and `xmax` returned by default, we can see that the x
range is about `(0.3, 9)`. To give us extra space to work with, we can
set `xlim` to `(1,10)`, forcing the range to be wider.

```julia
plot(net, xlim=(1, 10));
```

Knowing the coordinates, we can now add more information to the plot through
`RCall`. For this, I use the R functions `segments` and `text` to add side bars with
text on them.

```
using RCall # add (install) the RCall package prior to 'using' it
R"segments"([9, 9, 9], [0.8, 7.8, 9.8], [9, 9, 9], [7.2, 9.2, 10.2])
R"text"([9.5, 9.5, 9.5], [4, 8.5, 10], ["C", "B", "A"])
```
