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
`showedgenumber = true` option for this. (Use `shownodenumber = true` to see node numbers).

```@example adding_data
R"svg"(figname("adding_data1.svg"), width=6, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]); R"layout"([1 2]); # hide
net = readnewick("(A,((B,#H1),((C)#H1, D)));") # hide
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
| 1      | "edge number 1"  |
| 2      | "edge # 2" |

After including the DataFrames package, we can define it as so:
```@repl adding_data
using DataFrames
DataFrame(number=[1,2], label=["edge number 1","edge # 2"])
```
Using this data frame as input to the `edgelabel` option (`nodelabel` for nodes)
puts the text on the correct edges:
```@example adding_data
R"svg"(figname("edge_labels_example.svg"), width=4, height=3) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readnewick("(A,((B,#H1),(C,(D)#H1)));") # hide
plot(net, edgelabel=DataFrame(number = [1,2],
                              label = ["edge number 1", "edge # 2"]),
     edgelabelcolor="orangered", edgecex=[0.9,1.1], edgelabeladj=[.5,-.3]);
R"dev.off()" # hide
nothing # hide
```
![example2](../assets/figures/edge_labels_example.svg)

## Adding other annotations using R

We can use the return values of [`plot`](@ref) to get information on the coordinates of
different elements of the plot. Using this, we can add any other information we want.

The [`plot`](@ref) function returns the following named tuple:
```
(:xmin, :xmax, :ymin, :ymax,
 :node_x,    :node_y,    :node_y_lo, :node_y_hi,
 :edge_x_lo, :edge_x_hi, :edge_y_lo, :edge_y_hi,
 :node_data, :edge_data)
```
See the documentation for descriptions of these elements: [`plot`](@ref)

## Side clade bars example

Here's example code that adds bars to denote clades in the margin:

```@example adding_data
R"svg"(figname("side_bars.svg"), width=4, height=4) # hide
R"par"(mar=[.1,.1,.1,.1]) # hide
net = readnewick("(((((((t1,t2),t3),t4),t5),(t6,t7)),(t8,t9)),t10);");
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
res[[:xmin,:xmax]]
```

Looking at `xmin` and `xmax` returned by default, we can see that the x
range is about `(0.3, 9)`.
To give us extra space to work with, we can
set `xlim` to `(0.3,10)`, forcing the range to be wider on the right, for annotations.
We can also see that in this case there is some extra white space on the left,
so that we can increase `xmin` a little bit, finally settling on `xlim=(1,10)`.

```julia
plot(net, xlim=(1, 10));
```

Knowing the coordinates, we can now add more information to the plot through
`RCall`. For this, I use the R functions `segments` and `text` to add side bars with
text on them.

```julia
using RCall # add (install) the RCall package prior to 'using' it
R"segments"([9, 9, 9], [0.8, 7.8, 9.8], [9, 9, 9], [7.2, 9.2, 10.2])
R"text"([9.5, 9.5, 9.5], [4, 8.5, 10], ["C", "B", "A"])
```

# Beyond

To go beyond, we can access data on the node & edges to use them as we wish.
We can access the coordinates of points & segments and more data like this:

```@repl adding_data
res[:node_x] # x coordinate. similarly try res[:node_y]
hcat(res[:node_y_lo], res[:node_y_hi])
DataFrame(edge_x_lo=res[:edge_x_lo], edge_x_hi=res[:edge_x_hi],
          edge_y_lo=res[:edge_y_lo], edge_y_hi=res[:edge_y_hi])
res[:node_data]
res[:edge_data]
```
