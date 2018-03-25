# sexp uses helper functions defined in PhyloNetworks

doc"""
    function sexp(net::HybridNetwork)

Export a HybridNework object to the R language
as either `phylo` or `evonet` object (depending on degree of hybridization)
recognized by the R package `ape`.
Used by the `$object` syntax and by `@rput` to use a Julia object in R:
see the examples below. Makes it easy to plot a Julia tree or network
using plotting facilities in R.


code inspired from [Phylo.jl](https://github.com/richardreeve/Phylo.jl/blob/master/src/rcall.jl)

# Examples

```julia-repl
julia> using RCall
julia> using PhyloNetworks
julia> net = readTopology("(((A:.2,(B:.1)#H1:.1::0.9):.1,(C:.11,#H1:.01::0.1):.19):.1,D:.4);");
R> library(ape); # type $ to switch from julia to R
R> $net

Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "A" "B" "C" "D"

Rooted; includes branch lengths.
julia> @rput net # press the delete key to switch from R back to julia
R> net

Evolutionary network with 1 reticulation

--- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "A" "B" "C" "D"

Rooted; includes branch lengths.

R> str(net)
List of 7
 $ edge               : int [1:8, 1:2] 5 5 6 6 7 8 8 9 6 4 ...
 $ reticulation.length: num 0.01
 $ Nnode              : int 5
 $ edge.length        : num [1:8] 0.1 0.4 0.1 0.19 0.11 0.2 0.1 0.1
 $ reticulation       : int [1, 1:2] 7 9
 $ reticulation.gamma : num 0.1
 $ tip.label          : chr [1:4] "A" "B" "C" "D"
 - attr(*, "class")= chr [1:2] "evonet" "phylo"
NULL

R> plot(net)
```
""" #"

function sexp(net::HybridNetwork)
    PhyloNetworks.resetNodeNumbers!(net)
    ntips = length(net.leaf)
    totalnodes = length(net.node)
    Nnode = totalnodes - ntips
    o = sortperm([n.number for n in net.leaf])
    tipLabel = [net.leaf[i].name for i in o]
    edge = PhyloNetworks.generateMajorEdge(net)
    phy = Dict{Symbol, Any}() # dictionary exported with regular sexp at the end
    phy[:Nnode] = Nnode
    phy[Symbol("tip.label")] = tipLabel
    phy[:edge] = edge
    edgeLength = PhyloNetworks.generateMajorLength(net)
    nBL = sum([!isnull(e) for e in edgeLength]) # number of non-missing branch lengths
    if nBL>0
        phy[Symbol("edge.length")] = edgeLength
    end
    if net.numHybrids > 0
        reticulation = PhyloNetworks.generateMinorReticulation(net)
        reticulationGamma = PhyloNetworks.generateMinorReticulationGamma(net)
        reticulationLength = PhyloNetworks.generateMinorReticulationLength(net)
        phy[:reticulation] = reticulation
        if sum([!isnull(e) for e in reticulationGamma]) > 0
            phy[Symbol("reticulation.gamma")] = reticulationGamma
        end
        if sum([!isnull(e) for e in reticulationLength]) > 0
            phy[Symbol("reticulation.length")] = reticulationLength
        end
    end
    sobj = protect(sexp(phy)) # RObject
    if net.numHybrids == 0
        setclass!(sobj, sexp("phylo"))
    else
        setclass!(sobj, sexp(["evonet", "phylo"]))
    end
    unprotect(1)
    return(sobj)
end

doc"""
    rexport(net::HybridNetwork; mainTree=false, useEdgeLength=true)

Create an RObject of class `phylo` (and `evonet` depending on the number
of hybridizations) recognized by the `ape` library in R (S3 object). This
RObject can be evaluated using the tools available in the `ape` library in R.
For example, we can visualize the network using `ape`'s `plot` function.

not exported: [`sexp`](@ref) is the best way to go.

# Arguments

- useEdgeLength: if true, export edge lengths from `net`.
- mainTree: if true, minor hybrid edges are omitted.

# Examples

```julia-repl
julia> net = readTopology("(((A,(B)#H1:::0.9),(C,#H1:::0.1)),D);");
julia> phy = rexport(net)
RCall.RObject{RCall.VecSxp}
$Nnode
[1] 5

$edge
     [,1] [,2]
[1,]    5    6
[2,]    5    4
[3,]    6    8
[4,]    6    7
[5,]    7    3
[6,]    8    1
[7,]    8    9
[8,]    9    2

$tip.label
[1] "A" "B" "C" "D"

$reticulation
     [,1] [,2]
[1,]    7    9

$reticulation.gamma
[1] 0.1

attr(,"class")
[1] "evonet" "phylo"

julia> using RCall

julia> R"library(ape)";

julia> phy
RCall.RObject{RCall.VecSxp}

    Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "A" "B" "C" "D"

Rooted; no branch lengths.

R> phy

Evolutionary network with 1 reticulation

               --- Base tree ---
Phylogenetic tree with 4 tips and 5 internal nodes.

Tip labels:
[1] "A" "B" "C" "D"

Rooted; no branch lengths.

R> str(phy)
List of 5
$ Nnode             : int 5
$ edge              : int [1:8, 1:2] 5 5 6 6 7 8 8 9 6 4 ...
$ tip.label         : chr [1:4] "A" "B" "C" "D"
$ reticulation      : int [1, 1:2] 7 9
$ reticulation.gamma: num 0.1
- attr(*, "class")= chr [1:2] "evonet" "phylo"
```
""" #"
# worry about R object created within the function not accessible from outside:
# can it be garbage collected?
function rexport(net::HybridNetwork; mainTree::Bool=false, useEdgeLength::Bool=true)
    if mainTree == true && net.numHybrids > 0
        net = majorTree(net)
    end
    PhyloNetworks.resetNodeNumbers!(net)
    ntips = length(net.leaf)
    totalnodes = length(net.node)
    Nnode = totalnodes - ntips
    o = sortperm([n.number for n in net.leaf])
    tipLabel = [net.leaf[i].name for i in o]
    edge = PhyloNetworks.generateMajorEdge(net)
    R"""
    phy = list(Nnode = $Nnode, edge = $edge, tip.label = $tipLabel)
    """
    if useEdgeLength == true
        edgeLength = PhyloNetworks.generateMajorLength(net)
        if sum([!isnull(e) for e in edgeLength])>0
            R"""
            phy[['edge.length']] = $edgeLength
            """
        end
    end
    if net.numHybrids > 0
        reticulation = PhyloNetworks.generateMinorReticulation(net)
        reticulationGamma = PhyloNetworks.generateMinorReticulationGamma(net)
        R"""
        phy[['reticulation']] = $reticulation
        class(phy) <- c("evonet", "phylo")
        """
        if sum([!isnull(e) for e in reticulationGamma])>0
            R"phy[['reticulation.gamma']] = $reticulationGamma"
        end
        if useEdgeLength # extract minor edge lengths
            reticulationLength = PhyloNetworks.generateMinorReticulationLength(net)
            if sum([!isnull(e) for e in reticulationLength])>0
                R"""
                phy[['reticulation.length']] = $reticulationLength
                """
            end
        end
    elseif net.numHybrids == 0
        R"""
        class(phy) = "phylo"
        """
    end
    phy = reval("phy")
    return(phy)
end
