# notes to maintain documentation

## to make a local version of the website

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or interactively in `docs/`:
```shell
pkg> activate .
pkg> instantiate
pkg> dev ~/.julia/dev/PhyloPlots
pkg> add PhyloNetworks#master
julia> include("make.jl")
```

it will:
- tests the `jldoctest` blocks of examples in the docstrings
- creates or updates a `build/` directory with html files

## how it works: overview

- built with [Documenter](https://juliadocs.github.io/Documenter.jl), see
  its [doc](https://juliadocs.github.io/Documenter.jl/stable/man/syntax/)
  for more details on its md syntax.
- deployed [here](https://juliaphylo.github.io/PhyloPlots.jl/)
  (go to `dev/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

- `.github/workflows/ci.yml` asks to run `docs/make.jl` after a successful test & build.
  it installs `R`, since PhyloPlots depends on RCall which depends on R.
- the julia script `docs/make.jl` has these steps:
  1. install the master version of PhyloNetworks
  1. run `makedocs()` from `Documenter`: make the documentation;
     runs all `jldoctest` blocks in the source files,
     to check that the output in the blocks matches the actual output;
     translate the "Documenter md" documentation files into html files.
  2. run `deploydocs(...)` from Documenter: push files on `gh-pages` branch.

## what to update

- update the Documenter version in `docs/Project.toml`: check to see if a new
  version of Documenter was released. If so, up the dependency, check locally
  that it works, make any updates as needed.
- update Julia version in `.github/workflows/ci.yml`, `docs:` section.

## plots for the documentation

We chose to group all the output plots in the directory `assets/figures`.
Hence, the typical setup in a documentation page containing plots is:

    ```@setup name
    using PhyloNetworks, PhyloPlots, RCall
    mkpath("../assets/figures")
    figname(x) = joinpath("..", "assets", "figures", x)
    ```

The `mkpath` command is there to ensure that the target directory does indeed
exist. In theory, it only needs to be called once (by the first documentation
page being built). However, as this order might be subject to change over time,
it could be safer to include it on every such page.

After trial and error and discussions, we chose to use only the *SVG* format.
This format should ensure that when a plot is drawn again, it is identical,
making it more efficient for git to track changes.

**Warning**: If the same file name is re-used across documentation pages, only the
final version will be committed by git. Make sure to use different file names for
plots that are supposed to look different (across the whole site).
