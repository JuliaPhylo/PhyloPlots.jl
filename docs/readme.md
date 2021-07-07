# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl)
- deployed [here](https://cecileane.github.io/PhyloPlots.jl/)
  (go to `dev/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

## how it works: overview

- `.github/workflows/ci.yml` asks to run `docs/make.jl` after a successful test & build.
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
- update Julia version in `.github/workflows/ci.yml`, `matrix.version`
## The "Documenter md" format

### Note on the format

The documentation pages are all written in this format. It is a regular md, but
with extra blocks of codes (as `@example` and `@setup`) that contain Julia
commands. These lines will be executed during the `makedoc()` process. See the
`Documenter` [doc](https://juliadocs.github.io/Documenter.jl/stable/man/syntax/)
for more details on the syntax. For instance, @example blocks with the same "name"
are run in the same session. Otherwise, an @example blocks with no name
is run in its own anonymous Modules.

### Setting up the plot environment

Some of these blocs may contain plots, which are going to be drawn during the
process, requiring the use of `PhyloPlots` along with `RCall`. Hence,
before the doc is built, `.github/workflows/ci.yml` installs `R` on the server,
sets up the julia environment with dependencies like `PhyloPlots` before
starting the build in itself.

### Directory of the plots

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

### Format of the plots

After trial and error and discussions, we chose to use only the *SVG* format
for the plot. This format should ensure that when a plot is drawn again,
identical, in a new build, then Git will recognize that it has not change, and
hence not save a new version of it. This should ensure that the repository does
not grow unreasonably at each build of the doc, i.e. after each push to
master. The typical commands to save and display a plot should hence be:

    ```@example name
    R"svg"(figname("my_useful_name.svg"), width=4, height=3); # hide
    plot(net, :R);
    R"dev.off()" # hide
    nothing # hide
    ```
    ![my_useful_name](../assets/figures/my_useful_name.svg)

**Warning**: this is not like an interactive session. If the same file name
is re-used by some other documentation page for some other plot, only the
final version of the plot will be committed by git, with possible unintended
consequences. Make sure to use different file names for plots that are supposed
to look different (across the whole site).

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

## references

big difference:

    [blabla](@ref)
    [`blabla`](@ref)

The first version will look for a *section* header "blabla", to link to that section.
The secon version will look for a *function* named "blabla",
to link to the documentation for that function.
