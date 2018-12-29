# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl)
- deployed [here](https://cecileane.github.io/PhyloPlots.jl/)
  (go to `dev/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

## how it works: overview

- `.travis.yml` asks to run `docs/make.jl` after a successful test & build.
- the julia script `docs/make.jl` has 2 steps:
  1. run `makedocs()` from `Documenter`: make the documentation
  2. run `deploydocs(...)` also from Documenter. This step calls `mkdocs`,
     which turns the markdown files in `docs/.../*.md` into html files.

## what to update

- update Julia version in `.travis.yml`, Documentation section

## to make a local version of the website

```shell
julia --project=docs/ -e 'using Pkg; Pkg.instantiate(); Pkg.develop(PackageSpec(path=pwd()))'
julia --project=docs/ --color=yes docs/make.jl
```

or interactively in `docs/`:
```shell
pkg> activate .
pkg> instantiate
pkg> # dev PhyloNetworks # to get the j07 branch: master still at julia v0.6
pkg> dev ~/.julia/dev/PhyloPlots
julia> include("make.jl")
```

it will:
- tests the `jldoctest` blocks of examples in the docstrings
- creates or updates a `build/` directory with markdown files.
- does *not* convert the markdown files into html files.

To do this html conversion, use [MkDocs](http://www.mkdocs.org) directly,
and the `mkdocs-material` package (for the "material" theme).
First check/install `MkDocs`:

```shell
pip install --upgrade pip
pip install --upgrade mkdocs
pip install --upgrade mkdocs-material
pip install --upgrade python-markdown-math
pip install --upgrade Pygments
```
and check the installed versions
(in comments are versions that work okay together):
```shell
python --version # Python 3.5.5 :: Anaconda, Inc.
mkdocs --version              # v0.17.4  v1.0.4
pip show mkdocs-material      # v2.9.2   v3.2.0
pip show Pygments             # v2.2.0   v2.3.1
pip show pymdown-extensions   # v4.11    v4.11
pip show python-markdown-math # v0.6     v0.6
```

then use mkdocs to build the site.
this step creates a `site/` directory with html files.
they can be viewed at http://127.0.0.1:8000 (follow instructions)

```shell
mkdocs build
mkdocs serve
```
