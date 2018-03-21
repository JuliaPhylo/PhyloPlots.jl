# notes to maintain documentation

- built with [Documenter](https://juliadocs.github.io/Documenter.jl)
- deployed [here](http://cecile.github.io/PhyloPlots.jl/)
  (go to `latest/` or `stable/`)
  using github and files committed to the `gh-pages` branch.

## how it works: overview

- `.travis.yml` asks to run `docs/make.jl` after a successful test & build.
- the julia script `docs/make.jl` has 2 steps:
  1. run `makedocs()` from `Documenter`: make the documentation
  2. run `deploydocs(...)` also from Documenter. This step calls `mkdocs`,
     which turns the markdown files in `docs/.../*.md` into html files.

## what to update

- update Julia version in `docs/make.jl`

## to make a local version of the website

```shell
cd ~/.julia/v0.6/PhyloPlots/docs
julia --color=yes make.jl
```

first line: adapt to where the package lives  
second line:
- tests the `jldoctest` blocks of examples in the docstrings
- creates or updates a `build/` directory with markdown files.
- does *not* convert the markdown files into html files.

To do this html conversion, use MkDocs directly,
and the mkdocs-material package (for the "material" theme).
First check/install MkDocs:

```shell
pip install --upgrade pip
pip install --upgrade mkdocs
pip install --upgrade mkdocs-material
pip install --upgrade python-markdown-math
```
and check the installed versions:
```shell
python --version
mkdocs --version
pip show mkdocs-material
pip show Pygments
pip show pymdown-extensions
pip show python-markdown-math
```

then use mkdocs to build the site.
this step creates a `site/` directory with html files.
they can be viewed at http://127.0.0.1:8000 (follow instructions)

```shell
mkdocs build
mkdocs serve
```
