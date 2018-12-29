using Documenter, DocumenterMarkdown
using PhyloPlots

makedocs(sitename = "PhyloPlots.jl",
         modules = [PhyloPlots], # to list plot() methods from PhyloPlots only, not from Gadfly etc.
         format = Markdown(),
         Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true") # to make it easier when building locally, if using HTML (not markdown) format
        )

deploydocs(
    repo = "github.com/cecileane/PhyloPlots.jl.git",
    deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    make = () -> run(`mkdocs build`),
    target = "site" # which files get copied to gh-pages
)
