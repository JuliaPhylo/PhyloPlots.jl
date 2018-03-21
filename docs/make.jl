using Documenter, PhyloPlots

makedocs()

deploydocs(
    deps   = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/cecileane/PhyloPlots.jl.git",
    julia  = "0.6",
    osname = "linux"
)
