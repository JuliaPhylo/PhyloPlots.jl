using Documenter, DocumenterMarkdown

using Pkg
Pkg.add(PackageSpec(name="PhyloNetworks", rev="master"))

using PhyloPlots
DocMeta.setdocmeta!(PhyloPlots, :DocTestSetup, :(using PhyloPlots); recursive=true)

makedocs(
    sitename = "PhyloPlots.jl",
    authors = "Cécile Ané",
    modules = [PhyloPlots], # to list plot() methods from PhyloPlots only, not from Gadfly etc.
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    pages = [
        "home" => "index.md",
        "library" => [
            "public" => "lib/public.md",
            "internals" => "lib/internals.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/cecileane/PhyloPlots.jl.git",
    push_preview = true,
)
