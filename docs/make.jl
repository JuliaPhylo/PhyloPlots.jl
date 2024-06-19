using Documenter

using Pkg
Pkg.add(PackageSpec(name="PhyloNetworks", rev="master"))

using PhyloNetworks
using PhyloPlots
DocMeta.setdocmeta!(PhyloPlots, :DocTestSetup, :(using PhyloPlots); recursive=true)

makedocs(
    sitename = "PhyloPlots.jl",
    authors = "Cécile Ané and Guilhem Ané",
    modules = [PhyloPlots], # to list plot() methods from PhyloPlots only, not from Gadfly etc.
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    pages = [
        "home" => "index.md",
        "manual" => [
            "installation" => "man/installation.md",
            "getting started" => "man/getting_started.md",
            "untangling edges" => "man/untangling_edges.md",
            "better edges" => "man/better_edges.md",
            "adding data" => "man/adding_data.md",
        ],
        "library" => [
            "public" => "lib/public.md",
            "internals" => "lib/internals.md",
        ]
    ],
)

deploydocs(
    repo = "github.com/JuliaPhylo/PhyloPlots.jl.git",
    push_preview = true,
)
