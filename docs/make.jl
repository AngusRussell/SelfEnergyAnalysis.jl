using SelfEnergyAnalysis
using Documenter

DocMeta.setdocmeta!(SelfEnergyAnalysis, :DocTestSetup, :(using SelfEnergyAnalysis); recursive=true)

makedocs(;
    modules=[SelfEnergyAnalysis],
    authors="Angus Russell",
    repo="https://git.uwaterloo.ca//a5russel/SelfEnergyAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="SelfEnergyAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://a5russel.gitlab.io/SelfEnergyAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
