using SelfEnergyAnalysis
using Documenter

DocMeta.setdocmeta!(SelfEnergyAnalysis, :DocTestSetup, :(using SelfEnergyAnalysis); recursive=true)

makedocs(;
    modules=[SelfEnergyAnalysis],
    authors="Angus Russell",
    repo="https://git.uwaterloo.ca/QuINLab/Projects/ExcitonPolariton/exciton-polariton-self-energy-analysis/SelfEnergyAnalysis.jl",
    sitename="SelfEnergyAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://git.uwaterloo.ca/QuINLab/Projects/ExcitonPolariton/exciton-polariton-self-energy-analysis/SelfEnergyAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Guide" => Any[
        "Processing data " => "tut/process.md",
        "tut/mdc.md",
        "tut/edc.md",
        "tut/self.md",
        "tut/spec.md",
        ],
    "Library" => Any[
        "Public" => "Lib/public.md",
        "Internals" => "Lib/internals.md"
          
    ],
    "Further Reading" => "Lib/reading.md",
    ],
)
