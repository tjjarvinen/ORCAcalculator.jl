using ORCAcalculator
using Documenter

DocMeta.setdocmeta!(ORCAcalculator, :DocTestSetup, :(using ORCAcalculator); recursive=true)

makedocs(;
    modules=[ORCAcalculator],
    authors="Teemu Järvinen",
    sitename="ORCAcalculator.jl",
    format=Documenter.HTML(;
        canonical="https://tjjarvinen.github.io/ORCAcalculator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tjjarvinen/ORCAcalculator.jl",
    devbranch="main",
)
