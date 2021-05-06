using Photodynamics
using Documenter

DocMeta.setdocmeta!(Photodynamics, :DocTestSetup, :(using Photodynamics); recursive=true)

makedocs(;
    modules=[Photodynamics],
    authors="Zach Langford <langfzac@uw.edu> and contributors",
    repo="https://github.com/langfzac/Photodynamics.jl/blob/{commit}{path}#{line}",
    sitename="Photodynamics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://langfzac.github.io/Photodynamics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/langfzac/Photodynamics.jl",
)
