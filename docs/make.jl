using NLS_Fit
using Documenter

DocMeta.setdocmeta!(NLS_Fit, :DocTestSetup, :(using NLS_Fit); recursive=true)

makedocs(;
    modules=[NLS_Fit],
    authors="Vincent Picaud",
    repo="https://github.com/vincent-picaud/NLS_Fit.jl/blob/{commit}{path}#{line}",
    sitename="NLS_Fit.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vincent-picaud.github.io/NLS_Fit.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting started" => [
            "simple_fit.md",
            "recalibration.md",
            "mapped_parameters.md",
        ],
        "model2fit_bestiary.md",
        #            "getting_started.md",
        "API.md"
    ],
)

deploydocs(;
    repo="github.com/vincent-picaud/NLS_Fit.jl.git",
    devbranch="main",
)
