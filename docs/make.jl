using MDBM
using Documenter

makedocs(;
    modules=[MDBM],
    authors="Daniel Bachrathy",
    repo="https://github.com/bachrathyd/MDBM.jl/blob/{commit}{path}#L{line}",
    sitename="MDBM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
