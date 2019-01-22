include("MDBM__types.jl")


using LinearAlgebra

function b1(x,y)
    res=similar(x)
    res=(maximum(abs.(eigen(rand(5,5)).values)))
end
function b2(x,y)
    maximum(abs.(eigen(rand(5,5)).values)),x
end


function b3(x,y)
    [maximum(abs.(eigen(rand(5,5)).values)),x]
end

fa1=MemF(b1,Array{MDBMcontainer{Float64,Tuple{Float64, Float64}}}(undef, 0))
fa2=MemF(b2,Array{MDBMcontainer{Tuple{Float64, Float64},Tuple{Float64, Float64}}}(undef, 0))
fa3=MemF(b3,Array{MDBMcontainer{Vector{Float64},Tuple{Float64, Float64}}}(undef, 0))


# fa3=MemF((x,y)->(x^2.0)/y,Array{MDBMcontainer{Float64,Tuple{Float64, Float64}}}(undef, 0))
fa3(-1.40,-2.7584120)
fa3(-1.4,1.7)
# fa3([(-1.4,-2.7),(1.1,2.2)])
fa2(-1.40,-2.7584120)
fa2(-1.4,1.7)

# sok-sok futtatás az idő teszteléséhez
NN=2_000_000
N=150
@time for k in 1:NN
    fa3(ceil(rand(Float64)*N),ceil(rand(Float64)*N))
end
println(length(fa3.fvalarg))


@time for k in 1:NN
    fa3.f(ceil(rand(Float16)*N),ceil(rand(Float16)*N))
end
println("----------------")
