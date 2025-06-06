"""
    SortedCache{K,V}

A small wrapper that keeps:
  • `data::Dict{K,V}` for O(1) lookup/insert  
  • `keys::Vector{K}` (sorted) for position‐indexing/searchsortedfirst  

Insertion only does a binary‐search+`insert!` on the *key* array (cheap `K` copies,
not shifting huge containers), and storing the heavy `V` in a hash‐table.
"""
struct SortedCache{K,V}
    data::Dict{K,V}
    keys::Vector{K}
end

SortedCache{K,V}() where {K,V} = SortedCache(Dict{K,V}(), Vector{K}())

# insert or update
function insert!(sc::SortedCache{K,V}, key::K, val::V) where {K,V}
    if !haskey(sc.data, key)
     i = searchsortedfirst(sc.keys, key)
        Base.insert!(sc.keys, i, key)      # shifts only K’s, not whole containers
    end
    sc.data[key] = val
    return nothing
end

# lookup by key
get(sc::SortedCache{K,V}, key::K, default) where {K,V} = get(sc.data, key, default)
Base.getindex(sc::SortedCache{K,V}, key::K)  where {K,V}   = sc.data[key]

# “vector” interface: idx→ (funval,cval,callargs)
Base.length(sc::SortedCache) = length(sc.keys)
function Base.getindex(sc::SortedCache{K,V} , i::Integer) where {K,V}
    k = sc.keys[i]
    return sc.data[k]
end

# allow `searchsortedfirst(sc, key)` and iteration
Base.iterate(sc::SortedCache) = iterate(sc.keys)
Base.iterate(sc::SortedCache, st) = iterate(sc.keys, st)
Base.searchsortedfirst(sc::SortedCache{K,V}, key::K) where {K,V} =
    searchsortedfirst(sc.keys, key)


a=(1.2,4.4)
foo(aa)=sum(aa.^2)
foo(a)


SC=SortedCache{Tuple{Float64,Float64},Float64}()
SC
for c in 0:0.1:100
    for d in 0:0.1:100
    insert!(SC, (d,c), foo((d,c)))
    end
end
insert!(SC, (1.2,1.4), foo((1.2,100.4)))
@btime insert!(SC, (1.2,rand()), foo((1.2,rand())))

length(SC)
 using BenchmarkTools
 Base.summarysize(SC)/1024^2
 @btime SC[(1.2, 1.4)]

 A=zeros(Float64,1000,1000);
 Base.summarysize(A)*(2+1+2)/1024^2
@btime A[1,1]

SC[1]
SC.keys[[1,2]]


@time SC[(1.2, 1.41235)]