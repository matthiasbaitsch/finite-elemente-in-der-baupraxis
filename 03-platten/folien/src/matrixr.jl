# Matrix with redundant and negative redundant entries

struct MatrixR
    name::String
    indices::Vector{Set{Tuple{Int,Int}}}
    values::Vector{Any}
    MatrixR(name::String) = new(name, Vector{Set{Tuple{Int,Int}}}(), Vector{Any}())
end

function MatrixR(name, a::Symmetric)
    n = size(a, 1)
    ar = MatrixR(name)
    for i = 1:n, j = i:n
        ar[i, j] = a[i, j]
    end
    setnegatives!(ar)
    return ar
end

function setnegatives!(a::MatrixR)
    cnt = 0
    for k1 = length(a.indices):-1:2
        v1 = simplify(-a.values[k1])
        for k2 = (k1-1):-1:1
            i1, j1 = first(a.indices[k2])
            v2 = a.values[k2]
            if isequal(v1, v2)
                cnt += 1
                a.values[k1] = "-$(a.name)[$i1, $j1]"
            end
        end
    end
end

countnegatives(a::MatrixR) = count(e -> e isa String && startswith(e, "-"), a.values)

function Base.show(out::IO, a::MatrixR)
    nu = length(a.values) - countnegatives(a)
    println(
        out,
        "        # Matrix $(a.name) has $(nu) unique and negative unique entries")
    for k = eachindex(a.indices)
        print(out, "        ")
        for (i, j) ∈ a.indices[k]
            print(out, "$(a.name)[$i, $j] = ")
        end
        println(out, a.values[k])
    end
end

function Base.setindex!(a::MatrixR, value, i::Int, j::Int)
    idx = indexin(value, a.values)[1]
    if isnothing(idx)
        push!(a.indices, Set([(i, j)]))
        push!(a.values, value)
    else
        push!(a.indices[idx], (i, j))
    end
end

