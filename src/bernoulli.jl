function bernoulli_truncate(powers::Moment{N}, binary_vars)::Moment{N} where {N}
    iter_temp = collect(powers) # need to transform into array as Tuple is immutable
    for ind in binary_vars
        if powers[ind] > 1
            iter_temp[ind] = 1
        end
    end

    Tuple(iter_temp)
end

function bernoulli_iter_redundant(sys, binary_vars::AbstractVector{Int})
    redundant_iter_sub = Dict()

    # 0th and 1st order terms stay the same
    for iter in get_iter_M(sys)
        iter_temp = bernoulli_truncate(iter, binary_vars)

        if iter != iter_temp
            redundant_iter_sub[iter] = iter_temp
        end

    end

    redundant_iter_sub
end

function get_redundant_eqs(sys::Union{RawMomentEquations{N}, CentralMomentEquations{N}}, 
                           binary_vars) where {N}
    ret = Int[]

    for ind in binary_vars
        for (i, iter) in enumerate(get_iter_m(sys))
            if iter[ind] > 1
                push!(ret, i+N)
            end
        end
    end

    ret
end

function bernoulli_reduce(sys::RawMomentEquations{N}, binary_vars::AbstractVector{Int}) where {N}
    redundant_iter_sub = bernoulli_iter_redundant(sys, binary_vars)
    redundant_iter = keys(redundant_iter_sub)

    redundant_eqs = get_redundant_eqs(sys, binary_vars)

    iter_sub = [sys.μ[key] => sys.μ[val] for (key, val) in redundant_iter_sub]

    redundant_iter, redundant_eqs, iter_sub
end


function bernoulli_reduce(sys::CentralMomentEquations{N}, binary_vars::AbstractVector{Int}) where {N}
    redundant_iter_sub = bernoulli_iter_redundant(sys, binary_vars)
    redundant_iter = keys(redundant_iter_sub)

    redundant_eqs = get_redundant_eqs(sys, binary_vars)

    μ = define_μ(Moment{N}, sys.q_order, get_iter_all(sys))
    μ_redundant_sub = [μ[key] => μ[val] for (key, val) in redundant_iter_sub]

    clean_iter = setdiff(get_iter_M(sys), redundant_iter)
    central_to_raw = central_to_raw_moments(Moment{N}, sys.q_order, get_iter_all(sys))
    μ_clean_sub = Dict(μ[iter] => central_to_raw[iter] for iter in clean_iter)

    raw_to_central = raw_to_central_moments(Moment{N}, sys.q_order, get_iter_all(sys))
    iter_sub = Dict()
    for iter in redundant_iter
        M_temp = raw_to_central[iter]
        M_temp = mc_simplify(substitute(M_temp, μ_redundant_sub))
        M_temp = substitute(M_temp, μ_clean_sub)
        M_temp = mc_simplify(M_temp)
        iter_sub[sys.M[iter]] = M_temp
    end

    redundant_iter, redundant_eqs, iter_sub
end

"""
    bernoulli_moment_eqs(sys::MomentEquations, binary_vars::Array{Int,1})

Given `MomentEquations` and an array of indices specifying the species which molecule numbers are
binary variables (either 0 or 1), apply identities of Bernoulli variables to remove the redundant
ODEs and return the *cleaned up* `MomentEquations`. See [here](@ref geometric-and-conditional) for
example usage.
"""
function bernoulli_moment_eqs(sys::Union{RawMomentEquations{N}, CentralMomentEquations{N}},
                              binary_vars::AbstractVector{Int}) where {N}
    isempty(binary_vars) && return sys

    # Construct moment equations removing the redundant ones
    # noting the properties of the Bernoulli variables in the system

    redundant_iter, redundant_eqs, iter_sub = bernoulli_reduce(sys, binary_vars)

    # construct the cleaned moment equations
    clean_eqs = Equation[]
    for (i, eq) in enumerate(get_eqs(sys.odes))
        if !(i in redundant_eqs)
            clean_rhs = substitute(eq.rhs, iter_sub)
            clean_rhs = mc_expand(clean_rhs)
            push!(clean_eqs, Equation(eq.lhs, clean_rhs))
        end
    end

    # construct a new Raw/Central/MomentEquations system saving only the clean iterators
    iter_all = setdiff(get_iter_all(sys), redundant_iter)
    iter_m = setdiff(get_iter_m(sys), redundant_iter)
    iter_q = setdiff(get_iter_q(sys), redundant_iter)

    ## fixing ODE system to preserve consistent ordering of parameters
    iv = get_iv(sys.odes)
    ps = get_ps(sys.odes)

    vars = extract_variables(clean_eqs, Moment{N}, sys.q_order, iter_all)
    odes = ODESystem(clean_eqs, iv, vars, ps; name=Symbol(nameof(sys.odes),"_bernoulli"))

    if sys isa RawMomentEquations
        RawMomentEquations{N}(odes, sys.μ, sys.m_order, sys.q_order, iter_all, iter_m, iter_q)
    else
        CentralMomentEquations{N}(odes, sys.μ, sys.M, sys.m_order, sys.q_order, iter_all, iter_m, iter_q)
    end
end
