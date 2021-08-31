function zero_closure(sys::CentralMomentEquations{N}, binary_vars::AbstractVector{Int}=Int[]) where {N}
    closure = OrderedDict()
    closure_exp = OrderedDict()

    if !isempty(binary_vars)
        sys = bernoulli_moment_eqs(sys, binary_vars)
    end

    for i in get_iter_q(sys)
        closure[sys.M[i]] = 0
        closure_exp[sys.M[i]] = 0
    end

    close_eqs(sys, closure_exp, closure, true)
end

function zero_closure(sys::RawMomentEquations{N}, binary_vars::AbstractVector{Int}=Int[]) where {N}
    closure = OrderedDict()
    closure_exp = OrderedDict()

    if !isempty(binary_vars)
        sys = bernoulli_moment_eqs(sys, binary_vars)
    end

    raw_to_central = raw_to_central_moments(Moment{N}, sys.q_order, get_iter_all(sys))

    unique_iter_q = unique(sort(i) for i in get_iter_q(sys))
    sub = Dict()

    for i in unique_iter_q
        μi = sys.μ[i] - raw_to_central[i]
        closure[sys.μ[i]] = μi
        μi = mc_simplify(substitute(μi, closure_exp), expand=true)
        closure_exp[sys.μ[i]] = μi

        perms = multiset_permutations(i, length(i))

        for iter_perm in Iterators.drop(perms, 1)

            iter_perm_ind = sortperm(sortperm(iter_perm))
            for r in get_iter_all(sys)
                sub[sys.μ[r]] = sys.μ[r[iter_perm_ind]]
            end

            iter_perm = Tuple(iter_perm)
            closure[sys.μ[iter_perm]] = substitute(closure[sys.μ[i]], sub)
            closure_exp[sys.μ[iter_perm]] = substitute(closure_exp[sys.μ[i]], sub)
        end

    end

    close_eqs(sys, closure_exp, closure, true)
end
