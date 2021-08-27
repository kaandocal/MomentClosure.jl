function log_normal_closure(sys::MomentEquations{N}, binary_vars::AbstractVector{Int}=Int[]) where {N}
    closure = OrderedDict()
    closure_exp = OrderedDict() # here it does not play a role

    if isempty(binary_vars)
        isbernoulli = false
    else
        sys = bernoulli_moment_eqs(sys, binary_vars)
        isbernoulli = true
    end

    if sys isa CentralMomentEquations
        M = copy(sys.M)
        μ = central_to_raw_moments(N, sys.m_order, get_iter_all(sys), sys.μ, sys.M)
        μ_symbolic = sys.μ
    else
        μ = copy(sys.μ)
        μ_symbolic = copy(μ)
    end

    Σ = Dict()
    iter_1 = get_iter_1(sys)
    for j in 1:N
        for k in j:N
            eⱼ = iter_1[j]
            eₖ = iter_1[k]
            if sys isa CentralMomentEquations
                Σ[(j,k)] = 1. + M[eⱼ .+ eₖ] / μ[eⱼ] / μ[eₖ]
            else
                Σ[(j,k)] = μ[eⱼ .+ eₖ] / μ[eⱼ] / μ[eₖ]
            end
        end
    end

    unique_iter_q = unique(sort(i) for i in get_iter_q(sys))
    sub = Dict()

    for i in unique_iter_q
        term = prod(μ[iter_1[j]]^i[j] for j in 1:N)
        for j in 1:N
            for k in j+1:N
                term *= Σ[j,k]^(i[j]*i[k])
            end
            term *= Σ[j,j]^(i[j]*(i[j]-1)/2)
        end
        μ[i] = mc_simplify(term)
        closure[μ_symbolic[i]] = μ[i]

        perms = collect(multiset_permutations(i, length(i)))[2:end]

        for iter_perm in perms

            iter_perm_ind = sortperm(sortperm(iter_perm))
            for r in get_iter_all(sys)
                sub[μ_symbolic[r]] = μ_symbolic[r[iter_perm_ind]]
            end

            iter_perm = Tuple(iter_perm)
            μ[iter_perm] = substitute(μ[i], sub)
            closure[μ_symbolic[iter_perm]] = μ[iter_perm]
        end

    end

    if typeof(sys) == CentralMomentEquations
        # construct the corresponding truncated expressions of higher order
        # central moments from the obtained log-normal raw moment expressions
        raw_to_central = raw_to_central_moments(N, sys.q_order, 
                                                get_iter_all(sys), μ, sys.M,
                                                bernoulli=isbernoulli)
        central_to_raw = central_to_raw_moments(N, sys.q_order, get_iter_all(sys),
                                                sys.μ, sys.M)
        closure_M = OrderedDict()
        for i in get_iter_q(sys)
            closure_exp[M[i]] = raw_to_central[i]
            closure_M[M[i]] = closure[μ_symbolic[i]]-(central_to_raw[i]-M[i])
            closure_M[M[i]] = mc_simplify(closure_M[M[i]])
        end
        closure = closure_M
    else
        closure_exp = closure
    end

    close_eqs(sys, closure_exp, closure, false)

end
