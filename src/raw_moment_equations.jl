"""
    generate_raw_moment_eqs(rn::Union{ReactionSystem, ReactionSystemMod}, m_order::Int;
                            combinatoric_ratelaw=true, smap=speciesmap(rn))

Given a [`ReactionSystem`](https://catalyst.sciml.ai/stable/api/catalyst_api/#ModelingToolkit.ReactionSystem)
or [`ReactionSystemMod`](@ref), return the [`RawMomentEquations`](@ref) of the system generated up to `m_order`.

Notes:
- The expansion order ``q``, denoted by `q_order` throughout the docs, is automatically
  determined from the given polynomial form of the propensity functions, see the
  [tutorial](@ref main_tutorial) and the [theory section](@ref raw_moment_eqs) for
  more details on how `q_order` is obtained.
- `combinatoric_ratelaw=true` uses binomials in calculating the propensity functions
  of a `ReactionSystem`, see the notes for [`ModelingToolkit.jumpratelaw`]
  (https://mtk.sciml.ai/stable/systems/ReactionSystem/#ModelingToolkit.jumpratelaw).
  *Note* that this field is irrelevant using `ReactionSystemMod` as then the
  propensities are defined directly by the user.
- `smap` sets the variable ordering in the moment equations (which index corresponds to which species
  in the reaction network). By default, this is consistent with the internal system ordering
  accessible with [`speciesmap`](@ref).
"""
function generate_raw_moment_eqs(rn::Union{ReactionSystem,ReactionSystemMod}, m_order::Int;
                                 combinatoric_ratelaw=true, smap=speciesmap(rn))
    generate_raw_moment_eqs(Moment{numspecies(rn)}, rn, m_order; 
                            combinatoric_ratelaw=combinatoric_ratelaw, smap=smap)
end

function generate_raw_moment_eqs(::Type{Moment{N}},
                                 rn::Union{ReactionSystem,ReactionSystemMod}, m_order::Int;
                                 combinatoric_ratelaw=true, smap=speciesmap(rn)) where {N}
    S = get_S_mat(rn; smap)
    a = propensities(rn; combinatoric_ratelaw)

    term_factors, term_powers, poly_order = polynomial_propensities(a, rn; smap)

    q_order = poly_order + m_order - 1

    # iterator over all moments from lowest to highest moment order
    iter_all = construct_iter_all(Moment{N}, q_order)
    # iterator over the first order moments
    iter_1 = get_iter_1(iter_all)
    # iterator over raw moments up to order m
    iter_m = get_iter_m(iter_all, m_order)

    μ = define_μ(Moment{N}, q_order, iter_all)

    iv = get_iv(rn)
    D = Differential(iv)
    eqs = Equation[]

    for i in Iterators.flatten((iter_1, iter_m))
        push!(eqs, D(μ[i]) ~ get_dμ(i, iter_all, S, μ, term_factors, term_powers))
    end

    vars = extract_variables(eqs, Moment{N}, q_order, iter_all)
    odes = ODESystem(eqs, iv, vars, get_ps(rn); 
                     name=Symbol(nameof(rn),"_raw_moment_eqs_",m_order))

    RawMomentEquations{N}(
        odes,
        μ,
        m_order,
        q_order,
        iter_all,
        iter_m
    )
end

function get_dμ(i::Moment{N}, iter_all, S, μ, term_factors, term_powers) where {N}
    ret = 0
    for r = 1:size(S, 2)
        iter_j = filter(x -> all(x .<= i) && sum(x) <= sum(i) - 1, iter_all)
        for j in iter_j
            factor_j = 1
            for k = 1:N
                factor_j *= expected_coeff(S[k, r], i[k] - j[k]) * binomial(i[k], j[k])
            end
            suma = 0
            for k = 1:length(term_factors[r])
                suma += term_factors[r][k] * μ[j.+Tuple(term_powers[r][k])]
            end
            ret += factor_j * suma
        end
    end

    mc_expand(ret)
end
