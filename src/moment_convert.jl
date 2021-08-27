"""
    cumulants_to_raw_moments(N::Int, max_order::Int)

Express all `N`-variate cumulants ``κ`` up to order given by `max_order`
in terms of raw moments ``μ``. Return a Dictionary mapping from
vector ``\\mathbf{i}`` (that indicates the cumulant ``κ_{\\mathbf{i}}``)
to the corresponding raw moment expressions.

# Example

```julia
cumulants_to_raw_moments(2, 2)

Dict{Any,Any} with 5 entries:
  (1, 0) => μ₁₀(t)
  (2, 0) => μ₂₀(t) - ((μ₁₀(t))^2)
  (0, 1) => μ₀₁(t)
  (0, 2) => μ₀₂(t) - ((μ₀₁(t))^2)
  (1, 1) => μ₁₁(t) - ((μ₀₁(t))*(μ₁₀(t)))
```
"""
function cumulants_to_raw_moments(N::Int, max_order::Int, 
                                  iter_all=construct_iter_all(N, max_order),
                                  μ=define_μ(N, max_order, iter_all))
    # following Smith (1995)
    iter_1 = filter(x -> sum(x) == 1, iter_all)

    @assert length(μ) == length(iter_all) "passed arguments are inconsistent (μ vs N & order)"

    K = Dict{NTuple{N,Int},Any}()
    μ_star = Dict{NTuple{N,Int},Any}()

    μ_star[Tuple(zeros(N))] = 1.0
    for i in 1:N
        eᵢ = iter_1[i]
        K[eᵢ] = μ[eᵢ]
        μ_star[eᵢ] = -μ[eᵢ]
    end

    for order in 2:max_order

        iter_order = filter(x -> sum(x) == order, iter_all)
        for r in iter_order

            ind = findlast(!iszero, r)
            r_sub = r .- iter_1[ind]
            iter_i = filter(x -> all(x .<= r_sub), iter_all)

            suma = 0
            for i in iter_i
                factor = prod(binomial.(r_sub, i))
                suma += factor*μ[r.-i]*μ_star[i]
            end
            K[r] = mc_simplify(suma)

            suma = 0
            for i in iter_i
                factor = prod(binomial.(r_sub,i))
                suma += factor*(-K[r.-i])*μ_star[i]
            end
            μ_star[r] = mc_simplify(suma)
        end

    end

    K
end

"""
    cumulants_to_central_moments(N::Int, max_order::Int)

Express all `N`-variate cumulants ``κ`` up to order given by `max_order`
in terms of raw moments ``M``. Return a Dictionary mapping from
vector ``\\mathbf{i}`` (that indicates the cumulant ``κ_{\\mathbf{i}}``)
to the corresponding central moment expressions.
"""
function cumulants_to_central_moments(N::Int, max_order::Int,
                                      iter_all=construct_iter_all(N, max_order),
                                      μ=define_μ(N, 1, iter_all),
                                      M=define_M(N, max_order, iter_all))

    # obtain cumulants up to (m_order)^th order in terms of
    # central moments using formula from Balakrishan et al. (1998)

    K = Dict{NTuple{N,Int},Any}()
    M_star = Dict{NTuple{N,Int},Any}()

    iter_1 = filter(x -> sum(x) == 1, iter_all)

    M_star[Tuple(zeros(N))] = 1.0
    for i in 1:N
        eᵢ = iter_1[i]
        K[eᵢ] = μ[eᵢ]
        M_star[eᵢ] = 0.0
    end

    for order in 2:max_order
        # going up through the orders iteratively to build up the expressions (could do recursively?)

        iter_order = filter(x -> sum(x) == order, iter_all)
        for r in iter_order

            ind = findlast(!iszero, r)
            r_sub = r .- iter_1[ind]
            iter_i = filter(x -> all(x .<= r_sub), iter_all)
            # find the cumulant \kappa_{\bm{r}}}
            suma = 0.0
            for i in iter_i
                factor = prod(binomial.(r_sub, i))
                suma += factor*M[r.-i]*M_star[i]
            end
            K[r] = mc_simplify(suma)

            # Find the central moment \M^*_{\bm{r}}
            suma = 0.0
            for i in iter_i
                factor = prod(binomial.(r_sub, i))
                suma += factor*(-K[r.-i])*M_star[i]
            end
            suma -= -μ[iter_1[ind]]*M_star[r_sub]
            M_star[r] = mc_simplify(suma)

        end

    end

    K
end


function raw_to_central_moments(N::Int, order::Int,
                                iter_all=construct_iter_all(N, order),
                                μ=define_μ(N, order, iter_all),
                                M=define_M(N, order, iter_all); bernoulli=false)

    # Return a dictionary of central moments expressed in terms of raw moments
    # example use:
    # 1 raw_to_central = raw_to_central_moments(2, 3)
    # 2 M₁₂ = raw_to_central[(1,2)] = 2 μ₁₀ μ₀₁² - 2μ₁₁μ₀₁ - μ₀₂μ₁₀ + μ₁₂
    # note that μ is an optional argument which can be used to pass
    # arbitrary values/symbols for each raw moment (different from default μᵢ)

    iter_μ = filter(x -> sum(x) == 1, iter_all)
    if !(μ isa AbstractDict) || length(μ) != length(iter_all)
        if bernoulli
            iter_all = keys(μ)
        else
            error("passed arguments are inconsistent (μ vs N & order)")
        end
    end
    raw_to_central = Dict{NTuple{N,Int},Any}()

    for i in iter_all

        iter_j = Iterators.filter(x -> all(x .<= i), iter_all)
        suma = 0.0
        for j in iter_j
            term = μ[i.-j]
            for (k, e_k) in zip(1:N, iter_μ)
                term *= (-1)^(j[k])*binomial(i[k], j[k])*μ[e_k]^j[k]
            end
            suma += term
        end
        raw_to_central[i] = mc_simplify(suma)
    end

    raw_to_central

end

function central_to_raw_moments(N::Int, order::Int,
                                iter_all=construct_iter_all(N, order),
                                μ=define_μ(N, order, iter_all),
                                M=define_M(N, order, iter_all))

    # Return a dictionary of raw moments expressed in terms of central moments
    # example use:
    # 1 central_to_raw = central_to_raw_moments(2, 3)
    # 2 μ₁₂ = central_to_raw[(1,2)] = 2 M₁₁ μ₀₁ + M₀₂μ₁₀ + M₁₂ + μ₁₀ μ₀₁²

    iter_μ = filter(x -> sum(x) == 1, iter_all)

    central_to_raw = Dict{NTuple{N,Int},Any}()

    for i in iter_all

        iter_j = Iterators.filter(x -> all(x .<= i), iter_all)
        suma = 0.0
        for j in iter_j
            term = M[i.-j]
            for (k, e_k) in zip(1:N, iter_μ)
                term *= binomial(i[k], j[k])*μ[e_k]^j[k]
            end
            suma += term
        end
        central_to_raw[i] = mc_simplify(suma)
    end

    central_to_raw

end
