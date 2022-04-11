function gen_iter(n, d)
    # based on https://twitter.com/evalparse/status/1107964924024635392
    iter = NTuple{n, Int}[]
    for x in partitions(d + n, n)
        x = x .- 1
        if all(x .<= d)
            append!(iter, Tuple.(multiset_permutations(x, n)))
        end
    end
    iter
end

function construct_iter_all(N::Int, order::Int)
    mapreduce(vcat, 0:order) do d
        Base.sort(gen_iter(N, d), rev=true)
    end
end

# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))

# Expand a symbolic expression (no binomial expansion)
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys)) # apply distribution law 
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod))) # distributes terms until no longer possible
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x)) # 
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod)) # 
expand_expr = Fixpoint(PassThrough(Chain([expand_mod, flatten_mod]))) # apply flatten and distribution until no longer possible

function define_μ(iter::AbstractVector, iv::Union{Sym,Num})

    indices = map(trim_key, iter)

    μs = OrderedDict{eltype(iter), Any}()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            μs[idx] = 1
        else
            sym_name = Symbol('μ', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            term_raw = Term{Real}(sym_raw, [iv])
            μs[idx] = setmetadata(term_raw, Symbolics.VariableSource,
                                  (:momentclosure, sym_name))
        end
    end

    μs

end

define_μ(N::Int, order::Int, iv::Union{Sym,Num}) = define_μ(construct_iter_all(N, order), iv)

function define_μ(N::Int, order::Int, iter = construct_iter_all(N, order)) 
    @parameters t
    return define_μ(iter, value(t))
end


function define_M(iter::AbstractVector, iv::Union{Sym,Num})

    indices = map(trim_key, iter)

    Ms = OrderedDict{eltype(iter), Any}()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            Ms[idx] = 1
        elseif sum(idx) == 1
            Ms[idx] = 0
        else
            sym_name = Symbol('M', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            term_raw = Term{Real}(sym_raw, [iv])
            Ms[idx] = setmetadata(term_raw, Symbolics.VariableSource,
                                  (:momentclosure, sym_name))
        end
    end

    Ms

end

define_M(N::Int, order::Int, iv::Union{Sym,Num}) = define_M(construct_iter_all(N, order), iv)

function define_M(N::Int, order::Int, iter = construct_iter_all(N, order)) 
    @parameters t
    return define_M(iter, value(t))
end

#=
function extract_variables1(eqs::Array{Equation, 1}, μ, M=[])
    
    vars = vcat(values(μ)..., values(M)...)
    # extract variables from rhs of each equation
    eq_vars = unique(vcat(get_variables.(eqs)...))
    # need this as get_variables does not extract var from `Differential(t)(var(t))`
    diff_vars = [var_from_nested_derivative(eq.lhs)[1] for eq in eqs]
    # filter out the unique ones
    eq_vars = unique(vcat(eq_vars..., diff_vars...))
    # this should preserve the correct ordering

    intersect!(vars, eq_vars)

end
=#

function extract_variables(eqs::Array{Equation, 1}, μ, M=[])
    
    vars = vcat(values(μ)..., values(M)...)
    # extract variables from rhs of each equation
    eq_vars = unique(vcat(get_variables.(eqs)...))
    # get_variables changes the metadata so have to be careful here...
    intersect!(vars, eq_vars)
    # need this as get_variables does not extract var from `Differential(t)(var(t))`
    diff_vars = [var_from_nested_derivative(eq.lhs)[1] for eq in eqs]
    # filter out the unique ones
    unique(vcat(diff_vars, vars))
    # the correct ordering *should* be preserved

end

## Set of functions to deconstruct polynomial propensities ##

#=
    Consider a polynomial propensity function: a₁ = x^2 + c₁*y²*x,
    where x(t) and y(t) are the molecule numbers variables and c₁ is reaction parameter
    First we need to split a₁ into separate terms, i.e. x^2 and c₁*y²*x. Then we
    determine the independent multiplication factor in each term (1 and c₁). Finally,
    we obtain the power each variable is raised to in each term (x to power 2 in term 1;
    x to power 1 and y to power 2 in term 2). Having this information we can proceed
    in constructing raw moment equations. IF any propensity function is non-polynomial
    then the function `polynomial_propensities` will throw an error.
=#

isvar(x, vars) = any(isequal(x), vars)

"""
    Check that the given expression does NOT depend on the given variables `vars` (expr is constant wrt. vars)
    A variable here is defined as a function of the independent variable `iv`, e.g,  X(t) is variable, where t
    is the independent variable
"""

isconstant(expr, vars, iv) = !istree(expr) || (!isvar(expr, vars) && all(arg -> isconstant(arg, vars, iv), arguments(expr)))

function split_factor(expr::Pow, iv, smap, vars)
    base, exp = arguments(expr)
    (exp isa Int && exp >= 0) || error("Unexpected exponent: $expr")

    factor, powers = split_factor(base, iv, smap, vars)
    factor ^ exp, powers .* exp
end

function split_factor(expr::Mul, iv, smap, vars)
    powers = zeros(Int, length(smap))
    factor = 1

    for arg in arguments(expr)
        factor_arg, power_arg = split_factor(arg, iv, smap, vars)
        factor *= factor_arg
        powers .+= power_arg
    end

    factor, powers
end

function split_factor(expr::Div, iv, smap, vars)
    num, denom = arguments(expr)
    isconstant(denom, vars, iv) || error("The denominator $denom is not constant.")
    
    factor, powers = split_factor(num, iv, smap, vars)
    
    factor / denom, powers
end

function split_factor(expr, iv, smap, vars) 
    if isconstant(expr, vars, iv)
        expr, zeros(length(vars))
    elseif isvar(expr, vars)
        1, map(isequal(expr), vars)
    else
        error("Expression $expr could not be parsed correctly!")
    end
end

function polynomial_propensity(expr::Div, iv, smap, vars)
    num, denom = arguments(expr)
    isconstant(denom, vars, iv) || error("The denominator $denom is not constant.")
    factors, powers = polynomial_propensity(num, iv, smap, vars)
    factors ./ denom, powers
end

function polynomial_propensity(expr::Add, iv, smap, vars)
    factors = []
    powers = Vector{Int}[]
    
    for term in arguments(expr)
        factor_term, power_term = try
            split_factor(term, iv, smap, vars)
        catch e
            error("Propensity function $term is non-polynomial? \n" * string(e))
        end
        push!(factors, factor_term)
        push!(powers, power_term)
    end

    factors, powers
end

function polynomial_propensity(expr, iv, smap, vars)
    factor, powers = try
        split_factor(expr, iv, smap, vars)
    catch e
        error("Propensity function $expr is non-polynomial? \n" * string(e))
    end

    [ factor ], [ powers ]
end

function polynomial_propensities(arr::AbstractArray, iv::Sym, smap::AbstractDict)
    # Do we need sort here?
    vars = [ x for (x,_) in Base.sort(collect(smap), by=x->x[2]) ]

    all_factors = Array{Vector}(undef, size(arr))
    all_powers = Array{Vector{Vector{Int}}}(undef, size(arr))

    for (rind, expr) in enumerate(arr)
        expr = expand(expr)
        all_factors[rind], all_powers[rind] = polynomial_propensity(expr, iv, smap, vars)
    end

    max_power = maximum(sum.(vcat(all_powers...)))
    all_factors, all_powers, max_power
end

#=
function degree(p::Num, sym::AbstractVector) 
    p = value(p)
    sym = Set(value.(sym))
    if p isa Number 
        return 0
    elseif p ∈ sym
        return 1
    elseif p isa Symbolic
        return degree(p, sym)
    end
end
degree(p::Sym, sym::Set) = Int(p ∈ sym)
degree(p::Sym, sym::AbstractVector) = degree(p, Set(value.(sym)))
degree(p::Term, sym::Set) = Int(p ∈ sym)
degree(p::Term, sym::AbstractVector) = degree(p, Set(value.(sym)))
degree(p::Symbolics.Add, sym::AbstractVector) = degree(p, Set(value.(sym)))
degree(p::Symbolics.Mul, sym::AbstractVector) = degree(p, Set(value.(sym)))
degree(p::Symbolics.Pow, sym::AbstractVector) = degree(p, Set(value.(sym)))

poly_subs(ex::Union{Num, Number, Symbolics.Mul, Symbolics.Add}, subs::AbstractDict, ps::AbstractArray, flag::Bool = false) = poly_subs(ex, subs, Set(ps), flag) 
poly_subs(ex::Num, subs::AbstractDict, ps::Set = Set(), flag::Bool = false) = poly_subs(value(expand(ex)), subs, ps, flag) 
poly_subs(ex::Number, ::AbstractDict, ps::Set = Set(), ::Bool = false) = ex 
poly_subs(ex::Sym, subs::AbstractDict, ::Set, ::Bool = false) = haskey(subs, ex) ? subs[ex] : ex
poly_subs(ex::Term, subs::AbstractDict, ::Set, ::Bool = false) = haskey(subs, ex) ? subs[ex] : ex
poly_subs(ex::Symbolics.Pow, subs::AbstractDict, ::Set, flag::Bool = false) = substitute(ex, subs)
poly_subs(ex::Symbolics.Add, subs::AbstractDict, ps::Set = Set(), flag::Bool = false) = sum(poly_subs(term, subs, ps, flag) for term in arguments(ex))

function poly_subs(ex::Symbolics.Mul, subs::AbstractDict, ps = Set(), flag::Bool = false) 
    mono = 1
    coeff = ex.coeff
    if flag # flag controls whether there are variables not being substituted -> those are shifted to the coefficient
        for (base, exponent) in pairs(ex.dict)
            if base ∉ ps && haskey(subs, base) # very hacky but does the job ...
                mono *= base^exponent
            else
                coeff *= base^exponent
            end
        end
    else
        for (base, exponent) in pairs(ex.dict)
            if base ∉ ps 
                mono *= base^exponent
            else
                coeff *= base^exponent
            end
        end
    end
    return haskey(subs, mono) ? coeff*subs[mono] : ex
end


function taylor_expand(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, vars::AbstractVector, ref_point::AbstractVector, order::Int, aux_flag::Bool = true)
    @assert order >= 0 "Order of Taylor expansion must be non-negative"
    @assert length(ref_point) == length(vars) "Variables and reference point must have the same dimension" 
    N = length(vars)
    if aux_flag
        Δ = vars .- ref_point 
    else
        Δ = vars
    end
    sub_ref_point = Dict(vars[i] => ref_point[i] for i in 1:N)
    f_taylor = substitute(f, sub_ref_point)
    iter_all = filter(x -> sum(x) > 0, MomentClosure.construct_iter_all(N, order))
    for iter in iter_all
        f_taylor += substitute(derivative(f, vars, iter), sub_ref_point)*prod(Δ .^ iter)/MomentClosure.fact(iter)
    end
    return f_taylor
end

function derivative(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, vars::AbstractVector, order::Tuple)
    for i in 1:length(vars)
        if order[i] != 0
            f = nth_differential(f, vars[i], order[i])
        end
    end
    return f
end

function nth_differential(f::Union{Num, Term, Symbolics.Mul, Symbolics.Add, Symbolics.Pow}, var::Union{Num, Sym, Term}, n::Int) 
    for i in 1:n
        f = expand_derivatives(Differential(var)(f), true)
    end
    return f
end
=#
