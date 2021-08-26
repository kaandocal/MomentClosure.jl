import SymbolicUtils: <ₑ
<ₑ(a::SymbolicUtils.Symbolic, b::Num) = false
<ₑ(a::Num, b::SymbolicUtils.Symbolic) = true

Base.isequal(a::Symbolics.Add, b::Symbolics.Add) = isequal(a.coeff, b.coeff) && isequal(a.dict, b.dict)

function gen_iter(n::Int, d::Int)
    # based on https://twitter.com/evalparse/status/1107964924024635392
    iter = NTuple{n,Int}[]
    for x in partitions(d + n, n)
        x = x .- 1
        if all(x .<= d)
            ys = Set(multiset_permutations(x, n))
            for y in ys
                push!(iter, Tuple(y))
            end
        end
    end
    iter
end

function construct_iter_all(N::Int, order::Int)
    #Construct an ordered iterator going over all moments
    # sequentially in terms of order

    iters = NTuple{N,Int}[]
    for d in 0:order
        x = Base.sort(gen_iter(N, d), rev=true)
        append!(iters, x)
    end

    iters
end

# Trim a string of form "(a, b, c, d, ...)" to "abcd..."
trim_key(expr) = filter(x -> !(isspace(x) || x == ')' || x== '(' || x==','), string(expr))

# Expand a symbolic expression (no binomial expansion)
expansion_rule_mod = @acrule ~x * +(~~ys) => sum(map(y-> ~x * y, ~~ys))
expand_mod = Fixpoint(Prewalk(PassThrough(expansion_rule_mod)))
flatten_rule_mod = @rule(~x::isnotflat(+) => flatten_term(+, ~x))
flatten_mod = Fixpoint(PassThrough(flatten_rule_mod))
expand_expr = Fixpoint(PassThrough(Chain([expand_mod, flatten_mod])))

function define_μ(N::Int, order::Int, iter=construct_iter_all(N, order))
    indices = map(trim_key, iter)

    @parameters t

    μs = OrderedDict()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            μs[idx] = 1
        else
            sym_name = Symbol('μ', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            sym = SymbolicUtils.setmetadata(sym_raw, Symbolics.VariableSource, (:momentclosure, sym_name))
            μs[idx] = Term{Real}(sym, [t])
        end
    end

    μs
end


function define_M(N::Int, order::Int, iter=construct_iter_all(N, order))
    indices = map(trim_key, iter)

    @parameters t

    Ms = OrderedDict()
    for (i, idx) in enumerate(iter)
        if sum(idx) == 0
            Ms[idx] = 1
        elseif sum(idx) == 1
            Ms[idx] = 0
        else
            sym_name = Symbol('M', join(map_subscripts(indices[i])))
            sym_raw = Sym{FnType{Tuple{Any}, Real}}(sym_name)
            sym = SymbolicUtils.setmetadata(sym_raw, Symbolics.VariableSource, 
                                            (:momentclosure, sym_name))
            Ms[idx] = Term{Real}(sym, [t])
        end
    end

    Ms
end


function extract_variables(eqs::Array{Equation, 1}, N::Int, q_order::Int)
    iters = construct_iter_all(N, q_order)
    iter_μ = filter(x -> sum(x) > 0, iters)
    iter_M = filter(x -> sum(x) > 1, iters)

    μs = values(define_μ(N, q_order, iter_μ))
    Ms = values(define_M(N, q_order, iter_M))
    vars = vcat(μs..., Ms...)
    # extract variables from rhs of each equation
    eq_vars = unique(vcat(get_variables.(eqs)...))
    # need this as get_variables does not extract var from `Differential(t)(var(t))`
    diff_vars = [var_from_nested_derivative(eq.lhs)[1] for eq in eqs]
    # filter out the unique ones
    eq_vars = unique(vcat(eq_vars..., diff_vars...))
    # this should preserve the correct ordering

    vars = intersect!(vars, eq_vars)
    vars
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

isconstant(expr::Number, vars, iv) = false

function isconstant(expr::Symbolic, vars, iv)
    # Check that the given expression does NOT depend on the given variables `vars` (expr is constant wrt. vars)
    # A variable here is defined as a function of the independent variable `iv`, e.g,  X(t) is variable, where t
    # is the independent variable
    istree(expr) || return true

    args = arguments(expr)
    if length(args) == 1
        if isequal(args[1], iv)
            return !isvar(expr, vars)
        else
            return isconstant(args[1], vars, iv)
        end
    else
        return all(arg -> isconstant(arg, vars, iv), args)
    end
end

isvar(x::Symbolic, vars) = any(var -> isequal(x, var), vars)

function extract_pwr!(powers, expr::Symbolic, smap, vars)
    args = arguments(expr)
    if isvar(args[1], vars) && args[2] isa Int && args[2] >= 0
        idx = smap[args[1]]
        if powers[idx] != 0
            error(args[1], " occurring multiple times in a monomial is unexpected. In ", expr)
        else
            powers[idx] = args[2]
        end
    else
        error("Unexpected term: ", expr)
    end
end

function extract_mul(expr::Symbolic, smap, vars, iv)
    powers = zeros(Int, length(smap))
    factors = []

    for arg in arguments(expr)
        if isconstant(arg, vars, iv)
            push!(factors, arg)
        elseif operation(arg) == ^
            extract_pwr!(powers, arg, smap, vars)
        elseif isvar(arg, vars)
            powers[smap[arg]] = 1
        else
            error("Unexpected operation: ", operation(arg), " in ", expr)
        end
    end

    factor = isempty(factors) ? 1 : prod(factors)
    powers, factor
end


function polynomial_propensities(a::Vector, rn::Union{ReactionSystem, ReactionSystemMod}; smap=speciesmap(rn))
    R = length(a)
    N = numspecies(rn)
    # Can this replaced by states(rn)?
    vars = [x for (x,y) in Base.sort(collect(smap), by=x->x[2])]
    iv = rn.iv

    all_factors = [[] for i = 1:R]
    all_powers = [Vector{Int}[] for i = 1:R]

    for (rind, expr) in enumerate(a)

        expr = expand(expr)

        if isconstant(expr, vars, iv)

            push!(all_factors[rind], expr)
            push!(all_powers[rind], zeros(Int, N))

        elseif isvar(expr, vars)

            push!(all_factors[rind], 1)
            push!(all_powers[rind], map(v -> isequal(expr, v), vars))

        elseif operation(expr) == ^ #Symbolics.Pow

            try
                powers = zeros(Int, N)
                extract_pwr!(powers, expr, smap, vars)
                push!(all_factors[rind], 1)
                push!(all_powers[rind], powers)
            catch e
                error("Propensity function ", expr, " is non-polynomial? \n" * string(e))
            end

        elseif operation(expr) == * #Symbolics.Mul

            try
                powers, factor = extract_mul(expr, smap, vars, iv)
                push!(all_powers[rind], powers)
                push!(all_factors[rind], factor)
            catch e
                error("Propensity function ", expr, " is non-polynomial?\n" * string(e))
            end

        elseif operation(expr) == + #Symbolics.Add

            for term in arguments(expr)

                if isconstant(term, vars, iv)
                    push!(all_factors[rind], term)
                    push!(all_powers[rind], zeros(Int, N))
                else
                    try
                        op = operation(term)
                        if op == ^
                            powers = zeros(Int, N)
                            extract_pwr!(powers, term, smap, vars)
                            push!(all_factors[rind], 1)
                            push!(all_powers[rind], powers)
                        elseif op == *
                            powers, factor = extract_mul(term, smap, vars, iv)
                            push!(all_factors[rind], factor)
                            push!(all_powers[rind], powers)
                        elseif isvar(term, vars)
                            powers = zeros(Int, N)
                            powers[smap[term]] = 1
                            push!(all_factors[rind], 1)
                            push!(all_powers[rind], powers)
                        else
                            error("Unexpected operation: ", op, " in ", term)
                        end
                    catch e
                        error("Propensity function ", expr, " is non-polynomial?\n" * string(e))
                    end
                end
            end

        else

            error("Expression $expr could not be parsed correctly!")

        end

    end

    max_power = maximum(sum.(vcat(all_powers...)))
    all_factors, all_powers, max_power

end
