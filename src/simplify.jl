import SymbolicUtils: <ₑ
using AbstractAlgebra.Generic: MPoly
<ₑ(a::MPoly, b::MPoly) = false

let
	using SymbolicUtils: operation, arguments, symtype, similarterm, _promote_symtype, term
    using Symbolics: Sym, Add, Mul, Pow

    mc_similarterm(args...) = similarterm(args...)
    mc_similarterm(t, f, args, symtype) = f(args...)
    mc_similarterm(t, f, args) = mc_similarterm(t, f, args, _promote_symtype(f, args))
    mc_similarterm(::Term, f, args, symtype=nothing) = term(f, args...; type=symtype)

    function mc_similarterm(p::Union{Mul, Add, Pow}, f, args)
        if f === (+)
            T = _promote_symtype(f, args)
            Add(T, makeadd(1, 0, args...)...)
        elseif f == (*)
            T = _promote_symtype(f, args)
            Mul(T, makemul(1, args...)...)
        elseif f == (^) && length(args) == 2
            Pow(args...)
        else
            f(args...)
        end
    end

	function labels! end

	# Turn a Term into a multivariate polynomial
	function labels!(dicts, t::Sym)
		sym2term, term2sym = dicts
		if !haskey(term2sym, t)
			sym2term[t] = t
			term2sym[t] = t
		end
		return t
	end

	function labels!(dicts, t)
		if t isa Integer
			return t
		elseif istree(t) && (operation(t) == (*) || operation(t) == (+) || operation(t) == (-))
			tt = arguments(t)
			return mc_similarterm(t, operation(t), map(x->labels!(dicts, x), tt), symtype(t))
		elseif istree(t) && operation(t) == (^) && length(arguments(t)) > 1 && isnonnegint(arguments(t)[2])
			return mc_similarterm(t, operation(t), map(x->labels!(dicts, x), arguments(t)), symtype(t))
		else
			sym2term, term2sym = dicts
			if haskey(term2sym, t)
				return term2sym[t]
			end
			if istree(t)
				tt = arguments(t)
				sym = Sym{symtype(t)}(gensym(nameof(operation(t))))
				dicts2 = _dicts(dicts[2])
				sym2term[sym] = mc_similarterm(t, operation(t),
											map(x->to_mpoly(x, dicts)[1], arguments(t)),
											symtype(t))
			else
				sym = Sym{symtype(t)}(gensym("literal"))
				sym2term[sym] = t
			end

			term2sym[t] = sym

			return sym
		end
	end

	ismpoly(x) = x isa MPoly || x isa Integer
	isnonnegint(x) = x isa Integer && x >= 0

	_dicts(t2s=OrderedDict{Any, Sym}()) = (OrderedDict{Sym, Any}(), t2s)

	using SymbolicUtils.Rewriters
	using AbstractAlgebra.Generic: MPoly, ZZ, PolynomialRing, exponent_vector
	using AbstractAlgebra: ismonomial, symbols

	let
		mpoly_preprocess = [@rule(identity(~x) => ~x)
							@rule(zero(~x) => 0)
							@rule(one(~x) => 1)]

		mpoly_rules = [@rule(~x::ismpoly - ~y::ismpoly => ~x + -1 * (~y))
					   @rule(-(~x) => -1 * ~x)
					   @acrule(~x::ismpoly + ~y::ismpoly => ~x + ~y)
					   @rule(+(~x) => ~x)
					   @acrule(~x::ismpoly * ~y::ismpoly => ~x * ~y)
					   @rule(*(~x) => ~x)
					   @rule((~x::ismpoly)^(~a::isnonnegint) => (~x)^(~a))]
		global const MPOLY_CLEANUP = Fixpoint(Postwalk(PassThrough(RestartedChain(mpoly_preprocess))))
		MPOLY_MAKER = Fixpoint(Postwalk(PassThrough(RestartedChain(mpoly_rules))))

		global to_mpoly
		function to_mpoly(t, dicts=_dicts())
			# term2sym is only used to assign the same
			# symbol for the same term -- in other words,
			# it does common subexpression elimination
			t = MPOLY_CLEANUP(t)
			sym2term, term2sym = dicts
			labeled = labels!((sym2term, term2sym), t)

			if isempty(sym2term)
				return MPOLY_MAKER(labeled), Dict{Sym,Any}()
			end

			ks = Base.sort(collect(keys(sym2term)), lt=<ₑ)
			R, vars = PolynomialRing(ZZ, String.(nameof.(ks)))

			replace_with_poly = Dict{Sym,MPoly}(zip(ks, vars))
			t_poly = substitute(labeled, replace_with_poly, fold=false)
			MPOLY_MAKER(t_poly), sym2term
		end
	end

	function to_term(reference, x, dict)
		syms = Dict(zip(nameof.(keys(dict)), keys(dict)))
		dict = copy(dict)
		for (k, v) in dict
			dict[k] = _to_term(reference, v, dict, syms)
		end
		_to_term(reference, x, dict, syms)
	end

	function _to_term(reference, x::MPoly, dict, syms)

		function mul_coeffs(exps, ring)
			l = length(syms)
			ss = symbols(ring)
			monics = [e == 1 ? syms[ss[i]] : syms[ss[i]]^e for (i, e) in enumerate(exps) if !iszero(e)]
			if length(monics) == 1
				return monics[1]
			elseif length(monics) == 0
				return 1
			else
				return mc_similarterm(reference, *, monics, symtype(reference))
			end
		end

		monoms = [mul_coeffs(exponent_vector(x, i), x.parent) for i in 1:x.length]
		if length(monoms) == 0
			return 0
		elseif length(monoms) == 1
			t = !isone(x.coeffs[1]) ?  monoms[1] * Int(x.coeffs[1]) : monoms[1]
		else
			t = mc_similarterm(reference,
							+,
							map((x,y)->isone(y) ? x : Int(y)*x,
								monoms, x.coeffs[1:length(monoms)]),
							symtype(reference))
		end

		substitute(t, dict, fold=false)
	end

	function _to_term(reference, x, dict, vars)
		if istree(x)
			t=mc_similarterm(x, operation(x), _to_term.((reference,), arguments(x), (dict,), (vars,)), symtype(x))
		else
			if haskey(dict, x)
				return dict[x]
			else
				return x
			end
		end
	end


	global mc_expand
	"""
		mc_expand(expr)
	Expand expressions by distributing multiplication over addition.
	`a*(b+c)` becomes `ab+ac`. `expand` uses [AbstractAlgebra.jl](https://nemocas.github.io/AbstractAlgebra.jl/latest/) to construct
	dense Multi-variate polynomial to do this very fast.
	"""
	function mc_expand(x)
		to_term(x, to_mpoly(x)...)
	end
end
# Taken from SymbolicUtils:src/simplify_rules.jl

let
    using SymbolicUtils.Rewriters
    using SymbolicUtils: isnotflat, flatten_term, sort_args, is_literal_number,
                         hasrepeats, merge_repeats, _iszero, _isone, _isreal,
                         needs_sorting, _isinteger, istree, is_operation, pow,
                         @acrule
    
    PLUS_RULES = [
        @rule(~x::isnotflat(+) => flatten_term(+, ~x))
        @rule(~x::needs_sorting(+) => sort_args(+, ~x))
        SymbolicUtils.@ordered_acrule(~a::is_literal_number + ~b::is_literal_number => ~a + ~b)

        @acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
        @acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...))
        @acrule(*(~~x, ~α) + *(~~x, ~β) => *(~α + ~β, (~~x)...))

        @acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
        @acrule(*(~α::is_literal_number, ~x) + ~x => *(~α + 1, ~x))
        @rule(+(~~x::hasrepeats) => +(merge_repeats(*, ~~x)...))

        SymbolicUtils.@ordered_acrule((~z::_iszero + ~x) => ~x)
        @rule(+(~x) => ~x)
    ]

    TIMES_RULES = [
        @rule(~x::isnotflat(*) => flatten_term(*, ~x))
        @rule(~x::needs_sorting(*) => sort_args(*, ~x))

        SymbolicUtils.@ordered_acrule(~a::is_literal_number * ~b::is_literal_number => ~a * ~b)
        @rule(*(~~x::hasrepeats) => *(merge_repeats(^, ~~x)...))

        @acrule((~y)^(~n) * ~y => (~y)^(~n+1))
        SymbolicUtils.@ordered_acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

        SymbolicUtils.@ordered_acrule((~z::_isone  * ~x) => ~x)
        SymbolicUtils.@ordered_acrule((~z::_iszero *  ~x) => ~z)
        @rule(*(~x) => ~x)
    ]


    POW_RULES = [
        @rule(^(*(~~x), ~y::_isinteger) => *(map(a->pow(a, ~y), ~~x)...))
        @rule((((~x)^(~p::_isinteger))^(~q::_isinteger)) => (~x)^((~p)*(~q)))
        @rule(^(~x, ~z::_iszero) => 1)
        @rule(^(~x, ~z::_isone) => ~x)
        @rule(inv(~x) => ~x ^ -1)
    ]

    ASSORTED_RULES = [
        @rule(identity(~x) => ~x)
        @rule(-(~x) => -1*~x)
        @rule(-(~x, ~y) => ~x + -1(~y))
        @rule(~x::_isone \ ~y => ~y)
        @rule(~x \ ~y => ~y / (~x))
        @rule(~x / ~y => ~x * pow(~y, -1))
        @rule(one(~x) => one(symtype(~x)))
        @rule(zero(~x) => zero(symtype(~x)))
        @rule(conj(~x::_isreal) => ~x)
        @rule(real(~x::_isreal) => ~x)
        @rule(imag(~x::_isreal) => zero(symtype(~x)))
        @rule(ifelse(~x::is_literal_number, ~y, ~z) => ~x ? ~y : ~z)
    ]
   

 
    number_rules = [If(istree, Chain(ASSORTED_RULES)),
                    If(is_operation(+),
                       Chain(PLUS_RULES)),
                    If(is_operation(*),
                       Chain(TIMES_RULES)),
                    If(is_operation(^),
                       Chain(POW_RULES))] |> RestartedChain

    base_simplifier = Fixpoint(Postwalk(number_rules))
    expand_simplifier = If(istree, Chain((mc_expand, base_simplifier)))
    noexpand_simplifier = If(istree, base_simplifier)

    # Should be faster than SymbolicUtils.simplify as it skips
    # rules that are not really relevant (eg. trig identities)
    global function mc_simplify(expr; expand=false)
        simplify(expr; rewriter=expand ? expand_simplifier : noexpand_simplifier)
    end
end
