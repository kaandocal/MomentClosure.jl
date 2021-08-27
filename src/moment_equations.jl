abstract type MomentEquations{N} end

"""
$(TYPEDEF)

(N = number of species in the system)

Raw moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct RawMomentEquations{N} <: MomentEquations{N}
    """[`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/)
    consisting of the time-evolution equations of raw moments."""
    odes::ODESystem
    """Symbolic variables defining the raw moments."""
    μ::Dict
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector{NTuple{N,Int}}
    iter_m::Vector{NTuple{N,Int}}
    iter_q::Vector{NTuple{N,Int}}
end

function RawMomentEquations{N}(odes, μ, m_order, q_order, iter_all, 
                               iter_m=get_iter_m(iter_all, N, m_order)) where {N}
    RawMomentEquations{N}(odes, μ, m_order, q_order, iter_all, iter_m,
                          get_iter_q(iter_all, N, m_order, q_order))
end

"""
$(TYPEDEF)

Central moment equations generated for the given system plus a number of
helper parameters (used internally).

# Fields
$(FIELDS)
"""
struct CentralMomentEquations{N} <: MomentEquations{N}
    """[`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/)
    consisting of the time-evolution equations of central moments."""
    odes::ODESystem
    """Symbolic variables defining the means."""
    μ::Dict
    """Symbolic variables defining the central moments."""
    M::Dict
    """Order of moment equations."""
    m_order::Int
    """Expansion order."""
    q_order::Int
    """Vector of all index combinations up to `q_order`."""
    iter_all::Vector{NTuple{N,Int}}
    iter_m::Vector{NTuple{N,Int}}
    iter_q::Vector{NTuple{N,Int}}
end

function CentralMomentEquations{N}(odes, μ, M, m_order, q_order, iter_all,
                                   iter_m=get_iter_m(iter_all, N, m_order)) where {N}
    CentralMomentEquations{N}(odes, μ, M, m_order, q_order, iter_all, iter_m,
                              get_iter_q(iter_all, N, m_order, q_order))
end

"""
$(TYPEDEF)

Closed moment equations and the corresponding closure functions.

# Fields
$(FIELDS)
"""
struct ClosedMomentEquations{N,MET <: MomentEquations{N}} <: MomentEquations{N}
    """[`ModelingToolkit.ODESystem`](https://mtk.sciml.ai/stable/systems/ODESystem/)
    consisting of the time-evolution equations of *closed* moments."""
    odes::ODESystem
    """Dictionary of moment closure functions for each higher order moment."""
    closure::OrderedDict
    """Original raw or central moment equations (before closure was applied)."""
    open_eqs::MET
end

# a basic wrapper
function SciMLBase.ODEProblem(eqs::MomentEquations, args...; kwargs...)
    ODEProblem(eqs.odes, args...; kwargs...)
end

function Base.nameof(eqs::MomentEquations)
    nameof(eqs.odes)
end

get_iter_all(eqs::Union{RawMomentEquations,CentralMomentEquations}) = eqs.iter_all

# Iterate through all monomials of degree 1
function get_iter_1(eqs::Union{RawMomentEquations{N},CentralMomentEquations{N}}) where {N} 
    get_iter_1(get_iter_all(eqs), N)
end

# Iterate through all monomials of degree 1 < deg <= m
function get_iter_m(eqs::Union{RawMomentEquations,CentralMomentEquations}) 
    eqs.iter_m
end

# Iterate through all monomials of degree m < deg <= q
function get_iter_q(eqs::Union{RawMomentEquations,CentralMomentEquations})
    eqs.iter_q
end

function get_iter_M(eqs::Union{RawMomentEquations,CentralMomentEquations})
    vcat(get_iter_m(eqs), get_iter_q(eqs))
end

get_iter_1(iter_all, N) = @view iter_all[2:N+1]
get_iter_m(iter_all, N, m) = @view iter_all[N+2:sum(d->binomial(N-1+d,N-1), 0:m)]
get_iter_q(iter_all, N, m, q) = @view iter_all[sum(d->binomial(N-1+d,N-1), 0:m)+1:end]
get_iter_M(iter_all, N) = @view iter_all[N+2:end]

get_iter_all(eqs::ClosedMomentEquations) = get_iter_all(eqs.open_eqs)
get_iter_1(eqs::ClosedMomentEquations) = get_iter_1(eqs.open_eqs)
get_iter_m(eqs::ClosedMomentEquations) = get_iter_m(eqs.open_eqs)
get_iter_q(eqs::ClosedMomentEquations) = get_iter_q(eqs.open_eqs)
get_iter_M(eqs::ClosedMomentEquations) = get_iter_M(eqs.open_eqs)
