"""
Julia wrappers around adaptive multidimensional integration routines
from the [C Cubature Package](https://github.com/stevengj/cubature).

See the [README file](https://github.com/stevengj/Cubature.jl/blob/master/README.md) for documentation.   Note that we define
several variants of the functions in order to expose the full
functionality of the Cubature Package while still providing
simple interfaces to the more basic functionality (1d and >1d
integrals of scalar functions).
"""
module Cubature

using Cubature_jll

export hcubature, pcubature, hcubature_v, pcubature_v,
    hquadrature, pquadrature, hquadrature_v, pquadrature_v

# constants from cubature.h
const INDIVIDUAL = Int32(0)
const PAIRED = Int32(1)
const L2 = Int32(2)
const L1 = Int32(3)
const LINF = Int32(4)

const SUCCESS = Int32(0)
const FAILURE = Int32(1)

# type to distinguish cubature error codes from thrown exceptions
struct NoError <: Exception end # used for integrand_error when nothing thrown

mutable struct IntegrandData{F}
    integrand_func::F
    integrand_error::Any
    IntegrandData(f::F) where{F} = new{F}(f, NoError())
end

# C cubature code is not interrupt-safe (would leak memory), so
# use sigatomic_begin/end to defer ctrl-c handling until Julia code
import Base.sigatomic_begin, Base.sigatomic_end

# C-callable integrand wrappers around the user's Julia integrand
for fscalar in (false, true) # whether the integrand is a scalar
    for vectorized in (false, true) # whether multiple x are passed at once
        for xscalar in (false, true) # whether x is a scalar
            f = Symbol(string(fscalar ? :sintegrand : :integrand,
                              vectorized ? "_v" : ""))

            if xscalar
                f = Symbol(string("q",f))
                if vectorized
                    xex = :(unsafe_wrap(Array, x_, (convert(Int, npt),)))
                else
                    xex = :(unsafe_load(x_))
                end
            else
                if vectorized
                    xex = :(unsafe_wrap(Array, x_, (convert(Int, ndim),convert(Int, npt))))
                else
                    xex = :(unsafe_wrap(Array, x_, (convert(Int, ndim),)))
                end
            end

            if fscalar
                if vectorized
                    vex = :(unsafe_wrap(Array, fval_, (convert(Int, npt),)))
                    ex = :(func($xex, $vex))
                else
                    ex = :(unsafe_store!(fval_, func($xex)))
                end
            else
                if vectorized
                    vex = :(unsafe_wrap(Array, fval_, (convert(Int, fdim),convert(Int, npt))))
                else
                    vex = :(unsafe_wrap(Array, fval_, (convert(Int, fdim),)))
                end
                ex = :(func($xex, $vex))
            end

            body = quote
                func = d.integrand_func
                try
                    reenable_sigint() do
                        $ex
                    end
                    return SUCCESS
                catch e
                    d.integrand_error = e
                    return FAILURE
                end
            end

            if vectorized
                @eval function $f(ndim::UInt32, npt::UInt,
                                  x_::Ptr{Float64}, d,
                                  fdim::UInt32, fval_::Ptr{Float64})
                    $body
                end
            else
                @eval function $f(ndim::UInt32, x_::Ptr{Float64}, d,
                                  fdim::UInt32, fval_::Ptr{Float64})
                    $body
                end
            end
        end
    end
end

function integrands(d::D, xscalar, fscalar, vectorized) where {D}
    if xscalar
        if fscalar
            if vectorized
                return @cfunction(qsintegrand_v, Int32, (UInt32, UInt, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            else
                return @cfunction(qsintegrand, Int32, (UInt32, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            end
        else
            if vectorized
                return @cfunction(qintegrand_v, Int32, (UInt32, UInt, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            else
                return @cfunction(qintegrand, Int32, (UInt32, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            end
        end
    else
        if fscalar
            if vectorized
                return @cfunction(sintegrand_v, Int32, (UInt32, UInt, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            else
                return @cfunction(sintegrand, Int32, (UInt32, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            end
        else
            if vectorized
                return @cfunction(integrand_v, Int32, (UInt32, UInt, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            else
                return @cfunction(integrand, Int32, (UInt32, Ptr{Float64}, Ref{D}, UInt32, Ptr{Float64}))
            end
        end
    end
end

# low-level routine, not to be called directly by user
function cubature(xscalar::Bool, fscalar::Bool,
                  vectorized::Bool, padaptive::Bool,
                  fdim::Integer, f::F, # Force specialization on F
                  xmin_, xmax_,
                  reqRelError::Real, reqAbsError::Real, maxEval::Integer,
                  error_norm::Integer) where {F}
    dim = length(xmin_)
    if xscalar && dim != 1
        throw(ArgumentError("quadrature routines are for 1d only"))
    end
    if dim != length(xmax_)
        throw(ArgumentError("dimension mismatch between xmin and xmax"))
    end
    if fscalar && fdim != 1
        throw(ArgumentError("fdim = $fdim for scalar integrand"))
    end
    if fdim < 0
        throw(ArgumentError("fdim = $fdim < 0 is invalid"))
    end
    xmin = Float64[xmin_...]
    xmax = Float64[xmax_...]
    val = Vector{Float64}(undef, fdim)
    err = Vector{Float64}(undef, fdim)
    d = IntegrandData(f)
    fwrap = integrands(d, xscalar, fscalar, vectorized)
    # ccall's first arg needs to be a constant expression, so
    # we have to put the if statements outside the ccalls rather
    # than inside, unfortunately
    disable_sigint() do
        if padaptive
            if vectorized
                ret = ccall((:pcubature_v,libcubature), Int32,
                            (UInt32, Ptr{Cvoid}, Any,
                             UInt32, Ptr{Float64}, Ptr{Float64},
                             UInt, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, d, dim, xmin, xmax,
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            else
                ret = ccall((:pcubature,libcubature), Int32,
                            (UInt32, Ptr{Cvoid}, Any,
                             UInt32, Ptr{Float64}, Ptr{Float64},
                             UInt, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, d, dim, xmin, xmax,
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            end
        else
            if vectorized
                ret = ccall((:hcubature_v,libcubature), Int32,
                            (UInt32, Ptr{Cvoid}, Any,
                             UInt32, Ptr{Float64}, Ptr{Float64},
                             UInt, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, d, dim, xmin, xmax,
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            else
                ret = ccall((:hcubature,libcubature), Int32,
                            (UInt32, Ptr{Cvoid}, Any,
                             UInt32, Ptr{Float64}, Ptr{Float64},
                             UInt, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, d, dim, xmin, xmax,
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            end
        end
        if ret == SUCCESS
            return (val, err)
        else
            e = d.integrand_error
            throw(isa(e, NoError) ? ErrorException("hcubature error $ret") : e)
        end
    end
end

for f in (:hcubature, :pcubature, :hquadrature, :pquadrature)
    for vectorized in (true, false)
        g = Symbol(string(f, vectorized ? "_v" : ""))
        xscalar = f == :pquadrature || f == :hquadrature
        padaptive = f == :pcubature || f == :pquadrature
        @eval begin
            # vector integrands
            $g(fdim::Integer, f::Function, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0, error_norm=INDIVIDUAL) = cubature($xscalar, false, $vectorized, $padaptive, fdim, f, xmin, xmax, reltol, abstol, maxevals, error_norm)

            # scalar integrands
            $g(f::Function, xmin, xmax; reltol=1e-8, abstol=0, maxevals=0) = begin (val,err) = cubature($xscalar, true, $vectorized, $padaptive, 1, f, xmin, xmax, reltol, abstol, maxevals, INDIVIDUAL); (val[1], err[1]); end
        end
    end
end

end # module
