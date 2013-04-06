module Cubature

# Julia wrappers around adaptive multidimensional integration routines
# from the Cubature Package, at:
#              http://ab-initio.mit.edu/cubature
# by Steven G. Johnson
#
# See the README file for documentation.   Note that we define
# several variants of the functions in order to expose the full
# functionality of the Cubature Package while still providing
# simple interfaces to the more basic functionality (1d and >1d
# integrals of scalar functions).

export hcubature, pcubature, hcubature_v, pcubature_v, 
       hquadrature, pquadrature, hquadrature_v, pquadrature_v

const libcubature = Pkg.dir("Cubature", "deps", "libcubature")

# constants from cubature.h
const INDIVIDUAL = int32(0)
const PAIRED = int32(1)
const L2 = int32(2)
const L1 = int32(3)
const LINF = int32(4)

const SUCCESS = int32(0)
const FAILURE = int32(1)

# type to distinguish cubature error codes from thrown exceptions
type NoError <: Exception end # used for integrand_error when nothing thrown

type IntegrandData
    integrand_func::Function
    integrand_error::Any
    IntegrandData(f::Function) = new(f, NoError())
end

# C-callable integrand wrappers around the user's Julia integrand
for fscalar in (false, true) # whether the integrand is a scalar
    for vectorized in (false, true) # whether multiple x are passed at once
        for xscalar in (false, true) # whether x is a scalar
            f = symbol(string(fscalar ? :sintegrand : :integrand,
                              vectorized ? "_v" : ""))

            if xscalar
                f = symbol(string("q",f))
                if vectorized
                    xex = :(pointer_to_array(x_, (int(npt),)))
                else
                    xex = :(unsafe_ref(x_))
                end
            else
                if vectorized
                    xex = :(pointer_to_array(x_, (int(ndim),int(npt))))
                else
                    xex = :(pointer_to_array(x_, (int(ndim),)))
                end
            end

            if fscalar
                if vectorized
                    vex = :(pointer_to_array(fval_, (int(npt),)))
                    ex = :(func($xex, $vex))
                else
                    ex = :(unsafe_assign(fval_, func($xex)))
                end
            else
                if vectorized
                    vex = :(pointer_to_array(fval_, (int(fdim),int(npt))))
                else
                    vex = :(pointer_to_array(fval_, (int(fdim),)))
                end
                ex = :(func($xex, $vex))
            end

            body = quote
                d = unsafe_pointer_to_objref(d_)::IntegrandData
                func = d.integrand_func
                try
                    $ex
                    return SUCCESS
                catch e
                    d.integrand_error = e
                    return FAILURE
                end
            end

            if vectorized
                @eval function $f(ndim::Uint32, npt::Uint,
                                  x_::Ptr{Float64}, d_::Ptr{Void},
                                  fdim::Uint32, fval_::Ptr{Float64})
                    $body
                end
            else
                @eval function $f(ndim::Uint32, x_::Ptr{Float64},d_::Ptr{Void},
                                  fdim::Uint32, fval_::Ptr{Float64})
                    $body
                end
            end
        end
    end
end

cf(f,v) = cfunction(f, Int32, v ? (Uint32, Uint, Ptr{Float64}, Ptr{Void},
                                   Uint32, Ptr{Float64}) :
                                  (Uint32, Ptr{Float64}, Ptr{Void},
                                   Uint32, Ptr{Float64}))

# (xscalar, fscalar, vectorized) => function
const integrands = [ (false,false,false) => cf(integrand,false),
                     (false,false, true) => cf(integrand_v,true),
                     (false, true,false) => cf(sintegrand,false),
                     (false, true, true) => cf(sintegrand_v,true),
                     ( true,false,false) => cf(qintegrand,false),
                     ( true,false, true) => cf(qintegrand_v,true),
                     ( true, true,false) => cf(qsintegrand,false),
                     ( true, true, true) => cf(qsintegrand_v,true) ]

# low-level routine, not to be called directly by user
function cubature(xscalar::Bool, fscalar::Bool,
                  vectorized::Bool, padaptive::Bool,
                  fdim::Integer, f::Function, 
                  xmin_, xmax_, 
                  reqRelError::Real, reqAbsError::Real, maxEval::Integer,
                  error_norm::Integer)
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
    val = Array(Float64, fdim)
    err = Array(Float64, fdim)
    fwrap = integrands[(xscalar,fscalar,vectorized)]
    d = IntegrandData(f)
    # ccall's first arg needs to be a constant expression, so
    # we have to put the if statements outside the ccalls rather
    # than inside, unfortunately
    if padaptive
        if vectorized
            ret = ccall((:pcubature_v,libcubature), Int32,
                        (Uint32, Ptr{Void}, Any,
                         Uint32, Ptr{Float64}, Ptr{Float64},
                         Uint, Float64, Float64, Int32,
                         Ptr{Float64}, Ptr{Float64}),
                        fdim, fwrap, d, dim, xmin, xmax, 
                        maxEval, reqAbsError, reqRelError, error_norm,
                        val, err)
        else
            ret = ccall((:pcubature,libcubature), Int32,
                        (Uint32, Ptr{Void}, Any,
                         Uint32, Ptr{Float64}, Ptr{Float64},
                         Uint, Float64, Float64, Int32,
                         Ptr{Float64}, Ptr{Float64}),
                        fdim, fwrap, d, dim, xmin, xmax, 
                        maxEval, reqAbsError, reqRelError, error_norm,
                        val, err)
        end
    else
        if vectorized
            ret = ccall((:hcubature_v,libcubature), Int32,
                        (Uint32, Ptr{Void}, Any,
                         Uint32, Ptr{Float64}, Ptr{Float64},
                         Uint, Float64, Float64, Int32,
                         Ptr{Float64}, Ptr{Float64}),
                        fdim, fwrap, d, dim, xmin, xmax, 
                        maxEval, reqAbsError, reqRelError, error_norm,
                        val, err)
        else
            ret = ccall((:hcubature,libcubature), Int32,
                        (Uint32, Ptr{Void}, Any,
                         Uint32, Ptr{Float64}, Ptr{Float64},
                         Uint, Float64, Float64, Int32,
                         Ptr{Float64}, Ptr{Float64}),
                        fdim, fwrap, d, dim, xmin, xmax, 
                        maxEval, reqAbsError, reqRelError, error_norm,
                        val, err)
        end
    end
    e = d.integrand_error
    if ret == SUCCESS
        return (val, err)
    else
        throw(isa(e, NoError) ? ErrorException("hcubature error $ret") : e)
    end
end

for f in (:hcubature, :pcubature, :hquadrature, :pquadrature)
    for vectorized in (true, false)
        g = symbol(string(f, vectorized ? "_v" : ""))
        xscalar = f == :pquadrature || f == :hquadrature
        padaptive = f == :pcubature || f == :pquadrature
        @eval begin
            # vector integrands
            $g(fdim::Integer, f::Function, xmin, xmax, reqRelError::Real, reqAbsError::Real, maxEval::Integer, error_norm::Integer) = cubature($xscalar, false, $vectorized, $padaptive, fdim, f, xmin, xmax, reqRelError, reqAbsError, maxEval, error_norm)

            $g(fdim::Integer, f::Function, xmin, xmax, reqRelError::Real, reqAbsError::Real, maxEval::Integer) = $g(fdim, f, xmin, xmax, reqRelError, reqAbsError, maxEval, INDIVIDUAL)

            $g(fdim::Integer, f::Function, xmin, xmax, reqRelError::Real, reqAbsError::Real) = $g(fdim, f, xmin, xmax, reqRelError, reqAbsError, 0)

            $g(fdim::Integer, f::Function, xmin, xmax, reqRelError::Real) = $g(fdim, f, xmin, xmax, reqRelError, 0.0)
            
            # scalar integrands
            $g(f::Function, xmin, xmax, reqRelError::Real, reqAbsError::Real, maxEval::Integer) = begin (val,err) = cubature($xscalar, true, $vectorized, $padaptive, 1, f, xmin, xmax, reqRelError, reqAbsError, maxEval, INDIVIDUAL); (val[1], err[1]); end

            $g(f::Function, xmin, xmax, reqRelError::Real, reqAbsError::Real) = $g(f, xmin, xmax, reqRelError, reqAbsError, 0)

            $g(f::Function, xmin, xmax, reqRelError::Real) = $g(f, xmin, xmax, reqRelError, 0.0)
        end
    end
end

end # module
