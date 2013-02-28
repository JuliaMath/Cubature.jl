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

# global var used in place of a real closure (until closures are
# supported by cfunction).  Note that this makes cubature non-re-entrant
integrand_func = nothing

# global to save any exception returned by integrand_func
type NoError <: Exception end 
integrand_error = NoError()

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

            func = :((integrand_func::Function))

            if fscalar
                if vectorized
                    vex = :(pointer_to_array(fval_, (int(npt),)))
                    ex = :($func($xex, $vex))
                else
                    ex = :(unsafe_assign(fval_, $func($xex)))
                end
            else
                if vectorized
                    vex = :(pointer_to_array(fval_, (int(fdim),int(npt))))
                else
                    vex = :(pointer_to_array(fval_, (int(fdim),)))
                end
                ex = :($func($xex, $vex))
            end

            body = quote
                global integrand_func
                global integrand_error
                try
                    $ex
                    return SUCCESS
                catch e
                    integrand_error = e
                    return FAILURE
                end
            end

            if vectorized
                @eval function $f(ndim::Uint32, npt::Uint,
                                  x_::Ptr{Float64}, _::Ptr{Void},
                                  fdim::Uint32, fval_::Ptr{Float64})
                    $body
                end
            else
                @eval function $f(ndim::Uint32, x_::Ptr{Float64}, _::Ptr{Void},
                                  fdim::Uint32, fval_::Ptr{Float64})
                    $body
                end
            end
        end
    end
end

# (xscalar, fscalar, vectorized) => function
const integrands = [ (false,false,false) => integrand,
                     (false,false, true) => integrand_v,
                     (false, true,false) => sintegrand,
                     (false, true, true) => sintegrand_v,
                     ( true,false,false) => qintegrand,
                     ( true,false, true) => qintegrand_v,
                     ( true, true,false) => qsintegrand,
                     ( true, true, true) => qsintegrand_v ]

# low-level routine, not to be called directly by user
function cubature(xscalar::Bool, fscalar::Bool,
                  vectorized::Bool, padaptive::Bool,
                  fdim::Integer, f::Function, 
                  xmin_, xmax_, 
                  reqRelError::Real, reqAbsError::Real, maxEval::Integer,
                  error_norm::Integer)
    global integrand_func
    global integrand_error
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
    if integrand_func != nothing
        error("Nested calls to cubature are not supported")
        end
    fwrap = cfunction(integrands[(xscalar,fscalar,vectorized)],
                      Int32,
                      vectorized ? 
                      (Uint32, Uint, Ptr{Float64}, Ptr{Void},
                       Uint32, Ptr{Float64}) :
                      (Uint32, Ptr{Float64}, Ptr{Void},
                       Uint32, Ptr{Float64}))
    try
        integrand_func = f
        integrand_error = NoError()
        # ccall's first arg needs to be a constant expression, so
        # we have to put the if statements outside the ccalls rather
        # than inside, unfortunately
        if padaptive
            if vectorized
                ret = ccall((:pcubature_v,libcubature), Int32,
                            (Uint32, Ptr{Void}, Ptr{Void},
                             Uint32, Ptr{Float64}, Ptr{Float64},
                             Uint, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, C_NULL, dim, xmin, xmax, 
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            else
                ret = ccall((:pcubature,libcubature), Int32,
                            (Uint32, Ptr{Void}, Ptr{Void},
                             Uint32, Ptr{Float64}, Ptr{Float64},
                             Uint, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, C_NULL, dim, xmin, xmax, 
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            end
        else
            if vectorized
                ret = ccall((:hcubature_v,libcubature), Int32,
                            (Uint32, Ptr{Void}, Ptr{Void},
                             Uint32, Ptr{Float64}, Ptr{Float64},
                             Uint, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, C_NULL, dim, xmin, xmax, 
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            else
                ret = ccall((:hcubature,libcubature), Int32,
                            (Uint32, Ptr{Void}, Ptr{Void},
                             Uint32, Ptr{Float64}, Ptr{Float64},
                             Uint, Float64, Float64, Int32,
                             Ptr{Float64}, Ptr{Float64}),
                            fdim, fwrap, C_NULL, dim, xmin, xmax, 
                            maxEval, reqAbsError, reqRelError, error_norm,
                            val, err)
            end
        end
        e = integrand_error
        if ret == SUCCESS
            return (val, err)
        else
            throw(isa(e, NoError) ? ErrorException("hcubature error $ret") : e)
        end
    finally
        # clear global refs so as not to spoil garbage collection
        # and also to let subsequent calls know that we are not nested
        integrand_func = nothing
        integrand_error = NoError()
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
