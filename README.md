[![Build Status](https://travis-ci.org/JuliaMath/Cubature.jl.svg?branch=master)](https://travis-ci.org/JuliaMath/Cubature.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/0gsydibwwx59ruu7?svg=true)](https://ci.appveyor.com/project/StevenGJohnson/cubature-jl-ux9xb)


# The Cubature module for Julia

This module provides one- and multi-dimensional adaptive integration
routines for the [Julia language](http://julialang.org/), including
support for vector-valued integrands and facilitation of parallel
evaluation of integrands, based on the [Cubature
Package](https://github.com/stevengj/cubature) by Steven G. Johnson.

See also the [HCubature package] for a pure-Julia implementation of
h-adaptive cubature using the same algorithm (which is therefore much
more flexible in the types that it can integrate).

## h-adaptive versus p-adaptive integration

Adaptive integration works by evaluating the integrand at more and
more points until the integrand converges to a specified tolerance
(with the error estimated by comparing integral estimates with
different numbers of points).  The Cubature module implements two
schemes for this adaptation: *h-adaptivity* (routines `hquadrature`,
`hcubature`, `hquadrature_v`, and `hcubature_v`) and *p-adaptivity*
(routines `pquadrature`, `pcubature`, `pquadrature_v`, and
`pcubature_v`).  The h- and p-adaptive routines accept the same
parameters, so you can use them interchangeably, but they have
very different convergence characteristics.

**h-adaptive** integration works by recursively subdividing the
integration domain into smaller and smaller regions, applying the same
fixed-order (fixed number of points) integration rule within each
sub-region and subdividing a region if its error estimate is too
large.  (Technically, we use a Gauss-Kronrod rule in 1d and a
Genz-Malik rule in higher dimensions.)  This is well-suited for
functions that have localized sharp features (peaks, kinks, etcetera)
in a portion of the domain, because it will adaptively add more points
in this region while using a coarser set of points elsewhere.  The
h-adaptive routines should be *your default choice* if you know very
little about the function you are integrating.

**p-adaptive** integration works by repeatedly doubling the number of
points in the same domain, fitting to higher and higher *degree*
polynomials (in a stable way) until convergence is achieved to the
specified tolerance.  (Technically, we use Clenshaw-Curtis
quadrature rules.)  This is best-suited for integrating *smooth* functions
(infinitely differentiable, ideally analytic) in *low* dimensions
(ideally 1 or 2), especially when high accuracy is required.

One technical difference that is sometimes important for functions
with singularities at the edges of the integration domain: our
h-adaptive algorithm only evaluates the integrand at the *interior* of
the domain (never at the edges), whereas our p-adaptive algorithm also
evaluates the integrand at the edges.

(The names "h-adaptive" and "p-adaptive" refer to the fact that the
size of the subdomains is often denoted *h* while the degree of the
polynomial fitting is often called *p*.)

## Usage

Before using any of the routines below (and after installing, see above),
you should include `using Cubature` in your code to import the functions
from the Cubature module.

### One-dimensional integrals of real-valued integrands

The simplest case is to integrate a single real-valued integrand `f(x)`
from `xmin` to `xmax`, in which case you can call (similar to
Julia's built-in `quadgk` routine):

    (val,err) = hquadrature(f::Function, xmin::Real, xmax::Real;
                            reltol=1e-8, abstol=0, maxevals=0)

for h-adaptive integration, or `pquadrature` (with the same arguments)
for p-adaptive integration.  The return value is a tuple of `val` (the
estimated integral) and `err` (the estimated absolute error in `val`,
usually a conservative upper bound).  The required arguments are:

* `f` is the integrand, a function `f(x::Float64)` that accepts a real
  argument (in the integration domain) and returns a real value.

* `xmin` and `xmax` are the boundaries of the integration domain.  (That is,
  `f` is integrated from `xmin` to `xmax`.)  They must be *finite*; to
  compute integrals over infinite or semi-infinite domains, you can use
  a [change of variables](https://github.com/stevengj/cubature/blob/master/README.md#infinite-intervals).

There are also the following optional keyword arguments:

* `reltol` is the required *relative* error tolerance: the adaptive
  integration will terminate when `err` &le; `reltol`*|`val`|; the
  default is `1e-8`.

* The optional argument `abstol` is a required *absolute* error
  tolerance: the adaptive integration will terminate when `err` &le;
  `abstol`.  More precisely, the integration will terminate when
  *either* the relative- or the absolute-error tolerances are met.
  `abstol` defaults to 0, which means that it is ignored, but it
  can be useful to specify an absoute error tolerance for integrands
  that may integrate to zero (or nearly zero) because of large
  cancellations, in which case the problem is ill-conditioned and a
  small relative error tolerance may be unachievable.

* The optional argument `maxevals` specifies a (rough) maximum number
  of function evaluations: the integration will be terminated (and
  the current estimates returned) if this number is exceeded.  The
  default `maxevals` is 0, in which case `maxevals` is ignored (no
  maximum).

Here is an example that integrates f(x) = x^3 from 0 to 1, printing
the x coordinates that are evaluated:

    hquadrature(x -> begin println(x); x^3; end, 0,1)

and returning `(0.25,2.7755575615628914e-15)`, which is the correct
answer 0.25.  If we instead integrate from -1 to 1, the function may
never exit: the exact integral is zero, and it is nearly impossible to
satisfy the default `reltol` bound in floating-point arithmetic.  In
that case, you have to specify an `abstol` as explained above:

    hquadrature(x -> begin println(x); x^3; end, -1,1, abstol=1e-8)

in which case it quickly returns.

### Multi-dimensional integrals of real-valued integrands

The next simplest case is to integrate a single real-valued integrand `f(x)`
over a [multidimensional box](http://en.wikipedia.org/wiki/Hyperrectangle),
with each coordinate `x[i]` integrated from `xmin[i]` to `xmax[i]`.

    (val,err) = hcubature(f::Function, xmin, xmax;
                          reltol=1e-8, abstol=0, maxevals=0)

for h-adaptive integration, or `pcubature` (with the same arguments)
for p-adaptive integration.  The return value is a tuple of `val` (the
estimated integral) and `err` (the estimated absolute error in `val`,
usually a conservative upper bound).  The arguments are:

* `f` is the integrand, a function `f(x::Vector{Float64})` that accepts
  a vector `x` (in the integration domain) and returns a real value.

* `xmin` and `xmax` are arrays or tuples (or any iterable container)
  specifying the boundaries `xmin[i]` and `xmax[i]` of the integration
  domain in each coordinate.  They must have `length(xmin) == length(xmax)`.
  (As above, the components must be finite, but you can treat infinite
  domains via a change of variables).

* The optional keyword arguments `reltol`, `abstol`, and `maxevals`
  specify termination criteria as for `hquadrature` above.

Here is the same 1d example as above, integrating f(x) = x^3 from 0 to 1
while the x coordinates that are evaluated:

    hcubature(x -> begin println(x[1]); x[1]^3; end, 0,1)

which again returns the correct integral 0.25. The only difference from
before is that the argument `x` of our integrand is now an array, so
we must use `x[1]` to access its value.  If we have multiple coordinates, we use `x[1]`, `x[2]`, etcetera, as in this example integrating f(x,y) = x^3 y in the unit box [0,1]x[0,1] (the exact integral is 0.125):

    hcubature(x -> begin println(x[1],",",x[2]); x[1]^3*x[2]; end, [0,0],[1,1])

### Integrals of vector-valued integrands

In many applications, one wishes to compute integrals of several
different integrands over the same domain.  Of course, you could simply
call `hquadrature` or `hcubature` multiple times, once for each integrand.
However, in cases where the integrands are closely related functions, it
is sometimes much more efficient to compute them *together* for a given
point `x` than computing them separately.  For example, if you have
a complex-valued integrand, you could compute two separate integrals
of the real and imaginary parts, but it is often more efficient and
convenient to compute the real and imaginary parts at the same time.

The Cubature module supports this situation by allowing you to
integrate a *vector-valued* integrand, computing `fdim` real integrals
at once for any given dimension `fdim` (the dimension of the
*integrand*, which is independent of the dimensionality of the
integration *domain*).   This is achieved by calling one of:

    (val,err) = hquadrature(fdim::Integer, f::Function, xmin, xmax;
                            reltol=1e-8, abstol=0, maxevals=0,
                            error_norm = Cubature.INDIVIDUAL)
    (val,err) = hcubature(fdim::Integer, f::Function, xmin, xmax;
                          reltol=1e-8, abstol=0, maxevals=0,
                          error_norm = Cubature.INDIVIDUAL)

for h-adaptive integration, or `pquadrature`/`pcubature` (with the
same arguments) for p-adaptive integration.  The return value is a
tuple of two vectors of length `fdim`: `val` (the estimated integrals
`val[i]`) and `err` (the estimated absolute errors `err[i]` in
`val[i]`).  The arguments are:

* `fdim` the dimension (number of components) of the integrand,
  i.e. the number of real-valued integrals to perform simultaneously

* `f`, the integrand.  This is a function `f(x, v)` of two arguments:
  the point `x` in the integration domain (a `Float64` for
  `hquadrature` and a `Vector{Float64}` for `hcubature`), and the
  vector `v::Vector{Float64}` of length `fdim` which is used to output
  the integrand values.  That is, the function `f` should set `v[i]`
  to the value of the `i`-th integrand upon return.  (The return value
  of `f` is ignored.)  **Note**: the contents of `v` must be
  overwritten *in-place* by `f`.  If you are not setting `v[i]` individually,
  you should do `v[:] = ...` and *not* `v = ...`.

* `xmin` and `xmax` specify the boundaries of the integration domain,
  as for `hquadrature` and `hcubature` of scalar integrands above.

* The optional keyword arguments `reltol`, `abstol`, and `maxevals`
  specify termination criteria as for `hquadrature` above.

* The optional keyword argument `error_norm` specifies how the
  convergence criteria for the different integrands are combined.
  That is, given a vector `val` of integral estimates and a vector
  `err` of error estimates, how do we decide whether to stop?
  `error_norm` should be one of the following constants:

    * `Cubature.INDIVIDUAL`, the default.  This terminates the integration
      when all of the integrals, taken individually, converge.  That is,
      it checks `err[i]` &le; `reltol`*|`val[i]`| or `err[i]` &le;
      `abstol`, and only stops when one of these is true for *all* `i`.

    * `Cubature.PAIRED`.  This is like `Cubature.INDIVIDUAL`, but
      applies the convergence criteria to *consecutive pairs* of
      integrands, as if these integrands were real and imaginary parts
      of complex numbers.  (This is mainly useful for integrating
      complex functions in cases where you only care about error in
      the complex plane as opposed to error in the real and imaginary
      parts taken individually.)

    * `Cubature.L1`, `Cubature.L2`, or `Cubature.LINF`.  These
      terminate the integration when |`err`| &le;
      `reltol`*|`val`| or |`err`| &le; `abstol`, where |...|
      denotes a *norm* applied to the whole vector of errors or
      integrals.  In particular, the L1 norm (sum of absolute values),
      the L2 norm (the root-mean-square value), or the L-infinity norm
      (the maximum absolute value), respectively.  These are useful if
      you only care about the error in the vector of integrals taken
      as a whole in some norm, rather than the relative error in the
      components taken individually (which could be large if some of
      the components integrate almost to zero).   We provide three
      different norms for completeness, but probably the choice of
      norm doesn't matter too much; pick `Cubature.L1` if you aren't sure.

Here is an example, similar to above, which integrates a vector of
three integrands (x, x^2, x^3) from 0 to 1:

    hquadrature(3, (x,v) -> v[:] = x.^[1:3], 0,1)

returning `([0.5, 0.333333, 0.25],[5.55112e-15, 3.70074e-15,
2.77556e-15])`, which are of course the correct integrals.

### Parallelizing the integrand evaluation

These numerical integration algorithms actually call your integrand function
for batches of points at a time, not just point-by-point.   It is useful
to expose this information for parellelization: your code may be able
to evaluate the integrand in parallel for multiple points.

This is provided by a "vectorized" interface to the Cubature module:
functions `hquadrature_v`, `pquadrature_v`, `hcubature_v`, and
`pcubature_v`, which have *exactly* the same arguments as the
functions described in the previous sections, *except* that the
integrand function `f` must accept different arguments.

In particular, for the `_v` integration routines, the integrand must
be a function `f(x,v)` where `x` is an array of `n` points to evaluate
and `v` is an array in which to store the values of the integrands at
those points.  `n` is determined at runtime and varies between calls
to `f`. The shape of the arrays depends upon which routine is called:

* For `hquadrature_v` and `pquadrature_v` with real-valued integrands (no
  `fdim` argument), `x` and `v` are both 1d `Float64` arrays of length `n` of
  the points (input) and values (output), respectively.

* For `hcubature_v` and `pcubature_v` with real-valued integrands (no
  `fdim` argument) in `d` integration dimensions, `x` is a 2d `Float64`
  array of size `d`&times;`n` holding the points `x[:,i]` at which
  to evaluate the integrand, and `v` is a 1d `Float64` array of  length
  `n` in which to store the resulting integrand values.

* For `hquadrature_v` and `pquadrature_v` with vector-valued integrands (an
  `fdim` argument), `x` is a 1d `Float64` array of length `n` of points
  at which to evaluate the integrands, and `v` is a 2d `Float64` array
  of size `fdim`&times;`n` in which to store the values `v[:,i]` at these
  points.

* For `hcubature_v` and `pcubature_v` with vector-valued integrands
  (an `fdim` argument) in `d` integration dimensions, `x` is a 2d
  `Float64` array of length `d`&times;`n` of points `x[:,i]` at which
  to evaluate the integrands, and `v` is a 2d `Float64` array of size
  `fdim`&times;`n` in which to store the values `v[:,i]` at these
  points.

## Technical Algorithms and References

The h-adaptive integration routines are based on those described in:

* A. C. Genz and A. A. Malik, "An adaptive algorithm for numeric integration over an N-dimensional rectangular region,"  *J. Comput. Appl. Math.*, vol. 6 (no. 4), 295-302 (1980).
* J. Berntsen, T. O. Espelid, and A. Genz, "An adaptive algorithm for the approximate calculation of multiple integrals," *ACM Trans. Math. Soft.*, vol. 17 (no. 4), 437-451 (1991).

which we implemented in a C library, the [Cubature
Package](https://github.com/stevengj/cubature), that is called from Julia.

Note that we do ''not'' use any of the original DCUHRE code by Genz,
which is not under a free/open-source license.)  Our code is based in
part on code borrowed from the [HIntLib numeric-integration
library](http://mint.sbg.ac.at/HIntLib/) by Rudolf Schürer and from
code for Gauss-Kronrod quadrature (for 1d integrals) from the [GNU
Scientific Library](http://www.gnu.org/software/gsl/), both of which
are free software under the GNU GPL.  (Another free-software
multi-dimensional integration library, unrelated to our code here but
also implementing the Genz-Malik algorithm among other techniques, is
[Cuba](http://www.feynarts.de/cuba/).)

The `hcubature_v` technique is adapted from adapted from I. Gladwell,
"Vectorization of one dimensional quadrature codes," pp. 230--238 in
*Numerical Integration. Recent Developments, Software and
Applications*, G. Fairweather and P. M. Keast, eds., NATO ASI Series
C203, Dordrecht (1987), as described in J. M. Bull and T. L. Freeman,
"Parallel Globally Adaptive Algorithms for Multi-dimensional
Integration,"
<http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638>
(1994).

The p-adaptive integration algorithm is simply a tensor product of
nested [Clenshaw-Curtis
quadrature](http://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature)
rules for power-of-two sizes, using a pre-computed table of points and
weights up to order 2^20.

## Author

This module was written by [Steven G. Johnson](http://math.mit.edu/~stevenj/).
