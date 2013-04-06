using Base.Test
using Cubature

# The downside of having multiple interfaces to simplify 1d, scalar
# integrands, etcetera, is that we have lot of interfaces to test:

@test_approx_eq_eps hquadrature(cos, 0,1, 0, 1e-8)[1] sin(1) 1e-8
@test_approx_eq_eps pquadrature(cos, 0,1, 0, 1e-8)[1] sin(1) 1e-8

quad_tst1(x) = prod(cos(x))    
for dim = 0:3
    xmin = zeros(dim)
    xmax = Float64[1:dim]
    ans = prod(sin(xmax))
    @test_approx_eq_eps hcubature(quad_tst1, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    @test_approx_eq_eps pcubature(quad_tst1, xmin,xmax, 0, 1e-8)[1] ans 1e-8
end

quad_tst2(x, v) = for j = 1:length(v)
    v[j] = prod(cos(x)) * j
end
for fdim = 1:3
for dim = 0:3
    xmin = zeros(dim)
    xmax = Float64[1:dim]
    ans = prod(sin(xmax)) * [1:fdim]
    if dim == 1
        @test_approx_eq_eps hquadrature(fdim, quad_tst2, xmin,xmax, 0, 1e-8)[1] ans 1e-8
        @test_approx_eq_eps pquadrature(fdim, quad_tst2, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    end
    @test_approx_eq_eps hcubature(fdim, quad_tst2, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    @test_approx_eq_eps pcubature(fdim, quad_tst2, xmin,xmax, 0, 1e-8)[1] ans 1e-8
end
end

quad_tst0v(x,v) = for i = 1:length(x)
    v[i] = cos(x[i])
end
@test_approx_eq_eps hquadrature_v(quad_tst0v, 0,1, 0, 1e-8)[1] sin(1) 1e-8
@test_approx_eq_eps pquadrature_v(quad_tst0v, 0,1, 0, 1e-8)[1] sin(1) 1e-8

quad_tst1v(x,v) = for i = 1:length(v)
    v[i] = prod(cos(x[:,i]))
end
for dim = 0:3
    xmin = zeros(dim)
    xmax = Float64[1:dim]
    ans = prod(sin(xmax))
    @test_approx_eq_eps hcubature_v(quad_tst1v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    @test_approx_eq_eps pcubature_v(quad_tst1v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
end

quad_tst2v(x::Array{Float64,2}, v) = for i = 1:size(v,2)
    for j = 1:size(v,1)
        v[j,i] = prod(cos(x[:,i])) * j
    end
end
quad_tst2v(x::Array{Float64,1}, v) = for i = 1:size(v,2)
    for j = 1:size(v,1)
        v[j,i] = prod(cos(x[i])) * j
    end
end
for fdim = 1:3
for dim = 0:3
    xmin = zeros(dim)
    xmax = Float64[1:dim]
    ans = prod(sin(xmax)) * [1:fdim]
    if dim == 1
        @test_approx_eq_eps hquadrature_v(fdim, quad_tst2v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
        @test_approx_eq_eps pquadrature_v(fdim, quad_tst2v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    end
    @test_approx_eq_eps hcubature_v(fdim, quad_tst2v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
    @test_approx_eq_eps pcubature_v(fdim, quad_tst2v, xmin,xmax, 0, 1e-8)[1] ans 1e-8
end
end

# Check that nested calls work (although this is a suboptimal way to 
# compute multidimensional integrals):

@test_approx_eq_eps hquadrature(x -> hquadrature(cos, 0,2, 0, 1e-8)[1]*cos(x),
                                0,1, 0,1e-8)[1] sin(1)*sin(2) 1e-8
