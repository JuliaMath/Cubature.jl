using Cubature, Test

# The downside of having multiple interfaces to simplify 1d, scalar
# integrands, etcetera, is that we have lot of interfaces to test:

@testset "Cubature" begin
    @test isapprox(hquadrature(cos, 0,1, abstol=1e-8)[1], sin(1), atol=1e-8)
    @test isapprox(pquadrature(cos, 0,1, abstol=1e-8)[1], sin(1), atol=1e-8)

    quad_tst1(x) = prod(cos.(x))
    for dim = 0:3
        xmin = zeros(dim)
        xmax = 1:dim
        ans = prod(sin.(xmax))
        @test isapprox(hcubature(quad_tst1, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        @test isapprox(pcubature(quad_tst1, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
    end

    quad_tst2(x, v) = for j = 1:length(v)
        v[j] = prod(cos.(x)) * j
    end
    for fdim = 1:3
    for dim = 0:3
        xmin = zeros(dim)
        xmax = 1:dim
        ans = prod(sin.(xmax)) * (1:fdim)
        if dim == 1
            @test isapprox(hquadrature(fdim, quad_tst2, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
            @test isapprox(pquadrature(fdim, quad_tst2, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        end
        @test isapprox(hcubature(fdim, quad_tst2, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        @test isapprox(pcubature(fdim, quad_tst2, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
    end
    end

    quad_tst0v(x,v) = for i = 1:length(x)
        v[i] = cos(x[i])
    end
    @test isapprox(hquadrature_v(quad_tst0v, 0,1, abstol=1e-8)[1], sin(1), atol=1e-8)
    @test isapprox(pquadrature_v(quad_tst0v, 0,1, abstol=1e-8)[1], sin(1), atol=1e-8)

    quad_tst1v(x,v) = for i = 1:length(v)
        v[i] = prod(cos.(x[:,i]))
    end
    for dim = 0:3
        xmin = zeros(dim)
        xmax = 1:dim
        ans = prod(sin.(xmax))
        @test isapprox(hcubature_v(quad_tst1v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        @test isapprox(pcubature_v(quad_tst1v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
    end

    quad_tst2v(x::Array{Float64,2}, v) = for i = 1:size(v,2)
        for j = 1:size(v,1)
            v[j,i] = prod(cos.(x[:,i])) * j
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
        xmax = 1:dim
        ans = prod(sin.(xmax)) * (1:fdim)
        if dim == 1
            @test isapprox(hquadrature_v(fdim, quad_tst2v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
            @test isapprox(pquadrature_v(fdim, quad_tst2v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        end
        @test isapprox(hcubature_v(fdim, quad_tst2v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
        @test isapprox(pcubature_v(fdim, quad_tst2v, xmin,xmax, abstol=1e-8)[1], ans, atol=1e-8)
    end
    end

    # Check that nested calls work (although this is a suboptimal way to
    # compute multidimensional integrals):

    @test isapprox(hquadrature(x -> hquadrature(cos, 0,2, abstol=1e-8)[1]*cos(x),
                            0,1, abstol=1e-8)[1], sin(1)*sin(2), atol=1e-8)
end