using MDBM
using Test

function tests()
    @testset "Testing the Multi-Dimensional Bisection Method (MDBM) object" begin

        @test begin
            foo(x,y)=x^2.0+y^2.0-4.0^2.0
            c(x,y)=x-y

            ax1=Axis([-5,-2.5,0,2.5,5],"x")
            ax2=Axis(-5:2:5.0,"b")

            mymdbm=MDBM_Problem(foo,[ax1,ax2],constraint=c)
            true
        end

        @test begin
            mymdbm=MDBM_Problem((x,y)->x^2.0+y^2.0-4.0^2.0,[-5:5,-5:5])
            true
        end
    end

    @testset "Testing the MDBM solve!" begin
        @test begin
            mymdbm=MDBM_Problem((x,y)->x^2.0+y^2.0-4.0^2.0,[-5:5,-5:5])
            iteration=2 #number of refinements (resolution doubling)
            solve!(mymdbm,iteration)
            true
        end
    end
end

@testset "vectorized evaluation (one call per stage)" begin
    f(x, y) = x^2 + y^2 - 4.0^2
    calls = Ref(0)
    npts = Ref(0)
    fv(pts) = (calls[] += 1; npts[] += length(pts); [f(p...) for p in pts])
    m_scalar = MDBM_Problem(f, [-5.0:5.0, -5.0:5.0])
    m_vector = MDBM_Problem(f, [-5.0:5.0, -5.0:5.0]; vectorized = fv)
    solve!(m_scalar, 3; verbosity = 0)
    solve!(m_vector, 3; verbosity = 0)
    @test getinterpolatedsolution(m_vector) == getinterpolatedsolution(m_scalar)
    @test length(m_vector.fc.fvalarg.keys) == length(m_scalar.fc.fvalarg.keys)
    @test npts[] == length(m_vector.fc.fvalarg.keys)   # every point went through fv
    @test calls[] < 20                                  # a handful of batches, not one call per point
end
