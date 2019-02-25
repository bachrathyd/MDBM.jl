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
