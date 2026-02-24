using MDBM


using BenchmarkTools
#-----------------------------

#This type of type definitsion is an overkill here!
#function foo_par3_codim1(x::Float64, y::Float64, z::Float64)::Float64
function foo_par3_codim1(x, y, z)
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0,x+y
end

# constraint - calculate only the points where the constraint is satisfied (e.g.: on the positiev side)
#function c(x::Float64, y::Float64, z::Float64)::Float64
function c(x, y, z)
    x^2.0 + y^2.0 - 0.5^2.0
end


mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0])#, constraint=c)
@time solve!(mymdbm, 5, doThreadprecomp=false, verbosity=0);

@benchmark begin
    
mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:1.0, -1.0:3.0, -1.0:3.0])#, constraint=c)
solve!(mymdbm, 5, doThreadprecomp=false, verbosity=0);

end