using MDBM
using GLMakie

GLMakie.closeall()
GLMakie.activate!(; title="2 parameters, codimension 1")

# ------------ 2D ------------
function foo_par2_codim1(x, y,z)
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0
end
function coo_par2(x, y,z)
    x ,y
end


    mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.1, -3.0:3.1, -3.0:3.1], constraint=coo_par2)
    solve!(mymdbm, 2, doThreadprecomp=false, verbosity=0)


 f = Figure(size=(1000, 600))
 ax1 = GLMakie.Axis3(f[1, 1])
 
 xyz_sol = getinterpolatedsolution(mymdbm)
 scatter!(ax1, xyz_sol..., markersize=15, color=:red, marker='x', strokewidth=3, label="solution")


