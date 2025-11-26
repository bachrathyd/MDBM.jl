
#5+5
#using Revise
using MDBM
# using Plots - for testing the output

# ------------ 2D ------------
function foo_par2_codim1(x, y)
    x^2.0 + y^2.0 - 2.0^2.0
end
function coo_par2(x, y)
    x + y
end
# ------------ 3D ------------
function foo_par3_codim1(x, y, z)
    x^2.0 + y^2.0 + z^2.0 - 2.0^2.0
end
function foo_par3_codim2(x, y, z)
    return x^2.0 + y^2.0 + z^2.0 - 2.0^2.0, x - sin(z * 2)
end
function coo_par3(x, y, z)
    x + y + z
end
# ------------ 4D ------------
function foo_par4_codim2(x, y, z, w)
    return x^2.0 + y^2.0 + z^2.0 + w^2.0 - 2.0^2.0, x - sin(z * 2)
end
function foo_par4_codim3(x, y, z, w)
    return x^2.0 + y^2.0 + z^2.0 + w^2.0 - 2.0^2.0, x - sin(z * 2), y - cos(w * 2)
end
function coo_par4(x, y, z, w)
    x + y + z + w
end
# ------------ 2D ------------

mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xy_sol = getinterpolatedsolution(mymdbm)
# scatter(xy_sol...,)

mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0], constraint=coo_par2)
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par2_codim1, [-3.0:3.0, -3.0:3.0], constraint=coo_par2)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xy_sol = getinterpolatedsolution(mymdbm)
# scatter(xy_sol...,)
# ------------ 3D ------------
mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyz_sol = getinterpolatedsolution(mymdbm)
# scatter(xyz_sol...,)
mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par3_codim1, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyz_sol = getinterpolatedsolution(mymdbm)
# scatter(xyz_sol...,)


mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyz_sol = getinterpolatedsolution(mymdbm)
# scatter(xyz_sol...,)


mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par3_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par3)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyz_sol = getinterpolatedsolution(mymdbm)
# scatter(xyz_sol...,)

# ------------ 4D ------------

mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyzw_sol = getinterpolatedsolution(mymdbm)
# scatter(xyzw_sol[1:3]...,)

mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par4_codim2, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyzw_sol = getinterpolatedsolution(mymdbm)
# scatter(xyzw_sol[1:3]...,)


mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0])
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyzw_sol = getinterpolatedsolution(mymdbm)
# scatter(xyzw_sol[1:3]...,)

mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
solve!(mymdbm, 2, doThreadprecomp=false, verbosity=1)
mymdbm = MDBM_Problem(foo_par4_codim3, [-3.0:3.0, -3.0:3.0, -3.0:3.0, -3.0:3.0], constraint=coo_par4)
solve!(mymdbm, 2, doThreadprecomp=true, verbosity=1)
xyzw_sol = getinterpolatedsolution(mymdbm)
# scatter(xyzw_sol[1:3]...,)
