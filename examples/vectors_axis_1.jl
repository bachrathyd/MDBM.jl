using MDBM
using LinearAlgebra, Arpack
function fooeig(x)
       mymx=Matrix{Float64}(I,5,5);
       mymx[2,1]=x;
       mymx[1,2]=x+5;
       abs(eigs(mymx, nev = 2)[1][1])-3
end
eig_problem=MDBM_Problem(fooeig,[Axis(-10:10)])


println(" X solution  //   Function val")
for k=1:12
solve!(eig_problem,1)
Xsol=getinterpolatedsolution(eig_problem)
print(Xsol[1])
print("  //   ")
println(maximum(abs.(fooeig.(Xsol[1]))))
end
