using MDBM
using PyPlot
pygui(true)

#Solution of an uncertain implicit equation

function foo(x,y)
    x^2.0+y^2.0-2.0^2.0+randn()
end

for k = 1:5
mymdbm=MDBM_Problem(foo,[-3.0:3.0,-3.0:3.0])
solve!(mymdbm,5)

x_sol,y_sol=getinterpolatedsolution(mymdbm)
fig = figure(1)#;clf()
scatter(x_sol,y_sol,s=4);

fig = figure(2)#;clf()
hist(x_sol,100)
end
