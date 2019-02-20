path(path,'/home/bachrathy/Documents/Git_repo/Multi-Dimensional-Bisection-Method/matlab/')
%% Multi-Dimensional Bisection Method

ax=[];
ax(1).val=-2:2;
ax(2).val=-2:2;
ax(3).val=-2:2;
Niteration=4;%take care, the large values can easily lead to memory problem
mdbmoption=mdbmset('isconstrained',0,'interporder',1,'connections',0,'checkneighbour',1);
%% function for which the roots are detected
[c1] = clock
tic
bound_fuction_name='fval_nsphere';
mdbm_sol=mdbm(ax,bound_fuction_name,Niteration,mdbmoption);
toc
clock-c1
figure(4),clf
ghendle=plot_mdbm(mdbm_sol,'k',[],0);
set(ghendle,'MarkerSize',0.5)