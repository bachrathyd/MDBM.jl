path(path,'/home/bachrathy/Documents/Git_repo/Multi-Dimensional-Bisection-Method/matlab/')
%% Multi-Dimensional Bisection Method
% -4 basic example -
% parameter dimension : 3
%co-dimension (number of equations): 1

% definition of the parameter space
%the limits and the initial mesh
ax=[];
ax(1).val=[-5,0,3,5];
ax(2).val=-7:1:7;
ax(3).val=-3:3

% number of iteration
Niteration=5;%take care, the large values can easily lead to memory problem
mdbmoption=mdbmset('isconstrained',1,'interporder',1,'connections',1,'checkneighbour',0)
%% function for which the roots are detected
[c1] = clock
tic
bound_fuction_name='fval_basic_4_pardim3_codim1';
mdbm_sol=mdbm(ax,bound_fuction_name,Niteration,mdbmoption);
toc
clock-c1
figure(4),clf
ghendle=plot_mdbm(mdbm_sol,'k',[],0);
set(ghendle,'Marker','o','MarkerSize',10)
set(ghendle,'markerfacecolor','b')
view(2)
% set(ghendle,'Marker','*')
% set(ghendle,'LineWidth',3)
% figure(4),clf
% % hold on
% plot_mdbm(mdbm_sol);
% shading interp
% lighting flat
% light('Position',[1 -1 -15],'Style','infinite');
% light('Position',[1 -1 25],'Style','infinite');
% light('Position',[10 10 10],'Style','infinite');
% light('Position',[-10 -10 -10],'Style','infinite');
% view(3)
axis tight