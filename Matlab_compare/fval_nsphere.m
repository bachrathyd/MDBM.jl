function H=fval_nsphere(ax)
H=zeros(1,size(ax,2));
for k=1:size(ax,2)
    H(k)=norm(ax(:,k)-0.4,2.7)-1.9;
end
end