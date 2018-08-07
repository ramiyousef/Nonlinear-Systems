function [rcon]=RCON(matrix)
[m,n]=size(matrix);
if m~=n
    error('matrix is not square');
end
y=zeros(1,m);
% x=zeros(1,m);
zero=zeros(n,m);
tempmat=[matrix,zero]; %adding a column of zeros to the ...
%matrix to use it in the gsrp function to find the inverse
%%
for j=1:m
    y(j)=0;
    for i=1:m
        y(j)=y(j)+abs(matrix(i,j));
    end
end
maxabssum=max(y);
% inversematrix=inverseidentity(matrix);
[inversematrix,x]=gsrp(tempmat);
%%
for j=1:m
    x(j)=0;
    for i=1:m
        x(j)=x(j)+abs(inversematrix(i,j));
    end
end
%%
invmaxabssum=max(x);
rcon=1/(invmaxabssum*maxabssum);

