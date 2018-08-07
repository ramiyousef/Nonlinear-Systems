%General jacobian for any set of functions and variables
%This program uses Richardson extrapolation with 3 pt. stencil
function jac=jaco(F,xg)
xg=xg(:);
n=length(F);
m=length(xg);
jac=zeros(n,m);
h=0.00001;% step size chosen to be small for this project because the range of x1alpha and beta is small 0<x<1
hm=zeros(1,m);
hm=hm(:);
for j=1:n
    for i=1:m
        hm(i)=h;
        partial=richardsondifferentiationarrays(F{j},xg,15,hm);
        jac(j,i)=partial;
        hm(i)=0;
    end
end
%% RCON check 
reciprocalcondition=RCON(jac);
if reciprocalcondition<0.001
    error('Results may be inaccurate')
end
end