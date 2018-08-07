function [answer]=richardsondifferentiationarrays(fu,x,kmax,h)
h=h(:);
x=x(:);
N=zeros(kmax);
N(1,1)=(feval(fu,x+h)-feval(fu,x-h))/(2*max(h));
for k=2:1:kmax
    h=h/2;
    N(k,1)=(feval(fu,x+h)-feval(fu,x-h))/(2*max(h));
    for j=2:1:k
        N(k,j)=(4^(j-1)*N(k,j-1)-N(k-1,j-1))/(4^(j-1)-1);
    end
    check=abs((N(k,k)-N(k-1,k-1))/N(k,k));
    if check  <  0.0000001
        break
    end
end
answer=N(k,k);
N(k+1:end,:)=[];
N(:,k+1:end)=[];