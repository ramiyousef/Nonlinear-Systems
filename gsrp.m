%% gauss jordan with scaled row pivoting ++++ inverse of the A matrix
function [identity,x]=gsrp(augmented)
[m,n]=size(augmented);%dimensions of augmented matrix
s=zeros(m,1);
identity=eye(m);
%% finding the scale of every row
for i=1:m
    s(i)=max(abs(augmented(i,1:m)));
    if s(i)==0
        error('no unique solution exists')
    end
end

%% scaled row pivoting application
for k = 1:(m-1)
    r = abs(augmented(k,k)/s(k));
    kp = k;
    for i = k+1:m
        t = abs(augmented(i,k)/s(i));
        if t>r
            r=t; % setting new pivot
            kp=i;
        end
    end
    
%% interchange rows
    temp=augmented(kp,:);
    tempid=identity(kp,:);
    augmented(kp,:)=augmented(k,:);
    identity(kp,:)=identity(k,:);
    augmented(k,:)=temp;
    identity(k,:)=tempid;
    
%% Gaussian elimination
    for i=(k+1):m
        pivotelement=augmented(k,k);
        factor=augmented(i,k)/pivotelement;
        augmented(i,:)=augmented(i,:)-augmented(k,:)*factor;
        identity(i,:)=identity(i,:)-identity(k,:)*factor;
    end
end
%% Separating the augmented matrix from the column vector b
b=augmented(:,n);
augmented(:,n)=[];
%% performing Back Substitution
x=zeros(m,1);
x(m)=b(m)/augmented(m,m);
for j=m-1:-1:1
    x(j)=(b(j)-augmented(j,j+1:m)*x(j+1:m))/augmented(j,j);
end
%% The rest of operations on the augmented matrix and the identity...
% matrix until reaching the  reduced row echelon form and the inverse
%respectively
temp=augmented(i,i);
augmented(m,:)=augmented(m,:)/temp;
identity(m,:)=identity(m,:)/temp;
for i=m-1:-1:1
    temp=augmented(i,i);
    augmented(i,:)=augmented(i,:)/temp;
    identity(i,:)=identity(i,:)/temp;
    for j=i:-1:1
        factor=augmented(j,i+1);
        augmented(j,:)=augmented(j,:)-factor*augmented(i+1,:);
        identity(j,:)=identity(j,:)-factor*identity(i+1,:);
    end
end
end

