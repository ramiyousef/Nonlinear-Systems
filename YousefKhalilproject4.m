
function [xn]=YousefKhalilproject4(xg)

%% DECLARATIONS
n=length(xg);
ov=ones(1,n); %used for the check condition
check=1; % larger than zero to enter the while loop
xg=xg(:); % forcing xg to be a column vector
dyold=1e10;
i=0;%counter to break the while loop if too many iterations take place

%% LOOP
while check>1e-5
    i=i+1;
    l=binaryalphabeta;
    F=[l{1}(xg);l{2}(xg)];
    jacobian=jaco(binaryalphabeta,xg); % jaco has a built-in RCON checker
    augmentedm=zeros(2,3);
    %merging the jacobian with the F vector to be used in gauss-jordan
    augmentedm(:,1:2)=jacobian;
    augmentedm(:,3)=F;
    %% gauss-jordan with scaled row pivoting to find J^-1F
    [inverse,jinvtimesF]=gsrp(augmentedm); % inverse is not used but is written...
    %to extract the second output(jintimesF) only
    %% Newton Raphson
    xn=xg-jinvtimesF;
    %% using 1-norm for dy
    Fnew=[l{1}(xn);l{2}(xn)];
    dy=ov*abs(Fnew-F);
    if  i>100 || (dy>dyold && i>1) || isreal(xn)==0
        warning('Method failed after %1.0f iterations, solution may be incorrect!!',i-1);
        break
    end
    check=ov*abs(xg-xn);
    xg=xn;
end

%% PLOTTING
A=3;
B=2;
DELTAGRT=zeros(1,100);
x = 0.01:0.01:1;
for i=1:100 
DELTAGRT(i) = (x(i)*(A+2*(B-A)*x(i))*(1-x(i))^2+(1-x(i))*(B+2*(A-B)*(1-x(i)))*x(i)^2+x(i)*log(x(i))+(1-x(i))*log(1-x(i)));
end
% DELTAGRT(i)=[];
y = @(x) (x*(A+2*(B-A)*x)*(1-x)^2+(1-x)*(B+2*(A-B)*(1-x))*x^2+x*log(x)+(1-x)*log(1-x));

x1alpha = xn(1);
x1beta = xn(2);
y1alpha = feval(y,x1alpha);
y1beta = feval(y,x1beta);
slope = (y1beta-y1alpha)/(x1beta-x1alpha);
b = y1alpha - (slope*x1alpha);

tangent = (slope.*x+b);

fh=figure;
set(fh, 'color','w')
colordef white;
hold all
plot(x,DELTAGRT)
plot(x,tangent)
grid on

title('\Delta G^m^i^x/RT vs x_1','fontsize',14)
xlabel('x_1','fontsize',11)
ylabel('\Delta G^m^i^x/RT','fontsize',13,'fontangle','normal','fontweight','bold')
xlabel('x_1','fontsize',13,'fontangle','normal','fontweight','bold')
hlegend=legend('\Delta G^m^i^x/RT','Mutual Tangent');

set(hlegend,'fontsize',13,'box','off','units','normalized','position',[0.35 0.35 0.001 0.001],'fontangle','normal','orientation','horizontal')

