function DR_FP_1D
clear;
close all;
global Nx dx tau Nt N lambda M x Mas L k
%% Initialize parameter

Nt=1;
Nx=64;
N=(Nt+1)*(Nx+1);
L=6;
dx=2*L/Nx;
x=(-L:dx:L-dx)';
tau=1/20;
Nx=Nx-1;
f0 = 2*exp(-(x-2).^2)+0.5*exp(-(x+2.5).^2);
Mas=sum(f0)*dx;
M=exp(-(x).^2/2)/(2*pi)^0.5;
M=M/(sum(M)*dx)*Mas;
lambda=1;

%% Prepare A for constraints 
A_rho=eye(Nx+1);
A_m=zeros(Nx+1,Nx+1);
h=1/2/dx;
% Center difference 
 for i = 2:Nx
     A_m(i,i+1)=h;
     A_m(i,i-1)=-h;
 end
 A_m(end,end-1)=-h;
 A_m(1,2)=h;
 A_m(1,1)=h;
 A_m(end,end)=-h;
 
A=[A_rho A_m];




%% Initialize u0
f=f0;
num=1;
B=A*A';
maxiter=5e3;
b=f;
T=5;
t=0;
while t<T
    u_y=[f;zeros(Nx+1,1)];
    for k = 1:maxiter
        u_xp=proxF(u_y);
        u_yp=u_y+proxG(2*u_xp-u_y,A,B,b)-u_xp;
        if norm(u_y-u_yp)/norm(u_y)<5e-6
            u_y=u_yp;
            u_x=u_xp;
            break
        end 
        u_y=u_yp;
        u_x=u_xp;
        
    end
    t=t+tau;
    [f,m]=decomp(u_xp);
    plot(x,f,'r--',x,M,'g',x,f0,'b')
    drawnow
    num=num+1;
    b=f;
    filename=['FP_step_',num2str(num)];
    save(filename);
end
end

function [rho,m] = decomp(u)
global Nx
rho=u(1:Nx+1);
m = u(Nx+2:end);
end


function [rho,m] = proxF_biscalar(M,mstar,rhostar)
global lambda tau
tol=1e-12;
a=1e-15;
b=4;
res=1;
while abs(res)>tol
    c=(a+b)/2;
    ha= a+2*tau*lambda*log(a/M)+2*tau*lambda-lambda*mstar^2*(a+2*lambda).^(-2)-rhostar;
    hc= c+2*tau*lambda*log(c/M)+2*tau*lambda-lambda*mstar^2*(c+2*lambda).^(-2)-rhostar;
    %hb=(1+2*lambda*tau)*b+2*tau*lambda*log(b/M)-epsi*lambda*mstar^2/(b/lambda+2*epsi)^2-rhostar;
    sign(ha);
    if sign(hc) == sign(ha)
        a=c;
    else
        b=c;
    end
    res=hc;
end
rho=c;
m=(1+2*lambda/rho)^(-1)*mstar;
end


function p=proxF(u)
global M
[rhostar,mstar]=decomp(u);
n=length(mstar);
rho=0*rhostar;
m=0*mstar;
for i = 1:n
    [rho(i),m(i)] = proxF_biscalar(M(i),mstar(i),rhostar(i));
end
p=[rho;m];
end


function p=proxG(u, A,B, b)
mid=B\(A*u-b);
%error=sum(mid-A*u+b)
p=u-A'*mid;
end



