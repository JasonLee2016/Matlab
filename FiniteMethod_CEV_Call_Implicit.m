%%Implicit Method

clear all; clc;

%Data
K=[95 100 105]; %Strike price
S0=100;%Initial Underlying price
r=0.1; %Risk-free rate
s=0.5; %Annualized volatility
T=0.25; %time to expiration in years
alpha=0.5; %CEV positive constant
beta=alpha*2;

%CEV Implicit
M=50;  %M+1 price rows (ii)in grid
 
Smax=150;
DS=Smax/M;
 
N = 200; %N+1 time cols (jj)in grid
Dt=T/N;
SN=linspace(0,Smax,M+1)';
ii=0:M;
jj=0:N;
 
a=0.5*Dt*(r*ii-s^2*ii./DS);
b=1+Dt*(r+s^2*ii./DS);
c=-0.5*Dt*(s^2*ii./DS +r*ii);

coeff = diag(a(3:M),-1) + diag(b(2:M))+ diag(c(2:M-1),1);
[L,U]=lu(coeff);
 
aux=zeros(M-1,1);

%Boundary conditions for call
for kk=1:length(K)
    gridx=zeros(M+1,N+1);
    gridx(:,N+1)=max(SN-K(kk),0);   
    gridx(1,:)=0;
    gridx(M+1,:)=Smax-K(kk)*exp(-r*Dt*(N-jj));
    for j=N:-1:1
    aux(1)=-a(2)*gridx(1,j);
    gridx(2:M,j) = U\(L\(gridx(2:M,j+1)+aux));
    end
    cCI(kk)=interp1(SN,gridx(:,1),S0)  %explicit solution
    
    if (1+r*dt)^(-1)*(a(1)+b(1)+c(1))~= 1 % Numeric stability test
      disp('There is no numerical stability for the imiplicit. Consider changing the values of M and N.');
    end
end









