%%Explicit Method

clear all; clc;

%Import Data
K=[95 100 105]; %Strike price
S0=100;%Initial Underlying price
r=0.1; %Risk-free rate
s=0.5; %Annualized volatility
T=0.25; %time to expiration in years
alpha=0.5; %CEV positive constant
beta=alpha*2;

%CEV Explicit
M=50;  %M+1 price rows (ii)in grid
 
Smax=150;
DS=Smax/M;
 
N = 200; %N+1 time cols (jj)in grid
Dt=T/N;
SN=linspace(0,Smax,M+1)';
ii=0:M;
jj=0:N;
 
a=0.5*Dt*(s^2*ii./DS - r*ii);
b=1-Dt*(s^2*ii./DS + r);
c=0.5*Dt*(s^2*ii./DS +r*ii);

%Boundary conditions for put
for kk=1:length(K)
    gridx=zeros(M+1,N+1);
    gridx(:,N+1)=max(K(kk)-SN,0);   
    gridx(1,:)=K(kk)*exp(-r*Dt*(N-jj));
    gridx(M+1,:)=0;
    for j = N:-1:1
        for i = 2:M
            gridx(i,j)=a(i)*gridx(i-1,j+1)+b(i)*gridx(i,j+1)+c(i)*gridx(i+1,j+1);
        end
    end
    pCE(kk)=interp1(SN,gridx(:,1),S0)  %explicit solution
    
    if T/N*M^2*s^2 > 1 %Numeric Stability Test
        disp('Explicit method:There is no numerical stability. Consider changing the values of M and N');
    end
end










