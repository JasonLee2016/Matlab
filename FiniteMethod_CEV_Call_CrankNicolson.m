%% CRANK-NICOLSON METHOD

clear all; clc;

%Import Data
K=[95 100 105]; %Strike price
S0=100;%Initial Underlying price
r=0.1; %Risk-free rate
s=0.5; %Annualized volatility
T=0.25; %time to expiration in years
alpha=0.5; %CEV positive constant
beta=alpha*2;%CEV positive constan

M=1000; N=200; %Price steps and time steps
Smax=150;
DS=Smax/M;
Dt=T/N;

%Set up the the Coefficients Matrix
SN=linspace(0,Smax,M+1)';
ii=0:M;
veti=SN/DS;
jj=0:N;
alpha = 0.25*Dt*(s^2*veti/DS-r*veti);
beta = -0.5*Dt*(s^2*veti/DS+r);
gamma = 0.25*Dt*(s^2*veti/DS+r*veti);
M1 = -diag(alpha(3:M),-1)+diag(1-beta(2:M))-diag(gamma(2:M-1),1);
[L,U] = lu(M1);
M2 = diag(alpha(3:M),-1)+diag(1+beta(2:M))+diag(gamma(2:M-1),1);

%Boundary conditions for call
for kk=1:length(K)
    gridx=zeros(M+1,N+1);
    gridx(:,N+1)=max(SN-K(kk),0);   
    gridx(1,:)=0;
    gridx(M+1,:)=Smax-K(kk)*exp(-r*Dt*(N-jj));
    for j = N:-1:1
        gridx(2:M,j)=U\(L\(M2*gridx(2:M,j+1)));
    
    end
    cCCK(kk)=interp1(SN,gridx(:,1),S0);
end

