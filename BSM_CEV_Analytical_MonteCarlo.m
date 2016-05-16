clear all; clc;

%Import Data
K=[95 100 105]; %Strike price
S0=100;%Initial Underlying price
r=0.1; %Risk-free rate
s=0.5; %Annualized volatility
T=0.25; %time to expiration in years
alpha=0.5; %CEV positive constant
beta=alpha*2;%CEV positive constan
q=0; %Continuous underlying asset yield
GBM=zeros(5,6); %GBM Call&Put option prices
CEV=zeros(5,6); %CEV Call&Put option prices

%Put & Call: GBM Analytical Estimation
[c, p] = blsprice(S0,K,r,T,s,q)
GBM(1,1:3)=c;
GBM(1,4:6)=p;

%Put & Call: CEV Analytical Estimation
v=(s^2)/(2*(r-q)*(alpha-1))*(exp(2*(r-q)*(alpha-1)*T)-1);%Ch-squared distribution parameter.
aA=(K*exp(-(r-q)*T).^(2*(1-alpha)))/((1-alpha)^2*v);
bA=1/(1-alpha);
cA=S0^(2*(1-alpha))/((1-alpha)^2*v);

if (alpha<1 && alpha>0)
cCA=S0*exp(-q*T)*(1-ncx2cdf(aA,bA+2,cA))-K.*exp(-r*T).*ncx2cdf(cA,bA,aA);
pCA=K.*exp(-r*T).*(1-ncx2cdf(cA,bA,aA))-S0*exp(-q*T).*ncx2cdf(aA,bA+2,cA);
CEV(1,1:3)=cCA;
CEV(1,4:6)=pCA;

elseif alpha>1
cCA=S0*exp(-q*T)*(1-ncx2cdf(cA,-bA,aA))-K.*exp(-r*T).*ncx2cdf(aA,2-bA,cA);
pCA=K.*exp(-r*T).*(1-ncx2cdf(aA,2-bA,cA))-S0*exp(-q*T).*ncx2cdf(cA,-bA,aA);
CEV(1,1:3)=cCA;
CEV(1,4:6)=pCA;

else
   disp('Error! Wrong value! Alpha has to be a positive constant');
end

%Put & Call: GBM Monte Carlo Simulation
mu=r-q; %expected return
sig=s; %expected volatility
nu=mu-sig^2/2 %drift
steps=1000; %number of time steps
dt=T/steps; %size of time steps
nsims=100000; % Number of simulated paths

C=cumprod(exp(nu*dt+sig*sqrt(dt)*randn(steps,nsims))); %Simulation Matrix
S=S0*C(1000,:);
CallPayOffT95=max(S-K(1),0);
CallPayOffT100=max(S-K(2),0);
CallPayOffT105=max(S-K(3),0);
PutPayOffT95=max(K(1)-S,0);
PutPayOffT100=max(K(2)-S,0);
PutPayOffT105=max(K(3)-S,0);
cGM95=mean(CallPayOffT95)*exp(-r*T);
cGM100=mean(CallPayOffT100)*exp(-r*T);
cGM105=mean(CallPayOffT105)*exp(-r*T);
pGM95=mean(PutPayOffT95)*exp(-r*T);
pGM100=mean(PutPayOffT100)*exp(-r*T);
pGM105=mean(PutPayOffT105)*exp(-r*T);
GBM(2,1:3)=[cGM95 cGM100 cGM105]
GBM(2,4:6)=[pGM95 pGM100 pGM105]

%Put & Call: CEV Monte Carlo Simulation
SC=ones(steps,nsims);
Z=randn(steps,nsims);
SC(1,:)=S0*exp((mu-0.5*sig^2*(S0^(2*(alpha-1))))*dt+...
            sig*sqrt(dt)*(S0^(alpha-1)).*Z(1,nsims));

for i=1:nsims
    for j=2:steps
        SC(j,i)=SC(j-1,i)*exp((mu-0.5*sig^2*(SC(j-1,i)^(2*(alpha-1))))*dt+...
            sig*sqrt(dt)*(SC(j-1,i)^(alpha-1))*Z(j,i));
    end
end
STC=SC(1000,:);
CallPayOffT95C=max(STC-K(1),0);
CallPayOffT100C=max(STC-K(2),0);
CallPayOffT105C=max(STC-K(3),0);
PutPayOffT95C=max(K(1)-STC,0);
PutPayOffT100C=max(K(2)-STC,0);
PutPayOffT105C=max(K(3)-STC,0);
cCV95=mean(CallPayOffT95C)*exp(-r*T);
cCV100=mean(CallPayOffT100C)*exp(-r*T);
cCV105=mean(CallPayOffT105C)*exp(-r*T);
pCV95=mean(PutPayOffT95C)*exp(-r*T);
pCV100=mean(PutPayOffT100C)*exp(-r*T);
pCV105=mean(PutPayOffT105C)*exp(-r*T);
CEV(2,1:3)=[cCV95 cCV100 cCV105]
CEV(2,4:6)=[pCV95 pCV100 pCV105]






