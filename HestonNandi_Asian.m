clear all
global r N HN lret
%Monte Carlo Simulation
t1 = datenum('02/24/2016');
t2 = datenum('06/30/2016');
nT = t2 - t1;
N = 100000;
S0 = 0.949;
K = 0.949;
H0 = 6.0175e-04;
Z0 = -1.3163;
P = zeros(nT,N);
H = zeros(nT,N);
Z = randn(nT,N);

%Use callibrated Heston-Nandi Model for generating simulated price paths
H(1,:) = 1.2262e-05 +0.8257*6.0175e-04+ 4.6326e-05*(-1.3163-6.1407.*sqrt(6.0175e-04)).^2;
P(1,:) = 0.949*exp(0.00357/365 -0.5.*H(1,:)+sqrt(H(1,:)).*Z(1,:));
for j =2:nT
H(j,:) = 1.2262e-05 +0.8257.*H(j-1,:)+ 4.6326e-05.*(Z(j-1,:)-6.1407.*sqrt(H(j-1,:))).^2;
P(j,:) = P(j-1,:).*exp(0.00357/365 -0.5.*H(j,:)+sqrt(H(j,:)).*Z(j,:));
end

%Generate Expected Paths
Payoff = zeros(1,N);
for k = 1:N
Payoff(1,k) = mean(P(:,k));
end
PayoffC1 = Payoff-0.949;
PayoffP1 = 0.949-Payoff;
Payoff2C = zeros(1,N);
for k = 1:N
Payoff2C(1,k)= max(PayoffC1(1,k),0);
end
Payoff2P = zeros(1,N);
for k = 1:N
Payoff2P(1,k)= max(PayoffP1(1,k),0);
end
AsianOptionC=sum(Payoff2C)/N*exp(-0.00357/365*127);
AsianOptionP=sum(Payoff2P)/N*exp(-0.00357/365*127);