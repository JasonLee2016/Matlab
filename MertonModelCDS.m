clear all
global B T se r q ps d1 d2 pvf y
%Inputs
B = 900000; %face value of zero coupon credit risky debt
T = 1; %time to maturity of credit risky debt
r = 0.04 %continuous compounded risk free rate
q = 20000; %number of shares outstanding
ps = 45; %market price per share
se = 0.55; %annual volatility of log price change equity
pvf = exp(-r*T); %present value factor, G0
xstart = [ps*q 0.25]; %initial values for optimization
% x(1) = V0 % market value real assets
% x(2) = sv % annual vol. log change mkt value real assets

%type merton_74Two equations in two unknowns
%type merton_74c
% options=optimset('Display','final','MaxFunEvals', 1000,...
% 'TolX',1e-6,'TolFun',1e-6);
% check both equations are individually zedro
options = optimset('display','final');
lb = [1e-9 1e-2];   % Lower bound
ub = [1e+7 0.5];  % Upper bound
opts=optimoptions(@fmincon,'Algorithm','interior-point');
problem=createOptimProblem('fmincon','objective',@Merton_74c,...
    'x0',xstart,'lb',lb,'ub',ub,'nonlcon',@Merton_74_constraint,...
    'options',opts);

% Construct a GlobalSearch object
gs = GlobalSearch;

% Run GlobalSearch
tic;
[xgs,~,~,~,solsgs] = run(gs,problem);
toc
xgs

% Calculate recovery rate, yield to maturity, risk neutral 
% probability, and CDS premium.
V0=xgs(1);
sv=xgs(2);
B0=B*exp(-r*T)*normcdf(d2,0,1)+V0*(1-normcdf(d1,0,1));
G0=B*exp(-r*T);
p=G0-B0;
rr=exp(r*T)*(V0/B)*(normcdf(-d1,0,1)/normcdf(-d2,0,1));
y=-log(B0/B)/T;
pai= (B0/G0-1)/(rr-1);
CDS=pai*(1-rr);


