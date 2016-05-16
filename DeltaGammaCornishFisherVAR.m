clear all; clc;

%Import & Clean Data
df = xlsread('assign2.xlsx','Sheet1');
df = df(:,2);
df1 = df(1:end-4);
df4 = df(5:end);
logret4 = log(df4./df1);
lbqtest(logret4)
S = 1.1324; %Current underlying price
PutD = -0.29189; %Put Delta
PutG = 0.039464; %Put Gamma
PutP = 0.01375; %Put price
CallD = 0.374102; %Call Delta
CallG = 0.049241; %Call Gamma
CallP = 0.0191; %Call price

%iid Distribution VAR 99%
alpha = 0.01;
ust = std(logret4)*sqrt(4);
mu = mean(logret4);
DollarVAR = -(exp(mu+norminv(alpha)*ust)-1);
total = 10000000;
VAR = total*DollarVAR;

%GARCH(1,1)
Model = arima('ARLags',1,'Variance',garch(1,1))
EstMdl = estimate(Model,logret4);
cv1 = forecast(EstMdl,1);
cst = cv1*sqrt(4);
DollarVAR1 = -(exp(mu+norminv(alpha)*cst)-1);
VAR1 = total*DollarVAR1;


%Delta-Gamma Method
NumPut = 1/abs(PutD); %per 1 underlying contract
DollarVARDG = DollarVAR*abs(1+NumPut*PutD)...
              -1/2*PutG*DollarVAR^2;
VARDG = total*DollarVARDG;          
CostHedged = (S - PutP*NumPut)*total;


%Cornish-Fisher Method
Qp = 1.2322; %Number of Put
Qc = -1.7117; %Number of Call

NetD = (CallD*Qc+PutD*Qp);
NetG = (CallG*Qc+PutG*Qp);
PortD = NetD*total; %Costless Collar Delta
PortG = NetG*total; %Costless Collar Delta
Cost = (CallP*Qc + PutP*Qc)*total; %Total out of pocket cost


a = total*NetD*S;
b = 1/2*total*NetG*S^2;
sb = ust;

first = b*sb^2;
second = a^2*sb^2 + 3*b^2*sb^4;
third = 9*a^2*b*sb^4 + 15*b^3*sb^6;
fourth = 3*a^4*sb^4 + 90*a^2*b^2*sb^6 + 105*b^4*sb^8;

uv = first;
varv = second - uv.^2;
skew = (third - 3*second.*uv + 2*uv.^3)./((varv.^0.5).^3);
sigma = (varv.^0.5);

kurtosis = ((fourth - 4*uv*third + 6*second*uv^2 - 3*uv^4)/(varv.^0.5)^4) ...
           - 3;
z = norminv(0.01);

wq = z + 1/6*((z^2)-1)*skew + 1/24*((z^3)-3*z)*kurtosis...
     -1/36*(2*z^3 - 5*z)*skew.^2;
 
wq2 = z + 1/6*((z^2)-1)*skew; 
VARS = -(wq2.*varv.^0.5 + uv);
VARSK = -(wq.*varv.^0.5 + uv);
