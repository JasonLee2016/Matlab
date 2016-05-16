%%Gram Charlier Model

clear all; clc;

%Import SPX options 
[M, labelM] = xlsread('assign1.xlsx', 'OTM');
global data;
global Pgk;
global GC_Pgk;
global x;
global z;

%Form a Data Array
S0=1972
dy=0.01971193 %Continuous Underlying Asset Yield
data(:,1) = M(:,4); %Call/Put 
data(:,2) = S0; %Underlying Price
data(:,3) = M(:,2); %Strike Price
data(:,4) = M(:,5); %Bid and Ask Mid-point 
data(:,5) = M(:,6); %Risk-free Rate
data(:,6) = M(:,1); %Time to Maturity
data(:,7) = dy; %dy
%data(:,8) = a(:,8); %highest bid
%data(:,9) = a(:,9); %lowest ask
 
[call, labelcall] = xlsread('assign2.xlsx','Call');
plot(call(:,1), call(:,2))
xlabel('Strike Price')
ylabel('Implied Volatility')
 
[put, labelput] = xlsread('assign2.xlsx','Put');
plot(put(:,1), put(:,2))
xlabel('Strike Price')
ylabel('Implied Volatility')
 
%Garman-Kolhagen Model Calibration
x0=0.2; %Initial volatiliy guess
options=optimset('Display', 'final', 'MaxFunEvals', 500, 'MaxIter', 1000, 'TolX', 1e-4, 'TolFun', 1e-4);
[x, GK_fval]= fminsearch(@impvol,x0,options)  %finds the volatility & min of sum of sq price errors;
sqrt(GK_fval);
data(:,11) = Pgk;

%Calculating Delta for Garman-Kohlhagen Model
data(:,10) = x*ones(size(a,1),1);
diff = data(:,5)- data(:,7);
GK_d1 = (log(data(:,2)./data(:,3)))+(diff+(((data(:,10).^2)./2).*data(:,6))./(data(:,10).*sqrt(data(:,6))));
 
for j=1:size(data,1)
    if data(j,1)== 1 %For call
     GK_delta(j) = exp(-data(j,6).*data(j,7)).*normcdf(GK_d1(j,1),0,1); %Call Delta = e^-rt*N(d1)
    else %For Put
     GK_delta(j) = exp(-data(j,6).*data(j,7)).*(normcdf(GK_d1(j,1),0,1)- 1); %Put Delta = e^-rt*(N(d1) - 1)
    end
end
GK_delta;
 
%Gram-Charlier Model Calibration
y=[0.2, 0.25, 3]; %Initial guess values for annualized imp vol, skewness, and kurtosis
GC_options=optimset('Display', 'final', 'MaxFunEvals', 500, 'MaxIter', 1000, 'TolX', 1e-4, 'TolFun', 1e-4);
[z,GC_fval] = fminsearch(@GC_impvol,y, GC_options)  %Finds the volatility & min of sum of sq price errors
sqrt(GC_fval)
data(:,12) = GC_Pgk

%Estimating the Delta of Gram Chalier Model
GCdata(:,1) = z(1,1)*ones(size(a,1),1); %Sigma
GCdata(:,2) = z(1,2)*ones(size(a,1),1); %Skewness
GCdata(:,3) = z(1,3)*ones(size(a,1),1); %Kurtosis
GCsigma = z(1,1)*sqrt(data(:,6));
GC_d1 = (log(data(:,2)./data(:,3))+(diff+(GCdata(:,3).^1)./2).*data(:,6))./(GCdata(:,1).*(data(:,6).^0.5));
GCC_part1 = normpdf(GC_d1,0,1);
GC_skewness = (GCdata(:,2).*normpdf(GC_d1,0,1).*(1-GC_d1.^2+3.*GC_d1.*GCsigma)-2.*(GCsigma.^2))./(factorial(3));
GC_kurtosis = (GCdata(:,3).*normpdf(GC_d1,0,1)).*3.*GC_d1.*(1+2.*(GCsigma.^2))+4*(GC_d1.^2).*(GCsigma)-(GC_d1.^3)-.4*GCsigma+3.*(GCsigma.^3)/(factorial(4));
GCC_delta = exp(-data(:,6).*data(:,7)).*(GCC_part1 - GC_skewness+GC_kurtosis) %Call Delta
GCP_delta = GCC_delta-exp(-data(:,6).*data(:,7)) %Put Delta
GC_delta  = zeros(size(data,1),1);
for j=1:size(data,1)
    if data(j,7)== 1 %calls
        GC_delta(j) = GCC_delta(j);
    else
        GC_delta(j) = GCP_delta(j);
    end
end
GCC_delta % Call Delta
GCP_delta % Put Delta
 
%Garman-Kohlhagen Model function
function L=impvol(x,data)
global  Pgk data x
if(x<=0)
   L=realmax;
else
    Pgk=zeros(size(data,1),1);
    for j=1:size(data,1)
       b = data(j,5)- data(j,7)
       PVS = exp(-data(j,6)*data(j,7));
       PVK = exp(-data(j,6)*data(j,5));
       d1=(log(data(j,2)./data(j,3)) +(b+0.5*x^2).*data(j,6))./(x*data(j,6).^0.5);
       d2= d1-(x*data(j,6).^0.5);
            if data(j,1)==1
                Pgk(j) = data(j,2)*PVS*normcdf(d1,0,1)-data(j,3)*PVK*normcdf(d2,0,1);
            else
                % Put Call Parity
                Pgk(j)= data(j,2)*PVS*normcdf(d1,0,1)-data(j,3)*PVK*normcdf(d2,0,1)+data(j,3)*PVK-data(j,2)*PVS;
            end
    end
    resid=(Pgk(:,1)-data(:,4))./data(:,4);
    L=resid'*resid;
end
 
%% Gram-Charlier Model function
function L=GC_impvol(y,data)
global  data GC_Pgk
b=y(1);
gamma1=y(2);
gamma2=y(3);
if(b<=0)
    L=realmax;
else
    GC_Pgk=zeros(size(data,1),1);

    for j=1:size(data,1)

        d1=(log(data(j,2)./data(j,3)) +(data(j,5)-data(j,7)+0.5*b^2).*data(j,6))./(b*data(j,6).^0.5);
        d2= d1-(b*data(j,6).^0.5);
        call =  data(j,2)*exp(-data(j,6)*data(j,7))*normcdf(d1,0,1)-data(j,3)*exp(-data(j,5)*data(j,6))*normcdf(d2,0,1);
        part1 = data(j,2)*exp(-data(j,6)*data(j,7))*normpdf(d1,0,1)*b*sqrt(data(j,6));
        skewness = (gamma1*(2*b*sqrt(data(j,6))-d1))/(factorial(3)*sqrt(data(j,6))); %%skewness
        kurtosis = (gamma2*(1-(d1^2)+3*d1*b*sqrt(data(j,6))-3*(b^2)*data(j,6)))/(factorial(4)*data(j,6)); %% kurtosis
        if data(j,1)==1 %call
            GC_Pgk(j) = call + part1*(skewness-kurtosis);
        else    %Put
            GC_Pgk(j) =  call + part1*(skewness-kurtosis)+ data(j,3)*exp(-data(j,5)*data(j,6))-data(j,2)*exp(-data(j,6)*data(j,7));
        end
    end
    resid=GC_Pgk(:,1)-data(:,4);
    L=resid'*resid; %sum of squared errors
    end
end