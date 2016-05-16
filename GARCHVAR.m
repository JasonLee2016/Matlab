%Import and clean weekly gasoline and diesel prices from FRED
clear all
c = fred('https://research.stlouisfed.org/fred2/')
isconnection(c)

dsl=fetch(c, 'WDFUELLA'); %diesel price
idsl=find((isfinite(dsl.Data(:,2)))); %remove null
DSL=dsl.Data(idsl,:); %cleaned data
DSLw = fints(DSL(:,1), DSL(:,2)); %time series
DSLw = fts2mat(DSLw,1); %matrix

gas=fetch(c, 'DGASUSGULF'); %gasoline price
igas=find((isfinite(gas.Data(:,2)))); %remove null
GAS=gas.Data(igas,:); %cleaned data
GASw = fints(GAS(:,1), GAS(:,2)); %time series
GASw = fts2mat(GASw,1); %matrix

[cDSLGAS,iDSL,iGAS] = intersect(DSLw(:,1),GASw(:,1)); %common dates DSL GAS
data1w = [DSLw(iDSL,:) GASw(iGAS,2)]; %Combine filtered DSL data with GAS price

%Historical price pattern and risks.
hold on
plot(data1w(:,1),data1w(:,2)) %plot time and weekly diesel price
plot(data1w(:,1),data1w(:,3)) %plot time and weekly gas price
title('Weekly price of gasoline and diesel show highly corrleated market price risks.')
legend('WDFUELLA','DGASUSGULF','Location','northwest')
ylabel('$/gal')
datetick('x',12)
axis tight
hold off

%LEE Oil&Gas gross margin/Loss function
q1=25000000; %average monthly gallons of gasoline sold
q2=15000000; %average monthly gallons of diesel sold
m=0.02 %wholesale margin per gallon delivered
nL=length((data1w(3:4:end,2)-data1w(1:4:end-2,2)));

%L = -m*(q1+q2) + 0.5*q1*dS(:,2) + 0.5*q2*dS(:,3);%historic losses
%Assumed 0.5 inventory restocking at week1 beginning and
%week2 end with price set at week1 beginning
L = -m*(q1+q2)+ 0.5*q1*(data1w(3:4:end,2)-data1w(1:4:end-2,2))+...
0.5*q2*(data1w(3:4:end,3)-data1w(1:4:end-2,3));
sum(L>0) %data1w(:,2)=gasoline price, data1w(:,3)=diesel price
Lsort = sortrows(L); %Empirical Loss Distribution

%Historical Loss plot
subplot(2,2,1)
hold on
plot(data1w(1:4:end-2,1),L)
title('4-week Loss shows violations of iid normal.')
ylabel('$/gal')
datetick('x',12)
axis tight
hold off

%Serial dependency analysis
autocorr(L);
autocorr(abs(L));
autocorr(L.^2);
lbqtest(L)
N=length(L);
Lsort=sort(L);

%Empirical loss distribution risk measures
alpha=0.95;
VaR_H =Lsort(ceil(alpha*N))
ES_H = mean(Lsort(ceil(alpha*N)+1:end))
ecdf(L);
title('Empirical CDF realized losses')
format bank
%*ones(nL,1)

%Normality test
subplot(2,2,2)
hold on
histogram(L,20)
title('Historical Distribution 4-week Loss shows non-Gaussian properties.')
legend('skewness = -0.57, kurtosis = 6.75','Location','northwest')
ylabel('frequency')
xlabel('$')
axis tight
hold off

%Serial dependency test
subplot(2,2,3)
autocorr(L);
title('4-week Loss acf plot shows significant serial dependency.')
legend('lbq(indepence), jb(normality),arch(homoscedasticity)tests all failed.')


%Conitional variance model
model = garch('Offset',NaN,'GARCHLags',1,'ARCHLags',1)
fitL =estimate(model,L);
kappa=fitL.Constant;
G=cell2mat(fitL.GARCH(1));
A=cell2mat(fitL.ARCH(1));
LRvar = 250*kappa/(1-(G+A))
sqrt(LRvar)
V = infer(fitL,L); %conditional variances
res = (L-fitL.Offset)./sqrt(V); %standardized residuals

[lbqtest(res) jbtest(res) archtest(res)]

kurtosis(res);
skewness(res);
mean(res);
std(res);

subplot(2,2,4),autocorr(res)
title('Fitted GARCH(1,1) model residual acf plot does not show significant serial dependency.')
legend('lbq(indepence), jb(normality),arch(homoscedasticity)tests all passed.')

%Comparision between last observed conditional variance,
%first forecasted conditional variance, and long
%term variance (mean reverting to long term volatility)
Vplus1 = kappa + G*V(end)+A*(L(end)-fitL.Offset)^2;
[sqrt(V(end)) sqrt(kappa/(1-(G+A))) sqrt(Vplus1)]

%99 VaR assuming Gaussian innovations.
VaR99=norminv(0.99,0,1)*sqrt(Vplus1)+fitL.Offset
VaR95=norminv(0.95,0,1)*sqrt(Vplus1)+fitL.Offset
