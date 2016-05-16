%%Implied Volatility for SPX OTM options

%Variables
S0 = 1972;
dy = 0.01 * 1.971193;
dt = datenum('08/31/2015');
 
%SPX Options Data from WRDS
[opt labelopt] = xlsread('SPX.xlsx','Sheet1');
opt(:,1) = x2mdate(opt(:,1)); %Expiration date
opt(:,2) = x2mdate(opt(:,2)); %Last Trade Date
opt=opt(opt(:,1)>dt,:); %Live options
size(opt)
T = opt(:,14); %Time to expiration in years
max(T), min(T)

%Continuous compounded risk free rate from WRDS
%Interpolate risk free rate to time to expiration
[rf labelrf] = xlsread('TCurve.xlsx','Sheet1');
rf(:,1) = x2mdate(rf(:,1));
%rf=rf(rf(:,1)==dt,:); %Only 08/31/2015
rfi = interp1(rf(:,2), rf(:,3), T, 'spline')
 
%After performing ITM calls/puts validation we read in the data again
[a labela] = xlsread('SPX.xlsx','Sheet2');
 
% Calculate implied volatility
vol = blsimpv(S0, a(:,4), a(:,15), a(:,14), 0.5*(a(:,5)+a(:,6)),[],dy, [], a(:,3)==1);
sum(isfinite(vol))
sum((isnan(vol)))
check = [a(:,1) a(:,14) a(:,3) a(:,4) (0.5*(a(:,5)+a(:,6))) vol]
check_new = [a(:,1) a(:,14) a(:,15) a(:,3) a(:,4) (0.5*(a(:,5)+a(:,6))) vol]
clean = find(isfinite(check(:,6)))
clean_new = find(isfinite(check_new(:,6)))

%OTM options for further analysis
outc = find(and(check(:,3) == 1, check(:,4) > S0));
length(outc)
outp = find(and(check(:,3) == 0, check(:,4) < S0));
length(outp)
outa = cat(1, outc, outp);
out = check(outa, :);
 
