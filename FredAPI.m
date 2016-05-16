%Data Import from FRED
c=fred('https://research.stlouisfed.org/fred2/');%Connect to API
series1='DCOILWTICO';
series2='DCOILWTICO';
fromdate='01/07/2007';
todate='12/07/2010';
d1=fetch(c,series1,fromdate,todate);
d2=fetch(c,series2,fromdate,todate);
d1.Data;
d2.data;
close(c); %Disconnet API

%%Loop Series ID equivalent to R


