
clc; clear; close all
mxvals =  [39 79 159];
myvals =  [79 159 319];
ntest = length(mxvals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);   % to hold errors

for jtest=1:ntest
    mx = mxvals(jtest);
    my = myvals(jtest);
    [h,k,error] = ADI(mx,my);
    E(jtest) = error;
    hvals(jtest) = h;
    
end

error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit

