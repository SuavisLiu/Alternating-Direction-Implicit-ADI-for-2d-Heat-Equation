clc; clear all; 
[h,k,error] = ADI(79,79)

%%
clc; clear; close all
mxvals =  [39 79 159 319];
myvals =  [39 79 159 319];
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


%%
function [h,k,error] = ADI(mx, my)
    
    ax = 0;
    bx = 1;
    ay = 0;
    by = 1;
    tfinal = 0.1;
    hx = (bx-ax)/(mx+1);
    hy = (by-ay)/(my+1);
    x = linspace(ax,bx,mx+2);
    y = linspace(ay,by,my+2);
    [X,Y] = meshgrid(x,y);
    X = X';
    Y = Y';
    
    k = 0.0025 * (hx+hy)/2;    % time step
    nsteps = ceil(tfinal / k);     % number of time steps
    
    f = @(x,y,t) exp(-2*t*pi.^2).*sin(pi*x).*sin(pi*y);
    
    % initial condition at t = 0
    u0 = f(X,Y,0);
    
    % set up the matrices. 
    rx = (1/2) * k /(hx^2);
    ry = (1/2) * k /(hy^2);
    e = ones(mx,1);
    I = speye(my);
    S = spdiags([e e],[-1 1],my,my);
    T = spdiags([e -2*e e], [-1 0 1], mx, mx);
    Dx = rx*kron(I, T);  % This is actually (k/2) * Dx^2
    Dy = ry*((kron(I, -2*I) + kron(S, I)));  % This is actually (k/2) * Dy^2
    II = speye(my*mx);

    % initialize u and time 
    tn = 0;
    u = u0;
    
    % main time-stepping loop:
    
    for n = 1:nsteps
        gstar0 = zeros(mx,my);
        gn0 = zeros(mx,my);
        gnp0 = zeros(mx,my);
              
        tnp = tn + k;   % t_{n+1}
        
        fstar = f(X,Y,(tn + 0.5*k));
        
        % boundary conditions
%         gstar0(1,:) = fstar(1,2:(m+1));   % x = 0
%         gstar0(mx,:) = fstar(mx+2,2:(mx+1)) % x = 1
%         
%         %gstar0
%         
%         gn0(:,1) = u(2:(m+1),1);      % y = 0
%         gn0(:,m) = u(2:(m+1),m+2);    % y = 1
%         %gn0
%         
         unp = f(X,Y,tnp);
%         gnp0(:,1) = unp(2:(m+1),1);   % y = 0
%         gnp0(:,m) = unp(2:(m+1),m+2); % y = 1
%         %gnp0
%         
        uint = u(2:(mx+1),2:(my+1));  % interior points
        
        % reshape the interior pts and bcs to m*m vector
        uint = reshape(uint,mx*my,1);
        gstar0 = reshape(gstar0,mx*my,1);
        gn0 = reshape(gn0,mx*my,1);
        gnp0 = reshape(gnp0,mx*my,1);
        
        
        % solve for the first equation 
        rhs1 = (II + Dy)*uint + rx*gstar0 + ry*gn0;
        ustar = (II - Dx)\rhs1;
        
        % solve for the second equation, uint is U^{n+1} now
        rhs2 = (II + Dx)*ustar + rx*gstar0 + ry*gnp0;
        uint = (II - Dy)\rhs2;
        
        % add the boundary values, m*m vector unp
        uint = reshape(uint,mx,my);
        u = unp;
        u(2:(mx+1),2:(my+1)) = uint;

        
        tn = tnp;
        
    end      % end of the for loop
    
    h = (hx+hy)/2;
    
    ufinal = f(X,Y,tnp);
    error = max(max(abs(u-ufinal)));
    contour3(X,Y,u,50)
    title(sprintf('Numerical Soln by ADI with mx = %3d, my = %3d, tfinal = %2.2f',mx,my,tfinal))
    
    input('Hit <return> to continue  ');

    contour3(X,Y,ufinal,50)
    title('True solution to the heat equation')
end

%%
function error_table(h,E)
%
% Print out table of errors, ratios, and observed order of accuracy.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

ntest = length(h);
ratio = nan(size(h));   % initialize storage
order = nan(size(h));   % initialize storage

for j=2:ntest
   ratio(j) = E(j-1)/E(j);
   order(j) = log(abs(ratio(j))) / log(abs(h(j-1)/h(j)));
end


% print out table:

disp(' ')
disp('      h        error       ratio       observed order')
for j=1:ntest
   disp(sprintf(' %9.5f  %12.5e %9.5f %15.5f',h(j),E(j),ratio(j),order(j)));
end
disp(' ')
end
%%
function error_loglog(h,E)
%
% Produce log-log plot of E vs. h.
% Estimate order of accuracy by doing a linear least squares fit.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

h = h(:);            % make sure it's a column vector
E = E(:);            % make sure it's a column vector
ntest = length(h);
clf
loglog(h,E,'o-')
axis([.5*min(h) 1.5*max(h)  .5*min(E) 1.5*max(E)])
title('log-log plot of errors vs. h')

% Estimate order of accuracy from least squares fit:
Ap = ones(ntest,2);
Ap(:,2) = log(h);
bp = log(E);
Kp = Ap\bp;
K = Kp(1);
p = Kp(2);
disp(' ')
disp(sprintf('Least squares fit gives E(h) = %g * h^%g',exp(K),p))
disp(' ')

% add graph of this line to loglog plot:
hold on
err1 = exp(K)*h.^p;
loglog(h,err1,'r')
legend('errors', 'least squares fit','Location','SouthEast')
hold off
end