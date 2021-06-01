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

