function alpha = est_alpha(r,g,buffer)
%
% E.A.K. Cohen and N.M. Adams, Dept of Mathematics, Imperial College London
%
% This piece of code estimates the algorithmic resolution limit by
% searching for a changepoint in the pair correlation function. The change
% point is interpreted to be the radius of correlation rho, and the
% algorithmic resolution limit is defined as alpha := rho/2.
% This is followed by an F-test at 5% level to ensure it is a genuine radius of
% correlation/algorithmic resolution limit.
%
% INPUTS:
%   r              vector of radial distances at which pair correlation is evaluated
%   g              vector of pair correlation values at radial distances in r. Must be
%                  same dimension as r 
%   buffer         number of points at either end of g which is are not
%                  considered as a possible value of rho, the radius of correlatiojn
% OUTPUTS:
%   alpha          algorithmic resolution limit


% create vector for test statistics
T = zeros(1,length(g)-2*buffer);

% For each value of g within the search range (controlled by buffer)
for ii = buffer:(length(g)-buffer)
    % Take all values to the right of position ii
    g_r = g(ii+1:end);
    % Compute the test statistic as the standard error of all points in g_r
    s_r = sqrt(var(g_r));
    T(ii-buffer+1) = (sqrt(length(g_r))./s_r);
end
% the estimate of rho is where the maximum value of T occurs
[~,I] = max(T);
true_index = I+buffer-1;
rho = r(true_index);
% the algorithmic resolution limit alpha is half rho
alpha = rho/2;

% Perform F-test
g_l = g(1:true_index-1);
s_l = sqrt(var(g_l));
% Compute F statistics
Fstat = (s_l^2)/(s_r^2);
Fcritical = finv(0.95,true_index-2,length(g)-true_index);
if Fstat < Fcritical
    alpha = 0;
    disp('Warning: NOT a genuine change point')
else
    disp('genuine change point')
end