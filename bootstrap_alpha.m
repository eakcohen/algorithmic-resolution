function [alpha_min,alpha_max,alpha] = bootstrap_alpha(r,G,num_bs,interval,buffer)
%
% E.A.K. Cohen and N.M. Adams, Dept of Mathematics, Imperial College London
%
% This function estimates the algorithmic resolution limit using an averaged 
% pair correlation function and computes a bootstrap interval
% 
%
% INPUTS:
%   r              vector of radial distances at which pair correlation is evaluated
%   G              a matrix of dimension n x m where each row is an independent esimate
%                  of the pair correlation function. 
%   num_bs         number of bootstrap estimates e.g. 1000
%   interval_size  value between 80 and 100 indicating percentage size of bootstrap interval. e.g. 95% 
%   buffer         number of points at either end of g which is are not
%                  considered as a possible value of rho, the radius of correlatiojn
% OUTPUTS:
%   alpha_min      lower point of bootstrap interval
%   alpha_max      upper point of bootstrap interval
%   alpha          algorithmic resolution limit

% Take average of all estimates with which to make point estimate of alpha
g = mean(G,1);
% estimate alpha
alpha = est_alpha(r,g,buffer);

% create vector for bootstrap estimates
alphas = zeros(1,num_bs);


for ii = 1:num_bs
    % resample rows of G with replacement
    bs = randsample(size(G,1),size(G,1),1);
    BSg = mean(G(bs,:),1);
    % calculate and store bootstrap estimate
    alphas(ii) = est_alpha(r,BSg,buffer);
end

% sort bootstrap estimates
alphas = sort(alphas);

% find upper and lower point of interval% bootstrap interval
lower_point = floor((100-interval)*num_bs./100);
upper_point = floor(interval*num_bs./100);

alpha_min = alphas(lower_point);
alpha_max = alphas(upper_point);