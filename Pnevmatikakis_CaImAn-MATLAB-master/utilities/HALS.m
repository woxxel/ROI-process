function [A, C, b, f] = HALS(Y, A, C, b, f, params)
%% Hierarchical alternating least square method for solving NMF problem  
% Y = A*C + b*f 

%input: 
%   Y:      d1 X d2 X T, raw data. It will be reshaped to (d1*d2) X T in this
%   function 
%   A:      (d1*d2) X K, initial value of spatial components 
%   C:      K X T, initial value of temporal components 
%   b:      (d1*d2) X 1, initial value of background spatial component 
%   f:      1 X T, initial value of background temporal component
%   params: parameters used in this function. 
%       bSiz:   blur size. A box kernel (bSiz X bSiz) will be convolved
%       with each neuron's initial spatial component, then all nonzero
%       pixels will be picked as pixels to be updated, and the rest will be
%       forced to be 0. 
%       maxIter: maximum iteration of iterating HALS. 

% Author: Pengcheng Zhou, Carnegie Mellon University, based on a python
% implementation from Johannes Friedrich, Columbia University, 2015.

%% parameters
if isfield(params, 'maxIter'), maxIter = params.maxIter; else maxIter=5; end
if isfield(params, 'search_method'); method=params.search_method; else method='ellipse'; end
if and(isfield(params, 'bSiz'), strcmpi(method, 'dilate'))
    params.se = strel('disk', params.bSiz);
end
if isfield(params,'init_method'); init_method = params.init_method; else init_method = 'greedy'; end

Y = reshape(Y,size(A,1),[]);
% search locations
if strcmpi(init_method,'greedy') || strcmpi(init_method,'greedy_corr');
    IND = determine_search_location(A,method,params);
else
    IND = true(size(A));
end
%% update spatial and temporal components neuron by neurons
Yac = max(Y - b*f,0);
for miter=1:maxIter
    % update neurons
    % remove ROIs with no variance (= all same entries)
    ind_del = find(std(A,0,1)==0); 
    A(:, ind_del) = []; 
    C(ind_del, :) = []; 
    IND(:, ind_del) = []; 
    
    %   temporal update
    C = HALS_temporal(Yac, A, C, 5);
    
    % remove ROIs with no temporal variance (= static)
    ind_del = find(std(C,0,2)==0); 
    A(:, ind_del) = []; 
    C(ind_del, :) = []; 
    IND(:, ind_del) = []; 
    
    %   spatial
    A = HALS_spatial(Yac, A, C, IND, 5);
end    

% update background
Ybg = Y-A*C;
for miter = 1:maxIter
    % temporal
    f = HALS_temporal(Ybg, b, f, 5);
    % spatial 
    b = HALS_spatial(Ybg, b, f, [], 5);     
end