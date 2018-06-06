function [Ain, Cin, b_in, f_in, center, res] = greedyROI_fill(Y, A, params)
% component initialization using a greedy algorithm to identify neurons in 2d or 3d calcium imaging movies
%
% Usage:     [Ain, Cin, bin, fin, center, res] = greedyROI2d(data, K, params)
%
% Input:
% Y          d1 x d2 x (d3 x) T movie, raw data, each column is a vectorized movie
% K          number of neurons to extract (if K is a vector then the algorithm is run multiple times with different parameters)
% params     tuning parameter for fine-tuning the shape (optional)
%            params.nIter: number of iterations for shape tuning (default 5)
%            params.gSig: variance of Gaussian kernel to use (default 5) If params.gSig is a cell, then the algorithm is run with multiple times with different parameters
%            params.gSiz: size of kernel (default 2*gSig+1)
%            params.nb: rank of background component (default 1)
%            params.save_memory: flag for processing data in chunks to save memory (default 0)
%            params.windowSiz: size of spatial window when computing the median (default 32 x 32)
%            params.chunkSiz: number of timesteps to be processed simultaneously if on save_memory mode (default: 100)
%            params.med_app: number of timesteps to be interleaved for fast (approximate) median calculation (default: 1, no approximation)
%            params.rolling_sum: flag for using rolling sum to detect new components (default: True)
%            params.rolling_length: length of rolling window (default: 100)
% ROI_list   Kn x 2 (or 3) matrix with user specified centroids. If this are present, then the algorithm only finds components around these centroids

%Output:
% Ain        (d) x K matrix, location of each neuron
% Cin        T x K matrix, calcium activity of each neuron
% center     K x 2 (or 3) matrix, inferred center of each neuron
% bin        (d) X nb matrix, initialization of spatial background
% fin        nb X T matrix, initalization of temporal background
% res        d1 x d2 x (d3 x) T movie, residual
%
% Author: Yuanjun Gao with modifications from Eftychios A. Pnevmatikakis and Weijian Yang

  use_sum = false;
%  if nargin < 4 || isempty(ROI_list)
  user_ROIs = 0;
%  else
%      user_ROIs = 1;
%      K = size(ROI_list,1);
%  end

dimY = ndims(Y) - 1;  % dimensionality of imaged data (2d or 3d)
sizY = size(Y);
T = sizY(end);        % # of timesteps
dx = sizY(1:dimY);    % # of voxels in each axis
d = prod(dx);         % total # of voxels  

%  if ~exist('K', 'var'),  %find K neurons
%      K = 30;  
%      warning(['number of neurons are not specified, set to be the default value', num2str(K)]);
%  end

if ~exist('params', 'var') params = []; end

if ~isfield(params, 'gSig') || isempty(params.gSig); 
    if dimY == 2; params.gSig = [5, 5]; else params.gSig = [5,5,5]; end
elseif length(params.gSig) == 1, params.gSig = params.gSig + zeros(1,dimY);
    if dimY == 3; params.gSig(3) = params.gSig(3)/2; end
end

if ~isfield(params, 'gSiz') || isempty(params.gSiz); 
    if ~iscell(params.gSig)
        params.gSiz = ceil(2*params.gSig + 1);
    else
        for j = 1:length(params.gSig)
            params.gSiz{j,1} = 2*params.gSig{j}+1; %cellfun(@times,params.gSig{j},num2cell(ones(size(params.gSig{j}))*2));
        end
    end
elseif length(params.gSiz) == 1, params.gSiz = params.gSiz + zeros(1,dimY);
    if dimY == 3; params.gSiz(3) = ceil(params.gSiz(3)/2); end
end

if isfield(params,'ssub'); 
    if ~iscell(params.gSig); params.gSig(1:2) = params.gSig(1:2)/params.ssub; params.gSiz(1:2) = ceil(params.gSiz(1:2)/params.ssub);
    else
        for j = 1:length(params.gSig)
            params.gSig{j,1} = params.gSig{j}/params.ssub; %cellfun(@times,params.gSig{j},num2cell(ones(size(params.gSig{j}))/params.ssub));
            params.gSiz{j,1} = params.gSiz{j}/params.ssub; %cellfun(@times,params.gSiz{j},num2cell(ones(size(params.gSiz{j}))/params.ssub));
        end
    end
end

if ~isfield(params,'nb'), nb = 1; else nb = params.nb; end

if ~isfield(params, 'nIter'), nIter = 5; 
else nIter = params.nIter; end

if ~isfield(params, 'save_memory'), save_memory = 0;
else save_memory = params.save_memory; end
    
if ~isfield(params, 'chunkSiz'), chunkSiz = 100; else chunkSiz = params.chunkSiz; end
if ~isfield(params, 'windowSiz'), windowSiz = 32; else windowSiz = params.windowSiz; end
if ~isfield(params, 'med_app'), med_app = 1; else med_app = params.med_app; end

if ~isfield(params,'rem_prct') || isempty(params.rem_prct); params.rem_prct = 20; end

Tint = 1:med_app:T;
% if save_memory
%     med = zeros(M,N);
%     for ii = 1:ceil(M/windowSiz)
%         intx = (ii-1)*windowSiz+1:min(ii*windowSiz,M);
%         for jj = 1:ceil(N/windowSiz)
%             inty = (jj-1)*windowSiz+1:min(jj*windowSiz,N);
%             med(intx,inty) = median(data(intx,inty,Tint),3);
%         end
%     end 
% else
%if dimY == 2; med = median(Y(:,:,Tint), 3); else med = median(Y(:,:,:,Tint), 4); end

%  figure
%  for i=1:6
  params.rem_prct = 10;
  %%% obtain percentiles of pixel values along time axis
  if dimY == 2
    med = prctile(Y(:,:,Tint),params.rem_prct, 3);
  else
    med = prctile(Y(:,:,:,Tint),params.rem_prct,4);
  end
  %  size(med)
%    subplot(2,3,i)
%    imagesc(med)
%    set(gca,'YDir','normal')
%  end

% end
disp('(initialize) components to be found:')
K = A.ct;

Y = bsxfun(@minus, Y, med);

if iscell(params.gSig); params.gSig = cell2mat(params.gSig); params.gSiz = cell2mat(params.gSiz); end

Ain = sparse(d,K);

Cin = zeros(K,T);
center = zeros(K,dimY);


gSig = params.gSig(1,:);
gSiz = params.gSiz(1,:);

gHalf = floor(gSiz / 2); %half size of the kernel, used to calculate margin
gSiz = 2 * gHalf + 1; %actual size

Atemp = spalloc(d,K,K*ceil(prod(gSiz))); %zeros(M, N, K(r));
%basis = spalloc(M*N,K,K*prod(gSiz));
trace = zeros(T, K);
centers = zeros(K, dimY);
    
%scan the whole image (only need to do this at the first iteration)
%%% rho obtains filtered image - pixel value along one time can be used as calcium trace
rho = imblur(Y, gSig, gSiz, dimY, save_memory, chunkSiz); %covariance of data and basis

if ~params.rolling_sum
    v = sum(rho.^2, dimY+1); %variance explained
else
    % running maximum
    avg_fil = ones(1,params.rolling_length)/params.rolling_length;
    rho_s = filter(avg_fil,1,rho.^2,[],dimY+1);
    v = max(rho_s,[],dimY+1);
end

for k = 1:K
  
    %%% for each ROI, obtain guess of temporal trace (filter x data) and feed this as initial guess
    centers(k,:) = A.centroid(k,:);
  
    %%% obtain window to run pseudo-NMF
    [y_idx x_idx] = find(A.footprint(:,:,k));
    iSig(1,1) = max(1,min(y_idx));
    iSig(1,2) = min(dx(1),max(y_idx));
    iSig(2,1) = max(1,min(x_idx));
    iSig(2,2) = min(dx(2),max(x_idx));
    
    iSigLen = iSig(:,2) - iSig(:,1) + 1;
    
    
    c_hat = [round(centers(k,1)) - A.extents(1,1)+1, round(centers(k,2)) - A.extents(2,1)+1];
    c_upper = [round(centers(k,1)) - A.extents(1,2), round(centers(k,2)) - A.extents(2,2)];
    
    %%% obtain first guesses: temporal guess = value at center pixel, spatial guess = pixels, that fluctuate with center
    if dimY == 2
        dataTemp = Y(iSig(1,1):iSig(1,2), iSig(2,1):iSig(2,2), :);
        if all(c_hat > 0) && all(c_upper < 0)
          traceTemp = squeeze(rho(c_hat(1), c_hat(2), :));
        else
          traceTemp = (reshape(A.footprint(:,:,k),1,prod(dx))*reshape(Y,prod(dx),T))';
        end
    else
        dataTemp = Y(iSig(1,1):iSig(1,2), iSig(2,1):iSig(2,2), iSig(3,1):iSig(3,2), :);
%              traceTemp = squeeze(rho(iHat(1), iHat(2), iHat(3), :));
    end
    
%          [coef, score] = finetune(dataTemp, cin_guess, nIter);
    [coef, score] = finetune(dataTemp, traceTemp, nIter);               %%% run semi-NMF on first guess    
    
    %%% write spatial guess to basis
    basis = zeros(dx);
    basis(iSig(1,1):iSig(1,2), iSig(2,1):iSig(2,2)) = coef;
    Atemp(:,k) = basis(:);    %%% needed?
    
    plt = false;
    if plt && A.status(k)
      figure
      subplot(2,2,1)
      plot(1:8989,traceTemp)
        
      subplot(2,2,2)
      imagesc(A.footprint(:,:,k))
      
      subplot(2,2,3)
      plot(1:8989,score)
      
      subplot(2,2,4)
      imagesc(reshape(Atemp(:,k),dx(1),dx(2)))
      
      if A.status(k)
        suptitle('to be filled')
      else
        suptitle('remainder ROI')
      end
    end
    
    %%% write first guess of temporal trace
    trace(:, k) = score';
    
    %%% activate this, to remove already assigned data-to-ROI from data-to-be-considered
    
    %%% update residual
    dataSig = bsxfun(@times, coef, reshape(score, [ones(1,dimY),T]));     %%% reconstruct data from temp. and spat. guesses
    Y(iSig(1,1):iSig(1,2), iSig(2,1):iSig(2,2), :) = Y(iSig(1,1):iSig(1,2), iSig(2,1):iSig(2,2), :) - dataSig;
    
    %%% update rho for next Ca trace
%          if k < K
%              dataTemp = imblur(basis, gSig, gSiz, dimY, 0, chunkSiz);
%              rhoTemp = bsxfun(@times, dataTemp, reshape(score, [ones(1,dimY),T]));
%              rho = rho  - rhoTemp;
%          end
end

Ain(:,1:K) = Atemp; %sparse(reshape(basis,d,K(r)));  
Cin(1:K,:) = trace';
center(1:K,:) = centers;
%      pause(5)

res = reshape(Y,d,T) + repmat(med(:),1,T);

%% clear data matrix from local memory (avoid out-of-memory? see greedyROI_corr.m)
clear Y
%  clear Y_res

%[b_in,f_in] = nnmf(max(res,0),nb);
%[b_in,f_in] = nnsvd(max(res,0),nb);
%  f_in = [mean(res);rand(nb-1,T)];
f_in = [mean(res);rand(nb-1,T)*max(res(:))];
for nmfiter = 1:100
    b_in = max((res*f_in')/(f_in*f_in'),0);
    f_in = max((b_in'*b_in)\(b_in'*res),0);
end


function [ain,cin] = finetune(data,cin,nIter)
    
    %rank-1 semi-NMF to fine tune the inferred components
    %
    %Input:
    %data   d1 x d2 x (d3 x) T matrix, small patch containing one neuron
    %trace  initial value for trace
    %nIter  number of coordinate descent steps
    %
    %Output:
    %basis  d1 x d2 (x d3) matrix, result of the fine-tuned neuron shape
    %trace  1 x T matrix, result of the neuron
    
    if ~exist('nIter', 'var'), nIter = 1; end
    data = reshape(data,prod(iSigLen),T);
    for iter = 1:nIter
        %nc = norm(cin)^2;
        ain = max(data*cin,0)/norm(cin);
        an = norm(ain);
        if an > 0
            ain = ain/an;
        else
            fprintf('found degenerate component!\n')
            %break
        end
        cin = data'*ain;
    end
    ain = reshape(ain,iSigLen(:)');
    cin = cin(:)';
end

function data = imblur(data, sig, siz, nDimBlur, save_memory, chunkSiz)
%Gaussian blur for high dimensional data
%Input:
%data       original data
%sig        std of gaussian kernel
%siz        size of kernel
%nDimBlur   number of dims to blur (default: ndims(data) - 1)
%
%Output:
%data       result after the Gaussian blur

if ~exist('nDimBlur', 'var'), nDimBlur = ndims(data) - 1; 
else nDimBlur = min(nDimBlur, ndims(data)); end

if length(sig) == 1, sig = sig * ones(1,nDimBlur); end
assert(nDimBlur == length(sig));

if length(siz) == 1, siz = siz * ones(1,nDimBlur); end
assert(nDimBlur == length(siz));

for i = 1:nDimBlur,
    if sig(i) > 0,
        x = -floor(siz(i) / 2):floor(siz(i) / 2);
        H = exp(-(x.^2/ (2 * sig(i)^2)));
        H = H' / norm(H(:));
        if nDimBlur > 1,
            indH = 1:nDimBlur; indH(i) = 1; indH(1) = i;
            H = permute(H, indH);
        end
        if save_memory
            L = size(data,ndims(data));
            for ci = 1:ceil(L/chunkSiz)
                int = (ci-1)*chunkSiz+1:min(ci*chunkSiz,L);
                data(:,:,int) = imfilter(data(:,:,int), H, 'same', 0);
            end
        else
            data = imfilter(data, H, 'same', 0);
        end
    end
end
end

end
