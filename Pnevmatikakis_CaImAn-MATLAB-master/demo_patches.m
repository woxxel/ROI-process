% demo script for splitting the field of view in patches and processing in parallel
% with or without memory mapping. See also run_pipeline.m for the complete  
% pre-processing pipeline of large datasets

clear;
%% setup path to file and package
gcp;                                           % start local cluster
path_to_package = '../ca_source_extraction';   % path to the folder that contains the package
addpath(genpath(path_to_package));
tic             
%  filename = '/Users/epnevmatikakis/Documents/Ca_datasets/Neurofinder/neurofinder.02.00/images/neurofinder0200_rig.tif'; 
%  filename = '/media/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/Session01/median/images_00.tif';
filename = '/home/mizuta/AlexCode/test_data/246/ImDat.h5';
savefile = '/home/mizuta/AlexCode/test_data/246/resultsCNMF.mat'

        % path to file (assumed motion corrected)
        
is_memmaped = true;        % choose whether you want to load the file in memory or not

%% load file

if is_memmaped
    if exist([filename(1:end-3),'mat'],'file')
        data = matfile([filename(1:end-3),'mat'],'Writable',true);
    else
        sframe=1;					% user input: first frame to read (optional, default 1)
        num2read=[];					% user input: how many frames to read   (optional, default until the end)
        chunksize=[];                                   % user input: read and map input in chunks (optional, default read all at once)
        data = memmap_file(filename,sframe,num2read,chunksize);
        %data = memmap_file_sequence(foldername);
    end
    sizY = size(data,'Y');                    % size of data matrix
else
    T = 2000;                                 % load only a part of the file due to memory reasons
    data = read_file(filename,1,T);
    sizY = size(data);
end
    
%% Set parameters
patch_size = [272,272];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [32,32];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 1000;                % number of components to be found
tau = 10;                 % std of gaussian kernel (size of neuron) 
p = 2;                   % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;         % merging threshold

options = CNMFSetParms(...
    'd1',sizY(1),'d2',sizY(2),...
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'ssub',2,...
    'tsub',3,...
    'p',p,...                                   % order of AR dynamics
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,... 
    'spatial_method','regularized',...
    'cnn_thr',0.2,...
    'patch_space_thresh',0.25,...
    'min_SNR',2);

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% classify components 

rval_space = classify_comp_corr(data,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the space correlation test

try  % matlab 2017b or later is needed for the CNN classifier
    [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);
end

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std); % event exceptionality
ind_exc = (fitness < options.min_fitness);

keep = (ind_corr | ind_cnn) & ind_exc;

%% run GUI for modifying component selection (optional, close twice to save values)
Cn = correlation_image_max(data);  % background image for plotting
run_GUI = false;
if run_GUI
    Coor = plot_contours(A,Cn,options,1); close;
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};
end

%% re-estimate temporal components
A_throw = A(:,~keep);
C_throw = C(~keep,:);
A_keep = A(:,keep);
C_keep = C(keep,:);
options.p = 2;      % perform deconvolution
P.p = 2;
[A2,b2,C2] = update_spatial_components(data,C_keep,f,[A_keep,b],P,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components_fast(data,A2,b2,C2,f,P,options);

disp(fprintf('saving results to %',savefile))
save(savefile,'A2','C2','S2','b','b2','f','f2','Cn','YrA2','-v7.3');

toc
%% plot results
figure;
plot_contours(A2,Cn,options,0);
plot_components_GUI(data,A2,C2,b,f2,Cn,options);