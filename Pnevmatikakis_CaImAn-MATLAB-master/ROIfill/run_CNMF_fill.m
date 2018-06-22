
%%% input:
%%%     path        - path to h5-file
%%%     extents     
%%%     A           - structure with information about ROI cluster:
%%%             extents       - max/min positions of ROI-cluster boundaries (+margin)
%%%             footprint     - footprints of each ROI
%%%             total:        - projected footprints (needed?)
%%%             ct:           - number of ROIs (needed?)

%%% programm should load data from cropped h5 file and find positions of A in cropped region

function [A,C,keep,fitness] = run_CNMF_fill(path,A_in)
  
%    clear;
  tic
  
  %% load file
  %  gcp;                            % start cluster
  %  addpath(genpath('utilities'));
  %  addpath(genpath('deconvolution'));
  
  %  nam = 'demoMovie.tif';          % insert path to tiff stack here
  nam = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/Session01/images/ImagingData_MF1_LK1.h5';
  %  nam = '/media/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/Session01/median/images_00.tif';
  sframe=1;						% user input: first frame to read (optional, default 1)
  num2read=[];					% user input: how many frames to read   (optional, default until the end)
  Y = read_file_crop(nam,A_in.extents);
%    Y = Y./2^16;
  %Y = Y - min(Y(:)); 
  if ~isa(Y,'single');    Y = single(Y);  end         % convert to single
  
  [d1,d2,T] = size(Y);                                % dimensions of dataset
  d = d1*d2;                                          % total number of pixels
  
  %% Set parameters
  
%    K = 1000;                                           % number of components to be found
%    tau = 5;                                          % std of gaussian kernel (size of neuron) 
  p = 2;
  
  options = CNMFSetParms(...   
      'd1',d1,'d2',d2,...                         % dimensionality of the FOV
      'p',p,...                                   % order of AR dynamics
      'merge_thr',0.80,...                        % merging threshold  
      'nb',2,...                                  % number of background components    
      'min_SNR',3,...                             % minimum SNR threshold
      'space_thresh',0.5,...                      % space correlation threshold
      'cnn_thr',0.2...                            % threshold for CNN classifier    
      );
  
%    options.min_fitness
  options.min_fitness = -50;
  
  %% Data pre-processing
  [P,Y] = preprocess_data(Y,p);
  
  disp('preprocessing done')
  %% fast initialization of spatial components using greedyROI and HALS
  [Ain,Cin,bin,fin,center] = initialize_components_fill(Y,A_in,options,P);  % initialize
  
  % display centers of found components
  Cn =  correlation_image(Y); %reshape(P.sn,d1,d2);  %max(Y,[],3); %std(Y,[],3); % image statistic (only for display purposes)
  Cn_tmp = zeros(512,512);
  Cn_tmp(A_in.extents(1,1):A_in.extents(1,2),A_in.extents(2,1):A_in.extents(2,2)) = Cn;
%    figure;imagesc(Cn);
%        axis equal; axis tight; hold all;
%        scatter(center(:,2),center(:,1),'mo');
%        title('Center of ROIs found from initialization algorithm');
%        drawnow;
  
  
      %% manually refine components (optional)
  refine_components = false;  % flag for manual refinement
  if refine_components
      [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
  end
  
  A=Ain;
  b=bin;
  C=Cin;
  f=fin;
  Yr = reshape(Y,d,T);
%    figure('position',[100 100 1500 1200])
  for iter=1:1
%      subplot(4,3,iter)
    %% update spatial components
    [A,b,C] = update_spatial_components(Yr,C,f,[A,b],P,options);
    
    %% update temporal components
    P.p = 0;    % set AR temporarily to zero for speed
    [C,f,P,S,YrA] = update_temporal_components(Yr,A,b,C,f,P,options);
    
%      [A,C,K,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);
    
    %%% get fitness values, plot those and plot threshold
%      fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
%      fitness
%      hold on
%      
%      imagesc(Cn_tmp)
%      
%      for j=1:size(A,2)
%        A_tmp = zeros(512,512);
%        A_tmp(A_in.extents(1,1):A_in.extents(1,2),A_in.extents(2,1):A_in.extents(2,2)) = reshape(A(:,j),d1,d2);
%        A_tmp = A_tmp/sum(A_tmp(:));
%        centroid = [sum((1:512)*A_tmp),sum(A_tmp*(1:512)')];
%        contour(A_tmp,[0.01 0.01]*max(A_tmp(:)),'Color','r')
%        text(centroid(2),centroid(1),sprintf('fitness: %4.2g',fitness(j)))
%      end
%      hold off
%      xlim([A_in.extents(2,1),A_in.extents(2,2)])
%      ylim([A_in.extents(1,1),A_in.extents(1,2)])
%      
%      subplot(4,6,2*(iter+6)-1)
%      hold on
%      size(C)
%      XC = zeros(size(A,2));
%      for j=1:size(A,2)
%        plot(C(j,:))
%        for k=1:size(A,2)
%          XC(j,k) = xcorr(C(j,:),C(k,:),0,'coeff');
%        end
%      end
%      hold off
%      xlim([0 8989])
%      
%      subplot(4,6,2*(iter+6))
%      imagesc(XC)
%      set(gca,'CLim',[0,1])
%      title(sprintf('corr: %5.3g',mean(XC(:))))
%      colorbar
    
  end
  
  %% classify components
  
  rval_space = classify_comp_corr(Y,A,C,b,f,options);
  ind_corr = rval_space > options.space_thresh;           % components that pass the correlation test
                                          % this test will keep processes
  
  %% further classification with cnn_classifier
  try  % matlab 2017b or later is needed
      [ind_cnn,value] = cnn_classifier(A,[d1,d2],'cnn_model',options.cnn_thr);
  catch
      ind_cnn = true(size(A,2),1);                        % components that pass the CNN classifier
  end     
  
  %% event exceptionality
  fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std);
  ind_exc = (fitness < options.min_fitness);
  
  keep = [ind_corr,ind_cnn,ind_exc];
  %% select components
%    keep = (ind_corr | ind_cnn) & ind_exc;
%    
%    %% display kept and discarded components
%    A_keep = A(:,keep);
%    C_keep = C(keep,:);
%    
%  %    figure;
%  %        subplot(121); montage(extract_patch(A(:,keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15]);
%  %            title('Kept Components');
%  %        subplot(122); montage(extract_patch(A(:,~keep),[d1,d2],[30,30]),'DisplayRange',[0,0.15])
%  %            title('Discarded Components');
%    
%    %% merge found components
%  %    [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A_keep,b,C_keep,f,P,S,options);
%    [Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A,b,C,f,P,S,options);
%    
%    %%
%    display_merging = true; % flag for displaying merging example
%    if and(display_merging, ~isempty(merged_ROIs))
%        disp('plot')
%        i = 1; %randi(length(merged_ROIs));
%        ln = length(merged_ROIs{i});
%        figure;
%            set(gcf,'Position',[300,300,(ln+2)*300,300]);
%            for j = 1:ln
%  %                subplot(1,ln+2,j); imagesc(reshape(A_keep(:,merged_ROIs{i}(j)),d1,d2)); 
%  %                hold on
%                subplot(1,ln+2,j); imagesc(reshape(A(:,merged_ROIs{i}(j)),d1,d2)); 
%  %                subplot(1,ln+2,j); imagesc(reshape(A_in(:,merged_ROIs{i}(j)),d1,d2)); 
%                    title(sprintf('Component %i',j),'fontsize',16,'fontweight','bold'); axis equal; axis tight;
%            end
%            subplot(1,ln+2,ln+1); imagesc(reshape(Am(:,K_m-length(merged_ROIs)+i),d1,d2));
%                    title('Merged Component','fontsize',16,'fontweight','bold');axis equal; axis tight; 
%            subplot(1,ln+2,ln+2);
%  %                plot(1:T,(diag(max(C_keep(merged_ROIs{i},:),[],2))\C_keep(merged_ROIs{i},:))'); 
%                hold on
%                plot(1:T,(diag(max(C(merged_ROIs{i},:),[],2))\C(merged_ROIs{i},:))'); 
%                plot(1:T,(diag(max(C_in(merged_ROIs{i},:),[],2))\C_in(merged_ROIs{i},:))','r'); 
%                hold all; plot(1:T,Cm(K_m-length(merged_ROIs)+i,:)/max(Cm(K_m-length(merged_ROIs)+i,:)),'--k')
%                title('Temporal Components','fontsize',16,'fontweight','bold')
%            drawnow;
%    end
%  
%    %% refine estimates excluding rejected components
%  
%    Pm.p = p;    % restore AR value
%    [A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
%    [C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);
%  
%  
%  %    %% do some plotting
%  %    [A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
%  %    K_m = size(C_or,1);
%  %    [C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)
%  
%  %    figure;
%  %    [Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
%    %savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)
%  
%    %% display components
%  
%  %    plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options);
%  
%    %% make movie
%    if (0)  
%        make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)
%    end
%    toc
end