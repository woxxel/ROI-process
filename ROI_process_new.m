

function ROI_process_new(start_idx,end_idx,sz_median,do_dewarp,scattered,redo)
    
    pathMouse = uigetdir('Choose a mouse folder to process (completely!)');
%      pathMouse = '/media/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245';
%      pathMouse = '/home/mizuta/AlexCode/test_data';
%      pathMouse = sprintf('/media/mizuta/Analyze_AS1/%d',mouse);
%      pathMouse = sprintf('/media/mizuta/Analyze_AS1/1-photon/%d',mouse);
    [sessionList, nSessions] = getSessions(pathMouse);
    
    %% construct suffix for filenames if not specified
%      if (nargin < 5)
      suffix_MF = sprintf('_MF%d',sz_median);
      suffix_LK = sprintf('_LK%d',do_dewarp);
      suffix = sprintf('%s%s',suffix_MF,suffix_LK);
%      else
%        suffix_MF = '';
%        suffix_LK = '';
%      end
    
    parameter = set_parameter(sz_median);
    
    for s = start_idx:end_idx
	
      path = set_paths(sessionList{s},suffix_MF,suffix_LK,suffix)
      
      disp(sprintf('\t ### Now processing session %d ###',s))
      path.handover = path.images;
      
      if ~exist(path.H5,'file') || redo
        if scattered
          pathTmp = pathcat(path.session,'images');
          create_tiff_stacks(pathTmp,path.images,1000);
        end
        
        pathTmp = pathcat(path.session,'imageStacks');
        if isempty(dir(pathcat(pathTmp,'*.tif')))
%  	  disp('no tifs found')
          if ~isempty(dir(pathcat(pathTmp,'*.raw')))
	    pathRaw = dir(pathcat(pathTmp,'*.raw'));
	    pathRaw = pathcat(pathTmp,pathRaw.name);
            raw2tiffstacks(pathRaw);
          end
        end
        
        %% median filtering
        if sz_median
          if ~exist(path.median,'dir') || redo
              median_filter(path.images,path.median,parameter);
          else
              disp(sprintf('Path: %s already exists - skip median calculation',path.median))
          end
          path.handover = path.median;
        else
          disp('---- median filtering disabled ----')
        end
        
        % image alignment
        if do_dewarp
          if ~exist(path.LKalign,'dir') || redo
              tiff_align(path.handover,path.LKalign);
          else
              disp(sprintf('Path: %s already exists - skip image dewarping',path.LKalign))
          end
          path.handover = path.LKalign;
        else
          disp('---- LK-dewarping disabled ----')
        end
        
        tiff2h5(path.handover,path.H5);
      end
      
      if ~exist(path.reduced,'file') || redo
          reduce_data(path.H5,path);
      else
          disp(sprintf('Path: %s already exists - skip reduced image calculation',path.reduced))
      end
      
      
      if ~exist(path.CNMF,'file') || redo
          disp('do CNMF')
          CNMF_frame(path,parameter.npatches,parameter.K,parameter.tau,0);
      else
          disp(sprintf('Path: %s already exists - skip CNMF',path.CNMF))
      end
      
      
%        if ~exist(path.CNMF_post,'file')
%            ROI_post_procession_CNMF(path,parameter);
%        else
%            disp(sprintf('Path: %s already exists - skip CNMF post-procession',path.CNMF_post))
%        end
%        rmdir(path.images,'s');
    end
end



function [parameter] = set_parameter(sz_median)
    
    parameter = struct();
    
    %% parameter for preCellDeconv
    %% for tiff-image->stack conversion
    parameter.nsubFiles = 2000;
    
    %% median filtering
    parameter.filtersize = [sz_median,sz_median,1];
    
    %% parameter for CNMF
    parameter.npatches = 1;                      % how many patches are processed in parallel
    K = 100;                                    % first guess of the number of neurons to be found
    parameter.K = ceil(K/parameter.npatches);                           
    parameter.tau = 8;                           % guess of average neuron radius (in pixel)
    
    %% parameter for post-procession
    parameter.sd = 40;                 % multiple of STD for thresholding ROI images
    parameter.thr_size = [20 400];     % upper and lower threshold for ROI size (realistic pyramidal neuron size ~20mum length, wikipedia)
    parameter.thr_pos = [5 507];       % threshold for ROI-position (5 off the border)
    parameter.perc = 0.2;              % threshold for fraction of common pixels between ROIs
    
    %% parameter for session matching
    parameter.max_dist = 12;
    parameter.num_ses_thr = 3;
    parameter.SI_thr = 0.5;
end


function [path] = set_paths(pathIn,suffix_MF,suffix_LK,suffix)
    
    path = struct();
    
    path.session = pathIn;
    
    %% for pre-procession
    path.images = pathcat(path.session,'imageStacks');
    path.median = pathcat(path.images,sprintf('median%s',suffix_MF));
    path.LKalign = pathcat(path.images,sprintf('LKalign%s',suffix_MF));
    path.reduced = pathcat(path.session,sprintf('reduced%s.mat',suffix));
    
    %% for CNMF algorithm
    path.H5 = pathcat(path.session,sprintf('ImagingData%s.h5',suffix));
    path.CNMF = pathcat(path.session,sprintf('resultsCNMF%s.mat',suffix));
    path.plotContour = pathcat(path.session,sprintf('contour%s.png',suffix));
    path.CNMF_post = pathcat(path.session,sprintf('ROIs_CNMF%s.mat',suffix));    
    
end