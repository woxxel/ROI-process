

function ROI_post_procession_CNMF(path,parameter)
    
    tic
    
    ROIs_ini = load(path.CNMF,'A2');
    ROIs_ini = ROIs_ini.A2;
    
    %%% get dimensions of data
    nROI = size(ROIs_ini,2);
    ROIcross = zeros(0);
    
    ROIs = struct('filter',[],'area',[],'norm',[],'centroid',[]);
    ncells = 0;
    
    for n=1:nROI
	ROI_tmp = reshape(ROIs_ini(:,n),512,512);         % might need to transpose (col-major vs row-major storage)
	
	%% finding and separating blobs
	ROI_BWL = bwlabel(full(ROI_tmp),8);	% shouldn't be needed
	
	blobs = regionprops(ROI_BWL, 'PixelList');
	nBlobs = size(blobs, 1);
	
	for i = 1:nBlobs
	    ROI_cd = struct('filter',[],'area',[],'norm',[],'centroid',[]);
	    ROI_cd.filter = sparse((ROI_BWL == i).*ROI_tmp);
	    ROI_cd = ROI_update_CNMF(ROI_cd,'normalize','norm','centroid','area');
	    if nBlobs > 1
	      disp(sprintf('blob #%d size: %d',i,ROI_cd.area))
	    end
	    
	    %%% check the updated blob for area- and position threshold (if so, blob -> ROI)
	    check_area = (parameter.thr_size(1) < ROI_cd.area) && (ROI_cd.area < parameter.thr_size(2));
	    check_pos = (sum(parameter.thr_pos(1) < ROI_cd.centroid) == 2) && (sum(ROI_cd.centroid < parameter.thr_pos(2))==2);
	    
	    %% checks should be less about "thresholding". Rather get something non-parametric:
	    %%% assume gaussian distribution of areas
	    %%% assume certain level of activity (?, a bit tricky, since there are also low-active neurons that might be highly active in other sessions)
	    %%% assume certain shape of PSD (?, a bit tricky, as CNMF already kinda tries to characterize by this)
	    %%% assume independence from background (?, a bit tricky, as CNMF tries to maximize this anyhow)
	    %%% 
	    if check_area && check_pos  %% when all checks are passed, store ROI
		ncells = ncells + 1;
		ROIs(ncells) = ROI_cd;
	    end
	end
    end
    
    disp(sprintf('ROIs found: %d / %d',ncells,nROI))
    
    save(path.CNMF_post,'ROIs');
    toc
end


        