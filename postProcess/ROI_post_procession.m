

function ROI_post_procession(pathMouse,parameter)
    
    [sessionList, nSessions] = getSessions(pathMouse);
    %% read in ROI data from each session and translate to more efficient structure
    nSessions = 10;
    
    for s=1:nSessions
        tic
        pathH5 = dir(pathcat(sessionList{s},'*c.h5'));
        pathH5 = pathcat(sessionList{s},pathH5.name);
        
        pathResults = pathcat(sessionList{s},'ROIs.mat');
%          if exist(pathResults)==2
%              disp(strcat('Results file ',pathResults,' already exists'))
%              continue
%          end
        ROIs_ini = h5read(pathH5,'/MERGED/A');  % path in h5-file should be changed
        
        
%          pathH5 = dir(pathcat(sessionList{s},'*c.mat'));
%          pathH5 = pathcat(sessionList{s},pathH5.name);
%          pathResults = pathcat(sessionList{s},'ROIs_CNMF.mat');
        
        
        
        %%% get dimensions of data
        nROI = size(ROIs_ini,3);
        ROIcross = zeros(0);
%  %          
%  %  %          c = parcluster
%  %  %          p = gcp('nocreate'); % If no pool, do not create new one.
%  %  %          if isempty(p)
%  %  %              poolsize = 0;
%  %  %              parpool(c)
%  %  %          else
%  %  %              poolsize = p.NumWorkers
%  %  %          end
%  %  %          poolsize
        ROIs = struct('extent',[],'filter',[],'area',[],'norm',[],'centroid',[]);
        ncells = 0;
%         
        
        for n=1:nROI
%              [labindex, numlabs]
            %%% thresholding background noise
            ROI_tmp = ROIs_ini(:,:,n);         % might need to transpose (col-major vs row-major storage)
            thr = mean(ROI_tmp(:)) + parameter.sd*std(ROI_tmp(:));
            ROI_tmp(ROI_tmp<thr) = 0;
            
            %% finding and separating blobs
            ROI_BWL = bwlabel(ROI_tmp,8);
            
            blobs = regionprops(ROI_BWL, 'PixelList');
            nBlobs = size(blobs, 1);
                
            for i = 1:nBlobs
                blob_tmp = ROI_BWL == i;
                
                blobPix = blobs(i).PixelList;
                
                %%% define area around blob
                minBlob = [min(blobPix(:,2)),min(blobPix(:,1))]; % y,x
                maxBlob = [max(blobPix(:,2)),max(blobPix(:,1))];
                %%% and fill the holes
                blob_tmp(minBlob(1):maxBlob(1),minBlob(2):maxBlob(2)) = imfill(blob_tmp(minBlob(1):maxBlob(1),minBlob(2):maxBlob(2)),'holes');
                
                %%% check the updated blob for area- and position threshold (if so, blob -> ROI)
                blob = regionprops(blob_tmp, 'Area','Centroid');       %%% check, if only one blob comes out!
                
                check_area = (parameter.thr_size(1) < blob(1).Area) && (blob(1).Area < parameter.thr_size(2));
                check_pos = (sum(parameter.thr_pos(1) < blob(1).Centroid) == 2) && (sum(blob(1).Centroid < parameter.thr_pos(2))==2);
                
                if check_area && check_pos  %% when all checks are passed, store filter area and position
                    ncells = ncells + 1;
                    
%                      ROI_add = struct('extent',[],'filter',[],'area',[],'norm',[],'centroid',[]);
%                      
%                      ROI_add.extent = vertcat(minBlob, maxBlob);
%                      ROI_add.filter = ROI_tmp(minBlob(1):maxBlob(1),minBlob(2):maxBlob(2));
%                      ROI_add.area = blob(1).Area;
%                      ROI_add = ROI_update(ROI_add,'norm','centroid','normalize');
%                      
%                      ROIs = [ROIs, ROI_add];
                    
                    ROIs(ncells).extent = vertcat(minBlob, maxBlob);
                    ROIs(ncells).filter = ROI_tmp(minBlob(1):maxBlob(1),minBlob(2):maxBlob(2));
                    ROIs(ncells).area = blob(1).Area;
                    ROIs(ncells) = ROI_update(ROIs(ncells),'normalize','norm','centroid');
%                      
%                      %%% compute the crosstalk between already detected ROIs
                    for m = 1:ncells-1
                        cross_tmp = ROIs_dotproduct(ROIs(ncells),ROIs(m),'bw');
                        
                        ROIcross(ncells,m) = cross_tmp;
                        ROIcross(m,ncells) = cross_tmp;
                    end
                end
            end
        end
        
        %%% now discard the ROIs with large crosstalk (keep larger ones)    // or rather smaller ones??
        passfail = false(ncells);
        size_ratio = zeros(1,ncells);
        for n=1:ncells
            cross_scaled = ROIcross(n,:)/ROIs(n).area;

            for m = 1:ncells
                size_ratio(1,m) = ROIs(n).area/ROIs(m).area;
            end
            passfail(n,:) = ((cross_scaled <= parameter.perc) | ((cross_scaled > parameter.perc) & (size_ratio > 1)));
            passfail(n,n) = true;
        end
        passfail(ncells+1:end,:) = false;
        passfail(:,ncells+1:end) = false;
        pass = sum(passfail,2) == ncells;
        
        ROIs = ROIs(pass);
        disp(sprintf('Session %d, ROIs found: %d / %d',s,sum(pass),nROI))
        
        save(pathResults,'ROIs');
        toc
    end
end


        