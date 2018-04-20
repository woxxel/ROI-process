

function ROI_match_sessions_CNMF(pathMouse,parameter)
    
    tic
    [sessionList, nSessions] = getSessions(pathMouse);
    pathShift = dir(pathcat(pathMouse,'all_shift*.txt'));
    pathResults = pathcat(pathMouse,'All_ROIs_CNMF.mat');
    pathHisto = pathcat(pathMouse,sprintf('matchSessions_bin_d=%4.2f.mat',parameter.max_dist));
    
%      if exist(pathResults)==2
%          disp(strcat('Results file ',pathResults,' already exists'))
%          return
%      end
    
    shift = load(pathcat(pathMouse,pathShift.name));
    
    nSessions = 10;
    ncells = 0;
    cum_cell = 0;
    ROIs = [];

    %% load ROIs from all Sessions
    for s = 1:nSessions
        loadDat = load(pathcat(sessionList{s},'ROIs_CNMF.mat'),'ROIs');
        loadROIs = loadDat.ROIs;
        
        ncells_session = length(loadROIs);
        
        for n = ncells_session:-1:1 % start with last, to allow for preallocation of structures
            loadROIs(n).filter = imtranslate(full(loadROIs(n).filter),[shift(s,2),shift(s,1)]);
            loadROIs(n) = ROI_update(loadROIs(n),'truncate');
  %              loadROIs(n).extent = loadROIs(n).extent + shift(s,:);   % correct for shift in between sessions
            loadROIs(n) = ROI_update_CNMF(loadROIs(n),'normalize','norm','centroid','area');       % update ROIs
              %%% should one take into account additional shifted borders due to dewarping/alignment?
        end
        
        ROIs = [ROIs, loadROIs];
        ncells = ncells + ncells_session;
        cum_cell = [cum_cell ncells];
    end
    
    %% register the following to estimate probability of same ROI
    %% store values in histogram (how fine?)
    %% fit 2 functions to histogram and obtain 95% (or whatever) probability of match)
    
    nbins = 50;
    
    if exist(pathHisto)==2
        disp('load Histogram from previously processed data')
        load(pathHisto)
    else
        [prg_str1 prg_str2] = prepare_progress_report('ROIs processed: ',ncells);
        
        % set up histograms and according ranges for gathering data
        
%          c = parcluster
%          p = gcp('nocreate'); % If no pool, do not create new one.
%          if isempty(p)
%              poolsize = 0;
%              parpool(c)
%          else
%              poolsize = p.NumWorkers
%          end
        
        
        %% ziv gets bimodal histograms for all data! why? cause of restriction to <12mum?
        %% set up histograms
        
        max_dist = parameter.max_dist*512/530.684;
        histo.dist = zeros(nbins,2);
        
        max_corr = 1;
        histo.corr = zeros(nbins,2);
        
%          for n = 1:ncells
%              ROIs(n) = ROI_update(ROIs(n),'centroid');
%          end
        % set up variables
        n_border = 0;
        sn = 0;
        NN = zeros(ncells,nSessions,3);
        NN(:,:,2) = max_dist; %% need one NN per session! % or need ones(ncells,1)* ?
        tic
        for n = 1:ncells
            
    %          histo_dist_tmp = zeros(nbins,1);
    %          histo_corr_tmp = zeros(nbins,1);
            
            if n > n_border
                sn = sn+1;
                n_border = cum_cell(sn+1);
            end
            
            sm = sn;
            m_border = cum_cell(sm+1);
            %% here, define tmp-histograms to add to general one if not working in parfor?!
            for m = n+1:ncells
                
                if m > m_border
                    sm = sm+1;
                    m_border = cum_cell(sm+1);
                end
                if m > n_border   % remove iterating over same-session
                    
                    %% add centroid distance to histogram
                    tmp_dist = norm(ROIs(n).centroid - ROIs(m).centroid);
                    if tmp_dist <= max_dist
                        idx = max(1,ceil(nbins*tmp_dist/max_dist));
                        histo.dist(idx,2) = histo.dist(idx,2) + 1;
                        
                        %% this is not the "real" correlation, they are using in Ziv - try that one
                        %% add correlation value to histogram
                        tmp_corr = dot(ROIs(n).filter(:),ROIs(m).filter(:))/(ROIs(n).norm*ROIs(m).norm);
                        if tmp_corr > 0
                            idx = ceil(nbins*tmp_corr/max_corr);
                            histo.corr(idx,2) = histo.corr(idx,2) + 1;
                        end
                        
                        if tmp_dist < NN(n,sm,2)
                            NN(n,sm,:) = [m,tmp_dist,tmp_corr];
                        end
                        
                        if tmp_dist < NN(m,sn,2)
                            NN(m,sn,:) = [n,tmp_dist,tmp_corr];
                        end
                    end
                    if m==m_border
                        if NN(n,sm,2) < max_dist
                            idx = max(1,ceil(nbins*NN(n,sm,2)/max_dist));
                            histo.dist(idx,1) = histo.dist(idx,1) + 1;
                            histo.dist(idx,2) = histo.dist(idx,2) - 1;
                            if NN(n,sm,3) > 0
                                idx = max(1,ceil(nbins*NN(n,sm,3)/max_corr));
                                histo.corr(idx,1) = histo.corr(idx,1) + 1;
                                histo.corr(idx,2) = histo.corr(idx,2) - 1;
                            end
                        end
                    end
                    %% could add others?! but what...?
                    %% ratio of sizes?
                    %% shape similarity (?)
                    %% other stuff?
                end
            end
            
            if mod(n,100)==0
                now_time = toc;
                fprintf(1,prg_str1,n)
                fprintf(1,prg_str2,now_time)
            end
        end
        
        
        
        disp(strcat('saving results in ', pathHisto))
        save(pathHisto,'histo','NN','-v7.3')
    end
    
    %% normalize histograms
%      histo.dist(:,1) = histo.dist(:,1)/sum(histo.dist(:,1));
%      histo.dist(:,2) = histo.dist(:,2)/sum(histo.dist(:,2));
%      histo.corr(:,1) = histo.corr(:,1)/sum(histo.corr(:,1));
%      histo.corr(:,2) = histo.corr(:,2)/sum(histo.corr(:,2));
    
    histo.dist_x = linspace(0,parameter.max_dist,nbins);
    histo.corr_x = linspace(0,1,nbins);
    
    fig_dist = figure;
    b_dist = bar(histo.dist_x,histo.dist,1);
    b_dist(1).FaceColor = 'green';
    b_dist(2).FaceColor = 'red';
    xlabel('distance in \mu m')
    ylabel('number of pairs')
    legend('nearest neighbour','others')
    title(sprintf('max dist = %4.2f',parameter.max_dist))
    saveas(fig_dist,pathcat(pathMouse,sprintf('matchHistoDist_d=%4.2f_nSes=%d.jpg',parameter.max_dist,nSessions)),'jpg')
    
    fig_corr = figure;
    b_corr = bar(histo.corr_x,histo.corr,1);
    b_corr(1).FaceColor = 'green';
    b_corr(2).FaceColor = 'red';
    xlabel('spatial correlation')
    ylabel('number of pairs')
    legend('nearest neighbour','others')
    title(sprintf('max dist = %4.2f',parameter.max_dist))
    saveas(fig_corr,pathcat(pathMouse,sprintf('matchHistoCorr_d=%4.2f_nSes=%d.jpg',parameter.max_dist,nSessions)),'jpg')
    
    
    toc        
        
    
    
    
    
%      do this, according to Ziv, 2017:
%              build (bayesian) estimate to match cells
            
            %% if cells can change their shape - why do we merge all ROI-areas? shouldnt ROI areas then be dependent on session as well??? (otherwise recording background noise...)
            
%              delete unsure cell matches

    
end





function [prg_str1, prg_str2] = prepare_progress_report(dspl_string,num)%%% preparing progress output
    req_blanks1 = ceil(log10(num))+1;
    time_string = ', elapsed time (seconds): ';
    length_time = 6;
    req_blanks2 = length(time_string) + length_time;
    
    req_blanks = req_blanks1 + req_blanks2;
    fprintf(1,[dspl_string '(' num2str(num) ') ' blanks(req_blanks)]);
    backspace_string = '';
    for j = 1:req_blanks
        backspace_string = strcat(backspace_string,'\b');
    end
    prg_str1 = [backspace_string '%-' num2str(req_blanks1) 'd'];
    prg_str2 = [time_string '%' num2str(length_time) '.1f'];
end