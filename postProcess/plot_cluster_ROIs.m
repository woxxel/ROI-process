

function plot_cluster_ROIs(PC_fields,clusters,bh,c)
  
  nC = size(clusters,1);
  nSes = size(clusters,2);
  
  if nargin < 4
    c = randi(nC)
  end
  
  
  first = true;
  
  figure('position',[100 100 600 1200])
%    subplot(1,2,1)
  hold on
  for s = 1:nSes
    
    n = clusters(c,s).ROI_ID;
    
    if ~isnan(n)
      if first
        centr = round(clusters(c,s).centroid)
        x_lims = [max(1,centr(2)-15),min(512,centr(2)+15)];
        y_lims = [max(1,centr(1)-15),min(512,centr(1)+15)];
        
        [X,Y] = meshgrid(1:diff(x_lims)+1,1:diff(y_lims)+1);
        first = false;
      end
      A_tmp = reshape(clusters(c,s).A,512,512);
      A_tmp = A_tmp(y_lims(1):y_lims(2),x_lims(1):x_lims(2));
      A_tmp(A_tmp==0) = NaN;
      surf(X,Y,-2*A_tmp+s)
%        contour(A_tmp,[0.1,0.1]*max(clusters(c,s).A),'k')
      
    end
  end
  hold off
  view([14,11])
  rotate3d on
  zlim([0,nSes+1])
  set(gca,'ZDir','reverse')
%    xlim([centr(2)-15,centr(2)+15])
%    ylim([centr(1)-15,centr(1)+15])
  
  
%    sigma = 1.5;
%    sz = 5;    % length of gaussFilter vector
%    x = linspace(-sz / 2, sz / 2, sz);
%    gaussFilter = exp(-x .^ 2 / (2 * sigma ^ 2));
%    gaussFilter = gaussFilter / sum (gaussFilter); % normalize
%    
%    for s = 1:nSes
%      n = clusters(c,s).ROI_ID;
%      if ~isnan(n)
%        subplot(nSes,2,2*s)
%        hold on
%        bar(PC_fields(c).firingmap(s,:),'FaceColor','r')
%        normfields_filt = conv(PC_fields(c).firingmap(s,:),gaussFilter,'same');
%        plot(normfields_filt,'b-')
%        hold off
%      end
%    end
  
  
%    normfields = PC_fields(c).firingmap./max(PC_fields(c).firingmap,[],2);
%    filt = ones(1,5);
  
%    figure
%    plot(gaussFilter)
  
  
%    figure
%    subplot(1,4,1:3)
%    imagesc(normfields)
%    colormap('jet')
%    colorbar
%    
%    subplot(1,4,4)
%    hold on
%    plot(PC_fields(c).status,linspace(1,15,15),'kx')
%    plot(PC_fields(c).MI,linspace(1,15,15),'rx')
%    hold off
%    ylim([0.5,15.5])
%    set(gca,'YDir','reverse')
  
%    pause(1)
  figure('position',[100 100 2000 1200])
  for s=1:nSes
    
    nsd = 4;
    prc = 20;
    subplot(5,3,s);
    tight_axes(gca)
    
    hold on
    barh(-PC_fields(c).firingmap(s,:)*900/max(PC_fields(c).firingmap(s,:)),'r')
    plot(bh(s).location/20)
    
    modeS = prctile(clusters(c,s).S(clusters(c,s).S>0),prc);                    %% get mode from overall activity
    activity = floor(sqrt(clusters(c,s).S/(modeS*nsd)));         %% only activity from actual times
    if length(activity) > 0
      idx = activity & bh(s).longrunperiod;
      activity_lr = zeros(1,8989);
      activity_lr(idx) = activity(idx);
      scatter(find(activity_lr>0),bh(s).location(activity_lr>0)/20,2*activity_lr(activity_lr>0)+5,'r','fill')
      
      idx = activity & ~bh(s).longrunperiod;
      activity_nlr = zeros(1,8989);
      activity_nlr(idx) = activity(idx);
      scatter(find(activity_nlr>0),bh(s).location(activity_nlr>0)/20,2*activity_nlr(activity_nlr>0)+5,'k','fill')
    end
    
    hold off
    set(gca,'YAxisLocation','right')
    ylim([0 80])
    xlim([-1000,8989])
    
    if PC_fields(c).status(s)
      title(sprintf('PC, MI = %4.2g (%4.2g)',PC_fields(c).MI(s),PC_fields(c).MI_frac(s)))
    else
      title(sprintf('nPC, MI = %4.2g (%4.2g)',PC_fields(c).MI(s),PC_fields(c).MI_frac(s)))
    end
    
  end
  
end