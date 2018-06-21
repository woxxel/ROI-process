

function plot_cluster_ROIs(PC_fields,clusters,bh,para,c)
  
  nC = size(clusters,1);
  nSes = size(clusters,2);
  
  if nargin < 5
    c = randi(nC)
  end
  
  
%    first = true;
%    
%    figure('position',[100 100 600 1200])
%  %    subplot(1,2,1)
%    hold on
%    for s = 1:nSes
%      
%      n = clusters(c,s).ROI_ID;
%      
%      if ~isnan(n)
%        if first
%          centr = round(clusters(c,s).centroid)
%          x_lims = [max(1,centr(2)-15),min(512,centr(2)+15)];
%          y_lims = [max(1,centr(1)-15),min(512,centr(1)+15)];
%          
%          [X,Y] = meshgrid(1:diff(x_lims)+1,1:diff(y_lims)+1);
%          first = false;
%        end
%        A_tmp = reshape(clusters(c,s).A,512,512);
%        A_tmp = A_tmp(y_lims(1):y_lims(2),x_lims(1):x_lims(2));
%        A_tmp(A_tmp==0) = NaN;
%        surf(X,Y,-2*A_tmp+s)
%  %        contour(A_tmp,[0.1,0.1]*max(clusters(c,s).A),'k')
%        
%      end
%    end
%    hold off
%    view([14,11])
%    rotate3d on
%    zlim([0,nSes+1])
%    set(gca,'ZDir','reverse')
  
  
  
  
  
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
    
    subplot(5,3,s);
    tight_axes(gca)
    
    plot_placefield(gca,PC_fields,c,s,clusters,bh,para)
    
  end
  
end