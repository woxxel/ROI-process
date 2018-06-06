

function plot_matching(ROI_cluster)
  
  counts = [ROI_cluster.ct];
  scores = [ROI_cluster.score];
  
  nSes = max(counts);
  
  avg_score = zeros(nSes,1);
  
  for ct = 1:nSes
    avg_score(ct) = nanmean(scores(counts==ct));
  end
  
  figure('position',[100 100 1200 900])
%    subplot(2,4,[1:2])
%    plot(counts,'bx')
%    xlabel('cluster ID')
%    ylabel('ROI count')
  
%    subplot(2,4,3)
%    histogram(counts,'FaceColor','b')
%    xlabel('ROI count')
%    xlim([0,nSes])
%    ylim([0,500])
  
  subplot(2,4,4)
  hold on
  scatter(counts,scores,'kx')
  plot(1:nSes,avg_score,'r-','LineWidth',2)
  hold off
  xlabel('ROI count')
  ylabel('ROI score')
  
  subplot(2,4,[5:6])
  plot(scores,'rx')
  xlabel('cluster ID')
  ylabel('cluster score')
  
  subplot(2,4,7)
  histogram(scores,'FaceColor','r')
  xlabel('ROI count')
  xlim([0.5,1])
  
  plotPath = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/245/matching_results.png';
  print(plotPath,'-dpng','-r300')
  disp(sprintf('figure saved under %s',plotPath))
  
%    subplot(2,4,8)
  
end