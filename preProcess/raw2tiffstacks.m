
function pathH5 = raw2tiffstacks(pathRaw)
  
  tic
%    pathConf = '/home/wollex/Data/Documents/Uni/2016-XXXX_PhD/Japan/Work/Data/512x512x8989.conf';
  pathConf = '/home/mizuta/AlexCode/512x512x8989.conf';
%    pathRaw = dir(pathcat(pathIn,'*.raw'));
%    pathRaw = pathcat(pathIn,pathRaw(1).name);
  [pathFolder,fileName,ext] = fileparts(pathRaw);
%    pathRaw = pathcat(pathFolder,[fileName,'.raw']);
  pathH5 = pathcat(pathFolder,[fileName,'.h5']);
  
  system(sprintf('h5import %s -c %s -outfile %s',pathRaw,pathConf,pathH5))
  
  h52tiffstacks(pathFolder,2000)
  toc
end


function h52tiffstacks(pathFolder,nTiff)

  %%% path          - path to file directory
  %%% pathStacks    - path to created stacks
  %%% nTiff         - number of Tiffs per stack
  
  % ---------------------------------------------------------------------------------------------------

  %  path
  %  pathStacks
  %  if exist(pathStacks,'dir')
  %    rmdir(pathStacks,'s')
  %  end

  tiffs = struct;
  fileNames = dir(pathcat(pathFolder,'*.h5'));

  file = pathcat(pathFolder,fileNames(1).name);
  disp(sprintf('reading data from %s',file))

  I = h5read(file,'/DATA');

  height = size(I,1);
  width = size(I,2);
  t = size(I,3);
  nStacks = ceil(t/nTiff);

  %  img = zeros(height,width,nTiff,'uint16');

  %  c=1;
  for i = 1:nStacks
    
    idx_first = (i-1)*nTiff+1;
    idx_last = min(t,i*nTiff);
    
    img = I(:,:,idx_first:idx_last);
  %        tiffld = Tiff(pathcat(path,file),'r');
  %      img(:,:,c) = tiffld.read;
      
  %      if mod(c,nTiff)==0 || i==length(fileNames)
    svFile = sprintf('%s/stack%02d.tif',pathFolder,i);
    disp(sprintf('tiff#=%d, save stuff to %s',i,svFile))
    saveastiff(img,svFile);
  %        c=0;
  %      end
  %      c=c+1;
  %    end
  end
end