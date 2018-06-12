
function create_tiff_stacks(path,pathStacks,nTiff)

  %%% path          - path to file directory
  %%% pathStacks    - path to created stacks
  %%% nTiff         - number of Tiffs per stack
  
% ---------------------------------------------------------------------------------------------------

  %  path
  %  pathStacks
  if exist(pathStacks,'dir')
    rmdir(pathStacks,'s')
  end


  tiffs = struct;
  fileNames = dir(pathcat(path,'*.tif'));
  file = pathcat(path,fileNames(1).name);

  tiffs.InfoImage = imfinfo(file);
  width = tiffs.InfoImage.Width;
  height = tiffs.InfoImage.Height;

  if length(fileNames) == 1   %% only one tiff stack present -> burst
    nframes = 8989;
    [img, nframes] = imread_big(file,nframes);
    
    nStacks = ceil(nframes/nTiff);
    
    [~, stackName, ~] = fileparts(file);
    img(1:10,1:10,1)
    for n = 1:nStacks
      idx_start = (n-1)*nTiff+1;
      idx_end = min(nframes,n*nTiff);
      
      svFile = pathcat(pathStacks,sprintf('%s_%02d.tif',stackName,n));
      disp(sprintf('saving tiff-stack #%d to %s',n,svFile))
      saveastiff(img(:,:,idx_start:idx_end),svFile);
      
    end
    
  else                        %% 
    
    if length(tiffs.InfoImage) > 1    %% check if already in stacks
      disp('tiffs are already in stacks - still want to proceed?')
    else
      
      [~, stackName, ~] = fileparts(file);
      stackName = stackName(1:end-5);
      img = zeros(height,width,nTiff,'uint16');
      
      c=1;
      for i = 1:length(nframes)
          
          file = fileNames(i).name;
          tiffld = Tiff(pathcat(path,file),'r');
          img(:,:,c) = tiffld.read;
          
          if mod(c,nTiff)==0 || i==length(fileNames)
            svFile = pathcat(pathStacks,sprintf('%s_%02d.tif',stackName,n));
            disp(sprintf('tiff-stack #%d saved to %s',n,svFile))
            saveastiff(img(:,:,1:c),svFile);
            c=0;
          end
          c=c+1;
          
      end
    end
  end
end