
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
tiffs.stacksize = length(tiffs.InfoImage);
%  tiffs.tifflib = Tiff(tiffs.file_name, 'r');
width = tiffs.InfoImage.Width;
height = tiffs.InfoImage.Height;
%  bits = tiffs.BitDepth;

img = zeros(height,width,nTiff,'uint16');

c=1;
for i = 1:length(fileNames)
  file = fileNames(i).name;
  
%    while
    tiffld = Tiff(pathcat(path,file),'r');
    img(:,:,c) = tiffld.read;
    
    if mod(c,nTiff)==0 || i==length(fileNames)
      svFile = sprintf('%s/stack%02d.tif',pathStacks,floor(i/c));
      disp(sprintf('tiff#=%d, save stuff to %s',i,svFile))
      saveastiff(img(:,:,1:c),svFile);
      c=0;
    end
    c=c+1;
%    end
end