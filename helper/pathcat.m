
function [pathout] = pathcat(varargin)
    if nargin < 2
        error('Too few input arguments. Specify at least two to concatenate')
    end
    
    if ispc     % returns "true", if MATLAB version is for windows, else "false"
        path_sep = '\';
    else
        path_sep = '/';
    end
    
    pathout = varargin{1};
    for i = 2:nargin
        if strcmp(pathout(end),path_sep)
            pathout = [pathout,varargin{i}];
        else
            pathout = strcat(pathout,path_sep,varargin{i});
        end
    end
end