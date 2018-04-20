%%% truncates "ROIs" at regions where they overlap "borders" and calculates their norm

function [ROI] = ROI_update(ROI,varargin)
    
    for i = 1:nargin-1
        switch varargin{i}
            case 'truncate'
                ROI = truncate(ROI);
            case 'norm'
                ROI = get_norm(ROI);
            case 'centroid'
                ROI = get_centroid(ROI);
            case 'normalize'
                ROI = normalize(ROI);
            otherwise
                disp('option is not specified - please check and fix that!')
        end
    end
end

function [ROI] = truncate(ROI)

    borders = [1, 512];
    
    %% adjust ROI size, if shift crosses 0-pixel
    border_entry = find(ROI.extent < borders(1));
    if length(border_entry) > 0
        %% cutting filter accordingly
        for i = 1:length(border_entry)
            idx = border_entry(i);
            offset = borders(1) - ROI.extent(idx);
            ROI.extent(idx) = 1;
            if idx <= 2 % y position
                ROI.filter = ROI.filter(1+offset:end,:);
            else        % x position
                ROI.filter = ROI.filter(:,1+offset:end);
            end
        end
    end
    
    %% adjust ROI size, if shift crosses 512-pixel (end)           
    border_entry = find(ROI.extent > borders(2));
    if length(border_entry) > 0
        %% cutting filter accordingly
        for i = 1:length(border_entry)
            idx = border_entry(i);
            offset = ROI.extent(idx) - borders(2);
            ROI.extent(idx) = borders(2);
            if idx <= 2 % y position
                ROI.filter = ROI.filter(1:end-offset,:);
            else        % x position
                ROI.filter = ROI.filter(:,1:end-offset);
            end
        end
    end
    
    ROI = normalize(ROI);
    ROI = get_centroid(ROI);
    ROI = get_norm(ROI);
end

function [ROI] = normalize(ROI)
    ROI.filter = ROI.filter/sum(ROI.filter(:));
end


function [ROI] = get_centroid(ROI)
    
    range_x = ROI.extent(1,2):ROI.extent(2,2);
%      cent_x = sum(logical(ROI.filter())*range_x')/ROI.area;
    cent_x = sum(ROI.filter()*range_x');
    
    range_y = ROI.extent(1,1):ROI.extent(2,1);
%      cent_y = sum(range_y*logical(ROI.filter()))/ROI.area;
    cent_y = sum(range_y*ROI.filter());
    
    ROI.centroid = [cent_y, cent_x];
end


function [ROI] = get_norm(ROI)
    ROI.norm = norm(ROI.filter(:));
end