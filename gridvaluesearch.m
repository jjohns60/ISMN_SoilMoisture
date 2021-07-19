%% INPUTS
% data_hires: high resolution data input, of which values are extracted from
% lat_hires: higher resolution latitude values corresponding to input data
% lon_hires: higher resolution longitude values corresponding to input data
% lat_lowres: center grid cell latitudes to fit data to
% lon_lowres: center grid cell longitudes to fit data to
% method: 'list' returns cell array with lists containing all data values falling within each cell
%         'average' returns array with single average value
%         'nearest' returns the pixel value closest to the center point of the new grid
%         'linear' performs linear interpolation to regrid data inputs, requires consistent grid spacing
%% OUTPUTS
% DATA: a cell array with a list of all values that fell within the indicated grid
% lat_lowres: the original lower resolution grid values corresponding with the grid cell center latitudes
% lon_lowres: the original lower resolution grid values corresponding with the grid cell center longitudes

%% NOTE: Current bug, function assumes lat/lon grids are input as latitude changing by row and longitude by column
function [DATA,lat_lowres,lon_lowres] = gridvaluesearch(data_hires,lat_hires,lon_hires,lat_lowres,lon_lowres,method)
    
    if nargin == 5
        method = 'list';
    end
    
    %data_hires = val_rows;
    %lat_hires = val_lats;
    %lon_hires = val_lons;
    %lat_lowres = LAT;
    %lon_lowres = LON;
    
    %method is effective for data that is not necessarily higher resolution
    %potential for mismatch of data to true coordinates if irregular grid
    %method assumes data is input as such longitude is horizontal, latitude vertical
    if strcmp(method,'linear')
           
        [X,Y] = meshgrid(lon_hires(1,:),lat_hires(:,1));
        V = data_hires;
        
        [Xq,Yq] = meshgrid(unique(lon_lowres),flipud(unique(lat_lowres)));
        DATA = interp2(X,Y,V,Xq,Yq);
        
        %DATA = flipud(DATA);        
        return 
    end
    
    %method is effective for data that is not necessarily higher resolution
    if strcmp(method,'nearest')
        
        %hires data coordinates to search
        X = [lat_hires(:) lon_hires(:)];
        
        %query points for the new grid spacing
        Y = [lat_lowres(:) lon_lowres(:)];
        
        %identify nearest points
        idx = knnsearch(X,Y);
        
        %create array including indices of higher res data
        idx = reshape(idx,size(lat_lowres));
        idx = num2cell(idx);
        
        %apply function to all cells to extract corresponding datapoint
        wrapper = @(x) gridindex(x,data_hires);
        DATA = cellfun(wrapper,idx);
        
        return
    end
        
    %identify lower resolution cells that have data present
    Xedges = lon_lowres(~isnan(lon_lowres));
    Xedges = unique(Xedges);
    Xedges_i = zeros(length(Xedges)+1,1);
    Xedges_i(1) = Xedges(1) + ((Xedges(1) - Xedges(2))/2);
    Xedges_i(end) = Xedges(end) + ((Xedges(end) - Xedges(end-1))/2);
    Xedges_i(2:end-1) = .5*(Xedges(1:end-1) + Xedges(2:end));
    
    Yedges = lat_lowres(~isnan(lat_lowres));
    Yedges = unique(Yedges);
    Yedges_i = zeros(length(Yedges)+1,1);
    Yedges_i(1) = Yedges(1) + ((Yedges(1) - Yedges(2))/2);
    Yedges_i(end) = Yedges(end) + ((Yedges(end) - Yedges(end-1))/2);
    Yedges_i(2:end-1) = .5*(Yedges(1:end-1) + Yedges(2:end));
    
    %identify all cells that include points
    %h = histogram2();
    [~,~,~,Xidx,Yidx] = histcounts2(lon_hires,lat_hires,Xedges_i,Yedges_i);
    
    
    %pre-allocate data array
    if strcmp(method,'list')
        DATA = cell(size(lon_lowres));
    elseif strcmp(method,'average')
        DATA = zeros(size(lon_lowres));
        DATA_cnt = zeros(size(lon_lowres));
    end
    
    %prep data inputs to column vectors
    Xidx = repmat(Xidx(:),size(data_hires,3),1);
    Yidx = repmat(Yidx(:),size(data_hires,3),1);
    Aidx = (1:(length(Xidx)))'; %all positions in data input as col vector
    data_hires = data_hires(:); %convert to column vector
    
    %trim to locations with data
    idx = (Xidx > 0);
    Xidx = Xidx(idx);
    Yidx = Yidx(idx);
    Aidx = Aidx(idx);
        
    %use the identified index of cells with data to extract all of the higher resolution data in the specified cells
    %append these data lists to a cell array (DATA)
    if strcmp(method,'list')
        
        for i = 1:(length(Xidx))
            DATA{Yidx(i),Xidx(i)} = [DATA{Yidx(i),Xidx(i)} data_hires(Aidx(i))];
        end
        
        %correct data orientation
        DATA = flipud(DATA);
        
    elseif strcmp(method,'average')
        
        for i = 1:(length(Xidx))
            %get data value at each index
            d = data_hires(Aidx(i));
            if ~isnan(d)
                DATA(Yidx(i),Xidx(i)) = DATA(Yidx(i),Xidx(i)) + d; %sum
                DATA_cnt(Yidx(i),Xidx(i)) = DATA_cnt(Yidx(i),Xidx(i)) + 1; %total number of values at each index
            end
        end
        
        %ensure all locations with no values are not averaged (returns NaN instead)
        DATA_cnt(DATA_cnt == 0) = NaN;
        DATA = DATA./DATA_cnt;
        
        %correct data orientation
        DATA = flipud(DATA);
        
    end  
end

