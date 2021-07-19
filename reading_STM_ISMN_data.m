%% Input: set this equal to 1 if you are using a Windows machine, 0 for Mac
Windows = 0;

%% Input: update to your 'SCAN' path
folders_path ='';

%defines the timesteps to check if they exist within each file
time_window = [datetime(2015,4,1,0,0,0) datetime(2015,12,31,23,59,59)];

%identify all valid folders within the folder path
folders = dir([folders_path]);
idx = [folders.isdir];
folders = folders(idx);
folders = folders(3:end);

%changes slash type if windows or mac
if Windows == 0
    slash = '/';
elseif Windows == 1
    slash = '\';
end

%will return an array at the end that stores info from all the stations
%that have data within your time window
FID = fopen('valid_stations.csv','w');
fprintf(FID,['FileID,','n_obs,','NetworkID,','NetworkID2,','SiteID2,','latitude,','longitude,','elevation_m,','depth_m,','depth_m2,','SensorID,', '\n']);

%loop through each data folder
for i = 1:length(folders)
    
    %identify the folder name
    folder = folders(i).name;
    
    %get list of all files in the folder including soil moisture
    files = dir([folders_path folder slash '*_sm_*' ]);
    
    %loop through files in folder
    for ii = 1:length(files)
        
        %get file name
        file = files(ii).name;
        
        %check if name contains the target variable
        if contains(file,'sm_0.050800') || contains(file,'sm_0.101600') || contains(file,'sm_1.016000')
        %if contains(file,'sm_0.101600')
            
            %show the file name during processing
            disp(file)
            
            %read in file data as text
            T = fileread([folders_path folder slash file]);
            
            %replace all consecutive spaces/tabs with single ','
            T = regexprep(T,'\t',' ');
            T = regexprep(T,' +',',');
    
            %splits the file by lines
            T_lines = splitlines(T);
            
            %gets the first line of data (header)
            station_info = T_lines(1);
            
            %separates the data from the header
            data = T_lines(3:end);
            
            %recombine data into single array
            data = cellfun(@(x) strsplit(x,','),data,'UniformOutput',false);
            
            %determine maximum table width
            M = max(cellfun(@length,data));
            
            %preallocate storage for each variable           
            dt = zeros([length(data),1]) + NaT; %store datetimes
            sm = zeros([length(data),1]) + NaN; %store volumetric SM
            f1 = cell([length(data),1]); %store flag 1
            f2 = cell([length(data),1]); %store flag 2
            for iii = 1:length(data)
                data_i = data{iii};
                if length(data_i) >= 5
                    try
                        dt(iii,1) = datetime([data_i{1} ' ' data_i{2}],'InputFormat','yyyy/MM/dd HH:mm');
                    catch
                        dt(iii,1) = datetime([data_i{1} ' ' data_i{2}],'InputFormat','MM/dd/yyyy HH:mm');
                    end
                    sm(iii,1) = str2double(data_i{3});
                    f1(iii,1) = data_i(4);
                    f2(iii,1) = data_i(5);
                end
            end
            
            %combine variables into table
            T = table(dt,sm,f1,f2,'VariableNames',{'timestamp','SM_vv','Flag1','Flag2'});
            
            %trim table to your time period
            idx = (T.timestamp >= time_window(1) & T.timestamp <= time_window(2));
            T = T(idx,:);
            
            %checks if there is valid data within time period
            if sum(idx) > 0
                %append to list of valid stations
                fprintf(FID,[file,',',num2str(sum(idx)),',',station_info{1} '\n']);
                
                %save modified file as .csv
                writetable(T,[folders_path folder slash file(1:end-4) '.csv']);
                
            end
        end
    end
end

%saves/closes list of stations with valid data
fclose(FID);