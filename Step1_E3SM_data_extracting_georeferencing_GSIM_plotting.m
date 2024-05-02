%% Step 1: Setup path, casename, etc.
clc;clear;

% case name, usually embedded in the MOSART output file name
casename = 'extendedOutput.v3.LR.historical_0101';

% this is the path to the raw MOSART output directory
E3SMdir = ['P:\global\cfs\cdirs\e3sm\www\Tutorials\2024\simulations\' casename '\archive\rof\hist\'];

% path to directory for web publishing of the diagnostic visuals 
wwwdir = 'P:\global\cfs\cdirs\e3sm\www\tizhou\2024_Tutorial\';

% years for analysis
years = 2000:2014;

diagdir = [wwwdir casename '/'];

if ~exist(diagdir, 'dir')
    mkdir(diagdir);
end

numon = (years(end)-years(1)+1)*12;

E3SMflow = struct;

i = 1;
for year = years
    disp(year);
    for m = 1:12
        filename = [casename '.mosart.h0.' num2str(year) '-' sprintf('%02d',m) '.nc'];
        if i == 1
            mask   = ncread([E3SMdir filename],'mask');
            lat    = ncread([E3SMdir filename],'lat');
            lon    = ncread([E3SMdir filename],'lon');
            uparea = ncread([E3SMdir filename],'areatotal2');
        end
        wrmflow = ncread([E3SMdir filename],'RIVER_DISCHARGE_OVER_LAND_LIQ');
        E3SMflow.wrmflow(:,:,i) = wrmflow;
        i=i+1;
    end
end

E3SMflow.lat    = lat;
E3SMflow.lon    = lon;
E3SMflow.mask   = mask;
E3SMflow.uparea = uparea;

save([diagdir '/' casename '.mat'],'E3SMflow');
%% Step 2 Comparing simulated annual mean streamflow with GSIM streamflow

%%%%%%%%%%%% set path to the gauge metadata
%%%% metadata file of GSIM that has observed gauge lat lon and drainage area
%%%% This file includes 25765 gauges, which is a subset of the entire
%%%% dataset (30959 gauges). The removed gauges are associated with very
%%%% small drainage area (<1km2), which is not meaningful to be included.
gaugef = 'GSIM_catchment_characteristics_all_1km2.csv';
g = readtable(gaugef);

%%%%%%%%%%%% load the observed streamflow dataset (GSIM)
%the data has been reorganized to a 1380 * 30961 matrix. 1380 is the month
%number from 1901.1 to 2015.12. 30961 include two colums for year and month plus
%streamflow at 30959 gauge locations reported by GSIM
load ('GSIM.mat');

%%%%%%%%%%%% load E3SM simulated streamflow dataset
%the loaded "E3SMflow" struct includes 5 variables which are direct outputs of the
%river model (MOSART).
%wrmflow is a 3D matrix for the streamflow with dimension lon*lat*(month number)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LowRes
% load ([diagdir '/' casename '.mat']);
resolution = 0.5; % resolution of MOSART output
radius = 1; % search radius (number of grids around the center point)
error = 0.02; % drainage area bias allowed between sim and obs, 0.01 means 1%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat = E3SMflow.lat;
lon = E3SMflow.lon;
uparea = E3SMflow.uparea;

bins = floor(table2array(g(:,8:9))/resolution);
latlon = (bins+0.5)*resolution; %% move the observed lat long to grid center

%%%%%%%% define the export matrix
export = zeros(size(latlon,1),8); %annual mean of sim, annual mean of obs, error for area, lat, lon
kkkk = 0;

for i = 1:size(latlon,1)
    disp(i);
    lat_o = latlon(i,2);
    lon_o = latlon(i,1);
    area_o = table2array(g(i,14)); % % estimated drainage area (km2) from GSIM

    area_s = zeros(0,1);
    e_s = zeros(0,1);
    latlon_s = zeros(0,2);

    kk = 1;
    for ii = -radius:radius
        for jj = -radius:radius
            area_s(kk) = uparea((1+((lon_o+jj*resolution)-(-180+resolution/2))/resolution),(1+((lat_o+ii*resolution)-(-90+resolution/2))/resolution))/1000000;
            e_s(kk) = abs(area_s(kk)-area_o)/area_o;
            latlon_s(kk,1:2) = [lat_o+ii*resolution lon_o+jj*resolution];
            kk=kk+1;
        end
    end

    centerid = (length(e_s)-1)/2+1; %the id of the center grid in the searching area
    %%%%%%%%%%%%% searching the smallest error location within a box
    if e_s(centerid) < error
        ll = [lat_o lon_o];
        e = e_s(centerid);
        area_sim = area_s(centerid);
    else
        ll = latlon_s(e_s == min(e_s),:);
        area_sim = area_s(e_s == min(e_s));
        if size(ll,1)>1
            ll = ll(1,:);
            area_sim = area_sim(1);
        end
        e = min(e_s);
    end
    %%%%%%%%%%%%% use the center location
    %     ll = [lat_o lon_o];
    %     e = e_s(centerid);
    %%%%%%%%%%%%%%%%%%

    orignid = table2array(g(i,2));
    extract = GSIM(:,[1 2 orignid+2]); % extract the observed streamflow for this gauge

    monthly_mean = zeros(12,1)+nan;
    for m = 1:12
        if sum(extract(:,2) == m)>0
            monthly_mean(m) = mean(extract(extract(:,2)==m,3),'omitnan');
        end
    end

    annual_mean_obs = mean(monthly_mean); % this is observed annual mean streamflow
    %%%%%%%%%%%%%%%
    sim = squeeze(E3SMflow.wrmflow((1+(ll(2)-(-180+resolution/2))/resolution),(1+(ll(1)-(-90+resolution/2))/resolution),:));
    mmat = reshape(sim,12,[]);
    monthly_mean_sim = mean(mmat,2,'omitnan');
    annual_mean_sim = mean(monthly_mean_sim);

    %plot([monthly_mean,monthly_mean_sim]);


    export(i,1) = orignid; % original id of this gauge in GSIM (unique)
    export(i,2:3) = ll; % latlon of observation
    export(i,4) = area_o;
    export(i,5) = area_sim;
    export(i,6) = e * 100; % from fraction to percentage of the drainage area bias
    export(i,7) = annual_mean_obs;
    export(i,8) = annual_mean_sim;
end

%%% remove the gauges with nan flow
export(isnan(export(:,7)),:) = [];
export(isnan(export(:,8)),:) = [];


%% Step 3 Generate comparison plots for each selected gauge

mkdir([diagdir 'gauge_plots']);
cut = 10; % set the max area error (%) for all plots
sub = export;
sub(sub(:,6)>cut,:)=[];

sub(sub(:,7)<100,:)=[];
sub(sub(:,4)<2500*20,:)=[];

orignidlist = table2array(g(:,2));

% Extract river names and drainage obs from 'sub'
tempTable = table(sub(:,1), sub(:,4), 'VariableNames', {'Index', 'Drainage_obs'});
% Allocate space for the river names
tempTable.River = cell(size(sub, 1), 1);

% Loop through 'sub' to assign the correct river names from 'g'
for i = 1:size(sub, 1)
    % Find the corresponding row in 'g'
    rowIndex = find(orignidlist == sub(i, 1));
    
    % Assign the river name
    if ~isempty(rowIndex)
        riverName = g{rowIndex, 5};
        if iscell(riverName)
            tempTable.River{i} = riverName{1};  % Assuming river name is the first element in the cell
        else
            tempTable.River{i} = riverName;
        end
    end
end

% Find unique rivers
uniqueRivers = unique(tempTable.River);

% Initialize topTwoPerRiver with the same columns as tempTable but no rows
topTwoPerRiver = tempTable([], :);

for i = 1:length(uniqueRivers)
    riverData = tempTable(strcmp(tempTable.River, uniqueRivers{i}), :);
    [~, sortIdx] = sort(riverData.Drainage_obs, 'descend');
    % Selecting largest and smallest gauge for each river
    if length(sortIdx) > 1
        selectedIdx = [sortIdx(1); sortIdx(end)]; % Largest and smallest
    else
        selectedIdx = sortIdx; % Only one data point available
    end
    
    if ~isempty(selectedIdx)
        topTwoPerRiver = [topTwoPerRiver; riverData(selectedIdx, :)];
    end
end

% Use the indices to reconstruct the 'sub' matrix with only the relevant rows
finalSubIndices = topTwoPerRiver.Index;

subsub = zeros(length(finalSubIndices),8);
for i = 1:length(finalSubIndices)
 subsub(i,:) = sub(sub(:,1) == finalSubIndices(i),:);
end

sub = subsub;

[subnum,b] = size(sub);
rivers = cell(subnum, 1);
stations = cell(subnum, 1);

for i = 1:subnum
    extract_obs = GSIM(:,[1 2 sub(i,1)+2]); % extract the observed streamflow for this gauge
    extract_sim = squeeze(E3SMflow.wrmflow((1+(sub(i,3)-(-180+resolution/2))/resolution),(1+(sub(i,2)-(-90+resolution/2))/resolution),:));
    River = table2cell(g(orignidlist==sub(i,1),5)); rivers{i} = River;
    Station = table2cell(g(orignidlist==sub(i,1),6)); stations{i} = Station;

    monthly_mean_obs = zeros(12,1)+nan;
    for m = 1:12
        if sum(extract(:,2) == m)>0
            monthly_mean_obs(m) = mean(extract_obs(extract_obs(:,2)==m,3),'omitnan');
        end
    end

    monthly_mean_sim = mean(reshape(extract_sim,12,[]),2,'omitnan');

    R2 = 1 - sum((monthly_mean_obs - monthly_mean_sim).^2) / sum((monthly_mean_obs - mean(monthly_mean_obs)).^2);
    RMSE = sqrt(mean((monthly_mean_sim - monthly_mean_obs).^2));
    NSE = 1 - sum((monthly_mean_sim - monthly_mean_obs).^2) / sum((monthly_mean_obs - mean(monthly_mean_obs)).^2);
    PBIAS = (sum(monthly_mean_sim - monthly_mean_obs) / sum(monthly_mean_obs)) * 100;
    RMAE = mean(abs(monthly_mean_sim - monthly_mean_obs))/mean(monthly_mean_obs) * 100;
    MAPE = mean(abs((monthly_mean_obs - monthly_mean_sim) ./ monthly_mean_obs)) * 100;

    sub(i,9) = R2;
    sub(i,10) = RMSE;
    sub(i,11) = NSE;
    sub(i,12) = PBIAS;
    sub(i,13) = MAPE;

    fg = figure('Visible', 'off');
    set(fg, 'Position', [50 80 300 300]);
    plot(1:12,monthly_mean_sim,'linewidth',2)
    hold on
    plot(1:12,monthly_mean_obs,'linewidth',2)
    hold off
    metricsStr = sprintf('R²: %.2f\nRMSE: %.2f\nNSE: %.2f\nPBIAS: %.2f%%\nMAPE: %.2f%%', R2, RMSE, NSE, PBIAS, MAPE);

    % Add a text box with the metrics
    annotation(fg, 'textbox', [0.05, 0.05, 0.25, 0.25], 'String', metricsStr, 'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', 8);

    axis tight
    legend({'E3SM','GSIM'},'Location','southoutside','NumColumns',2)
    xlabel('Month')
    ylabel('Discharge (m^3/s)')
    titleStr = strjoin([string(Station), ' Station at ', string(River), ' River']);
    title(titleStr, 'FontSize', 8);
    set(gcf,'PaperPositionMode','auto')
    print(fg,'-dpng', [diagdir 'gauge_plots\' num2str(sub(i,1)) '.png'],'-r200');
    close (fg)
end

headers = {'ID','Lat (model)','Lon (model)', 'Drainage_obs','Drainage_model','A_error (%)','Annual meanQ obs','Annual meanQ sim', 'R2', 'RMSE', 'NSE', 'PBIAS', 'MAPE'};
matrixTable = array2table(sub, 'VariableNames', headers);

% Add the River and Station names to the table
matrixTable.River = rivers;
matrixTable.Station = stations;

% Write the table to a CSV file
writetable(matrixTable, [diagdir 'gauge_data.csv']);
