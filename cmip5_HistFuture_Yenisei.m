%% cmip5_HistFuture_Yenisei.m
% Summary: Watershed basin and lat-lon box analysis of CMIP5 data
% Specify historical and future RCP to create continuous extraction
% This program will scrape user defined variables from the CMIP5 archive on
% ingrid. It uses opendap.
% Inputs:  basin shapefile, CMIP5 data from LDEO server 
% Outputs: CMIP5 data for basin area of interest; Pulls down continuous 
% historical to future (RCP) runs in one go.

% Uses: 1) "m_map" scripts from https://www.eoas.ubc.ca/~rich/map.html
%       2) Basin shapefile from WRI 
%       3) Other shapefiles from http://www.naturalearthdata.com/
% Author: Michael J. Puma, 2017

%CHANGE Directory Stuff
%BIG CHANGE rewrote with mapping toolbox

%% Start clean
clear variables; close all; clc;

%% Directories
root_dir = '/Users/fullerjohnjef/climate/';

%% Basin specific information
basinID = 3; basin_name='Yenisey'; area_txt='Yenisey Basin';
file_name = 'wri_basins/wribasin.shp';

%% Lat-Lon limits for plotting
lonmin = 65; lonmax = 125; latmin = 40; latmax = 80;
lonlim=[lonmin lonmax];
latlim=[latmin latmax];

flag_basin_surrounding=false; % false=>use shpfile to define basin mask, true=> use surrounding lat/lon box

%% Folders
folder_feature = [root_dir 'MapFeatures/'];
folder_output = [root_dir 'CMIP5_output_continuous/'];
output_folder = [folder_output basin_name '/'];

%% Physical and cultural features for plots (download files from http://www.naturalearthdata.com/)
%using toolbox .shp for land; other one had inconsistent polygon topology
%(funky biz at poles)
country = shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp'], ...
    'UseGeoCoords', true);
provinces = shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces_lines_shp.shp'], ...
    'UseGeoCoords', true);
rivers_lakes = shaperead([folder_feature '10m_physical/ne_10m_rivers_lake_centerlines.shp'], ...
    'UseGeoCoords', true);
ocean = shaperead([folder_feature '10m_physical/ne_10m_ocean.shp'], ...
    'UseGeoCoords', true);
land = shaperead('landareas', 'UseGeoCoords', true);
coastline = shaperead([folder_feature '10m_physical/ne_10m_coastline.shp'], ...
    'UseGeoCoords', true);


%% Load shapefile for watershed
basin_shp_all = shaperead([folder_feature  file_name]);
basin_shp = basin_shp_all(basinID);

%% Map Proj
map_proj_name='Equidistant Cylindrical';
%% Plot river basin 
% figure;
% hold on;
% ax = worldmap(latlim, lonlim); shading flat;
% setm(ax, 'Grid', 'off', 'FontSize', 16, 'FontWeight', 'bold');
% geoshow([country.Lat],[country.Lon],'Color','black');
% patchm([basin_shp.Lat],[basin_shp.Lon],'g');
% goeshow([rivers_lakes.Lat],[rivers_lakes.Lon]);
% title([area_txt],'FontSize',16);
% hold off

%% User defined variables
% Scenario Name Options
% grid_yrs define the time periods to retain the spatial fields.
scen_name1='historical'; scen_name2='rcp85';
%scen_name1='historical'; scen_name2='rcp45';

window_ave=[1980 2099];

%% Variables
dataset_name='atmos'; var_name='tas'; unit_name='K'; scale_units=1;
% dataset_name='atmos'; var_name='tasmin'; unit_name='K'; scale_units=1;
%dataset_name='atmos'; var_name='tasmax'; unit_name='K'; scale_units=1;
%dataset_name='atmos'; var_name='pr'; unit_name='kg m^-^2 s^-^1'; scale_units=1;
%dataset_name='atmos'; var_name='pr'; unit_name='mm day^-^1'; scale_units=86400;
%dataset_name='atmos'; var_name='evspsbl'; unit_name='mm day^-^1'; scale_units=86400;
%dataset_name='atmos'; var_name='huss'; unit_name=''; scale_units=1; % Near Surface Specific Humidity
%dataset_name='atmos'; var_name='rlds'; unit_name='W m^-^2'; scale_units=1; % Surface Downwelling Longwave Radiation
%dataset_name='atmos'; var_name='rsds'; unit_name='W m^-^2'; scale_units=1; % Surface Downwelling Shortwave Radiation
%dataset_name='atmos'; var_name='rsus'; unit_name='W m^-^2'; scale_units=1; % Surface Upwelling Shortwave Radiation
%dataset_name='atmos'; var_name='clt'; unit_name='percent'; scale_units=1; %Total Cloud Fraction
%dataset_name='land'; var_name='mrso'; unit_name='kg m^-^2'; scale_units=1; %Total Soil Moisture Content
%dataset_name='land'; var_name='mrro'; unit_name='kg m^-^2'; scale_units=1; %Total runoff
%dataset_name='land'; var_name='mrsos'; unit_name='kg m^-^2'; scale_units=1; %Total surface runoff
%dataset_name='land'; var_name='lai'; unit_name='unitless'; scale_units=1; % Leaf Area Index
%dataset_name='landIce'; var_name='snc'; unit_name='percent'; scale_units=1; % Surface Snow Area Fraction
%dataset_name='landIce'; var_name='snw'; unit_name='kg m^-^2'; scale_units=1; % Surface Snow Amount

%% Static Variables
root_name='http://strega.ldeo.columbia.edu:81/home/.haibo/.CMIP5';
%root_name='http://stregone.ldeo.columbia.edu:81/home/.haibo/.CMIP5'; % alternative, works for CCSM4

% 1st Column: Model /// 2nd Column: land mask
model_names={
%     'bcc-csm1-1','bcc-csm1-1';
%     'CanAM4','CanAM4';
%     'CanCM4','CanCM4';
%     'CanESM2','CanESM2';
%     'CCSM4','CCSM4';
%     'CNRM-CM5','CNRM-CM5';
%     'CSIRO-Mk3-6-0','CSIRO-Mk3-6-0';
%     'GFDL-CM3','GFDL-ESM2M';
%     'GFDL-ESM2M','GFDL-ESM2M';
%     'GFDL-ESM2G','GFDL-ESM2M';
%     'GFDL-HIRAM-C180','GFDL-HIRAM-C180';
    'GISS-E2-R','GISS-E2-R';
%     'GISS-E2-H','GISS-E2-R';
%     'HadGEM2-A','HadGEM2-ES';
%     'HadGEM2-CC','HadGEM2-ES';
%     %       'HadGEM2-ES','HadGEM2-ES';
%     'inmcm4','inmcm4';
%     'IPSL-CM5A-LR','IPSL-CM5A-LR';
%     'MIROC-ESM','MIROC-ESM';
%     'MIROC-ESM-CHEM','MIROC-ESM-CHEM';
%     'MIROC4h','MIROC4h';
%     'MIROC5','MIROC5';
%     'MPI-ESM-LR','MPI-ESM-LR';
%     'MPI-ESM-P','MPI-ESM-P';
%     'MRI-CGCM3','MRI-CGCM3';
%     'NorESM1-M','NorESM1-M';
    };

% Possible Ensemble Names
ensemb_names={
    'r1i1p1';
    'r2i1p1';
    'r3i1p1';
    'r4i1p1';
    'r5i1p1';
    'r6i1p1';
    'r7i1p1';
    'r8i1p1';
    'r9i1p1';
    'r10i1p1';
    'r1i1p2';
    'r2i1p2';
    'r3i1p2';
    'r4i1p2';
    'r5i1p2';
    'r1i1p3';
    'r2i1p3';
    'r3i1p3';
    'r4i1p3';
    'r5i1p3';
    'r1i1p124';
    };

%% Output Structure

% We loop through model/ensemble/year
for i_model=1:size(model_names,1);          % MODEL LOOP
    
    % flag to save output
    out_flag=0;
    
    ar5data.('scenario')=[scen_name1 ' to ' scen_name2];
    ar5data.('dataset')=dataset_name;
    ar5data.('variable')=var_name;
    ar5data.('units')=unit_name;
    ar5data.('basin')=basin_name;
    n=model_names(i_model,1);
    ar5data.('model')=n{1};
    
    %% Load Land Fraction and lat/lon data
    nc_landfrac=ncgeodataset([root_name '/.fixed/.atmos/.' model_names{i_model,2} '/.sftlf.nc/.sftlf/dods']);
    landfrac=double(nc_landfrac.data('sftlf'))./100;
    lat_data=double(nc_landfrac.data('lat'));
    lon_data=double(nc_landfrac.data('lon'));
    
    % if plotted naively will look like east-west hemispheres are swapped,
    % done for -180-180 formatting (climate data hemispheres swapped later
    % for compatibility)
    if max(lon_data(:))>180+5
        lon_data=lon_data-180;
    end
    
    %% Create Basin Mask
    
    % Create New Dummy Mask of ones
    basin_mask=nan(length(lat_data),length(lon_data));

    % Clear variables from previous model
    clear lonlat lonlat_loc lonlat_loc_ind lat_data_cut lon_data_cut
    clear lonlat_loc_ind lat_data_cut lon_data_cut

    % Place lat-lon in array for inpolygon
    jj = 1;
    for i=1:length(lon_data)
        for j=1:length(lat_data)
            lonlat(jj,1) = lon_data(i);
            lonlat(jj,2) = lat_data(j);
            jj = jj+1;
        end
    end

    % Find grid cells within basin polygon
    if flag_basin_surrounding
       [IN,ON] = inpolygon(lonlat(:,1),lonlat(:,2),lonlim,latlim);
    else
       [IN,ON] = inpolygon(lonlat(:,1),lonlat(:,2),basin_shp.X,basin_shp.Y);
    end
    lonlat_loc = IN+ON; clear IN ON;
        lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
        lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
        lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
        for i = 1:size(lonlat_loc_ind,1)
            basin_mask(lat_data == lat_data_cut(i),lon_data == lon_data_cut(i))=1;
        end
        
    % Figure: Basin check 1
%     figure; shading flat; pcolor(lon_data,lat_data,landfrac.*basin_mask),title('model land fraction, basin mask'),colorbar, grid off;
%     plot(g);
    % Figure: Basin check 2
%     figure; hold on; shading flat;
%     worldmap('World');
%     pcolorm(lat_data, lon_data, landfrac);
%     geoshow([country.Lat], [country.Lon], 'Color', 'black');
%     patchm([basin_shp.Y], [basin_shp.X], 'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 2);
%     plot(g); % fail to make it stop :)


    % save lon/lat, basin mask, and land fraction mask data
    ar5data.('lon')=lon_data;
    ar5data.('lat')=lat_data;
    ar5data.('basinmask')=basin_mask;
    ar5data.('landfrac')=landfrac;
    
    % Repeat the basin mask for 12 months. This will be used for the
    % filtering
    clear basin_mask_month
    for b_mon=1:12
        basin_mask_month(b_mon,:,:)=basin_mask;
    end
    
    %%  ENSEMBLE LOOP
    for i_ensemb=1:length(ensemb_names)    % ENSEMBLE LOOP
        
        % For the model names, I need to use the dash to grab the data.
        % However, if I want to use a model name as the name of a structure
        % field, I have to use an underscore
        modname=model_names{i_model,1};
        TF = strfind(modname,'-');
        modname(TF)='_';
        
        % The way I have this setup: will download all data (12 months) at a time.
        % Then I will geographically subet the data in MATLAB.
        
        % Base Name Strings for historical and future scenario, dataset, variable, model, ensemble
        nc_name_base1=[root_name '/.' scen_name1 '/.' dataset_name  '/.mon/.' var_name '/.' ...
            model_names{i_model,1} '/.' ensemb_names{i_ensemb} '/.' var_name];
        
        % Base Name String for historical and future scenario, dataset, variable, model, ensemble
        nc_name_base2=[root_name '/.' scen_name2 '/.' dataset_name  '/.mon/.' var_name '/.' ...
            model_names{i_model,1} '/.' ensemb_names{i_ensemb} '/.' var_name];
        
        % Check to see if these runs/ensembles exists
        [~,STATUS1] = urlread(nc_name_base1);
        [~,STATUS2] = urlread(nc_name_base2);
        
        % If they both exist does, open a data structure
        if (STATUS1+STATUS2)==2
            
            % flag to save output
            out_flag=1;
            
            % Print current model to screen
            [model_names{i_model} '---' ensemb_names{i_ensemb}]
            
            % Determine full temporal range, including start and end year for
            % this run.
            nc_time_only=ncgeodataset([nc_name_base1 '/dods']);
            ser_date=nj_time(nc_time_only,'T');
            [yr_vect1,mon_vect1,~,~,~,~] = datevec(ser_date);
            min_yr1=min(yr_vect1);
            
            nc_time_only=ncgeodataset([nc_name_base2 '/dods']);
            ser_date=nj_time(nc_time_only,'T');
            [yr_vect2,mon_vect2,D,H,MI,S] = datevec(ser_date);
            min_yr2=min(yr_vect2);
            
            % Now, create new time vector, with 3rd column flag indicating
            % which run that date is associated with
            time_vect_sim1=[yr_vect1 mon_vect1 ones(size(yr_vect1))];
            time_vect_sim2=[yr_vect2 mon_vect2 ones(size(yr_vect2)).*2];
            
            time_vect_all=[time_vect_sim1; time_vect_sim2];
            
            % year vector
            min_yr=min(time_vect_all(:,1));  % start year for this simulation
            max_yr=max(time_vect_all(:,1));  % end year for this simulation
            yrs=min_yr:max_yr;    % all years for this simulation
             
            % Pull out lon/lat data;
            nc_time_only=ncgeodataset([nc_name_base2 '/dods']);
            
            %% YEAR LOOP
            % initialize output matrix
            clim_all_months=[];
            
            % Find closest match to selected window
            yrs=yrs(find(yrs>=window_ave(1) & yrs<=window_ave(2)));
            
            % counter for number of years; used for saving data
            yr_step=0;
            for i_yr=yrs
                yr_step=yr_step+1;
                
                % Index locations for this year
                yr_locs=find(time_vect_all(:,1)==i_yr);
                
                % flags indicating which run(s) have this data
                mod_data=unique(time_vect_all(yr_locs,3));
                
                clear clim_data
                clear clim_data1
                clear clim_data2
                
                if sum(mod_data)==1
                    % Load data from HISTORICAL RUN, if available
                    nc_name_time=['/T/' num2str(i_yr) '/' num2str(min_yr1) '/sub/12/mul/dup/12/add/RANGEEDGES/dods'];
                    
                    % Open a netcdf data structure and pull out the climate data for all months in the current year;
                    nc_data=ncgeodataset([nc_name_base1 nc_name_time]);
                    clim_data_or=double(nc_data.data(var_name)).*scale_units;
                    clim_data(:,:,1:floor(size(clim_data_or,3)/2))=clim_data_or(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3));
                    clim_data(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3))=clim_data_or(:,:,1:floor(size(clim_data_or,3)/2));
                    clear clim_data_or;
                    
%                     figure; hold on; shading flat;
%                     worldmap('World');
%                     pcolorm(lat_data, lon_data, squeeze(mean(clim_data(1:12,:,:),1)));
%                     geoshow([country.Lat], [country.Lon], 'Color', 'black');
                    
                    % Mask out non-basin cells
                    clim_data=clim_data.*basin_mask_month;
%                     
%                     figure; hold on; shading flat;
%                     worldmap('World');
%                     pcolorm(lat_data, lon_data, squeeze(mean(clim_data(1:12,:,:),1)));
%                     geoshow([country.Lat], [country.Lon], 'Color', 'black');
%                     plot(g); % fail to make it stop :)
                   
                elseif sum(mod_data)==2
                    % Load data from RCP RUN, if available
                    nc_name_time=['/T/' num2str(i_yr) '/' num2str(min_yr2) '/sub/12/mul/dup/12/add/RANGEEDGES/dods'];
                    
                    % Open a netcdf data structure and pull out the climate data for all months in the current year;
                    nc_data=ncgeodataset([nc_name_base2 nc_name_time]);
                    clim_data_or=double(nc_data.data(var_name)).*scale_units;
                    clim_data(:,:,1:floor(size(clim_data_or,3)/2))=clim_data_or(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3));
                    clim_data(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3))=clim_data_or(:,:,1:floor(size(clim_data_or,3)/2));
                    clear clim_data_or;
                    
                    % Mask out non-basin cells
                    clim_data=clim_data.*basin_mask_month;
                    
                elseif sum(mod_data)==3
                    
                    % Load data from HISTORICAL RUN, if available
                    nc_name_time=['/T/' num2str(i_yr) '/' num2str(min_yr1) '/sub/12/mul/dup/12/add/RANGEEDGES/dods'];
                    
                    % Open a netcdf data structure and pull out the climate data for all months in the current year;
                    nc_data=ncgeodataset([nc_name_base1 nc_name_time]);
                    clim_data_or=double(nc_data.data(var_name)).*scale_units;
                    clim_data1(:,:,1:floor(size(clim_data_or,3)/2))=clim_data_or(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3));
                    clim_data1(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3))=clim_data_or(:,:,1:floor(size(clim_data_or,3)/2));
                    clear clim_data_or;
                    
                    % Load data from RCP RUN, if available
                    nc_name_time=['/T/' num2str(i_yr) '/' num2str(min_yr2) '/sub/12/mul/dup/12/add/RANGEEDGES/dods'];
                    
                    % Open a netcdf data structure and pull out the climate data for all months in the current year;
                    nc_data=ncgeodataset([nc_name_base2 nc_name_time]);
                    clim_data_or=double(nc_data.data(var_name)).*scale_units;
                    clim_data2(:,:,1:floor(size(clim_data_or,3)/2))=clim_data_or(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3));
                    clim_data2(:,:,floor(size(clim_data_or,3)/2)+1:size(clim_data_or,3))=clim_data_or(:,:,1:floor(size(clim_data_or,3)/2));
                    clear clim_data_or;
                    
                    clim_data=[clim_data1; clim_data2];
                    
                    % Sometimes same year is available from both runs if this
                    % is the case, just take the first year from the historical
                    % run.
                    % Mask out non-basin cells
                    if size(clim_data,1)~=12
                        clim_data=clim_data1.*basin_mask_month;
                    else
                        clim_data=clim_data.*basin_mask_month;
                    end
                    
                end
                
                %% Month Loop (NEW)
                % each month separately
                for i_mon=1:12
                    
                    % load current month
                    curr_month=double(squeeze(clim_data(i_mon,:,:)));
                    
                    % Find Latitude and Longitudes for region
                    [i_lat]=find(isnan(nanmean(curr_month,2))==0); lat_basin_out=lat_data(i_lat);
                    [i_lon]=find(isnan(nanmean(curr_month,1))==0); lon_basin_out=lon_data(i_lon);
                    
                    % Average over basin
                    [cosmean,cosgrid,datagrid,new_lon,new_lat] = coswt(curr_month,lat_data,lon_data,...
                        min(lon_basin_out),max(lon_basin_out),min(lat_basin_out),max(lat_basin_out));
                    clim_var_ave(yr_step,i_mon)=cosmean;
                    
                    % Save spatial data
                    %clim_var_spatial(yr_step,i_mon,:,:)=curr_month(i_lat,i_lon);
                    % CHANGE
                    clim_var_spatial(yr_step,i_mon,:,:)=curr_month;
                    
                end % month loop
                
                % Clear Previous Years Data
                clear clim_data
                clear clim_data1
                clear clim_data2
                
            end
            
            ar5data.('lon_basin')=lon_basin_out;
            ar5data.('lat_basin')=lat_basin_out;
            ar5data.(ensemb_names{i_ensemb}).('yrs')=yrs;
            ar5data.(ensemb_names{i_ensemb}).('clim_var_ave')=clim_var_ave; % average across the watershed
            ar5data.(ensemb_names{i_ensemb}).('clim_grid')=clim_var_spatial; % spatial data, limited to the watershed
            
            
        end % status boolean
        
        clear clim_var_ave
        clear clim_var_spatial
        
    end % ensemble loop
    
    % Save output for current model
    % not enough room for grids in dropbox.........
    if out_flag==1
        if flag_basin_surrounding
            save([output_folder   var_name  '.surrounding.'  modname '.' scen_name1 '-' scen_name2 '.ar5.monmean.mat'],'-struct','ar5data','-v7.3');
        else
            save([output_folder   var_name  '.'  modname '.' scen_name1 '-' scen_name2 '.ar5.monmean.mat'],'-struct','ar5data','-v7.3');
        end
    end
    
    clear ar5data
    
end % model loop
