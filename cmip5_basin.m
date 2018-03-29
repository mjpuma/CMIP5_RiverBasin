
%% cmip5_basin_Yenisei.m
% Summary: Watershed basin, crop, and country-level analysis of CMIP5 data
% Inputs:  CMIP5 data processed by script cmip5_HistFuture_FullGrid.m
%          CRU/GPCP/CMAP observations previous formatted by matlab scripts

% Uses: 1) "m_map" scripts from https://www.eoas.ubc.ca/~rich/map.html
%       2) Basin shapefile from WRI 
%       3) Other shapefiles from http://www.naturalearthdata.com/
% Author: Michael J. Puma, 2017

%CHANGE directory stuff
%CHANGE rewrote m_map code for matlab mapping toolbox
%CHANGE land shapefile source from Natural Earth to mapping toolbox, 
%previous one threw error, 'does not form clean polygon topology'
%CHANGE program runs through data twice, 1st time with mask for lat/lon
%limit, second for basin mask
%CHANGE ended up downloading data for world then masking from there,
%accordingly renamed script cmip5_basin
%CHANGE rewrote most of actual plotting into a seperate function plot_basin
%for easier reading
%CHANGE added process for making a scatterplot of elevation and given
%var
%CHANGE added process for making tables of r and p values of elevation and
%given var
%CHANGE messed around a bit with file naming process to account for
%increased amount of stuff, also changed storage of maps/plots so that each
%basin folder has subfolders for variables which contain the maps/plots
%CHANGE now using "Yenisei" instead of "Yenisey"


root_dir='C:/Users/jefuller/climate/';

model_r_values=nan(50,50,2);
model_p_values=nan(50,50,2);

num_sites=10;
for model_res=1:2
    if model_res==1
        %% coarse models>= 2 degrees
        model_names={
%             'bcc-csm1-1','bcc-csm1-1';
%             'CanAM4','CanAM4';
%             'CanCM4','CanCM4';
%             'CanESM2','CanESM2';
%             'GFDL-CM3','GFDL-ESM2M';
%             'GFDL-ESM2M','GFDL-ESM2M';
%             'GFDL-ESM2G','GFDL-ESM2M';
%             'GISS-E2-R','GISS-E2-R';
%             'GISS-E2-H','GISS-E2-R';
%             'IPSL-CM5A-LR','IPSL-CM5A-LR';
%             'MIROC-ESM','MIROC-ESM';
%             'MIROC-ESM-CHEM','MIROC-ESM-CHEM';
%             'NorESM1-M','NorESM1-M'
             };
        res='Lower Resolution'; res_file ='low';
    elseif model_res==2
        %% fine models <2 degrees
        model_names={
            'CCSM4','CCSM4';
%             'CNRM-CM5','CNRM-CM5';
%             'CSIRO-Mk3-6-0','CSIRO-Mk3-6-0';
%             'HadGEM2-A','HadGEM2-ES';
%             'HadGEM2-CC','HadGEM2-ES';
%             'inmcm4','inmcm4';
%             'MIROC4h','MIROC4h';
%             'MIROC5','MIROC5';
%             'MPI-ESM-P','MPI-ESM-P';
%             'MPI-ESM-LR','MPI-ESM-LR';
%             'MRI-CGCM3','MRI-CGCM3'
               };
        res='Higher Resolution'; res_file ='hi';
    end
        
    ensemble_names={
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

    %% Start clean
    clearvars -except i_site num_sites variable_flag model_res res res_file model_names combined_model_names ensemble_names lat_box_names lat_boxes root_dir model_r_values model_p_values;
    close all; clc;
    
    %% Variables
    basin_name='Yenisei';
%     basin_name='Lena';
%     basin_name='Amur';
%     basin_name='Ob';
%     basin_name='Indigirka';
    
    variable='pr';
%     variable='tas'; 
%     variable='tasmin';
%     variable='snw'; % Surface Snow Amount
%     variable='snc'; % Surface Snow Area Fraction
%     variable='mrro'; % Total Runoff
%     variable='mrros'; % Total Surface Runoff
%     variable='mrsos'; % Moistuire in Upper Portion of Soil Column
%     variable='huss'; % Near-Surface Specific Humidity
    
    flag_obs=true;
    flag_basin_map=false; % whether or not to create a plain map of the basin
    plot_p1=false;
    plot_p2=false;
    plot_diff=false;
    plot_elevation=false;
    plot_diff_scatter=false;
    lat_box_table=false;
    save_results=false;
        
    flag_world=false; % false=>plot basin/latlim, lonlimit, else world 
    flag_land_mask=false; % if surrounding basin area intersects with ocean, only use over-land data
    flag_no_outliers=false;
    
    custom_save_string='';
%     custom_save_string=['1-deg-lat-box.'];
%     custom_save_string='no-outliers.';

    
    custom_var_lims=true; % defaults to range between rounded min/max values

    if custom_var_lims
        bas_var_p1_lims=[0 400]; bas_var_p2_lims=[0 400]; bas_var_diff_lims=[-1*10^-5 .5*10^-5]; 
        sur_var_p1_lims=[0 500]; sur_var_p2_lims=[0 500]; sur_var_diff_lims=[-1*10^-5 .5*10^-5];
        bas_orog_lims=[0 2000]; sur_orog_lims=[0 1500];
    end
    
    % Output folders
    base_figs_out=[root_dir 'figures/nonmaps/'];
    base_maps_out=[root_dir 'figures/maps/'];
    
    % Input folder
    folder_feature=[root_dir 'MapFeatures/'];
    cmip_folder=[root_dir 'CMIP5_FullGrid/'];
    
    %% User defined variables
    % Scenario Name Options
    scen_name1='historical'; window_cut=[1980 2004];
    %scen_name2='rcp45';  window_cut2=[2040 2064];
    scen_name2='rcp85';  window_cut2=[2075 2099];
    
    % window_cut define the time periods to retain the spatial fields.
    window_ave =[1980 2099]; % total year available
    
%     mon_comp=1:12;
    mon_comp=[1 2 12];
%     mon_comp=[6 7 8];
%     mon_comp=[3 4 5];
%     mon_comp=[9 10 11];
%     mon_comp=7;
%     mon_comp=[10 11 12 1 2 3];
%     mon_comp=[4 5 6 7 8 9];
    
    if flag_obs
        cru_years =1980:2004;
        cmap_years=1980:2004;
        gpcp_years=1980:2004;
    end
    
    num_months=12;
    num_ensembles=10;
    num_cells_cut=25; %% increase for large areas and high resolution models
    
    yrs_all=window_ave(1):window_ave(2);
    cut_yrs=window_cut(1):window_cut(2);
    cut_yrs2=window_cut2(1):window_cut2(2);
    cut_yrs=cut_yrs-1979;
    cut_yrs2=cut_yrs2-1979;
    
    if strcmp(variable,'pr')
        %dataset_name='atmos'; var_name='pr'; unit_name='mm day^-^1'; var_name_cru='pre'; unit_name_cru='mm/month'; var_name_cmap='prec'; title_name='Precipitation'; labelaxis='P'; L=[-55 -45 -35 -25 -15 -5  5 15 25 35 45 55]; c_flag=18;  precip_flag=1; scale =1;
        variable_flag=0; dataset_name='atmos'; var_name='pr'; unit_name='kg m^-^2 s^-^1'; var_name_cru='pre'; unit_name_cru='mm/month'; var_name_cmap='prec'; var_name_gpcp='prcp'; title_name='Precipitation'; labelaxis='P'; L=0:200:8*200; c_flag=19;  precip_flag=1; scale=1; Ldiff=-500:100:500/2; c_flagdiff=19;
%         bas_var_p1_lims=[4*10^-3 15*10^-3]; bas_var_p2_lims=[4*10^-3 15*10^-3]; bas_var_diff_lims=[.4*10^-3 2*10^-3]; sur_var_p1_lims=[2*10^-3 15*10^-3]; sur_var_p2_lims=[2*10^-3 15*10^-3]; sur_var_diff_lims=[0*10^-3 2*10^-3];
    elseif strcmp(variable,'tas')
        variable_flag=1; dataset_name='atmos'; var_name='tas'; unit_name='K'; var_name_cru='tmp'; unit_name_cru='^oC'; title_name='Temperature'; labelaxis='T'; c_flag=9;  scale =1; Ldiff=-5:1:5; c_flagdiff=11;L=-10:5:10;% L=15:5:35; % %L=10:5:30;%
        convert_F=0;
%         bas_var_p1_lims=[-14 0]; bas_var_p2_lims=[-14 0]; bas_var_diff_lims=[4.5 8]; sur_var_p1_lims=[-15 10]; sur_var_p2_lims=[-15 10]; sur_var_diff_lims=[4 10];
%         bas_var_p1_lims=[-33 -13]; bas_var_p2_lims=[-33 -13]; bas_var_diff_lims=[4.5 8]; sur_var_p1_lims=[-33 -5]; sur_var_p2_lims=[-33 -5]; sur_var_diff_lims=[2 15];
%         bas_var_p1_lims=[5 15]; bas_var_p2_lims=[5 15]; bas_var_diff_lims=[4.5 8]; sur_var_p1_lims=[1 23]; sur_var_p2_lims=[1 23]; sur_var_diff_lims=[2 15];   
    elseif strcmp(variable,'tasmin')
        variable_flag=2; dataset_name='atmos'; var_name='tasmin'; unit_name='K'; title_name='Temperature Min'; labelaxis='T'; c_flag=9;  scale =1; Ldiff=-5:1:5; c_flagdiff=11;L=-10:5:10;% L=15:5:35; % %L=10:5:30;%
%         bas_var_p1_lims=[0 80]; bas_var_p2_lims=[0 80]; bas_var_diff_lims=[-15 5]; sur_var_p1_lims=[0 80]; sur_var_p2_lims=[0 80]; sur_var_diff_lims=[-15 5];
    elseif strcmp(variable,'snw')
        variable_flag=3; dataset_name='landIce'; var_name='snw'; unit_name='kg m^-^2'; title_name='Surface Snow Amonut'; c_flag=50; c_flagdiff=26; scale_units=1; 
%         bas_var_p1_lims=[0 150]; bas_var_p2_lims=[0 150]; bas_var_diff_lims=[-30 5]; sur_var_p1_lims=[0 150]; sur_var_p2_lims=[0 150]; sur_var_diff_lims=[-30 5];
    elseif strcmp(variable,'snc')
        variable_flag=4; dataset_name='landIce'; var_name='snc'; unit_name='%'; title_name='Surface Snow Area Fraction'; c_flag=50; c_flagdiff=26; scale_units=1; 
%         bas_var_p1_lims=[0 90]; bas_var_p2_lims=[0 90]; bas_var_diff_lims=[-5 -3]; sur_var_p1_lims=[0 90]; sur_var_p2_lims=[0 90]; sur_var_diff_lims=[-20 5];
    elseif strcmp(variable,'mrro')
        variable_flag=5; dataset_name='land'; var_name='mrro'; unit_name='kg m^-^2 s^-^-1'; title_name='Total Runoff'; c_flag=19; c_flagdiff=19; scale_units=1;
%         bas_var_p1_lims=[0 25*10^-6]; bas_var_p2_lims=[0 25*10^-6]; bas_var_diff_lims=[-2*10^-6 5*10^-6]; sur_var_p1_lims=[0 25*10^-6]; sur_var_p2_lims=[0 25*10^-6]; sur_var_diff_lims=[-2*10^-6 5*10^-6];
    elseif strcmp(variable,'mrsos')
        variable_flag=6; dataset_name='land'; var_name='mrsos'; unit_name='kg m^-^2'; title_name='Moisture in Upper Portion of Soil Column'; c_flag=11; c_flagdiff=11; scale_units=1;
%         bas_var_p1_lims=[10 60]; bas_var_p2_lims=[10 60]; bas_var_diff_lims=[-1 1]; sur_var_p1_lims=[10 60]; sur_var_p2_lims=[10 60]; sur_var_diff_lims=[-1.5 1.5];
    elseif strcmp(variable,'mrros')
        variable_flag=7; dataset_name='land'; var_name='mrros'; unit_name='kg m^-^2 s^-^-1'; title_name='Total Surface Runoff'; c_flag=19; c_flagdiff=10; scale_units=1;
%         bas_var_p1_lims=[0 8*10^-6]; bas_var_p2_lims=[0 8*10^-6]; bas_var_diff_lims=[-2.5*10^-6 .5*10^-6]; sur_var_p1_lims=[0 9*10^-6]; sur_var_p2_lims=[0 9*10^-6]; sur_var_diff_lims=[-3.5*10^-6 1*10^-6];
    elseif strcmp(variable,'huss')
        variable_flag=8; dataset_name='atmos'; var_name='huss'; unit_name=''; title_name='Near-Surface Specific Humidity'; c_flag=19; c_flagdiff=10; scale_units=1;
    end
    
    %% Map Resolution and Projection
    % Can either use custom map projection or let worldmap() do the work
    % (for yenisey defaults to eqdconic)
    flag_high=0;
    flag_worldmap=true; % false=>no worldmap, true=> yes worldmap
    map_proj_name='eqdcylin';  % Equidistant Cylindrical, flat proj, usually use
    %map_proj_name='lambert';  % Lambert Conformal Conic, good for midlatitudes, about size of US
    %map_proj_name='mollweid';                % Mollweide, shrinks poles, good for global
    %map_proj_name='hammer'; latlim=[-90 90];    % shrinks poles, good for global
    
    %% Colors
    color_ocean=rgb('LightSlateGrey');
    color_land=rgb('CornSilk');
    color_basin=rgb('DarkGreen');
    color_rivers=rgb('LightBlue');
    
    %% Basin Shapefile
    if strcmp(basin_name,'Yenisei')
        basinID =3; area_txt='Yenisei River Basin'; lonmin=70; lonmax=115; latmin=45; latmax=85;
    elseif strcmp(basin_name,'Lena')
        basinID=252; area_txt='Lena River Basin'; lonmin=100; lonmax=145; latmin=45; latmax=80;
    elseif strcmp(basin_name,'Amur')
        basinID=23; area_txt='Amur River Basin'; lonmin=110; lonmax=155; latmin=40; latmax=70;
    elseif strcmp(basin_name,'Ob')
        basinID=253; area_txt='Ob River Basin'; lonmin=50; lonmax=80; latmin=45; latmax=85;
    elseif strcmp(basin_name,'Indigirka')
        basinID=1; area_txt='Indigirka River Basin'; lonmin=130; lonmax=160; latmin=60; latmax=75;
    end
    
    %% Place plotting limits into vector
    lonlim=[lonmin lonmax];
    latlim=[latmin latmax];
    
    %% Create latitude boxes latitude-separated analysis
    latstep=5;
    lat_boxes=latmin:latstep:latmax;
    lat_box_names=strings(size(lat_boxes,2)-1,1);
    for i=1:size(lat_boxes,2)-1
        lat_box_names(i)=[num2str(lat_boxes(i)) char(176) 'N-' num2str(lat_boxes(i+1)) char(176) 'N'];
    end
    
    %% Physical and cultural features for plots
    country=shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp'],...
        'UseGeoCoords', true);
    provinces=shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces_lines_shp.shp'],...
        'UseGeoCoords', true);
    rivers_lakes=shaperead([folder_feature '10m_physical/ne_10m_rivers_lake_centerlines.shp'],...
        'UseGeoCoords', true);
    ocean=shaperead([folder_feature '10m_physical/ne_10m_ocean.shp'],...
        'UseGeoCoords', true);
    land=shaperead('landareas', 'UseGeoCoords', true);
    coastline=shaperead([folder_feature '10m_physical/ne_10m_coastline.shp'],...
        'UseGeoCoords', true);
    
    %% Define Watershed
    file_name='wri_basins/wribasin.shp';
    basin_shp_all=shaperead([folder_feature  file_name], 'UseGeoCoords', true);
    basin_shp=basin_shp_all(basinID);
    
    %% Plot river basin
    if flag_basin_map
        figure;
        hold on;
        if flag_worldmap
            if flag_world
                ax=worldmap('World'); shading flat;
            else
                ax=worldmap(latlim, lonlim); shading flat;
            end
            setm(ax, 'FFaceColor', color_ocean);
            setm(ax, 'frame', 'on', 'Grid', 'off', 'FontSize', 16, ...
                'FontName', 'helvetica', 'FontWeight', 'bold', 'LabelFormat', 'none');
        else
            axesm('MapProjection', map_proj_name, 'MapLonLimit', lonlim, 'MapLatLimit', ...
                latlim, 'FFaceColor', color_ocean, 'frame', 'on', 'Grid', 'off', ...
                'FontSize', 16, 'FontName', 'helvetica', 'FontWeight', 'bold'); shading flat;
            tightmap;
        end
        geoshow(land, 'FaceColor', color_land, 'FaceAlpha', 1);
        geoshow([country.Lat], [country.Lon], 'Color', 'black');
        patchm([basin_shp.Lat], [basin_shp.Lon], color_basin);
        geoshow(rivers_lakes, 'Color', color_rivers);
        title(area_txt, 'fontsize', 16);
        if save_results
            print('-djpeg','-painters','-r600',[base_maps_out basin_name '/Basin_Map_' area_txt '.jpg'])
        end
    end
    %% Plot vectors
    %CS={'k','b','r','g','m','c'};%,'y'
    CS =[0.1059,0.6196,0.4667;...
        0.8510,0.3725,0.0078;...
        0.4588,0.4392,0.7020;...
        0.9059,0.1608,0.5412;...
        0.4000,0.6510,0.1176;...
        0.9020,0.6706,0.0078;...
        0.6510,0.4627,0.1137;...
        0.4000,0.4000,0.4000];
    MS={'o','s','d','+','*','s','x','o'};%'+','*','x','^','>','<','p','h'};
    LS={'-.','--','-',':'};
    
    %% Get colormap for plots
    cmap=colormap_function(c_flag);
    cmapdiff=colormap_function(c_flagdiff);
    
    %% Create a monthly name vector
    monthName={...
        'JAN'
        'FEB'
        'MAR'
        'APR'
        'MAY'
        'JUN'
        'JUL'
        'AUG'
        'SEP'
        'OCT'
        'NOV'
        'DEC'};
    
    %% Determine name of seasons for output file
    seas='';
%     for s=1:length(mon_comp)
%         tempMon=monthName{mon_comp(s)};
%         seas(s)=tempMon(1);
%     end
%     seasall=seas;
    if length(mon_comp)==12
        seas='ANNUAL';
    else
        for i=1:length(mon_comp)-1
            seas=strcat(seas,monthName{mon_comp(i)});
            seas=strcat(seas,'-');
        end
        seas=strcat(seas,monthName{mon_comp(end)});
    end
    
    %% Source folders
    if flag_obs
        
        cru_folder=[root_dir 'CRU_data/'];
        
        if variable_flag == 0
            cmap_folder=[root_dir 'CMAP_data/'];
            gpcp_folder=[root_dir 'GPCP_data/'];
        end
        
        %% Load CRU data
        load([cru_folder var_name_cru '.' basin_name '.cru3_10.mat']);
        
        %% Load CMAP and GPCP Precipitation data
        if variable_flag == 0
            load([cmap_folder var_name_cmap '.' basin_name '.cmap.mat']);
            load([gpcp_folder var_name_gpcp '.' basin_name '.gpcp.mat']);
        end
    end

    %% OBSERVATIONS
    if variable_flag == 0 && flag_obs % Precipitation
    %% Plot Maps
        %% CRU data
        
        basin_delta_plot_lon=nanmean(diff(crudata.lon_basin));
        basin_delta_plot_lat=nanmean(diff(crudata.lat_basin));
        longitude=double(crudata.lon_basin);
        latitude= double(crudata.lat_basin);
        
        obs_map=nanmean(crudata.basin_grid,1);
        obs_map=squeeze(nansum(obs_map));
        obs_map(obs_map==0)=NaN;
        
        %% Plot CRU data
        cru_lims=[min(L),max(L)];
        cru_label=[title_name '(' unit_name ')'];
        plot_basin(obs_map,cru_lims,cmap,cru_label,latitude-basin_delta_plot_lat/2, longitude-basin_delta_plot_lon/2,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
        title(['CRU TS 3.1, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
        if save_results
            print('-depsc2','-painters','-r600',[base_maps_out var_name '.map.CRU.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.eps' ])
            %print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '.map.CRU.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.jpg' ])
        end
        
        %% CMAP
        basin_delta_plot_lon=double(abs(nanmean(diff(cmapdata.lon_basin))));
        basin_delta_plot_lat=double(abs(nanmean(diff(cmapdata.lat_basin))));
        longitude=double(cmapdata.lon_basin);
        latitude= double(cmapdata.lat_basin);
        
        obs_map=nanmean(cmapdata.basin_grid,1);
        obs_map=squeeze(nansum(obs_map));
        obs_map(obs_map==0)=NaN;
        
        %% Plot CMAP
        cmap_lims=[min(L) max(L)];
        cmap_label=[title_name '(' unit_name ')'];
        plot_basin(obs_map,cmap_lims,cmap,cmap_label,latitude-basin_delta_plot_lat/2, longitude-basin_delta_plot_lon/2,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
        title(['CMAP V3.1, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
        if save_results
            print('-depsc2','-painters','-r600',[base_maps_out var_name '.map.CMAP.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.eps' ])
            %print('-djpeg', '-painters','-r600',[base_maps_out basin_name '/' var_name '.map.CMAP.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.jpg' ])
        end
        %% GPCP
        basin_delta_plot_lon=double(abs(nanmean(diff(gpcpdata.lon_basin))));
        basin_delta_plot_lat=double(abs(nanmean(diff(gpcpdata.lat_basin))));
        longitude=double(gpcpdata.lon_basin);
        latitude= double(gpcpdata.lat_basin);
        
        obs_map=nanmean(gpcpdata.basin_grid,1);
        obs_map=squeeze(nansum(obs_map));
        obs_map(obs_map==0)=NaN;
        
        %% Plot GPCP
        gpcp_lims=[min(L) max(L)];
        gpcp_label=[title_name '(' unit_name ')'];
        plot_basin(obs_map,gpcp_lims,cmap,gpcp_label,latitude-basin_delta_plot_lat/2, longitude-basin_delta_plot_lon/2,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
        title(['GPCP V2.2, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
        if save_results
            %print('-depsc2','-painters','-r600',[base_maps_out var_name '.map.GPCP.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.eps' ])
            %print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '.map.GPCP.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.jpg' ])
        end
    elseif variable_flag == 1 && flag_obs % Temperature
        if convert_F == 1
            var_plot=9*nanmean(crudata.basin_mean(:,mon_comp))/5+32;
        else
            var_plot=nanmean(crudata.basin_mean(:,mon_comp));
        end
        
        basin_delta_plot_lon=nanmean(diff(crudata.lon_basin));
        basin_delta_plot_lat=nanmean(diff(crudata.lat_basin));
        longitude=double(crudata.lon_basin);
        latitude= double(crudata.lat_basin);
        obs_map=nanmean(crudata.basin_grid,1);
        obs_map=squeeze(nanmean(obs_map));
        [longitude, latitude]=meshgrid(longitude, latitude);
                
        %% Plot CRU map
        cru_lims=[min(L) max(L)];
        cru_label=[title_name '(' unit_name ')'];
        plot_basin(obs_map,cru_lims,cmap,cru_label,latitude-basin_delta_plot_lat/2, longitude-basin_delta_plot_lon/2,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
        title(['CRU TS 3.1, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
        hold off
        if save_results
            print('-depsc2','-painters','-r600',[base_maps_out var_name '.map.CRU.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.eps' ])
    %         print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '.map.CRU.' seas '.'  basin_name '.' num2str(window_cut(1)) '.' num2str(window_cut(2)) '.jpg' ])
        end
    end

    %% Load model information
    for i_model=1:size(model_names,1)          % MODEL LOOP
        modname=model_names{i_model,1} %#ok<NOPTS> 

        TF=strfind(modname,'-'); modname(TF)='_';
        filename=[cmip_folder   var_name  '.'  modname '.' scen_name1 '-' scen_name2 '.ar5.monmean.mat'];
        file_TF=exist(filename, 'file');
        % Load maps of ensemble averaged variables
        if file_TF==2
            % Load file
            load(filename);
            
             % initialize grids for averaging
            grid1=nan(size(ensemble_names,1),size(lat,1),size(lon,1));
            grid2=nan(size(ensemble_names,1),size(lat,1),size(lon,1));
            % ENSEMBLE LOOP
            for i_ensemb=1:size(ensemble_names, 1)
                S=exist(ensemble_names{i_ensemb}, 'var');
                ensembname=ensemble_names{i_ensemb};
                if S==1
                    cur_ensemble=eval(ensembname);
                    model_period1.grid=squeeze(mean(cur_ensemble.clim_grid(cut_yrs,:,:,:),1));
                    model_period2.grid=squeeze(mean(cur_ensemble.clim_grid(cut_yrs2,:,:,:),1));
                    
                    % Figure: Basin check 2
%                     figure; hold on; shading flat;
%                     worldmap('World');
%                     pcolorm(lat, lon, squeeze(model_period1.grid(1,:,:)));
%                     geoshow([country.Lat], [country.Lon], 'Color', 'black');
%                     patchm([basin_shp.Lat], [basin_shp.Lon], 'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 2);
%                     plot(g);

                    if variable_flag == 0
                        ens_grid1=squeeze(nanmean(model_period1.grid(mon_comp,:,:),1))*365;
                        ens_grid2=squeeze(nanmean(model_period2.grid(mon_comp,:,:),1))*365;
                    elseif variable_flag == 1 || variable_flag==2
                        ens_grid1=squeeze(nanmean(model_period1.grid(mon_comp,:,:),1))-273.15;
                        ens_grid2=squeeze(nanmean(model_period2.grid(mon_comp,:,:),1))-273.15;
                    else
                        ens_grid1=squeeze(nanmean(model_period1.grid(mon_comp,:,:),1));
                        ens_grid2=squeeze(nanmean(model_period2.grid(mon_comp,:,:),1));
                    end
                    
                    grid1(i_ensemb,:,:)=ens_grid1;
                    grid2(i_ensemb,:,:)=ens_grid2;
                end
            end % ensemble loop
           
            grid1=squeeze(nanmean(grid1,1));
            grid2=squeeze(nanmean(grid2,1));
                        
            % load elevation data
            if plot_elevation
                e_scen_names={'historical';
                    'amip';
                    'decadal2000';
                    'historicalNat';
                    };
                for i=1:size(e_scen_names,1)
                    filename=[cmip_folder 'orog.' modname '.' e_scen_names{i,1} '.ar5.monmean.mat'];
                    file_TF=exist(filename, 'file');
                    if file_TF==2
                        break
                    end
                end
            else
                orog_grid=nan(size(grid1,1),size(grid2,2)); % awkward way to avoid problems where orog_grid later called
            end

            load(filename);

            orog_grid=r0i0p0.clim_grid;
            
            for i_basin_scope=1:2 % 1=> basin+surrounding area data, 2=>
                % just basin data
                if i_basin_scope==1
                    basin_scope='surrounding';
                else
                    basin_scope='basin';
                end
                if flag_world
                    basin_scope='world';
                end
                if flag_land_mask
                    basin_scope=[basin_scope '.land_mask'];
                end
                if custom_var_lims
                    basin_scope=[basin_scope '.custom_lims'];
                end
                
                % Create Basin Mask
                basin_mask=nan(length(lat),length(lon));
                land_mask=nan(length(lat),length(lon));

                % Clear variables from previous model
                clear lonlat lonlat_loc lonlat_loc_ind lat_data_cut lon_data_cut
                clear lonlat_loc_ind lat_data_cut lon_data_cut

                % Place lat-lon in array for inpolygon
                jj = 1;
                for i=1:length(lon)
                    for j=1:length(lat)
                        lonlat(jj,1) = lon(i);
                        lonlat(jj,2) = lat(j);
                        jj = jj+1;
                    end
                end

                if i_basin_scope==1
                    % Find grid cells within latlim and basin+surrounding polygon
                    [IN,ON] = inpolygon(lonlat(:,1),lonlat(:,2),lonlim,latlim);
                    lonlat_loc = IN+ON; clear IN ON;
                    lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
                    lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
                    lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
                    for i = 1:size(lonlat_loc_ind,1)
                        basin_mask(lat == lat_data_cut(i),lon == lon_data_cut(i))=1;
                    end
                else
                    % Find grid cells within latlim and basin polygon
                    [IN,ON] = inpolygon(lonlat(:,1),lonlat(:,2),basin_shp.Lon,basin_shp.Lat);
                    lonlat_loc = IN+ON; clear IN ON;
                    lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
                    lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
                    lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
                    for i = 1:size(lonlat_loc_ind,1)
                        basin_mask(lat == lat_data_cut(i),lon == lon_data_cut(i))=1;
                    end
                end
              
                if ~flag_world
                    grid1=grid1.*basin_mask;
                    grid2=grid2.*basin_mask;
                    orog_grid=orog_grid.*basin_mask;
                end
                
                if flag_land_mask
                    % Clear variables from previous model
                    clear lonlat_loc lonlat_loc_ind lat_data_cut lon_data_cut
                    clear lonlat_loc_ind lat_data_cut lon_data_cut
                
                    % Find grid cells within land mask
                    [IN,ON] = inpolygon(lonlat(:,1),lonlat(:,2),[land.Lon],[land.Lat]);
                    lonlat_loc = IN+ON; clear IN ON;
                    lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
                    lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
                    lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
                    for i = 1:size(lonlat_loc_ind,1)
                        land_mask(lat == lat_data_cut(i),lon == lon_data_cut(i))=1;
                    end
                    
                    grid1=grid1.*land_mask;
                    grid2=grid2.*land_mask;
                    orog_grid=orog_grid.*land_mask;
                end

                griddiff=grid2-grid1;
                
                % Figure: Basin check 1
    %             figure; shading flat; pcolor(lon,lat,landfrac.*basin_mask),title('model land fraction, basin mask'),colorbar, grid off;
    %             plot(g);
    %             % Figure: Basin check 2
    %             figure; hold on; shading flat;
    %             worldmap('World');
    %             pcolorm(lat, lon, orog_grid);
    %             geoshow([country.Lat], [country.Lon], 'Color', 'black');
    %             patchm([basin_shp.Lat], [basin_shp.Lon], 'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 2);
    % %             plot(g); % fail to make it stop :)
    
                % remove outliers
                if flag_no_outliers
                    orog_grid(isoutlier(griddiff))=nan;
                    grid1(isoutlier(grid1))=nan;
                    grid2(isoutlier(grid2))=nan;
                    griddiff(isoutlier(griddiff))=nan;
                end

                % offset grid points for plotting
                if i_basin_scope==1
                    if strcmp(modname,'GISS_E2_R')
                        delta_plot_lat=0;
                        delta_plot_lon=0;
                    elseif strcmp(modname, 'CNRM_CM5')
                        delta_plot_lat=nanmean(diff(lon))/1.5;
                        delta_plot_lon=nanmean(diff(lat))/1.5;
                    else
                        delta_plot_lat=nanmean(diff(lon))/2;
                        delta_plot_lon=nanmean(diff(lat))/2;
                    end
                elseif i_basin_scope==2
                    if strcmp(modname,'GISS_E2_R')
                        delta_plot_lat=nanmean(diff(lon))/2;
                        delta_plot_lon=nanmean(diff(lat))/2;
                    elseif strcmp(modname, 'CCSM4')
                        delta_plot_lat=nanmean(diff(lon))/2;
                        delta_plot_lon=nanmean(diff(lat))/2;
                    elseif strcmp(modname, 'CNRM_CM5')
                        delta_plot_lat=nanmean(diff(lon))/1.9;
                        delta_plot_lon=nanmean(diff(lat))/1.9;
                    end
                end
                
                % set color lims
                if custom_var_lims
                    if i_basin_scope==1
                        var_p1_lims=sur_var_p1_lims;
                        var_p2_lims=sur_var_p2_lims;
                        var_diff_lims=sur_var_diff_lims;
                        orog_lims=sur_orog_lims;
                    else
                        var_p1_lims=bas_var_p1_lims;
                        var_p2_lims=bas_var_p2_lims;
                        var_diff_lims=bas_var_diff_lims;
                        orog_lims=bas_orog_lims;
                    end
                else
                    var_p1_lims=[round(min(min(grid1)),2) round(max(max(grid1)),2)]; 
                    var_p2_lims=[round(min(min(grid2)),2) round(max(max(grid2)),2)]; 
                    var_diff_lims=[round(min(min(griddiff)),2) round(max(max(griddiff)),2)];
                    orog_lims=[round(min(min(orog_grid)),2) round(max(max(orog_grid)),2)];
                end
                
                p1_file_ext=[modname '.' scen_name1 '-' scen_name2 '-' num2str(window_cut(1)) '-' num2str(window_cut(2)) '.' seas '.jpg'];
                p2_file_ext=[modname '.' scen_name1 '-' scen_name2 '-' num2str(window_cut2(1)) '-' num2str(window_cut2(2))  '.' seas '.jpg'];
                dif_file_ext=[modname '.' scen_name1 '-' scen_name2 '-diff-' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) '-' num2str(window_cut(1)) '-' num2str(window_cut(2)) '.' seas '.jpg'];

                % plot and save figure of mean variable: period 1
                if plot_p1
                    p1_label=[title_name ' (' unit_name ')'];
                    plot_basin(grid1,var_p1_lims,cmap,p1_label,lat-delta_plot_lat,lon-delta_plot_lat,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
                    title({[area_txt ' ' title_name ', ' scen_name1 ' (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],[seas ' ' model_names{i_model,1}]});
                    hold off
                    if save_results
                       print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '/Basin_Map_' area_txt '.' var_name '.' basin_scope '.'  custom_save_string p1_file_ext]);
                    end
                end
                
                % plot and save figure of mean variable: period 2
                if plot_p2
                    p2_label=[title_name ' (' unit_name ')'];
                    plot_basin(grid2,var_p2_lims,cmap,p2_label,lat-delta_plot_lat,lon-delta_plot_lon,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
                    title({[model_names{i_model,1} ' ' title_name ', ' scen_name2 ' (' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) ')'],seas});
                    hold off
                    if save_results
                        print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '/Basin_Map_' area_txt '.' var_name '.' basin_scope '.' custom_save_string p2_file_ext])
                    end

                end

                % plot and save figure of mean variable difference (period2-period1)
                if plot_diff
                    diff_label=['Difference in ' title_name ' (' unit_name ')'];
                    plot_basin(griddiff,var_diff_lims,cmapdiff,diff_label,lat-delta_plot_lat,lon-delta_plot_lon,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
%                     title({[model_names{i_model,1} ' ' title_name ', ' 'future' ' (' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) ') minus ' scen_name1 ' (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],seas});
                    title({[title_name ', Future (' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) ') Minus Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],seas});
                    hold off
                    if save_results
                        print('-djpeg','-painters','-r600',[base_maps_out basin_name '/' var_name '/Basin_Map_' area_txt '.' var_name '.' basin_scope '.'  custom_save_string dif_file_ext])
                    end
                end

                % plot elevation
                if plot_elevation
                    orog_label='elevation (km)';
                    orog_grid_scaled=orog_grid/1000;
                    orog_lims=[0 2];
                    plot_basin(orog_grid_scaled,orog_lims,cmapdiff,orog_label,lat-delta_plot_lat,lon-delta_plot_lon,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world);
%                     title([area_txt ' elevation ' model_names{i_model,1}], 'fontsize',16);
                    title([area_txt ' Elevation ']);
                    hold off
                    if save_results
                        print('-djpeg','-painters','-r600',[base_maps_out basin_name '/orog/Basin_Map_' area_txt '.orog.' basin_scope '.'  modname '.jpg'])
                    end
                end

                % plot period difference as function of elevation
                if plot_diff_scatter
                    flag_cutoff=false;
                    flag_p1=false; % linear fit
                    flag_p2=false; % polynoimial (deg 2) fit
                    flag_log=false; % log fit
                    flag_lat_subsets=true; % show different in lat bracket

                    cutoff=0;

                    orog_scale=.001;
                    orog_grid_scaled=orog_grid.*orog_scale;

                    var_diff_vec=reshape(griddiff,[1 size(griddiff,1)*size(griddiff,2)]);
                    var_diff_vec(isnan(var_diff_vec))=[];

                    orog_vec=reshape(orog_grid_scaled,[1 size(orog_grid_scaled,1)*size(orog_grid_scaled,2)]);
                    orog_vec(isnan(orog_vec))=[];

                    if flag_cutoff
                        cutoff=cutoff*orog_scale;
                        var_diff_vec(orog_vec<=cutoff)=[];
                        orog_vec(orog_vec<=cutoff)=[];
                    end
                    
                    figure;
                    hold on;
                    if flag_lat_subsets
                         for i_lat=1:size(lat_boxes,2)-1
                            x=nan(size(orog_grid,1),size(orog_grid,2));
                            x(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:)=orog_grid_scaled(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:);
                            x_vect=reshape(x,[size(x,1)*size(x,2) 1]);
                            x_vect(isnan(x_vect))=[];

                            y=nan(size(griddiff,1),size(griddiff,2));
                            y(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:)=griddiff(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:);
                            y_vect=reshape(y,[size(y,1)*size(y,2) 1]);
                            y_vect(isnan(y_vect))=[];

                            if flag_cutoff
                                y_vect(x_vect<=cutoff)=[];
                                x_vect(x_vect<=cutoff)=[];
                            end

                            if size(x_vect,1)>5 && size(x_vect,1)==size(y_vect,1)
                                lat_box_name=[num2str(lat_boxes(i_lat)) char(176) 'N-' num2str(lat_boxes(i_lat+1)) char(176) 'N'];
                                scatter(x_vect,y_vect,'filled','DisplayName',lat_box_name);
                            end
                        end
                    else
                        ax=scatter(orog_vec,var_diff_vec,'filled','DisplayName','none');
                    end
                    xl=xlim;
                    x=linspace(0.001,xl(2));
                    if flag_p1
                        p1=polyfit(orog_vec,var_diff_vec,1);
                        y1=polyval(p1,x);
                        eq1=[num2str(p1(1)) 'x + ' num2str(p1(2))]';
                        plot(x,y1,'r','DisplayName',eq1);
                    elseif flag_p2
                        p2=polyfit(orog_vec,var_diff_vec,2);
                        y2=polyval(p2,x);
                        eq2=[num2str(p2(1)) 'x^2 + ' num2str(p2(2)) 'x + ' num2str(p2(3))];
                        plot(x,y2,'g','DisplayName',eq2);
                    elseif flag_log
                        logx=log(orog_vec(orog_vec~=0));
                        logy=log(var_diff_vec(var_diff_vec~=0));
                        pnl=polyfit(logx,var_diff_vec(orog_vec~=0),1);
                        ynl=polyval(pnl,log(x(x~=0)));
                        eqnl=[num2str(pnl(2)) ' + ' num2str(pnl(1)) 'ln(x)'];
                        plot(x(x~=0),ynl, 'DisplayName', eqnl);
                    end
                    legend('show', 'location', 'best');
                    xlabel('Elevation (km)');
                    ylabel([title_name ' Change (' unit_name ')']);
                    hold off;
%                     title({[model_names{i_model,1} ' ' title_name ', ' 'future' ' (' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) ') minus ' scen_name1 ' (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],seas});
                    title({[title_name ', Future (' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) ') Minus Historical' ' (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],seas});
                    if save_results
                        if flag_lat_subsets
                            print('-djpeg','-painters','-r600',[base_figs_out basin_name '/' var_name '/scatter_' area_txt '.elevation.' var_name '.' basin_scope custom_save_string '.' custom_save_string 'lat_subsets.'  dif_file_ext])
                        else
                            print('-djpeg','-painters','-r600',[base_figs_out basin_name '/' var_name '/scatter_' area_txt '.elevation.' var_name  '.' basin_scope '.' custom_save_string dif_file_ext])
                        end
                    end
                end % if plot dif scatter
                
                if lat_box_table
                    % find and save current model r and p values for lat boxes
                    flag_cutoff=false;

                    orog_scale=.001;
                    orog_grid_scaled=orog_grid*orog_scale;

                    % get lat-separated r and p values for basin
                    for i_lat=1:size(lat_boxes,2)-1
                        x=nan(size(orog_grid,1),size(orog_grid,2));
                        x(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:)=orog_grid_scaled(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:);
                        x_vect=reshape(x,[size(x,1)*size(x,2) 1]);
                        x_vect(isnan(x_vect))=[];

                        y=nan(size(griddiff,1),size(griddiff,2));
                        y(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:)=griddiff(lat>=lat_boxes(i_lat)&lat<lat_boxes(i_lat+1),:);
                        y_vect=reshape(y,[size(y,1)*size(y,2) 1]);
                        y_vect(isnan(y_vect))=[];

                        if flag_cutoff
                            y_vect(x_vect<=cutoff)=[];
                            x_vect(x_vect<=cutoff)=[];
                        end

                        if size(x_vect,1)>5
                            [R,P,RL,RU]=corrcoef(x_vect,y_vect);
                            if i_basin_scope==1
                                model_r_values(i_lat,i_model*2,model_res)=R(1,2);
                                model_p_values(i_lat,i_model*2,model_res)=P(1,2);
                            else
                                model_r_values(i_lat,i_model*2-1,model_res)=R(1,2);
                                model_p_values(i_lat,i_model*2-1,model_res)=P(1,2);
                            end
                        end
                    end
                end % if lat box table
            end % basin scope loop
            
            clear lat lon landfrac model
            clear(ensemble_names{:});
            
        end
%                     nan(g); % fail to make it stop :)
    end % model loop
end % model res loop

if lat_box_table
    % make table of p and r values for all models surrounding and no
    % surrounding
    
    mod_scen_names=cell(1,size(model_names,1)*2);
    for i_scen=1:size(model_names,1)
        modname=model_names{i_scen,1};
        TF=strfind(modname,'-'); modname(TF)='_';
        mod_scen_names{1,i_scen*2-1}=modname;
        mod_scen_names{1,i_scen*2}=[modname '_plus_surroundings'];
    end
    
    lat_box_names_cell=cell(1,size(lat_box_names,1));
    for i_lat_box=1:size(lat_box_names,1)
        lat_box_names_cell{1,i_lat_box}=lat_box_names{i_lat_box};
    end
    
    model_r_values=reshape(model_r_values,size(model_r_values,1),size(model_r_values,2)*2);
    model_r_values(~any(~isnan(model_r_values),2),:)=[];
    model_r_values(:,~any(~isnan(model_r_values),1))=[];
    
    model_p_values=reshape(model_p_values,size(model_p_values,1),size(model_p_values,2)*2);
    model_p_values(~any(~isnan(model_p_values),2),:)=[];
    model_p_values(:,~any(~isnan(model_p_values),1))=[];
    
    lat_box_names_cell(:,size((model_r_values),1)+1:end)=[];
         
    T_R=array2table(model_r_values,'VariableNames',mod_scen_names,'RowNames',lat_box_names_cell);
    T_P=array2table(model_p_values,'VariableNames',mod_scen_names,'RowNames',lat_box_names_cell);    
    
%     figure
%     uitable('Data',T_R{:,:},'ColumnName',T_R.Properties.VariableNames,...
%     'RowName',T_R.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% 
%     figure
%     uitable('Data',T_P{:,:},'ColumnName',T_P.Properties.VariableNames,...
%     'RowName',T_P.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    
    if save_results
        writetable(T_R, [base_figs_out basin_name '/' var_name '/lat_table_' area_txt '.elevation.' var_name '.' custom_save_string 'r_values.' scen_name1 '-' scen_name2 '-diff-' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) '-' num2str(window_cut(1)) '-' num2str(window_cut(2)) '.' seas '.txt']);
        writetable(T_P, [base_figs_out basin_name '/' var_name '/lat_table_' area_txt '.elevation.' var_name '.' custom_save_string 'p_values.' scen_name1 '-' scen_name2 '-diff-' num2str(window_cut2(1)) '-' num2str(window_cut2(2)) '-' num2str(window_cut(1)) '-' num2str(window_cut(2)) '.' seas '.txt']);
    end
end