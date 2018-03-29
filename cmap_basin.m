%% Basin (or rectangular box) average of CMAP data (precipitation ONLY)
% The CPC Merged Analysis of Precipitation ("CMAP") is a technique which
% produces pentad and monthly analyses of global precipitation in which
% observations from raingauges are merged with precipitation estimates from
% several satellite-based algorithms (infrared and microwave). The analyses
% are are on a 2.5 x 2.5 degree latitude/longitude grid and extend back to
% 1979.

%CHANGE directory stuff
%BIG CHANGE rewritting m_map stuff with mappping toolbox
%CHANGE land shapefile source from Natural Earth to mapping toolbox, 
%previous one did not form clean polygon topology
%CHANGE "Yenisey" to "Yenisei"

for i_site = 17:17%1:num_sites
%% Start clean
clearvars -except i_site num_sites; 
close all; clc;
root_dir = '/Users/fullerjohnjef/climate/';

basefolder = [root_dir 'CMAP_data/'];

%% Folders
folder_feature = [root_dir 'MapFeatures/'];
folder_output = [root_dir 'CMIP5_output_continuous/'];
flag_high = 1; % 0 = 50m resolution, 1 = 10m resolution

%% Map Projections
% Can either use custom map projection or let worldmap() do the work
flag_worldmap = true; %false=> no worldmap true=> worldmap
map_proj_name='eqdcylin';  % Flat proj, usually use
%map_proj_name='lambert';  % Good for midlatitudes, about size of US
%map_proj_name='mollweid';                % shrinks poles, good for global
%map_proj_name='hammer'; latlim=[-90 90];               % shrinks poles, good for global

%% Select cmap years
cmap_years    = 1980:2004;

%% Colors
color_ocean = rgb('LightSlateGrey');
color_land = rgb('Cornsilk');
color_basin = rgb('DarkGreen');
color_rivers = rgb('LightBlue');

%%  Set flag_basin_shpfile
% 0 - Grid-based dataset (UNH/GRDC basins)
% 1 - Shapefile (uses inpologon)
% 2 - Latitude-Longitude Box

%% UNDP Analysis
if i_site == 1
    % 1. Niger River
    area_txt = 'Niger River'; basin_name = 'NigerRiver'; flag_basin_shpfile = 1; lonmin = -17; lonmax = 17; latmin = 0; latmax = 30;
    basinID =115; flag_wri = 1;  file_name = 'wri_basins/wribasin.shp';
elseif i_site == 2
    % % 2. Lake Chad
    area_txt = 'Lake Chad'; basin_name = 'LakeChad'; file_name = 'chadbasin_gen/chadbasin_gen.shp'; flag_basin_shpfile = 1; lonmin = 0; lonmax = 30; latmin = 0; latmax = 30;
    basinID =244; flag_wri = 0;
elseif i_site == 3
    % 3. Nile River
    area_txt = 'Nile River'; basin_name = 'NileRiver'; file_name = 'wri_basins/wribasin.shp'; flag_basin_shpfile = 1; lonmin = 18; lonmax = 48; latmin = -9; latmax = 32;
    basinID =108; flag_wri = 1;
elseif i_site == 4
    % 4. Congo River
    area_txt = 'Congo River'; basin_name = 'CongoRiver'; file_name = 'wri_basins/wribasin.shp'; flag_basin_shpfile = 1; lonmin = 0; lonmax = 50; latmin = -18; latmax = 12;
    basinID =156; flag_wri = 1;
elseif i_site == 5
    % 6. Okavango River
    area_txt = 'Okavango River'; basin_name = 'OkavangoRiver'; file_name = 'Okavango_megabasin/Okavango_basin.shp'; flag_basin_shpfile = 1; lonmin = 6; lonmax = 36; latmin = -36; latmax = -11;
    basinID =209; flag_wri = 0;
elseif i_site == 6
    % 7. Orange-Sengu River
    area_txt = 'Orange-Sengu River'; basin_name = 'OrangeSenguRiver'; file_name = 'wri_basins/wribasin.shp'; flag_basin_shpfile = 1; lonmin = 6; lonmax = 36; latmin = -36; latmax = -11;
    basinID =227; flag_wri = 1;
elseif i_site == 7
    % 8. Dnipro River
    area_txt = 'Dnipro River'; basin_name = 'DniproRiver';file_name = 'wri_basins/wribasin.shp';  flag_basin_shpfile = 1; lonmin = 20; lonmax = 45; latmin = 40; latmax = 60;
    basinID =24; flag_wri = 1;
elseif i_site == 8
    % 11. Tigris-Euphrates River
    area_txt = 'Tigris-Euphrates River'; basin_name = 'TigrisEuphratesRiver'; file_name = 'wri_basins/wribasin.shp'; flag_basin_shpfile = 1; lonmin = 30; lonmax = 60; latmin = 24; latmax = 48;
    basinID =83; flag_wri = 1;
elseif i_site == 9
    % 12. Kura-Aras River
    area_txt = 'Kura-Aras River'; basin_name = 'KuraArasRiver'; flag_basin_shpfile = 1; lonmin = 40; lonmax = 55; latmin = 36; latmax = 48;
    basinID =69; flag_wri = 1; file_name = 'wri_basins/wribasin.shp';
elseif i_site == 10
    % 14. Lake Baikal
    area_txt = 'Lake Baikal'; basin_name = 'LakeBaikal'; flag_basin_shpfile = 1; lonmin = 90; lonmax = 120; latmin = 42; latmax = 62;
    basinID =22; flag_wri = 1; file_name = 'wri_basins/wribasin.shp';
elseif i_site == 11
    % 15. Lake Titicaca-Rio Mauri
    area_txt = 'Lake Titicaca-Rio Mauri'; basin_name = 'LakeTiticaca'; flag_basin_shpfile = 1; lonmin = -76; lonmax = -64; latmin = -22; latmax = -12;
    basinID =214; flag_wri = 1; file_name = 'wri_basins/wribasin.shp';
elseif i_site == 12
    % 5. Lake Tanganyika **lat-lon box
    area_txt = 'Lake Tanganyika'; basin_name = 'LakeTanganyika'; flag_basin_shpfile = 2;
    lonmin = 20; lonmax = 40; latmin = -14; latmax = 4; lon_region_min = 28; lon_region_max = 33; lat_region_min = -10; lat_region_max = -2; 
    flag_wri = 0; basinID =-99;%156; 
 
elseif i_site == 13
    % 9. Drin River **lat-lon box
    area_txt = 'Drin River'; basin_name = 'DrinRiver'; flag_basin_shpfile = 2; 
    lon_region_min = 18; lon_region_max = 23; lat_region_min = 40.5; lat_region_max = 43.5; 
    lonmin = 16; lonmax = 24; latmin = 40; latmax = 45;
    basinID =65; flag_wri = 0;
elseif i_site == 14
    % 10. Lake Prespa **lat-lon box
    area_txt = 'Lake Prespa'; basin_name = 'LakePrespa'; flag_basin_shpfile = 2; 
    lon_region_min = 19.5; lon_region_max = 22.5; lat_region_min = 39.5; lat_region_max = 42.5; 
    lonmin = 16; lonmax = 24; latmin = 38; latmax = 44;
    basinID =-99; flag_wri = 0;
elseif i_site == 15
    % 16. Artibonite River **lat-lon box
    area_txt = 'Artibonite River'; basin_name = 'ArtiboniteRiver'; flag_basin_shpfile = 2; 
    lon_region_min = -73.5; lon_region_max = -70.5; lat_region_min = 17.5; lat_region_max = 20.5; 
    lonmin = -76; lonmax = -68; latmin = 16; latmax = 22;
    basinID =-99; flag_wri = 0;
elseif i_site == 16
    % 13. Caspian Sea ** gridded dataset
    area_txt = 'Caspian Sea'; basin_name='CaspianSea'; flag_basin_shpfile = 0;  lonmin = 30; lonmax = 70; latmin = 25; latmax = 70; basinID = [17 54 69 85 201 215 230 258 307 319 355 359 365 453 479 647 803 867 948 951 971 1215 1287 1315 1443 1476 1482 1506 1521 1527 1529 1534 1538 1570 2098 2102 2111 2127 2140 2152 2164 2174 2183 2191 2227 2268 2314 3863 3864 3865 3891 3930 3949 3978 4063 4079 4086 4116 4145 4146 4174 4234 4315 4316 4340 4361 4404];  %%% CaspianSea: Volga - 17;  Ural - 69; Kura -> 85; Atrek -> 54
    flag_wri = 0;
elseif i_site == 17    
    area_txt = 'Yenisei River'; basin_name = 'Yenisei'; flag_basin_shpfile = 1; lonmin = 70; lonmax = 115; latmin = 45; latmax = 85;
    basinID =3; flag_wri = 1;  file_name = 'wri_basins/wribasin.shp'; 
end

%% Place plotting limits into vector
lonlim=[lonmin lonmax];
latlim=[latmin latmax];


%% Physical and cultural features for plots
if flag_high == 1
    country = shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp'], ...
        'UseGeoCoords', true);
    provinces = shaperead([folder_feature '10m_cultural/10m_cultural/ne_10m_admin_1_states_provinces_lines_shp.shp'], ...
        'UseGeoCoords', true);
    rivers_lakes = shaperead([folder_feature '10m_physical/ne_10m_rivers_lake_centerlines.shp'], ...
        'UseGeoCoords', true);
    ocean = shaperead([folder_feature '10m_physical/ne_10m_ocean.shp'], ...
        'UseGeoCoords', true);
    land = shaperead('landareas', 'UseGeoCoords', true);
    coastline = shaperead([folder_feature '10m_physical/ne_10m_coastline.shp']);
elseif flag_high == 0
    country = shaperead([folder_feature '50m_cultural/ne_50m_admin_0_countries.shp'], ...
        'UseGeoCoords', true);
    provinces = shaperead([folder_feature '50m_cultural/ne_50m_admin_1_states_provinces_lines.shp'], ...
        'UseGeoCoords', true);
    rivers_lakes = shaperead([folder_feature '50m_physical/ne_50m_rivers_lake_centerlines.shp'], ...
        'UseGeoCoords', true);
    ocean = shaperead([folder_feature '50m_physical/ne_50m_ocean.shp'], ...
        'UseGeoCoords', true);
    land = shaperead('landareas', 'UseGeoCoords', true);
    coastline = shaperead([folder_feature '50m_physical/ne_50m_coastline.shp'], ...
        'UseGeoCoords', true);
end

%% Define Watershed
if (flag_basin_shpfile == 0) % Grid-based dataset (UNH/GRDC basins)
    % Allowed offset from grid cell centers
    basin_off=0.25; % count model cell if center +/- one quarter of model resolution is in watershed
    %basin_off=0;    % center MUST be in watershed
    
    % load basin data and create
    basins=flipud(load([folder_feature 'watershed/basin.grd'])); basins(find(basins<=-900))=nan;
    lat_basin=[-55.25:0.5:82.75]; lon_basin=[-179.75:0.5:179.75];
    
    % DISABLE:Change map orientation to match CMIP5 models
    basin_new=basins;%[basins(:,361:720) basins(:,1:360)];
    lon_basin_new=lon_basin;%+180;
    
    % Mask anything outside chosen basin
    basin_temp = zeros(size(basin_new));
    for i_basin =1:size(basinID,2)
        basin_temp(find(basin_new==basinID(i_basin)))=1;
    end
    basin_new(find(basin_temp~=1))=nan;
    
    % Plot Basin to check
    figure
    hold on;
    if flag_worldmap
        ax = worldmap(latlim, lonlim); shading flat;
        setm(ax, 'FFaceColor', color_ocean);
        setm(ax, 'frame', 'on', 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold');
    else
        axesm('MapProjection', map_proj_name, 'MapLonLimit', lonlim, 'MapLatLimit', ...
            latlim, 'FFaceColor', color_ocean, 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold'); shading flat;
        tightmap;
    end
    h = pcolorm(lat_basin-nanmean(diff(lat_basin))/2, ...
        lon_basin_new-nanmean(diff(lon_basin_new))/2, basin_new);
    set(h, 'LineStyle', 'none');
    geoshow([country.Lat], [country.Lon], 'Color', [.2 .2 .2]);
    %linem(lon_data_cut,lat_data_cut,'line','none','marker','square',...
    %    'markersize',3,'color','k','MarkerFaceColor','r');
    %title(['CRU TS 3.1, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
    hold off
    
    
    %m_proj(map_proj_name, 'lon', lonlim, 'lat', latlim); shading flat;
    %h=m_pcolor(lon_basin_new-nanmean(diff(lon_basin_new))/2, lat_basin-nanmean(diff(lat_basin))/2, basin_new);
    %set(h,'LineStyle','none')
    %m_patch(ocean.x,ocean.y,color_ocean,'FaceAlpha',1);
    %m_plot(country.x,country.y,'color',[.2 .2 .2]);
    %m_line(lon_data_cut,lat_data_cut,'line','none','marker','square','markersize',3,'color','k','MarkerFaceColor','r');
    %m_grid('linestyle', 'none','fontsize',16,'fontname','helvetica','FontWeight','bold');    
    %title(['CRU TS 3.1, Historical (' num2str(window_cut(1)) '-' num2str(window_cut(2)) ')'],'fontsize',16);
    %hold off   
elseif(flag_basin_shpfile == 1) % Shapefile (uses inpologon)
    if flag_wri == 1
        basin_shp_all = shaperead([folder_feature  file_name], 'UseGeoCoords', true);
        basin_shp = basin_shp_all(basinID);
    elseif flag_wri == 0
        basin_shp = shaperead([folder_feature  file_name], 'UseGeoCoords', true);
    end
    %% Plot river basin
    figure;
    hold on;
    if flag_worldmap
        ax = worldmap(latlim, lonlim); shading flat;
        setm(ax, 'FFaceColor', color_ocean);
        setm(ax, 'frame', 'on', 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold');
    else
        axesm('MapProjection', map_proj_name, 'MapLonLimit', lonlim, 'MapLatLimit', ...
            latlim, 'FFaceColor', color_ocean, 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold'); shading flat;
        tightmap;
    end
    geoshow(land, 'FaceColor', color_land);
     if i_site ~=4 && i_site ~=8 && i_site ~= 7
        %geoshow([ocean.Lat], [ocean.Lon], 'Color', color_ocean);
    end
    geoshow([country.Lat],[country.Lon],'Color',[.2 .2 .2]);
    geoshow(basin_shp, 'FaceColor', color_basin, 'FaceAlpha', 1);
    geoshow(rivers_lakes, 'Color', color_rivers);
    title(area_txt, 'fontsize', 16);
    hold off
    
elseif(flag_basin_shpfile == 2) % Latitude-Longitude Box
    
    bndry_lat = [lat_region_min; lat_region_max;lat_region_max;lat_region_min;lat_region_min];
    bndry_lon = [lon_region_min;lon_region_min;lon_region_max;lon_region_max;lon_region_min];
    
    lonlim=[lonmin lonmax];
    latlim=[latmin latmax];
    
    figure;
    hold on;
    if flag_worldmap
        ax = worldmap(latlim, lonlim); shading flat;
        setm(ax, 'FFaceColor', color_ocean);
        setm(ax, 'frame', 'on', 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold');
    else
        axesm('MapProjection', map_proj_name, 'MapLonLimit', lonlim, 'MapLatLimit', ...
            latlim, 'FFaceColor', color_ocean, 'Grid', 'off', 'FontSize', 16, ...
            'FontName', 'helvetica', 'FontWeight', 'bold'); shading flat;
        tightmap;
    end
    geoshow(land, 'FaceColor', color_land);
    %geoshow([ocean.Lat],[ocean.Lon],'Color', color_ocean);
    geoshow([country.Lat], [country.Lon], 'Color', [.2 .2 .2]);
    geoshow(rivers_lakes, 'Color', color_rivers);
    patchm(bndry_lon,bndry_lat,[0.1294,0.4431,0.7098],'FaceAlpha',0,'LineWidth',4);
    title(area_txt,'fontsize',16);
    hold off
    
end


%% Select months/years to composite
mon_comp = 1:12;

%% Create a monthly name vector
monthName = {...
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

%% CMAP Data
dataset_name='atmos'; var_name='prec'; unit_name='mm/day'; L=0:1:9;L=L*1.5; c_flag=20; title_name = 'Precipitation'; 
cmapfilename='cmap.nc';

% Create Year and Month Vectors for cmap.
mon_vect=repmat([1:12]',32,1);
yr_vect=floor((1979:(1/12):2010.98)');


% Pull out longitude and latitude
lon_data=ncread([basefolder cmapfilename],'lon');
lat_data=ncread([basefolder cmapfilename],'lat');
lon_data = [lon_data(73:144)-360;lon_data(1:72)];

%% Create Basin Mask
if flag_basin_shpfile == 0 % Gridded basin dataset
    
    % Create New Dummy Mask of ones
    basin_mask=nan(length(lat_data),length(lon_data));
    
    % find lat/lons of basin cells at original basin resolution.
    [lat_b,lon_b]=find(isnan(basin_new)==0);
    
    % Allowed offset (if zero, no allowed offset of model grid cell center
    lat_half_mod=nanmean(diff(lat_data)).*basin_off;
    lon_half_mod=nanmean(diff(lon_data)).*basin_off;
    
    for i_basin_cell=1:length(lon_b)
        % Find overlapping model grid cells
        lon_loc=find((lon_data+lon_half_mod)>=(lon_basin_new(lon_b(i_basin_cell))-.25) & (lon_data-lon_half_mod)<=(lon_basin_new(lon_b(i_basin_cell))+.25));
        lat_loc=find((lat_data+lat_half_mod)>=(lat_basin(lat_b(i_basin_cell))-.25) & (lat_data-lat_half_mod)<=(lat_basin(lat_b(i_basin_cell))+.25));
        
        if (isempty(lon_loc)==0 & isempty(lat_loc)==0 )
            
            for i_lon_fill=1:length(lon_loc)
                for i_lat_fill=1:length(lat_loc)
                    basin_mask(lat_loc(i_lat_fill),lon_loc(i_lon_fill))=1;
                end
            end
        end
    end
    
    % test
    %figure; pcolor(lon_basin_new,lat_basin,basin_new),title('original basin data'),shading flat, caxis([0 50]),colorbar;
    %figure; pcolor(lon_data,lat_data,landfrac),title('model land fraction'),colorbar
    %figure; pcolor(lon_data,lat_data,landfrac.*basin_mask),title('model land fraction, basin mask'),colorbar
elseif (flag_basin_shpfile == 1) % Shapefile of basin
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
    [IN ON] = inpolygon(lonlat(:,1),lonlat(:,2),basin_shp.Lon,basin_shp.Lat);
    lonlat_loc = IN+ON; clear IN ON;
    lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
    lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
    lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
    for i = 1:size(lonlat_loc_ind,1)
        basin_mask(find(lat_data == lat_data_cut(i)),find(lon_data == lon_data_cut(i)))=1;
    end
    
    % Figure: Basin check 1
    %figure; pcolor(lon_data,lat_data,basin_mask),title('model land fraction, basin mask'),colorbar, grid off; shading flat
 
elseif (flag_basin_shpfile == 2) % lat-lon   
%         % Create New Dummy Mask of ones
%     basin_mask=zeros(length(lat_data),length(lon_data));
%     lon_loc=find(lon_data>=lon_region_min & lon_data<=lon_region_max);
%     lat_loc=find(lat_data>=lat_region_min & lat_data<=lat_region_max);
%     
%     if (isempty(lon_loc)==0 && isempty(lat_loc)==0 )
%         for i_lon_fill=1:length(lon_loc)
%             for i_lat_fill=1:length(lat_loc)
%                 basin_mask(lat_loc(i_lat_fill),lon_loc(i_lon_fill))=1;
%             end
%         end
%     end
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
    [IN ON] = inpolygon(lonlat(:,1),lonlat(:,2),bndry_lon,bndry_lat);
    lonlat_loc = IN+ON; clear IN ON;
    lonlat_loc_ind=find(lonlat_loc>0); lonlat_loc_ind(isnan(lonlat_loc_ind))=[];
    lon_data_cut(:,1) = lonlat(lonlat_loc_ind,1);
    lat_data_cut(:,1) = lonlat(lonlat_loc_ind,2);
    for i = 1:size(lonlat_loc_ind,1)
        basin_mask(find(lat_data == lat_data_cut(i)),find(lon_data == lon_data_cut(i)))=1;
    end

end



%% TIME LOOP
% Place holders
n=0; % place holder
i_year=0; % Year place holder

for curryear=cmap_years;
    curryear
    
    i_year=i_year+1;
    
    % Month place holder
    i_mon=0;
    
    for currmonth=mon_comp
        i_mon=i_mon+1;
        n=n+1;
        
        % Identify month
        yrcut=find(mon_vect==currmonth & yr_vect==curryear);
        
        % Open netcdf data structures for these files
        cmap_varold=ncread([basefolder cmapfilename],var_name,...
            [1 1 1 yrcut],[ 144 72 1 1]);
        
        % Replace missing values with NaNs
        missing_value = ncreadatt([basefolder cmapfilename],var_name,'_FillValue');
        cmap_varold(find(cmap_varold==missing_value))=NaN;
        
        % Realign and scale
        cmap_varold = rot90(cmap_varold,3);
        cmap_var=[fliplr(cmap_varold(:,1:72)) fliplr(cmap_varold(:,73:144))];
%         figure; pcolor(lon_data,lat_data,cmap_var),title(['CMAP ' monthName{i_mon} ' ' title_name]),colorbar, grid off; shading flat;
        
        if flag_basin_shpfile <= 1 % BASIN DATA

            var_basin=cmap_var.*basin_mask;
            % Find Latitude and Longitudes for region
            [i_lat]=find(isnan(nanmean(var_basin,2))==0); lat_basin_out=lat_data(i_lat);
            [i_lon]=find(isnan(nanmean(var_basin,1))==0); lon_basin_out=lon_data(i_lon);

            % Average over basin
            [cosmean,cosgrid,datagrid,new_lon,new_lat] = coswt(var_basin,lat_data,lon_data,...
                min(lon_basin_out),max(lon_basin_out),min(lat_basin_out),max(lat_basin_out));
%             figure; pcolor(new_lon,new_lat,datagrid),title(['cmap ' monthName{i_mon} ' ' var_name]),colorbar, grid off; shading flat;
%             figure; pcolor(new_lon,new_lat,double(cosgrid)),title(['cmap ' monthName{i_mon} ' ' var_name]),colorbar, grid off; shading flat;

        elseif flag_basin_shpfile ==2 % LAT-LON BOX
            % Cosine Weighted Average
            [cosmean,cosgrid,datagrid,new_lon,new_lat] = coswt(cmap_var,lat_data,lon_data,lonmin,lonmax,latmin,latmax);
            %figure; pcolor(new_lon,new_lat,datagrid),title(['cmap ' monthName{i_mon} ' ' var_name]),colorbar, grid off; shading flat;
            %figure; pcolor(new_lon,new_lat,cosgrid),title(['cmap ' monthName{i_mon} ' ' var_name]),colorbar, grid off; shading flat;
        end
        
        cosmean = double(cosmean);
        cosgrid = double(cosgrid);
        
        VAR_cmap(i_year,i_mon)=cosmean*eomday(cmap_years(i_year),i_mon);
        VAR_cmap_grid(i_year,i_mon,:,:)=cmap_var*eomday(cmap_years(i_year),i_mon);%datagrid*eomday(cmap_years(i_year),i_mon);
        
    end % month loop

end     % year loop


%% Save data to matlab file
cmapdata.('lon_basin')=lon_data;%lon_basin_out;
cmapdata.('lat_basin')=lat_data;%lat_basin_out;
cmapdata.('yrs')=cmap_years;
cmapdata.('months')=monthName{mon_comp};
cmapdata.('basin_mean')=VAR_cmap; % average across the watershed
cmapdata.('basin_grid')=VAR_cmap_grid; % spatial data, limited to the watershed
save([basefolder var_name '.' basin_name '.cmap.mat'],'cmapdata');

%% Plot
figure
hold on
plot(mon_comp,mean(VAR_cmap),'Color','b','LineStyle','-','LineWidth',3,'Marker','none');
set(get(gcf,'CurrentAxes'),'FontSize',18)
xlim([1 12])
xlabel('month','FontSize',18);
ylabel([var_name ' (' unit_name ')'],'FontSize',18);
title(['Monthly ' title_name ' , ' area_txt],'FontSize',18)
%set(get(gcf,'CurrentAxes'),'FontSize',10)
legend(['cmap 3.1'],'Location','NorthEast');
box on
hold off

set(gcf,'Renderer','painters')
set(gcf,'OuterPosition',[292   999   845   570])
set(gcf,'PaperPositionMode','auto')

end

% The CPC Merged Analysis of Precipitation ("CMAP") is a technique which
% produces pentad and monthly analyses of global precipitation in which
% observations from raingauges are merged with precipitation estimates from
% several satellite-based algorithms (infrared and microwave). The analyses
% are are on a 2.5 x 2.5 degree latitude/longitude grid and extend back to
% 1979. These data are comparable (but should not be confused with)
% similarly combined analyses by the Global Precipitation Climatology
% Project which are described in Huffman et al (1997).
% 
% It is important to note that the input data sources to make these
% analyses are not constant throughout the period of record. For example,
% SSM/I (passive microwave - scattering and emission) data became available
% in July of 1987; prior to that the only microwave-derived estimates
% available are from the MSU algorithm (Spencer 1993) which is
% emission-based thus precipitation estimates are avaialble only over
% oceanic areas. Furthermore, high temporal resolution IR data from
% geostationary satellites (every 3-hr) became available during 1986; prior
% to that, estimates from the OPI technique (Xie and Arkin 1997) are used
% based on OLR from polar orbiting satellites.
% 
% The merging technique is thoroughly described in Xie and Arkin (1997).
% Briefly, the methodology is a two-step process. First, the random error
% is reduced by linearly combining the satellite estimates using the
% maximum likelihood method, in which case the linear combination
% coefficients are inversely proportional to the square of the local random
% error of the individual data sources. Over global land areas the random
% error is defined for each time period and grid location by comparing the
% data source with the raingauge analysis over the surrounding area. Over
% oceans, the random error is defined by comparing the data sources with
% the raingauge observations over the Pacific atolls. Bias is reduced when
% the data sources are blended in the second step using the blending
% technique of Reynolds (1988). Here the data output from step 1 is used to
% define the "shape" of the precipitation field and the rain gauge data are
% used to constrain the amplitude.