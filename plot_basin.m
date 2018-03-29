function plot_basin(grid,grid_lims,cmap,cmap_label,lat,lon,latlim,lonlim,map_proj_name,country,basin_shp,flag_worldmap,flag_world)
% plots basin data
% grid=basin data grid, grid_lims=[min max] for coloring, cmap=cmap,
% latlim=[latmin latmax], lonlim=[lonmin lonmax], lat=landfrac lat,
% lon=landfrac lon, map_proj_name=map proj from matlab's mapping toolbox,
% country=world countries shp file, basin_shp=basin shp file, 
% flag_worldmap=bool, use matlab worldmap function or use custom map_proj,
% flag_world=bool, map of whole world or basin (mostly for testing purposes)

figure;
hold on;
colormap(cmap)
if flag_worldmap
    if flag_world
        ax=worldmap('World'); shading flat;
    else
        ax=worldmap(latlim, lonlim); shading flat;
    end
    setm(ax, 'frame', 'on', 'Grid', 'off', ...
        'FontName', 'helvetica', 'FontWeight', 'bold');
else
    axesm('MapProjection', map_proj_name, 'MapLonLimit', lonlim, ...
        'MapLatLimit', latlim, 'Grid', 'off', ...
        'FontName', 'helvetica', 'FontWeight', 'bold');
    tightmap; shading flat;
end
hold on;
h=pcolorm(lat, lon, grid);
set(h, 'LineStyle', 'none');
caxis(grid_lims);
geoshow([country.Lat], [country.Lon], 'Color', 'black');
patchm([basin_shp.Lat], [basin_shp.Lon], 'FaceColor', 'none', 'LineStyle', '-', 'LineWidth', 2);
h=colorbar('FontName','Times',...
    'FontWeight','bold',...
    'fontname','helvetica');
ylabel(h,cmap_label);
