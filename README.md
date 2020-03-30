# Identification of local water resource vulnerability to rapid deglaciation

This repository contains data/code required to reproduce the results and figures from the paper "Identification of local water resource vulnerability to rapid deglaciation".  It is structured in the following way:

Code: 
* main.m: loads data, performs PCA + SOM + regression, produces figures
* do_SOM.m: function which creates the SOM
* get_numberOfPatterns.m: function which determines the size of SOM to use
* load_data.m: function which loads the required data structure
* plot_AB.m: function which plots the border of Alberta, for making figures
* plot_rivers.m: function which plots rivers in Canada, for making figures
* plot_AB_glaciers.m: function which plots glaciers in Alberta, for making figures

Data:
* Canada_Borders: shapefiles of Canadian provinces
* Rivers_Lakes: shapefiles of Canadian rivers, from HYDAT database
* AB_city_intake_dam_vars.mat: data of municipalities/dams (e.g. locations, elevations, etc)
* AB_climate_vars.mat: climate data for Alberta and at specific locations of dams/communities/stream gauges
* AB_glacier_fed_rivers.mat: lat/lon of glacier-fed-rivers (already routed through DEM), and names of glacier-fed rivers
* AB_rgi.mat: glacier data from Randolph Glacier Inventory
* AB_topo_lowres.mat: low-resolution topography (processed from Shuttle Radar Topography Mission; lowres as just used for plotting)
* flow_vars.mat: processed flow variables from HYDAT database (e.g. normalized/non-normalized seasonal/yearly discharge for 194 stream gauges)
* stream_gauge_vars.mat: data of stream gauges (e.g. drainage area, location, etc.)

## Requirements

The code is provided in MATLAB.  Self-organizing maps are created using the [SOM Toolbox](http://www.cis.hut.fi/somtoolbox/) (Vesanto et al., 2000).
