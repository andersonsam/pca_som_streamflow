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
* Canada_Borders
* Rivers_Lakes
* AB_city_intake_dam_vars.mat
* AB_climate_vars.mat
* AB_glacier_fed_rivers.mat
* AB_rgi.mat
* AB_topo_lowres.mat
* flow_vars.mat
* stream_gauge_vars.mat

## Requirements

The code is provided in MATLAB.  Self-organizing maps are created using the [SOM Toolbox](http://www.cis.hut.fi/somtoolbox/) (Vesanto et al., 2000).
