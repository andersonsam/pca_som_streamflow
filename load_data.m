function data = load_data()

%load .mat files
data.rgi = load('./data/AB_rgi.mat'); %glacier data from RGI
data.flow_vars = load('./data/flow_vars.mat'); %preprocessed flow variables (e.g. normalized/non-normalized discharge for 194 rivers/24 years)
data.stream_gauge_vars = load('./data/stream_gauge_vars.mat'); %data of stream gauges (e.g. drainage area, location, etc)
data.city_intake_dam_vars = load('./data/AB_city_intake_dam_vars.mat'); %data of municipalities/dams (e.g. locations, elevations, etc)
data.climate_vars = load('./data/AB_climate_vars.mat'); %climate data for province and individual locations of dams/municipalities/stream gauges
data.glacier_fed_rivers = load('./data/AB_glacier_fed_rivers.mat'); %lat/lon of glacier-fed-rivers (already routed through DEM)

%load river names from headwaters file
filename = 'AB_River_Headwaters.xlsx';
fileDir = './data/';
d = readtable(strcat(fileDir,filename));
riverName = table2cell(d(:,1));
data.glacier_fed_rivers.riverName = riverName;


