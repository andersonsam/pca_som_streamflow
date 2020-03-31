# Overview of Data

Below is a description of the data provided.  Please see "Identification of local water resource vulnerability to rapid deglaciation" for a more detailed description of the processing steps and sources of these data.

AB_rgi.mat:
* GLIMS_ID: identifier assigned by the Global Land Ice Measurements from Space service at NSIDC
* glacierArea: area of each glacier in Alberta, in km^2
* glacierLat: latitude of midpoint of each glacier in Alberta, in degrees N
* glacierLon: longitude of midpoint of each glacier in Alberta, in degrees E
* glacierZ: mean elevation of each glacier in Alberta, in m

flow_vars.mat:
* all_flow: 8760 x 194 matrix; rows are days, columns are stream gauges; time series of discharge for each stream gauge for 24 years
* all_flowseason: 194 x 365 matrix; rows are stream gauges, columns are days; each row is the seasonal discharge for each stream gauge, averaged over the 24 years of observations
* all_flowseason_norm: 194 x 365 matrix; each row of all_flowseason has been normalized to have zero mean and unity variance
* flow_years: 4656 x 365 matrix; rows are station-years (194 stations x 24 years = 4656 rows), columns are days (e.g. the first 24 rows are the 24 years of discharge of the first stream gauge, rows 25 - 48 are the 24 years of discharge of the second stream gauge, etc.)
* flow_years_norm: 4656 x 194 matrix; each row of flow_years has been normalized to have zero mean and unity variance

stream_gauge_vars.mat: 
* stationDrainageArea: drainage area of each stream gauge, in km^2
* stationElevation: elevation of each stream gauge, in m
* stationID: ID of each stream gauge in HYDAT database
* stationLat: latitude of each stream gauge, in degrees N
* stationLon: longitude of each stream gauge, in degrees E
* stationName: name of each stream gauge
* stationPercentGlaciers: area fraction of glaciation in each stream gauge’s watershed (between 0 and 1), as calculated by meltwater routing through the DEM
* stationStreamOrder: stream order of rivers at each stream gauge
* stationsWithGlaciers: binary variable (1 if glaciers are present upstream of stream gauge) as calculated by meltwater routing through DEM
* stationWithoutGlaciers: binary variable (1 if glaciers are not present upstream of stream gauge) as calculated by meltwater routing through DEM

AB_city_intake_dam_vars.mat:
* cityLat: latitude of municipalities, in degrees N
* cityLon: longitude of municipalities, in degrees E
* cityName: name of municipalities
* damDrainageArea: drainage area upstream of dams, in km^2
* damElevtion: elevation of each dam, in m
* damLat: latitude of dams, in degrees N
* damLon: longitude of dams, in degrees E
* damName: name of dams
* damPercentGlaciers: area fraction of glaciation in each dam’s watershed (between 0 and 1), as calculated by meltwater routing through the DEM
* damStreamOrder: stream order of rivers at each dam
* intakeDrainageArea: drainage area upstream of intakes, in km^2
* intakeElevation: elevation of each intake, in m
* intakeLat: latitude of intakes, in degrees N
* intakeLon: longitude of intakes, in degrees E
* intakePercentGlaciers: area fraction of glaciation in each intake’s watershed (between 0 and 1), as calculated by meltwater routing through the DEM
* intakeStreamOrder: stream order of rivers at each intake
* source: string description of the water supply source for each municipality
* sourceType: numerical classification of each source (1: groundwater; 2: river; 3: lake; 4: missing (NaN); 5: combination/other)
* viaLat: latitude of via location for each municipality, in degrees N
* viaLon: longitude of via location for each municipality, in degrees E

AB_climate_vars.mat:
* EJJAmean: temporal mean of June, July, and August evaporation in Alberta, from ERA-Interim reanalysis
* PJJAmean: temporal mean total of June, July, and August precipitation in Alberta, from ERA-Interim reanalysis
* TJJAmean: temporal mean of June, July, and August temperature in Alberta, from ERA-Interim reanalysis
* Tmean: temporal mean temperature in Alberta, from ERA-Interim reanalysis
* damEJJAmean: temporal mean of June, July, and August evaporation at location of dams
* damPJJAmean: temporal mean total of June, July, and August precipitation at location of dams
* damTJJAmean: temporal mean of June, July, and August temperature at location of dams
* damTmean: temporal mean temperature at locations of dams
* damPmean: temporal mean total yearly precipitation at locations of dams
* intakeEJJAmean: temporal mean of June, July, and August evaporation at location of intakes
* intakePJJAmean: temporal mean total of June, July, and August precipitation at location of intakes
* intakeTJJAmean: temporal mean of June, July, and August temperature at location of intakes
* intakeTmean: temporal mean temperature at locations of intakes
* intakePmean: temporal mean total yearly precipitation at locations of intakes
* stationEJJAmean: temporal mean of June, July, and August evaporation at location of stream gauges
* stationPJJAmean: temporal mean total of June, July, and August precipitation at location of stream gauges
* stationTJJAmean: temporal mean of June, July, and August temperature at location of stream gauges
* stationTmean: temporal mean temperature at locations of stream gauges
* stationPmean: temporal mean total yearly precipitation at locations of stream gauges
* latERA: latitude of grid points in ERA-Interim fields
* lonERA: longitude of grid points in ERA-Interim fields

AB_glacier_fed_rivers.mat:
* flowRoute: 9 x 1 cell; each cell is latitude/longitude of a glacier-fed river
* riverName: 9 x 1 cell; each cell is the name of the glacier-fed river whose path is found in flowRoute
