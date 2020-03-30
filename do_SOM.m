function [sM, sT, sMap, bmus] = do_SOM(data, nx_som, ny_som)

% NOTE: Based off of scipts by Valentina Radic.
%
% Syntax: [sM, sT, sMap, bmus] = do_SOM(data, nx_som, ny_som)
% 
% sM: This is the self-organized map.  Note: sM.codebook contains the cluster patterns.
% sMap: This is the original map (before training) -- used in some optional figures below.  
% bmus: This is a vector containing the number of the best matching unit for each row in 'data'.  The pattern of bmu 'x' is sM.codebook(x,:).
% data: matrix of doubles; rows: stations; columns: 'time series'.
% nx_som: number of columns in SOM.
% ny_som: number of rows in SOM.

en = nx_som*ny_som;
msize=[ny_som nx_som];

% performing linear initialization of nodes (i.e. the nodes are distibuted in the space of
% first two eignevectors from PCA

display('SOM Initialization...')
data(isnan(data))=0; %SOM can't deal with NaN's -- set them to zero
sMap=som_lininit(data,'msize',msize,'hexa','sheet');
% sMap=som_randinit(data,'msize',msize,'hexa','sheet');

% training SOM
display('Training SOM...')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[2 1],'trainlen',200);

% calulating error
[q,t]=som_quality(sM,data);

% calulating hits (frequencies) of occurences of each pattern, for each seasn
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

bmus=som_bmus(sM,data,1);