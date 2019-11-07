function [] = plot_fields_biogem_2d(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PKMAX,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
% plot_fields_biogem_2d
%
%   ***********************************************************************
%   *** biogem 2-D (LON-LAT) DATA PLOTTING ********************************
%   ***********************************************************************
%
%   plot_fields_biogem_2d(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PKMAX,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   plots the BIOGEM 2-D netCDF data file:
%   'fields_biogem_2d.nc'
%
%   plot_fields_biogem_2d(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PKMAX,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   takes 11 arguments:
%   PATH .... is the relative path to the experiments
%             THIS PARAMETER MUST BE PASSED AS A STRING:
%             -> e.g., 'genie_output'
%   EXP1 .... is the (first) experiment name
%             THIS PARAMETER MUST BE PASSED AS A STRING:
%             -> e.g., 'preindustrial_spinup'
%   EXP2 .... is the experiment name of 2nd, optional, netCDF file
%             THIS PARAMETER MUST BE PASSED AS A STRING:
%             -> e.g., 'historical_transient'
%             -> leave EXP2 blank, i.e., '', for no second netCDF
%             file
%   PT1 ..... is the (first) time-slice year
%             -> all valid time-slice years will be listed if an invalid
%                year is entered
%   PT2 ..... is the year of the 2nd, optional, time-slice
%             -> set PT2 to -1 for no second time-slice
%   PVAR1 ... id the name of 1st variable to be plotted
%             -> all valid valiable names will be listed if an invalid name
%                is entered
%             THIS PARAMETER MUST BE PASSED AS A STRING;
%             -> e.g., 'ocn_PO4'
%   PVAR2 ... id the name of 2nd, optional, variable
%   MASK .... is the filename containing the mask
%             -> the mask file must live in working directory
%                (i.e., the same directory as the .m file)
%             THIS PARAMETER MUST BE PASSED AS A STRING
%             -> leave MASK blank, i.e., '', for no mask
%   PKMAX ... is the MAXIMUM (K) level in the ocean model to be plotted
%             -> passing no value will result in the full grid being
%             plotted by default
%   CSCALE .. is the scale factor for the plot
%             -> e.g., to plot in units of micro molar (umol kg-1), enter: 1e-6
%                      to plot in units of PgC, enter: 1e15
%                      to plot negative values, enter: -1
%             -> the plot is auto-scaled if a value of zero (0.0) is entered
%   CMIN .... is the minimum scale value
%   CMAX .... is the maxumum scale value
%   CN ...... is the number of (contor) intervals between minimum and maximum
%             scale values
%   PTIT .... is the string for the alternative plot title
%             THIS PARAMETER MUST BE PASSED AS A STRING:
%             -> e.g., 'distribution of bottom-water phosphate concentrations'
%             -> if an empty (i.e., '') value is passed to this parameter
%                then a title is automaticaly generated
%   PDATA ... is the filename containing any overlay data set,
%             which must be formatted as seperated columns of:
%             lon, lat, value
%             -> the full filename must be give, including any extensions
%             THIS PARAMETER MUST BE PASSED AS A STRING
%             -> leave PDATA blank, i.e., '', for no overlay data
%
%   EXAMPLE;
%           plot_fields_biogem_2d('genie_output','experiment_1','',1994.5,-1,'ocn_sur_PO4','','',14,1e-6,0,2,20,'PO_{4}','')
%           will plot the time-slice cenetered on a time of 1994.5,
%           of bottom-water [PO4] in units of umol kg-1,
%           between 0 and 2 umol kg-1 with 20 contour intervals
%           and omitting levels greater than 14 (i.e., 15 and 16 for a
%           16-level ocean model configuration)
%
%   Edit the 'm' file to change other user settings;
%           lon_min = -180;        % STARTING LONGITUDE FOR X-AXIS
%           delta_lon = 60;        % INCREMENT OF LONGITUDE ON X-AXIS
%           contour_plot = 'n';    % OVERLAY CONTOL PLOT?
%           contour_mod = 1;       % NUMBER OF COLOR INTERVALS PER CONTOR
%           contour_mod_label = 5; % NUMBER OF LABELED CONTOURS PER CONTOUR
%           contour_label = 'n';   % LABEL CONTOURS?
%           contour_noneg = 'n';   % RESTRICT DATA PLOTTED TO > 0.0?
%           plot_log10 = 'n';      % PLOT LOG10 OF THE DATA
%           contour_zero = 'y';    % PLOT ZERO CONTOUR
%           colorbar_old = 'n';    % PLOT 'OLD' COLORBAR
%           data_offset = 0.0;     % data offset (273.15 for K -> C)
%           data_ij = 'n';         % DATA as (i,j)?
%           data_ij_mean = 'y';    % average DATA by cell? 
%           data_size = 50.0;      % SIZE OF OVERLAY DATA POINTS
%           data_anomoly = 'n';    % PLOT AS MODEL-DATA ANOMOLY ONLY?
%           data_only = 'n';       % PLOT ONLY DATA (no model values)?
%           dscrsz = 0.60;         % FRACTIONAL FIGURE WINDOW SIZE
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   10/07/05: added taylor diagram plotting
%   10/07/06: adjustments to use calc_find_ij_v100706 (in place of find_ij)
%             sorted out confusion between (lon,lat) of the data and (j,i) of the model grid ...
%             added stats save
%   10/07/16: added option for inputting (i,j) data
%   11/05/15: bug-fixed experiment differencing
%             added date stamp
%   11/05/31: cosmetic changes
%   11/07/31: added additional contour option
%             changed (default) colorbar plotting
%             changed log10 plotting behavior
%             changed differencing calculation for log10(data)
%   12/01/21: altered 'help' text
%             changed subroutine name: calc_find_ij_v100706 -> calc_find_ij
%             added algorithm to average data per cell (if requested)
%   12/01/23: minor bug-fix to internal gridding
%   12/01/24: reorganized sequence of lon axis vs overlay data processing 
%             rationalized 'user settings'
%   12/02/09: added in options for: anomoly plotting; data-only
%
%   ***********************************************************************

% \/\/\/ USER SETTINGS \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ %
% PARAMATER             % DEFAULT DESCRIPTION
lon_min = -180;         % [-180]  STARTING LONGITUDE FOR X-AXIS
delta_lon = 90;         % [  90]  INCREMENT OF LONGITUDE ON X-AXIS
contour_plot = 'y';     % [ 'y']  OVERLAY CONTOL PLOT?
contour_mod = 2;        % [   2]  NUMBER OF COLOR INTERVALS PER CONTOR
contour_mod_label = 4;  % [   4]  NUMBER OF LABELED CONTOURS PER CONTOUR
contour_label = 'y';    % [ 'y']  LABEL CONTOURS?
contour_noneg = 'n';    % [ 'n']  RESTRICT DATA PLOTTED TO > 0.0?
plot_log10 = 'n';       % [ 'n']  PLOT LOG10 OF THE DATA
contour_zero = 'y';     % [ 'y']  PLOT ZERO CONTOUR
colorbar_old = 'n';     % [ 'n']  PLOT 'OLD' COLORBAR
data_offset = 0.0;      % [ 0.0]  data offset (273.15 for K -> C)
data_ij = 'n';          % [ 'n']  DATA as (i,j)?
data_ij_mean = 'n';     % [ 'n']  average DATA by cell?
data_size = 50.0;       % [50.0]  SIZE OF OVERLAY DATA POINTS
data_anomoly = 'n';     % [ 'n']  PLOT AS MODEL-DATA ANOMOLY ONLY?
data_only = 'n';        % [ 'n']  PLOT ONLY DATA (no model values)?
dscrsz = 0.60;          % [0.60]  FRACTIONAL FIGURE WINDOW SIZE
% /\/\/\ USER SETTINGS /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ %

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% close open windows
% NOTE: don't clear variable space here ...
close all;
% set axes
lat_min = -090;
lat_max = +090;
lon_max = lon_min+360;
lon_offset = 0;
% set passed parameters
path = PATH;
exp_1 = EXP1;
exp_2 = EXP2;
timesliceid_1 = PT1;
timesliceid_2 = PT2;
dataid_1 = PVAR1;
dataid_2 = PVAR2;
kplotmax = PKMAX;
data_scale = CSCALE;
con_min = CMIN;
con_max = CMAX;
con_n = CN;
graph_title = PTIT;
maskid = MASK;
overlaydataid = PDATA;
% set default data scaling
if data_scale == 0.0
    data_scale = 1.0;
    clear con_min;
    clear con_max;
    con_n = 10;
end
% define grey color
color_g = [0.75 0.75 0.75];
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% set function name
str_function = 'plot-fields-biogem-2d';
%
% *********************************************************************** %

% *********************************************************************** %
% *** INITIALIZE ARRAYS ************************************************* %
% *********************************************************************** %
%
xm = [];
ym = [];
zm = [];
lonm = [];
lone = [];
lonw = [];
latm = [];
latn = [];
lats = [];
laym = [];
layt = [];
layb = [];
rawdata=[];
data_1=[];
data_2=[];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA FILE ********************************************* %
% *********************************************************************** %
%
% open netCDF file
ncid_1=netcdf.open([path '/' exp_1 '/biogem/fields_biogem_2d.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'grid_level');
grid_k1(:,:) = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_k1 = grid_k1';
% load and calculate remaining grid information
% calculate grid dimensions
varid  = netcdf.inqVarID(ncid_1,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
jmax = dimlen;
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
imax = dimlen;
grid_lon = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% Non-uniform grid
lat_max = sin(pi*lat_max/180.0);
lat_min = sin(pi*lat_min/180.0);
latn = sin(pi*latn/180.0);
lats = sin(pi*lats/180.0);
%initialize dummy topography
topo = zeros(jmax,imax);
layb = zeros(jmax,imax);
%i=1 longitude grid origin
grid_lon_origin = grid_lon_edges(1);
% maximum k-level
if isempty(kplotmax)
    kplotmax = 99;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
%
% *** SET TIME-SLICE **************************************************** %
%
% check that the year exists
varid  = netcdf.inqVarID(ncid_1,'time');
timeslices = netcdf.getVar(ncid_1,varid);
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
clear time;
while exist('time','var') == 0
    for n = 1:dimlen,
        if double(int32(100*timeslices(n)))/100 == timesliceid_1
            time = timesliceid_1;
            tid = n;
        end
    end
    if exist('time','var') == 0
        disp('   > WARNING: Year must be one of the following;');
        format long g;
        double(int32(100*timeslices(:)))/100
        format;
        timesliceid_1 = input('   > Time-slice year: ');
    end
end
%
% *** SET DATA FIELD **************************************************** %
%
% check that the variable name exists
varid = [];
while isempty(varid)
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid = n;
        end
    end
    if isempty(varid)
        disp('   > WARNING: Variable name must be one of the following;');
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
            varname
        end
        dataid_1 = input('   > Variable name: ','s');
    end
end
%
% *** LOAD DATA ********************************************************* %
%
% NOTE: flip array around diagonal to give (j,i) array orientation
data_1(:,:) = zeros(jmax,imax);
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
rawdata = netcdf.getVar(ncid_1,varid);
if length(dimids) == 3
    rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax,tid);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
elseif length(dimids) == 2
    rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
else
    data_1 = NaN*data_1;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET ALTERNATIVE GRIDDED DATASET *********************************** %
% *********************************************************************** %
%
% *** ALT EXPERIMENT **************************************************** %
%
% open new netCDF file if necessary, else reuse 1st netCDF dataset
if ~isempty(exp_2)
    % open netCDF file
    ncid_2 = netcdf.open([path '/' exp_2 '/biogem/fields_biogem_2d.nc'],'nowrite');
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
else
    ncid_2 = ncid_1;
end
%
% *** SET ALT TIME-SLICE ************************************************ %
%
if timesliceid_2 > 0.0
    % check that the year exists
    varid  = netcdf.inqVarID(ncid_2,'time');
    timeslices = netcdf.getVar(ncid_2,varid);
    [dimname, dimlen] = netcdf.inqDim(ncid_2,varid);
    clear time;
    while exist('time','var') == 0
        for n = 1:dimlen,
            if double(int32(100*timeslices(n)))/100 == timesliceid_2
                time = timesliceid_2;
                tid = n;
            end
        end
        if exist('time','var') == 0
            disp('   > WARNING: Year must be one of the following;');
            format long g;
            double(int32(100*timeslices(:)))/100
            format;
            timesliceid_2 = input('   > Time-slice year: ');
        end
    end
end
%
% *** SET ALT DATA FIELD ************************************************ %
%
if ~isempty(dataid_2)
    % check that the variable name exists
    varid = [];
    while isempty(varid)
        for n = 0:nvars-1,
            [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
            if strcmp(varname,dataid_2)
                varid = n;
            end
        end
        if isempty(varid)
            disp('   > WARNING: Variable name must be one of the following;');
            for n = 0:nvars-1,
                [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,n);
                varname
            end
            dataid = input('   > Variable name: ','s');
        end
    end
else
    for n = 0:nvars-1,
        [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,n);
        if strcmp(varname,dataid_1)
            varid = n;
        end
    end
end
%
% *** LOAD ALT DATA ***************************************************** %
%
data_2(:,:) = zeros(jmax,imax);
if (~isempty(exp_2) || (timesliceid_2 > 0.0) || ~isempty(dataid_2)),
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid);
    if length(dimids) == 3
        rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax,tid);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_2 = NaN*data_2;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
% create an output filestring for data and plot saving
% NOTE: only incorporate a single 'dataid' string for simplicity
% *********************************************************************** %
%
if ~isempty(maskid)
    if ~isempty(exp_2)
        filename = [exp_1, '_MINUS_' exp_2, 'y', '.', num2str(timesliceid_1), '.', dataid_1, '.', maskid];
        if timesliceid_2 > 0.0
            filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1, '.', maskid];
        end
    elseif timesliceid_2 > 0.0
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1, '.', maskid];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', maskid];
    end
else
    if ~isempty(exp_2)
        filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '.', dataid_1];
        if timesliceid_2 > 0.0
            filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1];
        end
    elseif timesliceid_2 > 0.0
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1];
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% set ocean grid value to give white when plotted
if strcmp(dataid_1,'grid_mask')
    data_1 = NaN;
end
if strcmp(dataid_2,'grid_mask')
    data_2 = NaN;
end
%
xm = lonm;
ym = latm;
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) > 90
            data_1(j,i) = NaN;
            data_2(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif grid_k1(j,i) > kplotmax
            data_1(j,i) = NaN;
            data_2(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        elseif ((data_1(j,i) < -1.0E6) || (data_1(j,i) > 1.0E30) || (data_2(j,i) < -1.0E6) || (data_2(j,i) > 1.0E30))
            data_1(j,i) = NaN;
            data_2(j,i) = NaN;
            topo(j,i) = +1.0;
            layb(j,i) = -1.0;
        else
            if plot_log10 == 'y'
                if ((data_1(j,i) > 0.0) && (data_2(j,i) > 0.0))
                    data_1(j,i) = log10(data_1(j,i))/data_scale;
                    data_2(j,i) = log10(data_2(j,i))/data_scale;
                else
                    data_1(j,i) = NaN;
                    data_2(j,i) = NaN;
                end
            else
                data_1(j,i) = data_1(j,i)/data_scale;
                data_2(j,i) = data_2(j,i)/data_scale;
            end
            topo(j,i) = -1.0;
            layb(j,i) = +1.0;
            if contour_noneg == 'y'
                if ((data_1(j,i) - data_2(j,i) - data_offset) < 0.0)
                    data_1(j,i) = 0.0;
                    data_2(j,i) = 0.0;
                end
            end
        end
    end
end
data = data_1 - data_2;
data = data - data_offset;
zm = data;
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % load overlay datafile
    overlaydatafile = [overlaydataid];
    overlaydata_raw = load(overlaydatafile,'-ascii');
    overlaydata_size = size(overlaydata_raw(:,:));
    nmax=overlaydata_size(1);
    % create (i,j) from (lon,lat) and vice versa (depending on data input type)
    if (data_ij == 'n')
        % convert (lon,lat) overlay data to (i,j)
        % NOTE: function 'calc_find_ij' takes input in order: (lon,lat)
        %       i.e., the same as the raw overlay data, which is (lon,lat) (i.e., (i,j)) format
        % NOTE: !!! gridded data is (j,i) !!!
        overlaydata_ij(:,:) = zeros(size(overlaydata_raw(:,:)));
        for n = 1:nmax,
            overlaydata_ij(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
        end
        overlaydata_ij(:,3) = overlaydata_raw(:,3);
    else
        % convert (i,j) overlay data to (lon,lat)
        % NOTE: save (i,j) data first
        overlaydata_ij(:,:) = overlaydata_raw(:,:);
        overlaydata_raw(:,1) = grid_lon_origin + 360.0*(overlaydata_raw(:,1) - 0.5)/jmax;
        overlaydata_raw(:,2) = 180.0*asin(2.0*(overlaydata_raw(:,2) - 0.5)/jmax - 1.0)/pi;
    end
    % remove data in land cells
    for n = 1:nmax,
        if isnan(zm(overlaydata_ij(n,2),overlaydata_ij(n,1)))
            overlaydata_raw(n,3) = NaN;
            overlaydata_ij(n,3)  = NaN;
        end
    end
    % convert lat to sin(lat) for plotting
    overlaydata(:,:) = overlaydata_raw(:,:);
    overlaydata(:,2) = sin(pi*overlaydata_raw(:,2)/180.0);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ij_mean == 'y')
        overlaydata_ij_old(:,:) = overlaydata_ij(:,:);
        overlaydata_ij(:,:) = [];
        overlaydata(:,:)    = [];
        m=0;
        for i = 1:imax,
            for j = 1:jmax,
                if (~isnan(zm(j,i)))
                    samecell_locations = find((overlaydata_ij_old(:,1)==i)&(overlaydata_ij_old(:,2)==j));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ij_old(samecell_locations,3));
                        overlaydata_ij(m,:) = [i j samecell_mean];
                        overlaydata(m,1) = grid_lon_origin + 360.0*(overlaydata_ij(m,1) - 0.5)/jmax;
                        overlaydata(m,2) = 2.0*(overlaydata_ij(m,2) - 0.5)/jmax - 1.0;
                        overlaydata(m,3) = samecell_mean;
                    end
                end
            end
        end
        nmax=m;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TAYLOR DIAGRAM **************************************************** %
% *********************************************************************** %
% calculate stats needed for Taylor Diagram (and plot it!)
% *********************************************************************** %
%
% *** 3D (GRIDDED) DATA ************************************************* %
%
if (~isempty(dataid_2))
    % transform data sets in vectors
    data_vector_1 = reshape(data_1(:,:),imax*jmax,1);
    data_vector_2 = reshape(data_2(:,:),imax*jmax,1);
    % filter data
    data_vector_1(find(data_vector_1(:) < -1.0E6)) = NaN;
    data_vector_1(find(data_vector_1(:) > 0.9E36)) = NaN;
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    % calculate stats
    % NOTE: STATM = allstats(Cr,Cf)
    % 	    STATM(1,:) => Mean
    % 	    STATM(2,:) => Standard Deviation (scaled by N)
    % 	    STATM(3,:) => Centered Root Mean Square Difference (scaled by N)
    % 	    STATM(4,:) => Correlation
    %       STATM(5,:) => N
    STATM = calc_allstats(data_vector_1,data_vector_2);
    % plot Taylor diagrams
    % plot Taylor diagrams
    figure;
    plot_taylor(STATM(2,2),STATM(2,1),STATM(4,2),2,'model',filename,1.5);
    print('-f1','-depsc2','-painters', [filename, '_Taylor.', str_date, '.eps']);
    figure;
    taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
    print('-depsc2', [filename, '_TaylorDiagram.', str_date, '.eps']);
end
%
% *** DISCRETE DATA ***************************************************** %
%
% NOTE: no scale transformatoin has been appplied
%       to either gridded or % overlay data
% NOTE: valid only for data on a single depth level
if ~isempty(overlaydataid)
    % set overlay data vector
    data_vector_1 = overlaydata(:,3);
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! data is (j,i) !!! (=> swap i and j)
    for n = 1:nmax,
        data_vector_2(n) = data(overlaydata_ij(n,2),overlaydata_ij(n,1));
    end
    % filter data
    data_vector_2(find(data_vector_2(:) < -1.0E6)) = NaN;
    data_vector_2(find(data_vector_2(:) > 0.9E36)) = NaN;
    % calculate stats
    STATM = calc_allstats(data_vector_1,data_vector_2);
    % plot Taylor diagrams
    figure;
    plot_taylor(STATM(2,2),STATM(2,1),STATM(4,2),2,'model',filename,1.5);
    print('-f1','-depsc2','-painters', [filename, '_Taylor.', str_date, '.eps']);
    figure;
    taylordiag_vargout = plot_taylordiag(STATM(2,1:2),STATM(3,1:2),STATM(4,1:2));
    print('-depsc2', [filename, '_TaylorDiagram.', str_date, '.eps']);
end
%
% *** SAVE STATS DATA *************************************************** %
%
if (~isempty(dataid_2) | ~isempty(overlaydataid))
    fid = fopen([filename '_STATS.' str_date '.dat'], 'wt');
    fprintf(fid, 'Stats summary: reference data');
    fprintf(fid, '\n');
    fprintf(fid, 'Mean                                               : %8.6e \n', STATM(1,1));
    fprintf(fid, 'Standard Deviation (scaled by N)                   : %8.6e \n', STATM(2,1));
    fprintf(fid, 'Centered Root Mean Square Difference (scaled by N) : %8.6e \n', STATM(3,1));
    fprintf(fid, 'Correlation                                        : %8.6e \n', STATM(4,1));
    fprintf(fid, '\n');
    fprintf(fid, 'Stats summary: model data');
    fprintf(fid, '\n');
    fprintf(fid, 'Mean                                               : %8.6e \n', STATM(1,2));
    fprintf(fid, 'Standard Deviation (scaled by N)                   : %8.6e \n', STATM(2,2));
    fprintf(fid, 'Centered Root Mean Square Difference (scaled by N) : %8.6e \n', STATM(3,2));
    fprintf(fid, 'Correlation                                        : %8.6e \n', STATM(4,2));
    fclose(fid);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** ANOMOLY PLOTTING DATA ADJUSTMENTS ********************************* %
% *********************************************************************** %
% 
if ~isempty(overlaydataid)
    % calculate molde-data anomoly
    if (data_anomoly == 'y')
        overlaydata(:,3) = data_vector_2(:) - data_vector_1(:);
    end
    % redefine model grid location values so as to all appear white
    if (data_only == 'y')
        zm = zeros(size(zm(:,:)));
        zm(find(zm(:,:) == 0)) = NaN;
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** TRANSFORM LON GRID ************************************************ %
% *********************************************************************** %
%
% extend gridded data in +/- longitude
xm_ex = [xm - 360.0 xm + 000.0 xm + 360.0];
ym_ex = [ym + 000.0 ym + 000.0 ym + 000.0];
zm_ex = [zm zm zm];
topo_ex = [topo - 360.0 topo + 000.0 topo + 360.0];
lonm_ex = [lonm - 360.0 lonm + 000.0 lonm + 360.0];
lone_ex = [lone - 360.0 lone + 000.0 lone + 360.0];
lonw_ex = [lonw - 360.0 lonw + 000.0 lonw + 360.0];
layb_ex = [layb - 360.0 layb + 000.0 layb + 360.0];
% shorten to conform to desired lon start value
lon_start = min(min(lonw));
i_start = round((lon_min-(lon_start-360.0))/(360.0/jmax)) + 1;
xm = xm_ex(:,i_start:i_start+imax-1);
ym = ym_ex(:,i_start:i_start+imax-1);
zm = zm_ex(:,i_start:i_start+imax-1);
topo = topo_ex(:,i_start:i_start+imax-1);
lonm = lonm_ex(:,i_start:i_start+imax-1);
lone = lone_ex(:,i_start:i_start+imax-1);
lonw = lonw_ex(:,i_start:i_start+imax-1);
layb = layb_ex(:,i_start:i_start+imax-1);
if ~isempty(overlaydataid)
    % force discrete data to lie within longitude plotting axis
    % (lon_min to lon_min + 360)
    for n = 1:nmax,
        if (overlaydata(n,1) < lon_min)
            overlaydata(n,1) = overlaydata(n,1) + 360.0;
        end
        if (overlaydata(n,1) > (lon_min + 360.0))
            overlaydata(n,1) = overlaydata(n,1) - 360.0;
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** PLOT MAIN FIGURE ************************************************** %
% *********************************************************************** %
%
% *** CONFIGURE AND CREATE PLOTTING WINDOW ****************************** %
%
% create figure
scrsz = get(0,'ScreenSize');
figure('Position',[((1 - dscrsz)/2)*dscrsz*scrsz(3) (1 - dscrsz)*dscrsz*scrsz(4) dscrsz*scrsz(3) 0.60*dscrsz*scrsz(4)])
clf;
% define plotting regions
fh(1) = axes('Position',[0 0 1 1],'Visible','off');
fh(2) = axes('Position',[0.10 0.05 0.65 0.90]);
fh(3) = axes('Position',[0.80 0.27 0.20 0.46],'Visible','off');
% define colormap
cmap = colormap(jet((2*(con_n+1))+1));
%cmap = colormap(hot((2*(con_n+1))+1));
% date-stamp plot
set(gcf,'CurrentAxes',fh(1));
text(0.95,0.50,[str_function, ' / ', 'on: ', str_date],'FontName','Arial','FontSize',8,'Rotation',90.0,'HorizontalAlignment','center','VerticalAlignment','top');
%
% *** SET PLOT SCALE **************************************************** %
%
% set minimum contour value
if exist('con_min','var') == 0
    con_min = min(min(zm));
end
% set maximum contour value
if exist('con_max','var') == 0
    con_max = max(max(zm));
end
% ensure min and max are not identical ...
if con_min == con_max
    if con_max == 0.0
        con_max = 1.0;
    else
        con_min = (1.00/1.01)*con_min;
        con_max = (1.01)*con_max;
    end
end
% if min > max, then reverse min and max
if con_min > con_max
    con_min_TEMP = con_min;
    con_max_TEMP = con_max;
    con_min = con_max_TEMP;
    con_max = con_min_TEMP;
end
%
% *** CREATE MAIN PLOT ************************************************** %
%
set(gcf,'CurrentAxes',fh(2));
hold on;
% set color and lat/lon axes and labels
caxis([con_min-(con_max-con_min)/con_n con_max]);
set(gca,'PlotBoxAspectRatio',[1.0 0.5 1.0]);
axis([lon_min lon_max lat_min lat_max]);
set(gca,'TickDir','out');
set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',lon_min:delta_lon:lon_max);
set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
if ~isempty(graph_title)
    title(graph_title,'FontSize',18);
else
    title(['Year: ',strrep(num2str(time),'_',' '),' / ','Data ID: ',strrep(dataid_1,'_',' ')],'FontSize',12);
end
% draw filled rectangles
for i = 1:imax,
    for j = 1:jmax,
        if topo(j,i) > layb(j,i)
            h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(zm(j,i)))
                h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            else
                col = round((1/2)+con_n*(zm(j,i)-con_min)/(con_max-con_min));
                if col < 1
                    col = 1;
                elseif col > con_n
                    col = 2*(con_n+1)+1;
                else
                    col = 2*col+1;
                end
                h = patch([lonw(j,i) lonw(j,i) lone(j,i) lone(j,i)],[lats(j,i) latn(j,i) latn(j,i) lats(j,i)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
%
% *** PLOT CONTINENTAL OUTLINE ****************************************** %
%
% draw continental outline
for j = 1:jmax,
    for i = 1:imax-1,
        if topo(j,i) > layb(j,i)
            if topo(j,i+1) <= layb(j,i+1)
                h = plot([lone(j,i) lone(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for i = 2:imax,
        if topo(j,i) > layb(j,i)
            if topo(j,i-1) <= layb(j,i-1)
                h = plot([lonw(j,i) lonw(j,i)],[lats(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
for i = 1:imax,
    for j = 1:jmax-1,
        if topo(j,i) > layb(j,i)
            if topo(j+1,i) <= layb(j+1,i)
                h = plot([lonw(j,i) lone(j,i)],[latn(j,i) latn(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for j = 2:jmax,
        if topo(j,i) > layb(j,i)
            if topo(j-1,i) <= layb(j-1,i)
                h = plot([lonw(j,i) lone(j,i)],[lats(j,i) lats(j,i)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
%
% *** OVERLAY CONTOURS ************************************************** %
%
% plot contours
if (contour_plot == 'y') && (data_only == 'n'),
    if ((con_min < 0.0) && (con_max > 0.0)),
        v = [con_min:(con_max-con_min)/(con_n/contour_mod):0.0];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
        set(h,'LineWidth',0.33);
        v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):0.0];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
        set(h,'LineWidth',0.66);
        if contour_label == 'y', clabel(C,h); end
        v = [0.0:(con_max-con_min)/(con_n/contour_mod):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
        set(h,'LineWidth',0.33);
        v = [0.0:(con_max-con_min)/(con_n/contour_mod_label):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
        set(h,'LineWidth',0.66);
        if contour_label == 'y', clabel(C,h); end
    elseif (con_min < 0.0),
        v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
        set(h,'LineWidth',0.33);
        v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-.');
        set(h,'LineWidth',0.66);
        if contour_label == 'y', clabel(C,h); end
    else
        v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
        set(h,'LineWidth',0.33);
        v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
        set(h,'LineWidth',0.66);
        if contour_label == 'y', clabel(C,h); end
    end
    if contour_zero == 'y',
        v = [-1.0e19:1.0e19:1.0E19];
        [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
        set(h,'LineWidth',1.0);
        if contour_label == 'y', clabel(C,h); end
    end
end
%
% *** OVERLAY DATA ****************************************************** %
%
% plot overlay data
if ~isempty(overlaydataid)
    scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata(:,3)/data_scale,'o','Filled','Sizedata',data_size,'MarkerEdgeColor','k');
end
%
% *** PLOT BORDER ******************************************************* %
%
% draw plot border
h = plot([lon_min lon_max],[lat_min lat_min],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_max],[lat_max lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_min lon_min],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lon_max lon_max],[lat_min lat_max],'k-');
set(h,'LineWidth',1.0);
%
hold off;
%
% *** CREATE COLOR BAR ************************************************** %
%
set(gcf,'CurrentAxes',fh(3));
hold on;
%
set(gca,'XTick',[],'YTick',[]);
axis([0 1 0 con_n+2]);
% draw and label color bar rectangles
c = 1;
h = fill([0.1 0.2 0.3],[c c-1.0 c],cmap(2*c-1,:));
if (colorbar_old == 'y'),
    str = ['below ',num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
else
    str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
end
textsize = 2+round(80/con_n);
if textsize > 10,
    textsize = 10;
end
if (colorbar_old == 'y'),
    text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
else
    text(0.40,c,str,'FontName','Arial','FontSize',textsize);
end
set(h,'LineWidth',0.5);
set(h,'EdgeColor','k');
for c = 2:con_n+1,
    h = fill([0.1 0.1 0.3 0.3],[c-1.0 c c c-1.0],cmap(2*c-1,:));
    if (colorbar_old == 'y'),
        str = [num2str(con_min + (c-2)*(con_max-con_min)/con_n),' - ',num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
    else
        str = [num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
    end
    textsize = 2+round(80/con_n);
    if textsize > 10,
        textsize = 10;
    end
    if (colorbar_old == 'y'),
        text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
    else
        text(0.40,c,str,'FontName','Arial','FontSize',textsize);
    end
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
end
c = con_n+2;
h = fill([0.1 0.2 0.3],[c-1.0 c c-1.0],cmap(2*c-1,:));
if (colorbar_old == 'y'),
    str = ['above ',num2str(con_min + (c-2)*(con_max-con_min)/con_n)];
    textsize = 2+round(80/con_n);
    if textsize > 10,
        textsize = 10;
    end
    text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
end
set(h,'LineWidth',0.5);
set(h,'EdgeColor','k');
%
hold off;
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
print('-dpsc2', [filename, '.', str_date, '.ps']);

% *** END *************************************************************** %
% close netCDF files
% *********************************************************************** %
%
netcdf.close(ncid_1);
