function [] = plot_fields_biogem_3d_k(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PK,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
% plot_fields_biogem_3d_k
%
%   *******************************************************************   %
%   *** biogem k-SECTION (LON-LAT) DIFFERENCE PLOTTING ****************   %
%   *******************************************************************   %
%
%   plot_fields_biogem_3d_k(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PK,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   plots a k-section through the BIOGEM 3-D netCDF data file;
%   'fields_biogem_3d.nc'
%   Additional difference plotting capabilities (experiment or time-slice).
%   Also water column integral capabilities.
%
%   plot_fields_biogem_3d_k(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PK,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   takes 14 arguments:
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
%   PT1 ..... is the year of the 1st time-slice
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
%   PK ...... is the level in the ocean model to be plotted (the 'k' slice)
%             -> a zero will result in a water column integral being
%             plotted
%   CSCALE .. is the scale factor for the plot
%             -> e.g., to plot micro molar (umol kg-1), enter; 1e-6
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
%   PDATA ... is the filename containing the an overlay data set,
%             which must be formatted as seperated columns of:
%             lon, lat, value
%             -> the full filename must be give, including any extensions
%             THIS PARAMETER MUST BE PASSED AS A STRING
%             -> leave PDATA blank, i.e., '', for no overlay data
%
%   EXAMPLE;
%   plot_fields_biogem_3d_k(PATH,EXP1,EXP2,PT1,PT2,PVAR1,PVAR2,MASK,PK,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%           plot_fields_biogem_3d_k('cgenie_output','experiment_1','',1994.5,-1,'ocn_PO4','','',8,1e-6,0,2,20,'PO_{4}','')
%           will plot the time-slice cenetered on a time of 1994.5,
%           of PO4 concentrations at the ocean surface (k = 8),
%           between 0 and 2 umol kg-1 in 20 contour intervals
%           and with a title of 'PO4'
%           experiment is called 'experiment_1'
%           and lives in the 'cgenie_output' subdirectory
%
%   Edit the 'm' file to change other parameters ('USER SETTINGS')
%           lon_min = -260;        % STARTING LONGITUDE FOR X-AXIS
%           delta_lon = 60;        % INCREMENT OF LONGITUDE ON X-AXIS
%           contour_plot = 'n';    % OVERLAY CONTOL PLOT?
%           contour_mod = 1;       % NUMBER OF COLOR INTERVALS PER CONTOR
%           contour_mod_label = 5; % NUMBER OF LABELED CONTOURS PER CONTOUR
%           contour_label = 'n';   % LABEL CONTOURS?
%           contour_noneg = 'n';   % RESTRICT DATA PLOTTED TO > 0.0?
%           plot_log10 = 'n';      % PLOT LOG10 OF THE DATA
%           data_offset = 0.0;     % data offset (273.15 for K -> C)
%           data_ijk = 'n';        % [ 'n']  DATA as (i,j,k)?
%           data_ijk_mean = 'y';   % average DATA by cell?
%           data_size = 50.0;      % SIZE OF OVERLAY DATA POINTS
%           data_anomoly = 'n';    % PLOT AS MODEL-DATA ANOMOLY ONLY?
%           data_only = 'n';       % PLOT ONLY DATA (no model values)?
%           data_uv = 'y';         % overlay (u,v) velocity data?
%           data_uv_scale = 1.0;   % scaling factor for vector length
%           dscrsz = 0.60;         % FRACTIONAL FIGURE WINDOW SIZE
%
%   *******************************************************************   %
%   *** HISTORY *******************************************************   %
%   *******************************************************************   %
%
%   10/06/19: CREATED ...
%   10/07/01: Passed reference data stats to plot_taylordiag
%   10/07/05: cosmetic changes ...
%   10/07/05: name changed: allstats -> calc_allstats
%   10/07/06: adjustments to use calc_find_ij_v100706 (in place of find_ij)
%   10/07/06: sorted out confusion between (lon,lat) of the data and (j,i) of the model grid ...
%   10/07/06: added stats save
%   10/07/08: filtered overlay data to avoid plotting on land!
%   10/07/16: added option for inputting (i,j) data
%   11/01/30: re-formatting
%             addition of 3D data Taylor Diagram analysis
%   11/01/31: 3D data Taylor Diagram analysis debugging ...
%   11/05/31: Added time-stamping to Taylor plots
%   11/05/31: cosmetic changes
%   12/01/21: changed subroutine name: calc_find_ij_v100706 -> calc_find_ij
%   12/01/23: fixed bug in plot save string
%             added quiver plot option (for velocity field)
%             added cell averaging + depth layer filtering of overlay data
%   12/01/24: reorganized sequence of lon axis vs overlay data processing 
%             added parameter to control vector length
%             rationalized 'user settings'
%   12/02/09: added in options for: anomoly plotting; data-only
%   12/05/23: changed order of streamfunction netCDF loading and grid setup
%             to prevent parameter value conflicts (with data netCDF
%             parameters)
%
%   *******************************************************************   %

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
data_offset = 0.0;      % [ 0.0]  data offset (273.15 for K -> C)
data_ijk = 'n';         % [ 'n']  DATA as (i,j,k)?
data_ijk_mean = 'n';    % [ 'n']  average DATA by cell?
data_size = 50.0;       % [50.0]  SIZE OF OVERLAY DATA POINTS
data_anomoly = 'n';     % [ 'n']  PLOT AS MODEL-DATA ANOMOLY ONLY?
data_only = 'n';        % [ 'n']  PLOT ONLY DATA (no model values)?
data_uv = 'n';          % [ 'n']  overlay (u,v) velocity data?
data_uv_scale = 1.0;    % [ 1.0]  scaling factor for vector length
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
kplot = PK;
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
str_function = 'plot-fields-biogem-3d(k)';
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
% additional
rawdata=[];
zm_count = [];
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPEN netCDF DATA & LOAD (OPTIONAL) MASK FILE ********************** %
% *********************************************************************** %
%
% open netCDF file
ncid_1=netcdf.open([path '/' exp_1 '/biogem/fields_biogem_3d.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_1);
% load mask data
% NOTE: flip in j-direction to make consistent with netCDF grid
maskfile = maskid;
if ~isempty(maskid)
    mask = load(maskfile,'-ascii');
    mask = flipdim(mask,1);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% load grid data
varid  = netcdf.inqVarID(ncid_1,'grid_level');
grid_k1 = netcdf.getVar(ncid_1,varid);
% flip array around diagonal to give (j,i) array orientation
grid_k1 = grid_k1';
% calculate grid dimensions
varid  = netcdf.inqVarID(ncid_1,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
jmax = dimlen;
varid  = netcdf.inqVarID(ncid_1,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
imax = dimlen;
varid  = netcdf.inqVarID(ncid_1,'zt');
[dimname, dimlen] = netcdf.inqDim(ncid_1,varid);
kmax = dimlen;
% load remaining grid information
varid  = netcdf.inqVarID(ncid_1,'zt');
grid_zt = netcdf.getVar(ncid_1,varid);
grid_zt = flipud(grid_zt);
varid  = netcdf.inqVarID(ncid_1,'zt_edges');
grid_zt_edges = netcdf.getVar(ncid_1,varid);
grid_zt_edges = flipud(grid_zt_edges);
% calculate topography
for i = 1:imax,
    for j = 1:jmax,
        if grid_k1(j,i) <= kmax
            topo(j,i) = -grid_zt_edges(grid_k1(j,i));
        else
            topo(j,i) = 0.0;
        end
        if kplot > 0
            laym(j,i) = -grid_zt(kplot);
            layb(j,i) = -grid_zt_edges(kplot);
            layt(j,i) = -grid_zt_edges(kplot+1);
        else
            laym(j,i) = 0.0;
            layb(j,i) = -grid_zt_edges(1);
            layb(j,i) = -grid_zt_edges(kmax);
        end
    end
end
if ~isempty(maskid)
    topo = mask.*topo;
end
% load and calculate remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% calculate cell masses
% NOTE: assume equal area grid, normaalized area
data_M = zeros(kmax,jmax,imax);
for k = 1:kmax,
    data_M(k,:,:) = 1027.649*1.0*(grid_zt_edges(k) - grid_zt_edges(k+1));
end
% Non-uniform grid
lat_max = sin(pi*lat_max/180.0);
lat_min = sin(pi*lat_min/180.0);
latn = sin(pi*latn/180.0);
lats = sin(pi*lats/180.0);
%i=1 longitude grid origin
grid_lon_origin = grid_lon_edges(1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET PRIMARY GRIDDED DATASET *************************************** %
% *********************************************************************** %
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
        dataid = input('   > Variable name: ','s');
    end
end
% load data
% flip array around diagonal to give (j,i) array orientation
[varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
rawdata = netcdf.getVar(ncid_1,varid);
if length(dimids) == 4
    rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid);
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 3
    rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
    for n = 1:kmax,
        data_1(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
    end
elseif length(dimids) == 2
    rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
    data_1(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
else
    data_1(:,:,:) = NaN;
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
    ncid_2 = netcdf.open([path '/' exp_2 '/biogem/fields_biogem_3d.nc'],'nowrite');
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_2);
else
    ncid_2 = ncid_1;
end
%
% *** ALT TIME-SLICE **************************************************** %
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
% *** ALT DATA FIELD **************************************************** %
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
end
%
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (timesliceid_2 > 0.0) || (~isempty(dataid_2))
    % load data
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid);
    if length(dimids) == 4
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid);
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 3
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
        for n = 1:kmax,
            data_2(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
        data_2(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_2(:,:,:) = NaN;
    end
else
    data_2(:,:,:) = 0.0;
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** OPTIONAL (u,v) VELOCITY DATA OPERLAY ****************************** %
% *********************************************************************** %
%
if (data_uv == 'y'),
    varid  = netcdf.inqVarID(ncid_1,'phys_u');
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
    rawdata = netcdf.getVar(ncid_1,varid);
    if length(dimids) == 4
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid);
        for n = 1:kmax,
            data_u(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 3
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
        for n = 1:kmax,
            data_u(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
        data_u(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_u(:,:,:) = NaN;
    end
    varid  = netcdf.inqVarID(ncid_1,'phys_v');
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_1,varid);
    rawdata = netcdf.getVar(ncid_1,varid);
    if length(dimids) == 4
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax,tid);
        for n = 1:kmax,
            data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 3
        rawdata(1:imax,1:jmax,1:kmax) = rawdata(1:imax,1:jmax,1:kmax);
        for n = 1:kmax,
            data_v(kmax - n + 1,1:jmax,1:imax) = rawdata(1:imax,1:jmax,n)';
        end
    elseif length(dimids) == 2
        rawdata(1:imax,1:jmax) = rawdata(1:imax,1:jmax);
        data_v(1:jmax,1:imax) = rawdata(1:imax,1:jmax)';
    else
        data_v(:,:,:) = NaN;
    end
else
    data_u(:,:,:) = zeros(size(data_1));
    data_v(:,:,:) = zeros(size(data_1));
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET OUTPUT FILESTRING ********************************************* %
% *********************************************************************** %
%
% create an output filestring for data and plot saving
% NOTE: only incorporate a single 'dataid' string for simplicity
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
        filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', 'k', num2str(kplot)];
        if timesliceid_2 > 0.0
            filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1, '.', 'k', num2str(kplot)];
        end
    elseif timesliceid_2 > 0.0
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid_1, '.', 'k', num2str(kplot)];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid_1, '.', 'k', num2str(kplot)];
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
data = data_1 - data_2;
data = data - data_offset;
if kplot > 0
    % process single depth layer
    zm(:,:) = data(kplot,:,:);
    z_u(:,:) = data_u(kplot,:,:);
    z_v(:,:) = data_v(kplot,:,:);
    speed = NaN(size(zm));
    for i = 1:imax,
        for j = 1:jmax,
            if topo(j,i) > layb(j,i)
                zm(j,i) = NaN;
                xm(j,i) = NaN;
                ym(j,i) = NaN;
                z_u(j,i) = NaN;
                z_v(j,i) = NaN;
            elseif (zm(j,i) < -1.0E6) || (zm(j,i) > 1.0E30)
                zm(j,i) = NaN;
                xm(j,i) = NaN;
                ym(j,i) = NaN;
                z_u(j,i) = NaN;
                z_v(j,i) = NaN;
            else
                if plot_log10 == 'y'
                    if (zm(j,i) > 0.0)
                        zm(j,i) = log10(zm(j,i)/data_scale);
                    else
                        if contour_noneg == 'y'
                            zm(j,i) = 0.0;
                        else
                            zm(j,i) = NaN;
                        end
                    end
                else
                    zm(j,i) = zm(j,i)/data_scale;
                end
                if (data_uv == 'y'), speed(j,i) = data_scale*(z_u(j,i)^2.0 + z_v(j,i)^2.0)^0.5; end
            end
        end
    end
else
    % create water column integral
    zm = zeros(jmax,imax);
    zm_count = zeros(jmax,imax);
    if (data_uv == 'y')
        z_u(:,:) = data_u(kmax,:,:);
        z_v(:,:) = data_v(kmax,:,:);
        speed = NaN(size(zm));
    end
    for j = 1:jmax,
        for i = 1:imax,
            if topo(j,i) == 0.0
                xm(j,i) = NaN;
                ym(j,i) = NaN;
                zm(j,i) = NaN;
                z_u(j,i) = NaN;
                z_v(j,i) = NaN;
            else
                for k = 1:kmax,
                    if ((data(k,j,i) > -1.0E6) && (data(k,j,i) < 0.9E36))
                        zm(j,i) = zm(j,i) + data_M(k,j,i)*data(k,j,i);
                        zm_count(j,i) = zm_count(j,i) + 1;
                    end
                end
                if (zm_count(j,i) > 0)
                    if plot_log10 == 'y'
                        if (zm(j,i) > 0.0)
                            zm(j,i) = log10(zm(j,i)/data_scale);
                        else
                            if contour_noneg == 'y'
                                zm(j,i) = 0.0;
                            else
                                zm(j,i) = NaN;
                            end
                        end
                    else
                        zm(j,i) = zm(j,i)/data_scale;
                    end
                else
                    zm(j,i) = NaN;
                end
            end
            if (data_uv == 'y'), speed(j,i) = data_scale*(z_u(j,i)^2.0 + z_v(j,i)^2.0)^0.5; end
        end
    end
end
% copy zm before it gets transformed ...
overlaydata_zm(:,:) = zm(:,:);
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
    if (data_ijk == 'n')
        % convert (lon,lat) overlay data to (i,j)
        % NOTE: function 'calc_find_ij' takes input in order: (lon,lat)
        %       i.e., the same as the raw overlay data, which is (lon,lat) (i.e., (i,j)) format
        % NOTE: !!! gridded data is (j,i) !!!
        overlaydata_ijk(:,:) = zeros(size(overlaydata_raw(:,:)));
        for n = 1:nmax,
            overlaydata_ijk(n,1:2) = calc_find_ij(overlaydata_raw(n,1),overlaydata_raw(n,2),grid_lon_origin,imax,jmax);
            overlaydata_ijk(n,3)   = calc_find_k(overlaydata_raw(n,3),kmax);
        end
        overlaydata_ijk(:,4) = overlaydata_raw(:,4);
        % delete data lines with depth levels not equal to kplot
        wronglayer_locations = find(overlaydata_ijk(:,3)~=kplot);
        wronglayer_n = size(wronglayer_locations);
        overlaydata_ijk(wronglayer_locations,:) = [];
        overlaydata_raw(wronglayer_locations,:) = [];
        nmax = nmax-wronglayer_n(1);
    else
        % convert (i,j) overlay data to (lon,lat)
        % NOTE: save (i,j,k) data first
        overlaydata_ijk(:,:) = overlaydata_raw(:,:);
        overlaydata_raw(:,1) = grid_lon_origin + 360.0*(overlaydata_raw(:,1) - 0.5)/jmax;
        overlaydata_raw(:,2) = 180.0*asin(2.0*(overlaydata_raw(:,2) - 0.5)/jmax - 1.0)/pi;
    end
    % remove data in land cells (or for k == 0 selection)
    for n = 1:nmax,
        if (isnan(overlaydata_zm(overlaydata_ijk(n,2),overlaydata_ijk(n,1))) || (kplot == 0))
            overlaydata_raw(n,4) = NaN;
            overlaydata_ijk(n,4) = NaN;
        end
    end
    % convert lat to sin(lat) for plotting
    overlaydata(:,:) = overlaydata_raw(:,:);
    overlaydata(:,2) = sin(pi*overlaydata_raw(:,2)/180.0);
    % grid (and average per cell) data if requested
    % NOTE: data vector length is re-calculated and the value of nmax reset
    if (data_ijk_mean == 'y')
        overlaydata_ijk_old(:,:) = overlaydata_ijk(:,:);
        overlaydata_ijk(:,:) = [];
        overlaydata(:,:)    = [];
        m=0;
        for i = 1:imax,
            for j = 1:jmax,
                if (~isnan(overlaydata_zm(j,i)))
                    samecell_locations = find((overlaydata_ijk_old(:,1)==i)&(overlaydata_ijk_old(:,2)==j)&(overlaydata_ijk_old(:,3)==kplot));
                    samecell_n = size(samecell_locations);
                    if (samecell_n(1) > 0)
                        m=m+1;
                        samecell_mean = mean(overlaydata_ijk_old(samecell_locations,4));
                        overlaydata_ijk(m,:) = [i j kplot samecell_mean];
                        overlaydata(m,1) = grid_lon_origin + 360.0*(overlaydata_ijk(m,1) - 0.5)/jmax;
                        overlaydata(m,2) = 2.0*(overlaydata_ijk(m,2) - 0.5)/jmax - 1.0;
                        overlaydata(m,3) = laym(j,i);
                        overlaydata(m,4) = samecell_mean;
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
%
% calculate stats needed for Taylor Diagram (and plot it!)
%
% *** 3D (GRIDDED) DATA ************************************************* %
%
if (~isempty(dataid_2))
    % transform data sets in vectors
    if kplot > 0
        data_vector_1 = reshape(data_1(kplot,:,:),imax*jmax,1);
        data_vector_2 = reshape(data_2(kplot,:,:),imax*jmax,1);
    else
        data_vector_1 = reshape(data_1(:,:,:),imax*jmax*kmax,1);
        data_vector_2 = reshape(data_2(:,:,:),imax*jmax*kmax,1);
    end
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
    data_vector_1 = overlaydata(:,4);
    % populate the gridded dataset vector with values corresponding to
    % the overlay data locations
    % NOTE: !!! data is (k,j,i) !!! (=> swap i and j)
    for n = 1:nmax,
        data_vector_2(n) = data(kplot,overlaydata_ijk(n,2),overlaydata_ijk(n,1));
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
    fid = fopen([filename '_STATS' '.dat'], 'wt');
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
        overlaydata(:,4) = data_vector_2(:) - data_vector_1(:);
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
z_u_ex = [z_u z_u z_u];
z_v_ex = [z_v z_v z_v];
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
z_u = z_u_ex(:,i_start:i_start+imax-1);
z_v = z_v_ex(:,i_start:i_start+imax-1);
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
fh(3) = axes('Position',[0.80 0.15 0.20 0.70],'Visible','off');
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
set(gca,'XLabel',text('String','Longitude','FontSize',15),'XTick',[lon_min:delta_lon:lon_max]);
set(gca,'YLabel',text('String','Latitude','FontSize',15),'YTick',[-1 -0.866 -0.5 0 0.5 0.866 1], 'YTickLabel',{'-90';'-60';'-30';'0';'30';'60';'90'});
if ~isempty(graph_title)
    title(graph_title,'FontSize',18);
else
    title(['Data: ',strrep(dataid_1,'_',' '),' / Level (k) = ', num2str(kplot)],'FontSize',12);
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
if (contour_plot == 'y') && (data_only == 'n'),
    v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k-');
    set(h,'LineWidth',0.25);
    v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
    [C,h] = contour(xm_ex,sin(pi*ym_ex/180.0),zm_ex,v,'k');
    set(h,'LineWidth',1.0);
    if plot_log10 == 'y'
        %%%%%%%%
    elseif contour_label == 'y'
        clabel(C,h);
    end
end
%
% *** OVERLAY VELOCITY FIELD ******************************************** %
%
% plot velocity field if requested
if (data_uv == 'y'),
    [h] = scatter(reshape(xm(:,:),jmax*imax,1),sin(pi*reshape(ym(:,:),jmax*imax,1)/180.0),2.5,'filled','k');
    [h] = quiver(xm,sin(pi*ym/180.0),z_u,sin(pi*z_v/180.0),data_uv_scale,'k','MaxHeadSize',0.0);
    set(h,'LineWidth',0.75);
end
%
% *** OVERLAY DATA ****************************************************** %
%
% plot overlay data
if ~isempty(overlaydataid)
    scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata(:,4)/data_scale,'o','Filled','Sizedata',data_size,'MarkerEdgeColor','k');
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
h = fill([0.1 0.2 0.3],[c-0.1 c-0.9 c-0.1],cmap(2*c-1,:));
if plot_log10 == 'y'
    str = ['below ',num2str(10^(con_min + (c-1)*(con_max-con_min)/con_n))];
else
    str = ['below ',num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
end
textsize = 2+round(80/con_n);
if textsize > 10,
    textsize = 10;
end
text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
set(h,'LineWidth',0.5);
set(h,'EdgeColor','k');
for c = 2:con_n+1,
    h = fill([0.1 0.1 0.3 0.3],[c-0.9 c-0.1 c-0.1 c-0.9],cmap(2*c-1,:));
    if plot_log10 == 'y'
        str = [num2str(10^(con_min + (c-2)*(con_max-con_min)/con_n)),' - ',num2str(10^(con_min + (c-1)*(con_max-con_min)/con_n))];
    else
        str = [num2str(con_min + (c-2)*(con_max-con_min)/con_n),' - ',num2str(con_min + (c-1)*(con_max-con_min)/con_n)];
    end
    textsize = 2+round(80/con_n);
    if textsize > 10,
        textsize = 10;
    end
    text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
    set(h,'LineWidth',0.5);
    set(h,'EdgeColor','k');
end
c = con_n+2;
h = fill([0.1 0.2 0.3],[c-0.9 c-0.1 c-0.9],cmap(2*c-1,:));
if plot_log10 == 'y'
    str = ['above ',num2str(10^(con_min + (c-2)*(con_max-con_min)/con_n))];
else
    str = ['above ',num2str(con_min + (c-2)*(con_max-con_min)/con_n)];
end
textsize = 2+round(80/con_n);
if textsize > 10,
    textsize = 10;
end
text(0.40,c-0.5,str,'FontName','Arial','FontSize',textsize);
set(h,'LineWidth',0.5);
set(h,'EdgeColor','k');
%
hold off;
%
% *** PRINT PLOT ******************************************************** %
%
set(gcf,'CurrentAxes',fh(1));
print('-dpsc2', [filename '.' str_date '.ps']);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% close netCDF files
netcdf.close(ncid_1);
if ~isempty(exp_2)
    netcdf.close(ncid_2);
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** OLD *************************************************************** %
% *********************************************************************** %

% % SAVE DIFFERENCE netCDF
%
% if (length(exp_2) > 0) | (timesliceid_2 > 0.0)
%     %
%     data = flipdim(data,1);
%     %
%     if length(exp_2) > 0
%         nc_name = [exp_1, '_MINUS_' exp_2, '.', num2str(timesliceid_1), '.nc'];
%         if timesliceid_2 > 0.0
%             nc_name = [exp_1, '_MINUS_' exp_2, '.', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.nc'];
%         end
%     elseif timesliceid_2 > 0.0
%         nc_name = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.nc'];
%     end
%     % Create NetCDF file
%     nc = netcdf(nc_name, 'clobber');
%     % Global attributes
%     nc.description = ['Difference data: '];
%     nc.author = 'Andy Rigdwell';
%     nc.date = datestr(now);
%     % Define dimensions
%     nc('lat') = jmax;
%     nc('lat_edges') = jmax + 1;
%     nc('lon') = imax;
%     nc('lon_edges') = imax + 1;
%     nc('zt') = kmax;
%     nc('zt_edges') = kmax + 1;
%     nc('time') = 1;
%     % Define variables - grid
%     nc{'lat'} = {'lat'};
%     nc{'lat_edges'} = {'lat_edges'};
%     nc{'lon'} = {'lon'};
%     nc{'lon_edges'} = {'lon_edges'};
%     nc{'zt'} = {'zt'};
%     nc{'zt_edges'} = {'zt_edges'};
%     nc{'time'} = {'time'};
%     % Define variables - data
%     nc{'grid_level'} = {'lat', 'lon'};
%     nc{'grid_topo'} = {'lat', 'lon'};
%     nc{dataid} = {'time' 'zt', 'lat', 'lon'};
%     % Attributes - grid
%     nc{'lat'}.units = 'degrees';
%     nc{'lat_edges'}.units = 'degrees';
%     nc{'lon'}.units = 'degrees';
%     nc{'lon_edges'}.units = 'degrees';
%     nc{'zt'}.units = 'm';
%     nc{'zt_edges'}.units = 'm';
%     nc{'time'}.units = 'yr';
%     % Attributes - data
%     nc{'grid_level'}.units = '1';
%     nc{'grid_topo'}.units = '1';
%     nc{dataid}.units = 'n/a';
%     % Put - grid
%     nc{'lat'}(:) = grid_lat;
%     nc{'lat_edges'}(:) = grid_lat_edges;
%     nc{'lon'}(:) = grid_lon;
%     nc{'lon_edges'}(:) = grid_lon_edges;
%     nc{'zt'}(:) = zt;
%     nc{'zt_edges'}(:) = zt_edges;
%     nc{'time'}(:) = 0.0;
%     % Put - data
%     nc{'grid_level'}(:,:) = grid_level;
%     nc{'grid_topo'}(:,:) = grid_topo;
%     nc{dataid}(1,:,:,:) = data(:,:,:);
%     % Close the file
%     nc = close(nc);
% end
