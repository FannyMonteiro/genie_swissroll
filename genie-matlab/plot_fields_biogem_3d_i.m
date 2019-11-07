function [] = plot_fields_biogem_3d_i(PATH,EXP1,EXP2,PT1,PT2,PVAR,MASK,PI,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
% plot_fields_biogem_3d_i
%
%   *********************************************************
%   *** biogem i-SECTION (LAT-LAY) + INTEGRATED PLOTTING  ***
%   *********************************************************
%
%   plot_fields_biogem_3d_i(PATH,EXP1,EXP2,PT1,PT2,PVAR,MASK,PI,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   plots slices and zonally averaged vertical sections from the BIOGEM 3-D netCDF data file:
%   'fields_biogem_3d.nc'
%
%   plot_fields_biogem_3d_i(PATH,EXP1,EXP2,PT1,PT2,PVAR,MASK,PI,CSCALE,CMIN,CMAX,CN,PTIT,PDATA)
%   takes 14 arguments;
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
%   PVAR .... id the name of the variable to be plotted
%             -> all valid valiable names will be listed if an invalid name
%                is entered
%             THIS PARAMETER MUST BE PASSED AS A STRING;
%             -> e.g., 'ocn_PO4'
%   MASK .... is the filename containing the meridional mask
%             to be used to construct the zonal average
%             -> the full filename must eb give, including any extensions
%             THIS PARAMETER MUST BE PASSED AS A STRING
%             -> leave MASK blank, i.e., '', for no mask
%   PI ...... is the meridional section to be plotted (the 'i' slice)
%             (in the absence of a mask being specified)
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
%           plot_fields_biogem_3d_i('genie_output','experiment_1','',1994.5,-1,'ocn_PO4','mask_Pacific.dat',1,1e-6,0,2,20,'PO_4','')
%           will plot the time-slice cenetered on a time of 1994.5,
%           of PO4 concentrations zonally averaged according to
%           the mask file mask_Pacific.dat,
%           between 0 and 2 umol kg-1 in 20 contour intervals
%           experiment is called 'experiment_1'
%           and lives in the 'genie_output' subdirectory
%
%   Edit the 'm' file to change other user settings;
%           contour_plot = 'n';    % OVERLAY CONTOL PLOT?
%           contour_mod = 1;       % NUMBER OF COLOR INTERVALS PER CONTOR
%           contour_mod_label = 5; % NUMBER OF LABELED CONTOURS PER CONTOUR
%           contour_label = 'n';   % LABEL CONTOURS?
%           contour_noneg = 'y';   % RESTRICT DATA PLOTTED TO > 0.0?
%           plot_log10 = 'n';      % PLOT LOG10 OF THE DATA
%           overlay_size = 50.0;   % SIZE OF OVERLAY DATA POINTS 
%           data_anomoly = 'n';    % PLOT AS MODEL-DATA ANOMOLY ONLY?
%           data_only = 'n';       % PLOT ONLY DATA (no model values)?
%           plot_opsi = 'g';       % PLOT OVERTURNING STREAMFUNCTION (basin)?
%           dscrsz = 0.60;         % FRACTIONAL FIGURE WINDOW SIZE
%
%   *******************************************************************   %
%   *** HISTORY *******************************************************   %
%   *******************************************************************   %
%
%   11/05/30: Added time-stamping
%   12/02/10: added in options for: anomoly plotting; data-only plotting
%             some code reorganisation / rationalization
%             added overturning streamfubnction contour plotting
%
%   *******************************************************************   %

% \/\/\/ USER SETTINGS \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ %
contour_plot = 'y';     % [ 'y']  OVERLAY CONTOUR PLOT?
contour_mod = 1;        % [   1]  NUMBER OF COLOR INTERVALS PER CONTOR
contour_mod_label = 5;  % [   5]  NUMBER OF LABELED CONTOURS PER CONTOUR
contour_label = 'y';    % [ 'y']  LABEL CONTOURS?
contour_noneg = 'n';    % [ 'n']  RESTRICT DATA PLOTTED TO > 0.0?
plot_log10 = 'n';       % [ 'n']  PLOT LOG10 OF THE DATA
data_offset = 0.0;      % [ 0.0]  data offset (273.15 for K -> C)
data_size = 50.0;       % [50.0]  SIZE OF OVERLAY DATA POINTS
data_anomoly = 'n';     % [ 'n']  PLOT AS MODEL-DATA ANOMOLY ONLY?
data_only = 'n';        % [ 'n']  PLOT ONLY DATA (no model values)?
plot_opsi = '';         % [  '']  PLOT OVERTURNING STREAMFUNCTION (basin)?
plot_opsi_min = -15;    % [ -15] 
plot_opsi_max = +15;    % [ +15] 
plot_opsi_dminor = 1;   % [   1] 
plot_opsi_dmajor = 5;   % [   5] 
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
D_min   = 0000;
D_max   = 5000;
zt_min = 0;
zt_max = 5000;
% set passed parameters
path = PATH;
exp_1 = EXP1;
exp_2 = EXP2;
timesliceid_1 = PT1;
timesliceid_2 = PT2;
dataid = PVAR;
iplot = PI;
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
str_function = 'plot-fields-biogem-3d(i)';
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
data_0=[];
data_1=[];
data_2=[];
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
% open streamfunction data (if selected)
if ~isempty(plot_opsi)
    % open netCDF file
    ncid_0=netcdf.open([path '/' exp_1 '/biogem/fields_biogem_2d.nc'],'nowrite');
    % read netCDf information
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid_0);
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
    end
end
if ~isempty(maskid)
    topo = mask.*topo;
else
    mask = zeros(jmax,imax);
    mask(:,iplot) = 1.0;
    topo = mask.*topo;
end
% load and calculate remaining grid information
varid  = netcdf.inqVarID(ncid_1,'lat');
grid_lat = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon');
grid_lon = netcdf.getVar(ncid_1,varid);
[latm laym] = meshgrid(grid_lat,-grid_zt);
varid  = netcdf.inqVarID(ncid_1,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid_1,varid);
varid  = netcdf.inqVarID(ncid_1,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid_1,varid);
[lats layb] = meshgrid(grid_lat_edges(1:jmax),-grid_zt_edges(1:kmax));
[latn layt] = meshgrid(grid_lat_edges(2:jmax+1),-grid_zt_edges(2:kmax+1));
% calculate cell volumes
% NOTE: assume equal area grid, normaalized area
data_V = zeros(kmax,jmax,imax);
for k = 1:kmax,
    data_V(k,:,:) = 1.0*(grid_zt_edges(k) - grid_zt_edges(k+1));
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID: STREAMFUNCTION *************************************** %
% *********************************************************************** %
%
if ~isempty(plot_opsi)
    % load grid information
    varid  = netcdf.inqVarID(ncid_0,'zt');
    [dimname, dimlen] = netcdf.inqDim(ncid_0,varid);
    zmax = dimlen;
    opsigrid_zt = netcdf.getVar(ncid_0,varid);
    varid  = netcdf.inqVarID(ncid_0,'lat');
    opsigrid_lat = netcdf.getVar(ncid_0,varid);
    [latm ztm] = meshgrid(grid_lat,-grid_zt);
    varid  = netcdf.inqVarID(ncid_0,'zt_edges');
    opsigrid_zt_edges = netcdf.getVar(ncid_0,varid);
    varid  = netcdf.inqVarID(ncid_0,'lat_edges');
    opsigrid_lat_edges = netcdf.getVar(ncid_0,varid);
    [opsilats zts] = meshgrid(opsigrid_lat_edges(1:jmax),opsigrid_zt_edges(1:zmax));
    [opsilatn ztn] = meshgrid(opsigrid_lat_edges(2:jmax+1),opsigrid_zt_edges(2:zmax+1));
%     %do something for zt
%     dz = zeros(size(ztn));
%     zts = dz;
%     for i=1:zmax
%         dz(i,:) = opsigrid_zt_edges(i+1)-opsigrid_zt_edges(i);
%     end
%     ztn = ztn./dz;
%     ztn = ztn.*zt_max./ztn(zmax,1);
%     zts(2:zmax,:) = ztn(1:zmax-1,:);
end
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
        if strcmp(varname,dataid)
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
    % load data
    varid = netcdf.inqVarID(ncid_2,dataid);
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_2,varid);
    rawdata = netcdf.getVar(ncid_2,varid);
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
% *** SET DATA ********************************************************** %
%
if (~isempty(exp_2)) || (timesliceid_2 > 0.0)
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
% *** SET STREAMFUNCTION DATASET **************************************** %
% *********************************************************************** %
%
if ~isempty(plot_opsi)
    % set variable name
    switch plot_opsi
        case {'g'}
            varid  = netcdf.inqVarID(ncid_0,'phys_opsi');
        case 'a'
            varid  = netcdf.inqVarID(ncid_0,'phys_opsia');
        case 'p'
            varid  = netcdf.inqVarID(ncid_0,'phys_opsip');
        otherwise
            disp('Unknown opsi definition.')
    end
    % open data
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid_0,varid);
    rawdata = netcdf.getVar(ncid_0,varid);
    if length(dimids) == 3
        rawdata = rawdata(1:jmax,1:zmax,tid);
        data_0 = rawdata(1:jmax,1:zmax)';
    elseif length(dimids) == 2
        rawdata = rawdata(1:jmax,1:zmax);
        data_0 = rawdata(1:jmax,1:zmax)';
    else
        data_0 = NaN*data_0;
    end
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
        filename = [exp_1, '_MINUS_' exp_2, 'y', '.', num2str(timesliceid_1), '.', dataid, '.', maskid];
        if timesliceid_2 > 0.0
            filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid, '.', maskid];
        end
    elseif timesliceid_2 > 0.0
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid, '.', maskid];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid, '.', maskid];
    end
else
    if ~isempty(exp_2)
        filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '.', dataid, '.', 'i', num2str(iplot)];
        if timesliceid_2 > 0.0
            filename = [exp_1, '_MINUS_' exp_2, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid, '.', 'i', num2str(iplot)];
        end
    elseif timesliceid_2 > 0.0
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '_MINUS_', 'y', num2str(timesliceid_2), '.', dataid, '.', 'i', num2str(iplot)];
    else
        filename = [exp_1, '.', 'y', num2str(timesliceid_1), '.', dataid, '.', 'i', num2str(iplot)];
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** FILTER & PROCESS RAW DATA ***************************************** %
% *********************************************************************** %
%
% *** INITALIZE ********************************************************* %
%
xm = latm;
ym = laym;
data = data_1 - data_2;
data = data - data_offset;
if ~isempty(plot_opsi)
    opsiym = ztm;
    opsidata = data_0;
end
% define initial array sizes
zm = zeros(kmax,jmax);
zm_count = zm;
zl = zeros(kmax,1);
zl_V = zl;
zz = zeros(kmax,jmax);
zz_V = zz;
%
% *** PROCESS MAIN DATASET ********************************************** %
%
for k = 1:kmax,
    for j = 1:jmax,
        for i = 1:imax,
            if (topo(j,i) < -grid_zt(k))
                if ((data(k,j,i) > -1.0E6) && (data(k,j,i) < 0.9E36))
                    zm(k,j) = zm(k,j) + data(k,j,i);
                    zm_count(k,j) = zm_count(k,j) + 1.0;
                    zl(k) = zl(k) + data_V(k,j,i)*data(k,j,i);
                    zl_V(k) = zl_V(k) + data_V(k,j,i);
                    zz(k,j) = zz(k,j) + data_V(k,j,i)*data(k,j,i);
                    zz_V(k,j) = zz_V(k,j) + data_V(k,j,i);
                end
            end
        end
        if (zm_count(k,j) > 0.0)
            if plot_log10 == 'y'
                if (zm(k,j) > 0.0)
                    zm(k,j) = log10(zm(k,j)/data_scale/zm_count(k,j));
                else
                    zm(k,j) = NaN;
                end
            else
                zm(k,j) = zm(k,j)/data_scale/zm_count(k,j);
                if contour_noneg == 'y'
                    if (zm(k,j) < 0.0)
                        zm(k,j) = 0.0;
                    end
                end
            end
        else
            zm(k,j) = NaN;
        end
    if (zz_V(k,j) > 0.0)
        zz(k,j) = zz(k,j)/data_scale/zz_V(k,j);
    else
        zz(k,j) = NaN;
    end
    end
    if (zl_V(k) > 0.0)
        zl(k) = zl(k)/data_scale/zl_V(k);
    else
        zl(k) = NaN;
    end
end
% set topography uniform in i-direction and equal to deepest point found
for j = 1:jmax,
    for i = 2:imax,
        if (topo(j,i) < topo(j,i-1))
            topo(j,1:i-1) = topo(j,i);
        else
            topo(j,i) = topo(j,i-1);
        end
    end
end
%
% *** PROCESS OVERTURNING STREAMFUNCTION (IF SELECTED) ****************** %
%
if ~isempty(plot_opsi)
    opsizm = zeros(zmax,jmax);
    for j = 1:jmax,
        for z = 1:zmax,
            if (opsidata(z,j) < -1.0E6) || (opsidata(z,j) > 1.0E30)
                opsizm(z,j) = NaN;
            else
                opsizm(z,j) = opsidata(z,j);
            end
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD (OPTIONAL) OVERLAY DATA ************************************** %
% *********************************************************************** %
%
if ~isempty(overlaydataid)
    % load overlay data
    overlaydatafile = [overlaydataid];
    overlaydata = load(overlaydatafile,'-ascii');
    overlaydata(:,2) = sin(pi*overlaydata(:,2)/180.0);
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
%         overlaydata(:,4) = data_vector_2(:) - data_vector_1(:);
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
axis([lat_min lat_max -D_max -D_min]);
set(gca,'TickDir','out');
set(gca,'XLabel',text('String','Latitude','FontSize',15),'XTick',[-90 -60 -30 0 30 60 90]);
set(gca,'YLabel',text('String','Depth (km)','FontSize',15),'YTick',[-D_max:1000:-D_min],'YTickLabel',{'5';'4';'3';'2';'1';'0'});
if ~isempty(graph_title)
    title(graph_title,'FontSize',18);
else
    if ~isempty(maskid)
        title(['Data ID: ',strrep(dataid,'_','-'),' / mask = ', strrep(maskid,'_','-')],'FontSize',12);
    else
        title(['Data ID: ',strrep(dataid,'_','-'),' / i = ', num2str(iplot)],'FontSize',12);
    end
end
% draw filled rectangles
for j = jmax:-1:1,
    for k = 1:kmax,
        if topo(j,iplot) > ym(k,j)
            h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],color_g);
            set(h,'EdgeColor',color_g);
        else
            if (isnan(zm(k,j)))
                h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],[1 1 1]);
                set(h,'EdgeColor',[1 1 1]);
            else
                col = round((1/2)+con_n*(zm(k,j)-con_min)/(con_max-con_min));
                if col < 1
                    col = 1;
                elseif col > con_n
                    col = 2*(con_n+1)+1;
                else
                    col = 2*col+1;
                end
                h = patch([lats(k,j) lats(k,j) latn(k,j) latn(k,j)],[layb(k,j) layt(k,j) layt(k,j) layb(k,j)],cmap(col,:));
                set(h,'EdgeColor',cmap(col,:));
            end
        end
    end
end
%
% *** PLOT CONTINENTAL OUTLINE ****************************************** %
%
for k = 1:kmax,
    for j = 1:jmax-1,
        if topo(j,iplot) > ym(k,j)
            if topo(j+1,iplot) <= ym(k,j+1)
                h = plot([latn(k,j) latn(k,j)],[layb(k,j) layt(k,j)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
    for j = 2:jmax,
        if topo(j,iplot) > ym(k,j)
            if topo(j-1,iplot) <= ym(k,j-1)
                h = plot([lats(k,j) lats(k,j)],[layb(k,j) layt(k,j)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
for j = 1:jmax,
    for k = 2: kmax,
        if topo(j,iplot) < ym(k,j)
            if topo(j,iplot) > ym(k-1,j)
                h = plot([lats(k,j) latn(k,j)],[layb(k,j) layb(k,j)],'k-');
                set(h,'LineWidth',1.0);
            end
        end
    end
end
%
% *** OVERLAY CONTOURS ************************************************** %
%
if (contour_plot == 'y') && (isempty(plot_opsi))
    v = [con_min:(con_max-con_min)/(con_n/contour_mod):con_max];
    [C,h] = contour(xm,ym,zm,v,'k-');
    set(h,'LineWidth',0.25);
    v = [con_min:(con_max-con_min)/(con_n/contour_mod_label):con_max];
    [C,h] = contour(xm,ym,zm,v,'k');
    set(h,'LineWidth',1.0);
    if plot_log10 == 'y'
        %%%%%%%%
    elseif contour_label == 'y'
        clabel(C,h);
    end
end
%
% *** OVERLAY CONTOURS -- OVERTURNING STREAMFUNCTION (IF SELECTED) ****** %
%
if ~isempty(plot_opsi)
    v = [0.0:plot_opsi_dminor:plot_opsi_max];
    [C,h] = contour(xm,opsiym,opsizm,v,'k-');
    set(h,'LineWidth',0.25);
    v = [plot_opsi_min:plot_opsi_dminor:0.0];
    [C,h] = contour(xm,opsiym,opsizm,v,'k:');
    set(h,'LineWidth',0.25);
    v = [0.0:plot_opsi_dmajor:plot_opsi_max];
    [C,h] = contour(xm,opsiym,opsizm,v,'k');
    set(h,'LineWidth',1.0);
    v = [plot_opsi_min:plot_opsi_dmajor:0.0];
    [C,h] = contour(xm,opsiym,opsizm,v,'k:');
    set(h,'LineWidth',1.0);
    if contour_label == 'y'
        clabel(C,h);
    end
end
%
% *** OVERLAY DATA ****************************************************** %
%
if ~isempty(overlaydataid)
    scatter(overlaydata(:,1),overlaydata(:,2),4,overlaydata(:,3)/data_scale,'o','Filled','Sizedata',overlay_size,'MarkerEdgeColor','k');
end
%
% *** PLOT BORDER ******************************************************* %
%
h = plot([lat_min lat_max],[-D_max -D_max],'k-');
set(h,'LineWidth',1.0);
h = plot([lat_min lat_max],[-D_min -D_min],'k-');
set(h,'LineWidth',1.0);
h = plot([lat_min lat_min],[-D_max -D_min],'k-');
set(h,'LineWidth',1.0);
h = plot([lat_max lat_max],[-D_max -D_min],'k-');
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
print('-dpsc2', [filename '.ps']);
%
% *** PLOT FIGURE (profile) ********************************************* %
%
figure
plot(zl(:),-grid_zt(:));
hold on;
scatter(zl(:),-grid_zt(:),25,'r');
axis([con_min con_max -grid_zt_edges(1) -grid_zt_edges(kmax+1)]);
xlabel(strrep(dataid,'_','-'));
ylabel('Elevation (m)');
if ~isempty(graph_title)
    title(graph_title,'FontSize',18);
else
    if ~isempty(maskid)
        title(['Data ID: ',strrep(dataid,'_','-'),' / i = ', strrep(maskid,'_','-')],'FontSize',12);
    else
        title(['Data ID: ',strrep(dataid,'_','-'),' / i = ', num2str(iplot)],'FontSize',12);
    end
end
print('-dpsc2', [filename '.PROFILE', str_date, '.ps']);
%
% *** PRINT PLOT (profile) ********************************************** %
%
fprint_1D2_d([flipud(grid_zt(:)) flipud(zl(:))],[filename '.PROFILE', str_date, '.res'])
%
% *** PLOT FIGURE (surface zonal mean) ********************************** %
%
figure
plot(grid_lat,zz(kmax,:));
hold on;
scatter(grid_lat,zz(kmax,:),25,'r');
axis([-90.0 90.0 con_min con_max ]);
xlabel('Latitude');
ylabel(strrep(dataid,'_','-'));
print('-dpsc2', [filename '.ZONAL.', str_date, '.ps']);
%
% *** SAVE DATA (surface zonal mean) ************************************ %
%
fprint_1Dn_d([flipud(grid_lat) rot90(zz(kmax,:),1)],[filename '.ZONAL.', str_date, '.res'])
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
if ~isempty(plot_opsi)
    netcdf.close(ncid_0);
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
