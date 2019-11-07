function [] = plot_extract_CCD(PATH,EXPID,MASK)
% plot_extract_CCD
%
%   ***********************************************************************
%   *** Extract CCD data **************************************************
%   ***********************************************************************
%
%
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   10/07/20: cleanup
%   10/08/27: ***
%   11/08/20: copied to new filename
%
%   ***********************************************************************

% \/\/\/ USER SETTINGS \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ %
ctrl_weight_CCD     = 1;     % weight the mean CCD position?
par_CCD_CaCO3       =  20.0; % (wt%)
par_CCD_CaCO3_min   =   1.0; % (wt%)
par_CCD_CaCO3_max   =  50.0; % (wt%)
par_CCD_depth_min  = 1000.0; % (m)
par_CCD_Ddepth_max = 2000.0; % (m)
lon_offset = -0;       % LONGITUDE OFFSET
lon_min = -180;        % STARTING LONGITUDE FOR X-AXIS
% /\/\/\ USER SETTINGS /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\ %

% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
%
% *********************************************************************** %
% 
% close open windows
% NOTE: don't clear variable space here ...
close all;
% set passed parameters
path = PATH;
expid = EXPID;
basinmaskid = MASK;

% open netCDF file
ncid=netcdf.open([path '/' expid '/sedgem/fields_sedgem_2d.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% SET UP GRID
% load grid dimension data
varid  = netcdf.inqVarID(ncid,'lat');
[dimname, dimlen] = netcdf.inqDim(ncid,varid);
jmax = dimlen;
grid_lat = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'lon');
[dimname, dimlen] = netcdf.inqDim(ncid,varid);
imax = dimlen;
grid_lon = netcdf.getVar(ncid,varid) + lon_offset;
[lonm latm] = meshgrid(grid_lon,grid_lat);
varid  = netcdf.inqVarID(ncid,'lat_edges');
grid_lat_edges = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'lon_edges');
grid_lon_edges = netcdf.getVar(ncid,varid) + lon_offset;
[lonw lats] = meshgrid(grid_lon_edges(1:imax),grid_lat_edges(1:jmax));
[lone latn] = meshgrid(grid_lon_edges(2:imax+1),grid_lat_edges(2:jmax+1));
% load grid data
varid  = netcdf.inqVarID(ncid,'grid_mask');
grid_mask = flipud(netcdf.getVar(ncid,varid));
grid_mask(grid_mask > 1.0) = 0.0;
varid  = netcdf.inqVarID(ncid,'grid_topo');
grid_depth = -flipud(netcdf.getVar(ncid,varid));
grid_depth(grid_depth < 0.0) = 0.0;
% load mask data
% NOTE: flip in j-direction and rotate to make consistent with netCDF grid
maskfile = [basinmaskid '.dat'];
if ~isempty(basinmaskid)
    grid_mask_basin = load(maskfile,'-ascii');
    grid_mask_basin = flipdim(grid_mask_basin,1);
    grid_mask_basin = rot90(grid_mask_basin(:,:));
    grid_mask = grid_mask_basin(:,:).*grid_mask(:,:);
end
if isempty(basinmaskid)
    basinmaskid = 'GLOBAL';
end

% LOAD DATA
% wt% CaCO3
varid = netcdf.inqVarID(ncid,'sed_CaCO3');
grid_CaCO3 = flipud(netcdf.getVar(ncid,varid));
grid_CaCO3(grid_CaCO3 > 100.0) = 0.0;
grid_CaCO3 = double(grid_CaCO3);
% various geochemsitry
varid = netcdf.inqVarID(ncid,'ocn_temp');
grid_geochem(:,:,1)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'ocn_sal');
grid_geochem(:,:,2)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'ocn_DIC');
grid_geochem(:,:,3)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'ocn_ALK');
grid_geochem(:,:,4)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'carb_conc_CO2');
grid_geochem(:,:,5)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'carb_conc_HCO3');
grid_geochem(:,:,6)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'carb_conc_CO3');
grid_geochem(:,:,7)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'carb_ohm_cal');
grid_geochem(:,:,8)  = flipud(netcdf.getVar(ncid,varid));
varid = netcdf.inqVarID(ncid,'carb_H');
grid_geochem(:,:,9)  = flipud(netcdf.getVar(ncid,varid));

% initialize search grid
grid_search = zeros(jmax,imax);
grid_CCD = zeros(jmax,imax);
grid_CCD_b = zeros(jmax,imax);
grid_CCD_t = zeros(jmax,imax);
grid_mask_slope = zeros(jmax,imax);

% find slope
for i = 1:imax,
    for j = 1:jmax,
        if grid_mask(j,i) ~= 1.0
            for di = -1:1,
                for dj = -1:1,
                    i_neigh = i + di;
                    j_neigh = j + dj;
                    if (i_neigh > 0) & (i_neigh <= imax) & (j_neigh > 0) & (j_neigh <= jmax) & (grid_mask(j_neigh,i_neigh) == 1.0)
                        grid_mask_slope(j_neigh,i_neigh) = 1.0;
                    end
                end
            end
        end
    end
end

%   *********************************************************

% find CCD
n = 0;
for i = 1:imax,
    for j = 1:jmax,
        if grid_mask(j,i) ~= 1.0
            grid_search(j,i) = 1.0;
            grid_depth(j,i) = -1.0;
            grid_CaCO3(j,i) = -100.0;
            grid_CCD(j,i)   = -1.0;
            grid_CCD_b(j,i) = -1.0;
            grid_CCD_t(j,i) = -1.0;
        else
            %disp(['> ' '(' num2str(i) ',' num2str(j) ')'])
            if ((grid_CaCO3(j,i) < par_CCD_CaCO3) & (grid_CaCO3(j,i) > par_CCD_CaCO3_min) & (grid_depth(j,i) > par_CCD_depth_min) & (grid_mask_slope(j,i) ~= 1.0))
%                 disp([' > (' num2str(i) ',' num2str(j) '); wt% CaCO3 = ' num2str(grid_CaCO3(j,i)) '; depth = ' num2str(grid_depth(j,i))])
                for di = -1:1,
                    for dj = -1:1,
                        i_neigh = i + di;
                        j_neigh = j + dj;
                        %disp(['   (' num2str(i_neigh) ',' num2str(j_neigh) ')'])
                        if ((i_neigh > 0) & (i_neigh <= imax) & (j_neigh > 0) & (j_neigh <= jmax) & (grid_mask(j_neigh,i_neigh) == 1.0) & (grid_mask_slope(j_neigh,i_neigh)) ~= 1.0)
                            if grid_depth(j,i) > (grid_depth(j_neigh,i_neigh))
                                if (grid_CaCO3(j_neigh,i_neigh) > par_CCD_CaCO3) & (grid_depth(j,i) < (grid_depth(j_neigh,i_neigh) + par_CCD_Ddepth_max))
                                    if (grid_CaCO3(j_neigh,i_neigh) < par_CCD_CaCO3_max)
                                        n = n + 1;
                                        grid_CCD(j,i)               = 1.0;
                                        grid_CCD(j_neigh,i_neigh)   = 2.0;
                                        grid_CCD_b(j,i)             = 1.0;
                                        grid_CCD_t(j_neigh,i_neigh) = 1.0;
                                        cell_CCD_b(n,1)    = i;
                                        cell_CCD_b(n,2)    = j;
                                        cell_CCD_b(n,3)    = grid_depth(j,i);
                                        cell_CCD_b(n,4)    = grid_CaCO3(j,i);
                                        cell_CCD_b(n,5:13) = grid_geochem(j,i,1:9);
                                        cell_CCD_t(n,1)    = i_neigh;
                                        cell_CCD_t(n,2)    = j_neigh;
                                        cell_CCD_t(n,3)    = grid_depth(j_neigh,i_neigh);
                                        cell_CCD_t(n,4)    = grid_CaCO3(j_neigh,i_neigh);
                                        cell_CCD_t(n,5:13) = grid_geochem(j_neigh,i_neigh,1:9);
%                                         disp(['   ' '(' num2str(i_neigh) ',' num2str(j_neigh) '); wt% CaCO3 = ' num2str(grid_CaCO3(j_neigh,i_neigh)) '; depth = ' num2str(grid_depth(j_neigh,i_neigh))])
                                    end
                                end
                            end
                        end
                    end
                end
                grid_search(j,i) = 1.0;
            end
        end
    end
end
%
nmax = n;
% disp(['>> n = ' num2str(n)])

%   *********************************************************

% set maskID
if ~isempty(basinmaskid)
    
end

% average data
cell_CCD_mean = zeros(1,13);
cell_CCD = zeros(nmax,13);
if (ctrl_weight_CCD)
% normalize average CCD depth (and properties) to par_CCD_CaCO3
    for n = 1:nmax,
        tmp_DCaCO3_b = par_CCD_CaCO3 - cell_CCD_b(n,4);
        tmp_DCaCO3_t = cell_CCD_t(n,4) - par_CCD_CaCO3;
        tmp_weight_b = tmp_DCaCO3_t/(cell_CCD_t(n,4) - cell_CCD_b(n,4));
        tmp_weight_t = tmp_DCaCO3_b/(cell_CCD_t(n,4) - cell_CCD_b(n,4));
        cell_CCD(n,:) = tmp_weight_b*cell_CCD_b(n,:) + tmp_weight_t*cell_CCD_t(n,:);
        cell_CCD_mean(1,:) = cell_CCD_mean(1,:) + cell_CCD(n,:);
    end
else
    % simple average of CCD depth (and properties)
    for n = 1:nmax,
        cell_CCD(n,:) = (cell_CCD_b(n,:) + cell_CCD_t(n,:))/2.0;
        cell_CCD_mean(1,:) = cell_CCD_mean(1,:) + cell_CCD(n,:);
    end
end
cell_CCD_mean(1,:) = cell_CCD_mean(1,:)/nmax;

% save global data
fid = fopen([expid '_data_summary_' basinmaskid '.txt'], 'wt');
fprintf(fid, 'Mean CCD properties \n');
fprintf(fid, '\n');
fprintf(fid, 'Number of CCD-spanning grid point pairs = %6.0f \n', nmax);
fprintf(fid, '\n');
fprintf(fid, 'Mean CCD depth      : %6.1f m \n', cell_CCD_mean(1,3));
fprintf(fid, 'Mean CaCO3 content  : %6.2f wt%% \n', cell_CCD_mean(1,4));
fprintf(fid, 'Mean temperature    : %6.2f K \n', cell_CCD_mean(1,5));
fprintf(fid, 'Mean salinity       : %6.2f o/oo \n', cell_CCD_mean(1,6));
fprintf(fid, 'Mean alkalinity     : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,7));
fprintf(fid, 'Mean DIC            : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,8));
fprintf(fid, 'Mean [CO2(aq)]      : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,9));
fprintf(fid, 'Mean [HCO3-]        : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,10));
fprintf(fid, 'Mean [CO3-]         : %8.3f umol kg-1 \n', 1e6*cell_CCD_mean(1,11));
fprintf(fid, 'Mean calcite ohmega : %6.3f \n', cell_CCD_mean(1,12));
fprintf(fid, 'Mean pH             : %6.3f \n', -log10(cell_CCD_mean(1,13)));
fclose(fid);

%   *********************************************************

% plot CaCO3-depth relationship
figure;
hold on;
scatter(reshape(-grid_depth(:,:),jmax*imax,1),reshape(grid_CaCO3(:,:),jmax*imax,1),30,'k');
scatter(-cell_CCD_b(:,3),cell_CCD_b(:,4),20,[1.0 0.5 0.0],'filled');
scatter(-cell_CCD_t(:,3),cell_CCD_t(:,4),20,[0.6 1.0 0.1],'filled');
axis([-6000.0 0.0 0 100.0]);
filename = [expid '_plot_scatter_CaCO3vsdepth_' basinmaskid '.ps'];
print('-dpsc2', filename);
% axis([-6000.0 -3000.0 0 25.0]);
% filename = [EXP '_plot_scatter2_CaCO3vsdepth.ps'];
% print('-dpsc2', filename);

% plot depth
figure;
grid_plot = [grid_depth, zeros(jmax,1)];
grid_plot = [zeros(1,imax+1); grid_plot];
pcolor(rot90(grid_plot(:,:),3));
% pcolor(grid_plot(:,:));
axis([1 imax+1 1 jmax+1]);
filename = [expid '_plot_grid_depth_' basinmaskid '.ps'];
print('-dpsc2', filename);

% plot mask
figure;
grid_plot = [grid_mask, zeros(jmax,1)];
grid_plot = [zeros(1,imax+1); grid_plot];
pcolor(rot90(grid_plot(:,:),3));
axis([1 imax+1 1 jmax+1]);
filename = [expid '_plot_grid_mask_' basinmaskid '.ps'];
print('-dpsc2', filename);

% plot slope
figure;
grid_plot = [grid_mask_slope, zeros(jmax,1)];
grid_plot = [zeros(1,imax+1); grid_plot];
pcolor(rot90(grid_plot(:,:),3));
axis([1 imax+1 1 jmax+1]);
filename = [expid '_plot_grid_slopeth_' basinmaskid '.ps'];
print('-dpsc2', filename);

% plot CaCO3
figure;
grid_plot = [grid_CaCO3, zeros(jmax,1)];
grid_plot = [zeros(1,imax+1); grid_plot];
pcolor(rot90(grid_plot(:,:),3));
axis([1 imax+1 1 jmax+1]);
filename = [expid '_plot_grid_CaCO3_' basinmaskid '.ps'];
print('-dpsc2', filename);

% plot CCD
figure;
grid_plot = [grid_CCD, zeros(jmax,1) + 3.0];
grid_plot = [zeros(1,imax+1) + 3; grid_plot];
pcolor(rot90(grid_plot(:,:),3));
axis([1 imax+1 1 jmax+1]);
filename = [expid '_plot_grid_CCD_' basinmaskid '.ps'];
print('-dpsc2', filename);

% % plot CCD - bottom
% figure;
% grid_plot = [grid_CCD_b, zeros(jmax,1)];
% grid_plot = [grid_plot; zeros(1,imax+1)];
% pcolor(grid_plot(:,:));
% axis([1 imax+1 1 jmax+1]);
% filename = [RID '_plot_grid_CCD_b.ps'];
% print('-dpsc2', filename);

% % plot CCD - top
% figure;
% grid_plot = [grid_CCD_t, zeros(jmax,1)];
% grid_plot = [grid_plot; zeros(1,imax+1)];
% pcolor(grid_plot(:,:));
% axis([1 imax+1 1 jmax+1]);
% filename = [RID '_plot_grid_CCD_t.ps'];
% print('-dpsc2', filename);

% save data
savedata = fliplr(rot90(grid_CCD(1:imax,1:jmax),1));
save([expid '_data_CCD_' basinmaskid '.dat'],'savedata','-ascii');
fprintf(' - CCD locations saved\n');
% savedata = flipdim(grid_CCD_b(1:imax,1:jmax),1);
% save([RID '_data_CCD_bottom.dat'],'savedata','-ascii');
% fprintf(' - CCD bottom locations saved\n');
% savedata = flipdim(grid_CCD_t(1:imax,1:jmax),1);
% save([RID '_data_CCD_top.dat'],'savedata','-ascii');
% fprintf(' - CCD top locations saved\n');
savedata = fliplr(rot90(grid_depth(1:imax,1:jmax),1));
save([expid '_data_depth_' basinmaskid '.dat'],'savedata','-ascii');
fprintf(' - depth grid saved\n');
savedata = fliplr(rot90(grid_mask(1:imax,1:jmax),1));
save([expid '_data_mask_' basinmaskid '.dat'],'savedata','-ascii');
fprintf(' - mask saved\n');
savedata = fliplr(rot90(grid_CaCO3(1:imax,1:jmax),1));
save([expid '_data_CaCO3_' basinmaskid '.dat'],'savedata','-ascii');
fprintf(' - wt%% CaCO3 saved\n');

% CLOSE
netcdf.close(ncid);

% END
disp(['END ...'])

%   *********************************************************
