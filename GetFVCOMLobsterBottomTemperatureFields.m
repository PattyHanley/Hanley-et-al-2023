%%%  This code was used to extract information from the FVCOM-GOM 30 year
%%%  hindcast, and to save these data to a netcdf file for future
%%%  processing



setup_nctoolbox

url='http://www.smast.umassd.edu:8080/thredds/dodsC/fvcom/hindcasts/30yr_gom3'; 
nc=ncgeodataset(url);

lat = nc.data('lat');
lon = nc.data('lon');

westLon = -69.5;
eastLon = -63.15;
northLat = 46.10;
southLat = 43.40;
lobsterbox = intersect( find(lat<=northLat & lat>=southLat), ...
    find(lon<=eastLon & lon>=westLon)  );
lat = double(lat(lobsterbox));
lon = double(lon(lobsterbox));

%%%%% BEGIN trim lobsterbox
% trim off narrow harbors on north shore
drop = find(lat > 45.25 & lon < -65.95);
lat(drop) = [];
lon(drop) = [];
lobsterbox(drop) = [];

drop = find(lat > 45.08 & lon < -67.1);
lat(drop) = [];
lon(drop) = [];
lobsterbox(drop) = [];

drop = find(lat > 44.5 & lon > -68.86 & lon < -68.72);
lat(drop) = [];
lon(drop) = [];
lobsterbox(drop) = [];

drop = find(lat > 44.37 & lon > -68.8 & lon < -68.67);
lat(drop) = [];
lon(drop) = [];
lobsterbox(drop) = [];

% trim off waters SE of Nova Scotia
drop = find(lat < 44.8 & lon > -65.75);
lat(drop) = [];
lon(drop) = [];
lobsterbox(drop) = [];

%%%%% END trim lobsterbox

 

% INTERPOLANT = scatteredInterpolant();
% INTERPOLANT.Points = double([lon lat]);
% INTERPOLANT.ExtrapolationMethod = 'none';

bathymetry = nc.data('h');
bathymetry = bathymetry(lobsterbox);
% INTERPOLANT.Values = double(bathymetry);
% gridBathymetry = INTERPOLANT(gridLon,gridLat);


basedate = datenum([1858 11 17 0 0 0]);
jd = datenum(datevec(double(nc.data('time') + basedate)));
tempDays=floor(jd);
UtempDays=unique(tempDays);

tempSize = nc.size('temp');

dimLat = length(lat);
dimLon = length(lon);

clear schema
schema.Name = '/';
schema.Format = 'netcdf4';
schema.Dimensions(1).Name = 'time';
schema.Dimensions(1).Length = Inf;
schema.Dimensions(2).Name = 'gridcellnumber';
schema.Dimensions(2).Length = dimLat;
schema.Dimensions(3).Name = 'date_str_length';
schema.Dimensions(3).Length = 11;

% if exist('FVCOMBayOfFundyBottomTemperature.nc') == 2
%     eval('!del FVCOMBayOfFundyBottomTemperature.nc')
% end

ncwriteschema('FVCOMBayOfFundyBottomTemperature.nc',schema);

nccreate('FVCOMBayOfFundyBottomTemperature.nc','temperature','Dimensions',...
    {'gridcellnumber',dimLat,'time',1});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','depth','Dimensions',...
    {'gridcellnumber',dimLat});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','envmntWindow','Dimensions',...
    {'gridcellnumber',dimLat,'time',1});
% nccreate('FVCOMBayOfFundyBottomTemperature.nc','mask','Dimensions',...
%     {'lat',dimLat,'lon',dimLon});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','lat','Dimensions',...
    {'gridcellnumber',dimLat});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','lon','Dimensions',...
    {'gridcellnumber',dimLon});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','date','Dimensions',...
    {'time',1});
nccreate('FVCOMBayOfFundyBottomTemperature.nc','date_str','Dimensions',...
    {'date_str_length',11,'time',1},'datatype','char');


% ncwrite('FVCOMBayOfFundyBottomTemperature.nc','mask',mask)
ncwrite('FVCOMBayOfFundyBottomTemperature.nc','lat',lat)
ncwrite('FVCOMBayOfFundyBottomTemperature.nc','lon',lon)

ncwrite('FVCOMBayOfFundyBottomTemperature.nc','depth',...
    bathymetry)

% temperatureFields = NaN*zeros(length(vecLat),length(vecLon));%,length(UtempDays));
% temperatureFields = NaN*zeros(length(lobsterbox),length(UtempDays));

%%% Dates for Brent W.
% UtempDays = UtempDays(UtempDays > datenum(2009,10,1));
%%%

nDays = length(UtempDays);

start = 1;
for d = 1:(nDays - 1)
    tic
    nc_idx = d - start + 1;
    
    % BEGIN calc daily mean; read in the hourly data; compress to a daily mean
%     itime = date_index(jd,datevec(UtempDays(d)),datevec(UtempDays(d+1)-0.01));
    itime = datefind(jd,datevec(UtempDays(d)),datevec(UtempDays(d+1)-0.01));
    itime = [min(itime) max(itime)];
    data = mean(squeeze(nc.data('temp',[itime(1) 45 1],[itime(2) 45 tempSize(3)])));
    % END calc daily mean
    
    data = data(lobsterbox)';
    ncwrite('FVCOMBayOfFundyBottomTemperature.nc','temperature',...
        data,[1 nc_idx]);
    ncwrite('FVCOMBayOfFundyBottomTemperature.nc','date',...
        UtempDays(d),nc_idx);
    ncwrite('FVCOMBayOfFundyBottomTemperature.nc','date_str',...
        datestr(UtempDays(d))',[1 nc_idx],[1 1]);    clc
    disp(sprintf('chill bro, i am %f percent done with region 1',d/length(UtempDays)));
    toc
    
    %     INTERPOLANT.Values = double(nanmean(data,1))';
    %     temperatureFields(:,:,d) = INTERPOLANT(gridLon,gridLat);
end
