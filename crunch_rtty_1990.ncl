
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/glade/u/home/capehart/heat_stress.ncl"



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("avgrunlen")
function avgrunlen( runs )
; find max run length
; Example:
; runs = dim_numrun_n(p.gt.5, 1, 0) ; where p(ntim) or p(ntim,jlat,ilon)
;                                   ; runs(ntim) or runs(ntim,jlat,ilon)
; maxrun = maxrunlen( runs )
;
; http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_numrun_n.shtml
;

local dimr, rankr, avgrun, i, j,length_of_run,runs_per_length,total_runs_per_length
begin

  dimr   = dimsizes(runs)         
  rankr  = dimsizes(dimr)

  if (.not.(rankr.eq.3)) then
      print("maxrunlen: rankr="+rankr+" not yet supported")
      exit
  end if


  if (rankr.eq.3) then
      avgrun = new ( dimr(1:rankr-1), "float")
      length_of_run    = tofloat(ispan(1,dimr(0),1)) ; time
      do j=0,dimr(1)-1
        do i=0,dimr(2)-1
           runs_per_length  = tofloat((/ runs(:,j,i) /))
           total_runs_per_length = sum(runs_per_length)

           if ( total_runs_per_length .gt. 0) then
               avgrun(j,i) =  (/ sum(length_of_run*runs_per_length) / total_runs_per_length /)
           else
               avgrun(j,i) = 0
           end if
        end do
      end do
      copy_VarCoords(runs(0,:,:), avgrun)
  end if

  return(avgrun)
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




undef("maxrunlen")
function maxrunlen( runs )                     
; find max run length
; Example:
; runs = dim_numrun_n(p.gt.5, 0, 0) ; where p(ntim) or p(ntim,jlat,ilon)
;                                   ; runs(ntim) or runs(ntim,jlat,ilon)
; maxrun = maxrunlen( runs )
;
; http://www.ncl.ucar.edu/Document/Functions/Built-in/dim_numrun_n.shtml
;

local dimr, rankr, maxrun, i, j, indrun
begin

  dimr   = dimsizes(runs)         
  rankr  = dimsizes(dimr)

  if (.not.(rankr.eq.1 .or. rankr.eq.3)) then
      print("maxrunlen: rankr="+rankr+" not yet supported")
      exit
  end if

  if (rankr.eq.1) then
      maxrun = new( 1, "integer", -999)
      indrun = ind( runs.ne.0 )   
      if (.not.ismissing(indrun(0))) then
          maxrun = indrun(dimsizes(indrun)-1) + 1
      else
          maxrun = 0
      end if
  end if

  if (rankr.eq.3) then
      maxrun = new ( dimr(1:rankr-1), "integer", -999)
      do j=0,dimr(1)-1
        do i=0,dimr(2)-1
           indrun = ind( runs(:,j,i).ne.0 )   
           if (.not.ismissing(indrun(0))) then
               maxrun(j,i) = indrun(dimsizes(indrun)-1) + 1
           end if
           delete(indrun)       ; may change size next iteration
        end do
      end do
      copy_VarCoords(runs(0,:,:), maxrun)
  end if

  return(maxrun)
end





begin

 ;   setfileoption("nc", "FileStructure", "Advanced")
 ;   setfileoption("nc", "Format",  "NetCDF4")


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  User Modification Area
   ;
   ;   User selects Repositories, Years and File Locations
   ;

   ;
   ; Repository, and Date Runs
   ;


   scenario = (/  "rtty_1990", \
                  "rtty_2020", \
                  "rtty_2030", \
                  "rtty_2050", \
                  "rtty_2080" /)

   ;
   ; Physics Options
   ;

   radiation      = (/ "RRTM", \
                       "RRTM", \
                       "RRTM", \
                       "RRTM", \
                       "RRTM" /)

   convection     = (/ "Tiedtke", \
                       "Tiedtke", \
                       "Tiedtke", \
                       "Tiedtke", \
                       "Tiedtke" /)

   microphysics   = (/ "Thompson", \
                       "Thompson", \
                       "Thompson", \
                       "Thompson", \
                       "Thompson" /)

   boundary_layer = (/ "YSU", \
                       "YSU", \
                       "YSU", \
                       "YSU", \
                       "YSU" /)
   ;
   ; Start and End Years
   ;

   year_starts = (/ 1990, 2020, 2030, 2050, 2080 /)
   year_ends   = (/ 2000, 2030, 2040, 2060, 2090 /)

   nyr = 11   ; # of years to run... currently each has 11 yr

   n_scenarios = dimsizes(scenario)

   ;
   ; Clipbox Latitude and Longitudes
   ;


   far_south =   23.00 ; degrees_north
   far_north =   55.00 ; degrees_north

   far_west  = -129.00 ; degrees_east
   far_east  =  -52.00 ; degrees_east


   ;
   ; Location of similations on Glade
   ;

   root_dir = "/glade/p/p66770001/DNV_Ensembles/ERAI/"



   ;
   ; Mass and Energy Buckets
   ;

   bucket_mm = 100.0 ; mm
   bucket_J  = 1.0e9 ; J

   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  Time Control Units and Calendar Declaration
   ;
   ;   These runs only use a 365 day calendar (no leap years)
   ;
   ;   We only want ONE reference time for our monthly time stamping
   ;


   time_units =  "days since 1990-01-01 00:00:00"
   year_units = "years since 0000-01-01 00:00:00"

   calendar = 0
      calendar@calendar = "365_day"


   hh_aux_plevl = (/ (/  0,  1,  2,  3,  4,  5 /),   \
                     (/  6,  7,  8,  9, 10, 11 /),   \
                     (/ 12, 13, 14, 15, 16, 17 /),   \
                     (/ 18, 19, 20, 21, 22, 23 /)  /)

   hh_aux_plevl@description = "hourly indicies for wrf_plev->aux files"
   hh_aux_plevl!0           = "6_hr_WRF_Plev_Time_Steps"
   hh_aux_plevl!1           = "1_hr_AUX_Time_Steps"


   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  Load Static Geospatial Fields and creating requested subset
   ;     Lat and Lon
   ;     Terrain Elevation and Orography Detail
   ;     Land Cover, Soil Type, Land Mask
   ;

   wrf_model_file = root_dir         + \
                    "rnty_1990-2000" + \
                    "/plevs/wrfout_d01_1998-08-01_00:00:00_PLEV.nc"

   wf = addfile(wrf_model_file, "r")

   lon0 = (/ wf->XLONG(0,0,:) /)
      lon0!0           = "lon0"
      lon0&lon0         =  lon0
      lon0@description = "longitude"


   lat0 = (/  wf->XLAT(0,:,0) /)
      lat0!0           = "lat0"
      lat0&lat0         =  lat0
      lat0@description = "latitude"

   lon = lon0({far_west:far_east})
      lon!0             = "lon"
      lon@description   = "longitude"
      lon@long_name     = "longitude"
      lon@standard_name = "longitude"
      lon@units         = "degrees_east"

   lat = lat0({far_south:far_north})
      lat!0             = "lat"
      lat@description   = "latitude"
      lat@long_name      = "latitude"
      lat@standard_name = "latitude"
      lat@units         = "degrees_north"



   soil_depth = (/ wf->ZS(0,:) /)
      soil_depth!0           = "soil_depth"
      soil_depth@description = "Center of Soil Layer"
      soil_depth@long_name   = soil_depth@description
      soil_depth@units       = "m"

   n_soil = dimsizes(soil_depth)

   soil_thickness = (/ wf->DZS(0,:) /)
      soil_thickness!0           = "soil_depth"
      soil_thickness@description = "Center of Soil Layer"
      soil_thickness@long_name = "Center of Soil Layer"
      soil_thickness@units       = "m"

   lat_s = ind( lat0 .eq. min(lat) )
   lat_n = ind( lat0 .eq. max(lat) )

   lon_w = ind( lon0 .eq. min(lon) )
   lon_e = ind( lon0 .eq. max(lon) )

   ny = dimsizes(lat)
   nx = dimsizes(lon)


   terrain_height = (/  wf->HGT(0, lat_s:lat_n, lon_w:lon_e)  /)
      terrain_height!0             = "lat"
      terrain_height!1             = "lon"
      terrain_height@description   = "Terrain Elevation"
      terrain_height@long_name   = "Terrain Elevation"
      terrain_height@standard_name = "surface_altitude"
      terrain_height@units         = "m"

   terrain_height_variance = (/ wf->VAR(0, lat_s:lat_n, lon_w:lon_e) /)
      terrain_height_variance!0             = "lat"
      terrain_height_variance!1             = "lon"
      terrain_height_variance@description   = "Orography Variance"
      terrain_height_variance@long_name   = "Orography Variance"
      terrain_height_variance@standard_name = "surface_altitude standard_error"
      terrain_height_variance@units         = "m2"

   soil_texture = (/  wf->ISLTYP(0, lat_s:lat_n, lon_w:lon_e)  /)
      soil_texture!0             = "lat"
      soil_texture!1             = "lon"
      soil_texture@description   = "Soil Texture"
      soil_texture@long_name   = "Soil Texture"
      soil_texture@standard_name = "soil_type"
      soil_texture@lookup_table_variable = "labels_for_variable_soil_texture"


   land_classification = toint( (/  wf->IVGTYP(0, lat_s:lat_n, lon_w:lon_e)  /) )
      land_classification!0             = "lat"
      land_classification!1             = "lon"
      land_classification@description   = "Land Use Category"
      land_classification@long_name   = "Land Use Category"
      land_classification@lookup_table_variable = "labels_for_variable_land_classification"

   land_mask = toint(  (/  wf->LANDMASK(0, lat_s:lat_n, lon_w:lon_e)  /) )
      land_mask!0             = "lat"
      land_mask!1             = "lon"
      land_mask@description   = "Land Mask"
      land_mask@long_name     = "Land Mask"
      land_mask@standard_name = "land_binary_mask"
      land_mask@comment1      = "0 = water;  1=land"



    temp_mask = land_mask
    temp_mask = where(land_classification.eq.16,0,land_mask)

   delete(wf)

   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  Initialize Output Arrays
   ;
   ;
   ;

   nm  = 12 * nyr ; # of total months in each ensemble battery

   ;
   ;  Monthly Floating Point Values
   ;

   t2m_daily_mean_monthly_mean = new((/ nm, ny, nx /), float)
      t2m_daily_mean_monthly_mean!0             = "month"
      t2m_daily_mean_monthly_mean!1             = "lat"
      t2m_daily_mean_monthly_mean!2             = "lon"
      t2m_daily_mean_monthly_mean@description   = "Mean Monthly Avg Daily 2-m Temp"
      t2m_daily_mean_monthly_mean@long_name     = "Mean Monthly Avg Daily 2-m Temp"
      t2m_daily_mean_monthly_mean@units         = "K"
      t2m_daily_mean_monthly_mean@standard_name = "air_temperature"
      t2m_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   t2m_daily_max_monthly_mean = t2m_daily_mean_monthly_mean
      delete(t2m_daily_max_monthly_mean@standard_name)
      t2m_daily_max_monthly_mean@description   = "Mean Monthly Max Daily 2-m Temp"
      t2m_daily_max_monthly_mean@long_name     = "Mean Monthly Max Daily 2-m Temp"
      t2m_daily_max_monthly_mean@units         = "K"
      t2m_daily_max_monthly_mean@standard_name = "air_temperature"
      t2m_daily_max_monthly_mean@cell_methods  = "hour: max; dom: mean"

   t2m_daily_min_monthly_mean = t2m_daily_mean_monthly_mean
      delete(t2m_daily_min_monthly_mean@standard_name)
      t2m_daily_min_monthly_mean@description   = "Mean Monthly Min Daily 2-m Temp"
      t2m_daily_min_monthly_mean@long_name     = "Mean Monthly Min Daily 2-m Temp"
      t2m_daily_min_monthly_mean@units         = "K"
      t2m_daily_min_monthly_mean@standard_name = "air_temperature"
      t2m_daily_min_monthly_mean@cell_methods  = "hour: min; dom: mean"

   t2m_daily_min_monthly_min = t2m_daily_mean_monthly_mean
      delete(t2m_daily_min_monthly_min@standard_name)
      t2m_daily_min_monthly_min@description   = "Minimum Monthly Daily 2-m Temp"
      t2m_daily_min_monthly_min@long_data     = "Minimum Monthly Daily 2-m Temp"
      t2m_daily_min_monthly_min@units         = "K"
      t2m_daily_min_monthly_min@standard_name = "air_temperature"
      t2m_daily_min_monthly_min@cell_methods  = "hour: min; dom: min"
      t2m_daily_min_monthly_min@crd_climate_change_index  = "TNn"

   t2m_daily_min_monthly_max = t2m_daily_mean_monthly_mean
      delete(t2m_daily_min_monthly_max@standard_name)
      t2m_daily_min_monthly_max@description    = "Maximum Monthly Min Daily 2-m Temp"
      t2m_daily_min_monthly_max@long_data      = "Maximum Monthly Min Daily 2-m Temp"
      t2m_daily_min_monthly_max@units         = "K"
      t2m_daily_min_monthly_max@standard_name = "air_temperature"
      t2m_daily_min_monthly_max@cell_methods  = "hour: min; dom: max"
      t2m_daily_min_monthly_max@crd_climate_change_index  = "TNx"

   t2m_daily_max_monthly_min = t2m_daily_mean_monthly_mean
      delete(t2m_daily_max_monthly_min@standard_name)
      t2m_daily_max_monthly_min@description   = "Minimum Monthly Max Daily 2-m Temp"
      t2m_daily_max_monthly_min@long_name     = "Minimum Monthly Max Daily 2-m Temp"
      t2m_daily_max_monthly_min@units         = "K"
      t2m_daily_max_monthly_min@standard_name = "air_temperature"
      t2m_daily_max_monthly_min@cell_methods  = "hour: max; dom: min"
      t2m_daily_max_monthly_min@crd_climate_change_index  = "TXn"

   t2m_daily_max_monthly_max = t2m_daily_mean_monthly_mean
      delete(t2m_daily_max_monthly_max@standard_name)
      t2m_daily_max_monthly_max@description   = "Maximum Monthly 2-m Temp"
      t2m_daily_max_monthly_max@long_name     = "Maximum Monthly 2-m Temp"
      t2m_daily_max_monthly_max@units         = "K"
      t2m_daily_max_monthly_max@standard_name = "air_temperature"
      t2m_daily_max_monthly_max@cell_methods  = "hour: max; dom: max"
      t2m_daily_max_monthly_max@crd_climate_change_index  = "TXx"

   t2m_daily_mean_daily_range = t2m_daily_mean_monthly_mean
      delete(t2m_daily_mean_daily_range@standard_name)
      t2m_daily_mean_daily_range@description   = "Maximum Monthly 2-m Temp Range"
      t2m_daily_mean_daily_range@long_name     = "Maximum Monthly 2-m Temp Range"
      t2m_daily_mean_daily_range@units         = "K"
      t2m_daily_mean_daily_range@standard_name = "air_temperature"
      t2m_daily_mean_daily_range@cell_methods  = "hour: max,min; dom: mean(max-min)"
      t2m_daily_mean_daily_range@crd_climate_change_index  = "DTR"



   nws_heat_index_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(nws_heat_index_daily_mean_monthly_mean@standard_name)
      nws_heat_index_daily_mean_monthly_mean@description   = "Mean Monthly Daily Avg NWS Heat Index"
      nws_heat_index_daily_mean_monthly_mean@long_name     = "Mean Monthly Daily Avg NWS Heat Index"
      nws_heat_index_daily_mean_monthly_mean@units         = "degF"
      nws_heat_index_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   nws_heat_index_daily_mean_monthly_max = t2m_daily_mean_monthly_mean
      delete(nws_heat_index_daily_mean_monthly_max@standard_name)
      nws_heat_index_daily_mean_monthly_max@description   = "Max Monthly Daily Mean NWS Heat Index"
      nws_heat_index_daily_mean_monthly_max@long_name     = "Max Monthly Daily Mean NWS Heat Index"
      nws_heat_index_daily_mean_monthly_max@units         = "degF"
      nws_heat_index_daily_mean_monthly_max@cell_methods  = "hour: mean; dom: mean"


   nws_heat_index_daily_max_monthly_max = t2m_daily_mean_monthly_mean
      delete(nws_heat_index_daily_max_monthly_max@standard_name)
      nws_heat_index_daily_max_monthly_max@description   = "Max Monthly Daily Max NWS Heat Index"
      nws_heat_index_daily_max_monthly_max@long_name     = "Max Monthly Daily Max NWS Heat Index"
      nws_heat_index_daily_max_monthly_max@units         = "degF"
      nws_heat_index_daily_max_monthly_max@cell_methods  = "hour: max; dom: max"


   e2m_deficit_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(e2m_deficit_daily_mean_monthly_mean@standard_name)
      e2m_deficit_daily_mean_monthly_mean!0             = "month"
      e2m_deficit_daily_mean_monthly_mean!1             = "lat"
      e2m_deficit_daily_mean_monthly_mean!2             = "lon"
      e2m_deficit_daily_mean_monthly_mean@description   = "Mean Monthly Avg 2-m Vap. Pres. Deficit"
      e2m_deficit_daily_mean_monthly_mean@long_name     = "Mean Monthly Avg 2-m Vap. Pres. Deficit"
      e2m_deficit_daily_mean_monthly_mean@units         = "Pa"
      e2m_deficit_daily_mean_monthly_mean@cell_methods  =  "hour: mean; dom: mean"

   q2m_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(q2m_daily_mean_monthly_mean@standard_name)
      q2m_daily_mean_monthly_mean!0             = "month"
      q2m_daily_mean_monthly_mean!1             = "lat"
      q2m_daily_mean_monthly_mean!2             = "lon"
      q2m_daily_mean_monthly_mean@description   = "Mean Monthly Avg 2-m Spec H."
      q2m_daily_mean_monthly_mean@long_name     = "Mean Monthly Avg 2-m Spec H."
      q2m_daily_mean_monthly_mean@units         = "kg3 kg-3"
      q2m_daily_mean_monthly_mean@standard_name = "specific_humidity"
      q2m_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   precip_daily_total_monthly_total = t2m_daily_mean_monthly_mean
      delete(precip_daily_total_monthly_total@standard_name)
      precip_daily_total_monthly_total@description   = "Total Monthly Precipitation"
      precip_daily_total_monthly_total@long_name   = "Total Monthly Precipitation"
      precip_daily_total_monthly_total@units         = "kg m-2"
      precip_daily_total_monthly_total@standard_name = "precipitation_amount"
      precip_daily_total_monthly_total@cell_methods  = "hour: sum; dom: sum"

   sfroff_daily_total_monthly_total = t2m_daily_mean_monthly_mean
      delete(sfroff_daily_total_monthly_total@standard_name)
      sfroff_daily_total_monthly_total@description   = "Total Monthly Surface Runoff"
      sfroff_daily_total_monthly_total@long_name     = "Total Monthly Surface Runoff"
      sfroff_daily_total_monthly_total@units         = "kg m-2"
      sfroff_daily_total_monthly_total@cell_methods  = "hour: sum; dom: sum"

   ugroff_daily_total_monthly_total = t2m_daily_mean_monthly_mean
      delete(ugroff_daily_total_monthly_total@standard_name)
      ugroff_daily_total_monthly_total@description   = "Total Monthly Soil Zone Drainage"
      ugroff_daily_total_monthly_total@long_name     = "Total Monthly Soil Zone Drainage"
      ugroff_daily_total_monthly_total@units         = "kg m-2"
      ugroff_daily_total_monthly_total@cell_methods  = "hour: sum; dom: sum"

   precip_daily_total_monthly_max = t2m_daily_mean_monthly_mean
      delete(precip_daily_total_monthly_max@standard_name)
      precip_daily_total_monthly_max@description   = "Monthly Max Daily Precipitation"
      precip_daily_total_monthly_max@long_name     = "Monthly Max Daily Precipitation"
      precip_daily_total_monthly_max@units         = "kg m-2"
      precip_daily_total_monthly_max@standard_name = "precipitation_amount"
      precip_daily_total_monthly_max@cell_methods  = "hour: sum; dom: max"
      precip_daily_total_monthly_max@crd_climate_change_index  = "Rx1day"

   precip_hourly_max_monthly_max = t2m_daily_mean_monthly_mean
      delete(precip_hourly_max_monthly_max@standard_name)
      precip_hourly_max_monthly_max@description   = "Monthly Max Hourly Precipitation"
      precip_hourly_max_monthly_max@long_name     = "Monthly Max Hourly Precipitation"
      precip_hourly_max_monthly_max@units         = "kg m-2"
      precip_hourly_max_monthly_max@standard_name = "precipitation_amount"
      precip_hourly_max_monthly_max@cell_methods  = "hour: max; dom: max"


    precip_nonzero_daily_total_monthly_mean = t2m_daily_mean_monthly_mean
      delete(precip_nonzero_daily_total_monthly_mean@standard_name)
      precip_nonzero_daily_total_monthly_mean@description   = "Monthly Mean Daily Non-Zero Precipitation"
      precip_nonzero_daily_total_monthly_mean@long_name     = "Monthly Mean Daily Non-Zero Precipitation"
      precip_nonzero_daily_total_monthly_mean@units         = "kg m-2"
      precip_nonzero_daily_total_monthly_mean@standard_name = "precipitation_amount"
      precip_nonzero_daily_total_monthly_mean@cell_methods  = "hour: total; dom: mean(P>1mm)"



    standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value = t2m_daily_mean_monthly_mean
      delete(standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@standard_name)
      standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@description   = "3-mo Standardized Precip-Evap Index (SPEI)"
      standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@long_name     = "3-mo Standardized Precip-Evap Index (SPEI)"
      delete(standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@units)
      standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@cell_methods  = "hour: mean/total; dom: mean"
      standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value@references = "http://sac.csic.es/spei/"


    standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value = t2m_daily_mean_monthly_mean
      delete(standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@standard_name)
      standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@description   = "6-mo Standardized Precip-Evap Index (SPEI)"
      standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@long_name     = "6-mo Standardized Precip-Evap Index (SPEI)"
      delete(standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@units)
      standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@cell_methods  = "hour: mean/total; dom: mean"
      standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value@references = "http://sac.csic.es/spei/"



    standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value = t2m_daily_mean_monthly_mean
      delete(standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@standard_name)
      standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@description   = "12-mo Standardized Precip-Evap Index (SPEI)"
      standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@long_name     = "12-mo Standardized Precip-Evap Index (SPEI)"
      delete(standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@units)
      standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@cell_methods  = "hour: mean/total; dom: mean"
      standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value@references = "http://sac.csic.es/spei/"



    standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value = t2m_daily_mean_monthly_mean
      delete(standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@standard_name)
      standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@description   = "24-mo Standardized Precip-Evap Index (SPEI)"
      standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@long_name     = "24-mo Standardized Precip-Evap Index (SPEI)"
      delete(standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@units)
      standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@cell_methods  = "hour: mean/total; dom: mean"
      standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value@references = "http://sac.csic.es/spei/"



    standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value = t2m_daily_mean_monthly_mean
      delete(standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@standard_name)
      standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@description   = "36-mo Standardized Precip-Evap Index (SPEI)"
      standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@long_name     = "36-mo Standardized Precip-Evap Index (SPEI)"
      delete(standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@units)
      standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@cell_methods  = "hour: mean/total; dom: mean"
      standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value@references = "http://sac.csic.es/spei/"


   length_of_dry_periods_daily_sum_monthly_mean = new((/ nm, ny, nx /), float)
      length_of_dry_periods_daily_sum_monthly_mean!0             = "month"
      length_of_dry_periods_daily_sum_monthly_mean!1             = "lat"
      length_of_dry_periods_daily_sum_monthly_mean!2             = "lon"
      length_of_dry_periods_daily_sum_monthly_mean@description   = "Monthly Mean Consecutive Dry Days"
      length_of_dry_periods_daily_sum_monthly_mean@long_name     = "Monthly Mean Consecutive Dry Days"
      length_of_dry_periods_daily_sum_monthly_mean@units         = "days"
      length_of_dry_periods_daily_sum_monthly_mean@cell_methods  = "hour: sum; dom: avgrun(P<1mm)"


   length_of_wet_periods_daily_sum_monthly_mean = new((/ nm, ny, nx /), float)
      length_of_wet_periods_daily_sum_monthly_mean!0             = "month"
      length_of_wet_periods_daily_sum_monthly_mean!1             = "lat"
      length_of_wet_periods_daily_sum_monthly_mean!2             = "lon"
      length_of_wet_periods_daily_sum_monthly_mean@description   = "Monthly Mean Consecutive Wet Days"
      length_of_wet_periods_daily_sum_monthly_mean@long_name     = "Monthly Mean Consecutive Wet Days"
      length_of_wet_periods_daily_sum_monthly_mean@units         = "days"
      length_of_wet_periods_daily_sum_monthly_mean@cell_methods  = "hour: sum; dom: avgrun(P>1mm)"



   simple_precipitation_intensity_index_monthly_value = t2m_daily_mean_monthly_mean
      delete(simple_precipitation_intensity_index_monthly_value@standard_name)
      simple_precipitation_intensity_index_monthly_value@description   = "Simple Precipitation Intensity Index"
      simple_precipitation_intensity_index_monthly_value@long_name     = "Simple Precipitation Intensity Index"
      simple_precipitation_intensity_index_monthly_value@units         = "kg m-2 dy-1"
      simple_precipitation_intensity_index_monthly_value@cell_methods  = "hour: max; dom: mean"
      simple_precipitation_intensity_index_monthly_value@comment1      = "eq: sum(prec)/count(prec_daily>1mm)"
      simple_precipitation_intensity_index_monthly_value@crd_climate_change_index  = "SDII"

   m10m_wind_speed_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(m10m_wind_speed_daily_mean_monthly_mean@standard_name)
      m10m_wind_speed_daily_mean_monthly_mean@description   = "Mean Monthly 10-m Wind Speed"
      m10m_wind_speed_daily_mean_monthly_mean@long_name     = "Mean Monthly 10-m Wind Speed"
      m10m_wind_speed_daily_mean_monthly_mean@units         = "m s-1"
      m10m_wind_speed_daily_mean_monthly_mean@standard_name = "wind_speed"
      m10m_wind_speed_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   u10m_wind_speed_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(u10m_wind_speed_daily_mean_monthly_mean@standard_name)
      u10m_wind_speed_daily_mean_monthly_mean@description   = "Westerly Mean Monthly 10-m Wind Speed"
      u10m_wind_speed_daily_mean_monthly_mean@long_name     = "Westerly Mean Monthly 10-m Wind Speed"
      u10m_wind_speed_daily_mean_monthly_mean@units         = "m s-1"
      u10m_wind_speed_daily_mean_monthly_mean@standard_name = "wind_speed"
      u10m_wind_speed_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   v10m_wind_speed_daily_mean_monthly_mean = t2m_daily_mean_monthly_mean
      delete(v10m_wind_speed_daily_mean_monthly_mean@standard_name)
      v10m_wind_speed_daily_mean_monthly_mean@description   = "Southerly Mean Monthly 10-m Wind Speed"
      v10m_wind_speed_daily_mean_monthly_mean@long_name     = "Southerly Mean Monthly 10-m Wind Speed"
      v10m_wind_speed_daily_mean_monthly_mean@units         = "m s-1"
      v10m_wind_speed_daily_mean_monthly_mean@standard_name = "wind_speed"
      v10m_wind_speed_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"




   ;
   ;  Monthly Depth Values
   ;

   swe_daily_mean_monthly_mean = new((/ nm, ny, nx /), float)
      swe_daily_mean_monthly_mean!0             = "month"
      swe_daily_mean_monthly_mean!1             = "lat"
      swe_daily_mean_monthly_mean!2             = "lon"
      swe_daily_mean_monthly_mean@description   = "Mean Monthly Snow Water Equiv."
      swe_daily_mean_monthly_mean@long_name    = "Mean Monthly Snow Water Equiv."
      swe_daily_mean_monthly_mean@units         = "kg m-2"
      swe_daily_mean_monthly_mean@standard_name = "liquid_water_content_of_surface_snow"
      swe_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   soil_mois_daily_mean_monthly_mean = new((/ nm, n_soil, ny, nx /), float)
      soil_mois_daily_mean_monthly_mean!0             = "month"
      soil_mois_daily_mean_monthly_mean!1             = "soil_depth"
      soil_mois_daily_mean_monthly_mean!2             = "lat"
      soil_mois_daily_mean_monthly_mean!3             = "lon"
      soil_mois_daily_mean_monthly_mean@description   = "Mean Monthly Percent Soil Saturation"
      soil_mois_daily_mean_monthly_mean@long_name     = "Mean Monthly Percent Soil Saturation"
      soil_mois_daily_mean_monthly_mean@units         = "fraction"
      soil_mois_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   vol_soil_wat_content_daily_mean_monthly_mean = soil_mois_daily_mean_monthly_mean
      vol_soil_wat_content_daily_mean_monthly_mean@description   = "Mean Monthly Vol. Soil Water Content"
      vol_soil_wat_content_daily_mean_monthly_mean@long_name     = "Mean Monthly Vol. Soil Water Content"
      vol_soil_wat_content_daily_mean_monthly_mean@units         = "m3 m-3"
      vol_soil_wat_content_daily_mean_monthly_mean@standard_name = "liquid_water_content_of_soil_layer"
      vol_soil_wat_content_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   soil_temp_daily_mean_monthly_mean = soil_mois_daily_mean_monthly_mean
      soil_temp_daily_mean_monthly_mean@description   = "Mean Monthly Soil Temperature"
      soil_temp_daily_mean_monthly_mean@long_name     = "Mean Monthly Soil Temperature"
      soil_temp_daily_mean_monthly_mean@units         = "K"
      soil_temp_daily_mean_monthly_mean@standard_name = "soil_temperature"
      soil_temp_daily_mean_monthly_mean@cell_methods  = "hour: mean; dom: mean"

   ;
   ;  Monthly Integer Values
   ;

   heating_degree_days_monthly_total =  new((/ nm, ny, nx /), integer)
      heating_degree_days_monthly_total!0             = "month"
      heating_degree_days_monthly_total!1             = "lat"
      heating_degree_days_monthly_total!2             = "lon"
      heating_degree_days_monthly_total@description   = "Monthly Heating Degree Days"
      heating_degree_days_monthly_total@long_name     = "Monthly Heating Degree Days"
      heating_degree_days_monthly_total@units         = "HDD"
      heating_degree_days_monthly_total@cell_methods  = "hour: mean; dom: accrued_excedence(avg_daily_t<291.48K)"

   cooling_degree_days_monthly_total = heating_degree_days_monthly_total
      cooling_degree_days_monthly_total@description   = "Monthly Cooling Degree Days"
      cooling_degree_days_monthly_total@long_name     = "Monthly Cooling Degree Days"
      cooling_degree_days_monthly_total@units         = "CDD"
      cooling_degree_days_monthly_total@cell_methods  = "hour: mean; dom: accrued_excedence(avg_daily_t>291.48K)"

   growing_degree_days_monthly_total = heating_degree_days_monthly_total
      growing_degree_days_monthly_total@description   = "Monthly Growing Degree Days"
      growing_degree_days_monthly_total@long_name   = "Monthly Growing Degree Days"
      growing_degree_days_monthly_total@units         = "GDD"
      growing_degree_days_monthly_total@cell_methods  = "hour: mean; dom: accrued_excedence(avg_daily_t>283.15K)"

   maize_growing_degree_days_monthly_total = heating_degree_days_monthly_total
      maize_growing_degree_days_monthly_total@description   = "Monthly Maize Growing Degree Days"
      maize_growing_degree_days_monthly_total@long_name     = "Monthly Maize Growing Degree Days"
      maize_growing_degree_days_monthly_total@units         = "MGDD"
      maize_growing_degree_days_monthly_total@cell_methods  = "hour: mean; dom: accrued_excedence(283.15K<avg_daily_t<303.15K)"

   number_of_frost_days_monthly_total =  heating_degree_days_monthly_total
      number_of_frost_days_monthly_total@description   = "Monthly Number of Frost Days"
      number_of_frost_days_monthly_total@long_name     = "Monthly Number of Frost Days"
      number_of_frost_days_monthly_total@units         = "days"
      number_of_frost_days_monthly_total@cell_methods  = "dom: count(TN<273.15)"
      number_of_frost_days_monthly_total@crd_climate_change_index  = "FD"

   number_of_summer_days_monthly_total =  heating_degree_days_monthly_total
      number_of_summer_days_monthly_total@description   = "Monthly Number of Summer Days"
      number_of_summer_days_monthly_total@long_name     = "Monthly Number of Summer Days"
      number_of_summer_days_monthly_total@units         = "days"
      number_of_summer_days_monthly_total@cell_methods  = "doy: count(TX>298.15K)"
      number_of_summer_days_monthly_total@crd_climate_change_index  = "SD"

   number_of_icing_days_monthly_total =  heating_degree_days_monthly_total
      number_of_icing_days_monthly_total@description   = "Monthly Number of Icing Days"
      number_of_icing_days_monthly_total@long_name     = "Monthly Number of Icing Days"
      number_of_icing_days_monthly_total@units         = "days"
      number_of_icing_days_monthly_total@cell_methods  = "dom: count(TX<273.15)"
      number_of_icing_days_monthly_total@crd_climate_change_index  = "ID"

   number_of_tropical_nights_monthly_total =  heating_degree_days_monthly_total
      number_of_tropical_nights_monthly_total@description   = "Monthly Number of Tropical Nights"
      number_of_tropical_nights_monthly_total@long_name     = "Monthly Number of Tropical Nights"
      number_of_tropical_nights_monthly_total@units         = "days"
      number_of_tropical_nights_monthly_total@cell_methods  = "dom: count(TN>293.15L)"
      number_of_tropical_nights_monthly_total@crd_climate_change_index  = "TN"

   number_of_days_precip_ge_010mm_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_precip_ge_010mm_monthly_total@description   = "Monthly Number of Days with Precip over 10 mm"
      number_of_days_precip_ge_010mm_monthly_total@long_name     = "Monthly Number of Days with Precip over 10 mm"
      number_of_days_precip_ge_010mm_monthly_total@units         = "days"
      number_of_days_precip_ge_010mm_monthly_total@cell_methods  = "dom: count(P>10mm)"

   number_of_days_precip_ge_020mm_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_precip_ge_020mm_monthly_total@description   = "Monthly Number of Days with Precip over 20 mm"
      number_of_days_precip_ge_020mm_monthly_total@long_name     = "Monthly Number of Days with Precip over 20 mm"
      number_of_days_precip_ge_020mm_monthly_total@units         = "days"
      number_of_days_precip_ge_020mm_monthly_total@cell_methods  = "doy: count(P>20mm)"

   number_of_days_precip_ge_050mm_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_precip_ge_050mm_monthly_total@description   = "Monthly Number of Days with Precip over 50 mm"
      number_of_days_precip_ge_050mm_monthly_total@long_name     = "Monthly Number of Days with Precip over 50 mm"
      number_of_days_precip_ge_050mm_monthly_total@units         = "days"
      number_of_days_precip_ge_050mm_monthly_total@cell_methods  = "dom: count(P>50mm)"
        
   number_of_days_precip_ge_075mm_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_precip_ge_075mm_monthly_total@description   = "Monthly Number of Days with Precip over 75 mm"
      number_of_days_precip_ge_075mm_monthly_total@long_name     = "Monthly Number of Days with Precip over 75 mm"
      number_of_days_precip_ge_075mm_monthly_total@units         = "days"
      number_of_days_precip_ge_075mm_monthly_total@cell_methods  = "dom: count(P>75mm)"

   number_of_days_precip_ge_100mm_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_precip_ge_100mm_monthly_total@description   = "Monthly Number of Days with Precip over 100 mm"
      number_of_days_precip_ge_100mm_monthly_total@long_name     = "Monthly Number of Days with Precip over 100 mm"
      number_of_days_precip_ge_100mm_monthly_total@units         = "days"
      number_of_days_precip_ge_100mm_monthly_total@cell_methods  = "dom: count(P>100mm)"


   length_of_dry_periods_daily_sum_monthly_max = new((/ nm, ny, nx /), integer)
      length_of_dry_periods_daily_sum_monthly_max!0             = "month"
      length_of_dry_periods_daily_sum_monthly_max!1             = "lat"
      length_of_dry_periods_daily_sum_monthly_max!2             = "lon"
      length_of_dry_periods_daily_sum_monthly_max@description   = "Monthly Max Consecutive Dry Days"
      length_of_dry_periods_daily_sum_monthly_max@long_name     = "Monthly Max Consecutive Dry Days"
      length_of_dry_periods_daily_sum_monthly_max@units         = "days"
      length_of_dry_periods_daily_sum_monthly_max@cell_methods  = "hour: sum; dom: maxrun(P<1mm)"
      length_of_dry_periods_daily_sum_monthly_max@crd_climate_change_index  = "CDD"


   length_of_wet_periods_daily_sum_monthly_max = new((/ nm, ny, nx /), integer)
      length_of_wet_periods_daily_sum_monthly_max!0             = "month"
      length_of_wet_periods_daily_sum_monthly_max!1             = "lat"
      length_of_wet_periods_daily_sum_monthly_max!2             = "lon"
      length_of_wet_periods_daily_sum_monthly_max@description   = "Monthly Max Consecutive Wet Days"
      length_of_wet_periods_daily_sum_monthly_max@long_name     = "Monthly Max Consecutive Wet Days"
      length_of_wet_periods_daily_sum_monthly_max@units         = "days"
      length_of_wet_periods_daily_sum_monthly_max@cell_methods  = "hour: sum; dom: maxrun(P>1mm)"
      length_of_wet_periods_daily_sum_monthly_max@crd_climate_change_index  = "CWD"

   number_of_days_with_snow_cover_monthly_total  =  heating_degree_days_monthly_total
      number_of_days_with_snow_cover_monthly_total@description   = "Monthly Number of Days with Snow Depth over 10 mm"
      number_of_days_with_snow_cover_monthly_total@long_name     = "Monthly Number of Days with Snow Depth over 10 mm"
      number_of_days_with_snow_cover_monthly_total@units         = "days"
      number_of_days_with_snow_cover_monthly_total@cell_methods  = "dom: count(SNOWH>100mm)"




   number_of_days_max_heat_index_080_to_090_monthly_total =  heating_degree_days_monthly_total
      number_of_days_max_heat_index_080_to_090_monthly_total@description   = "Monthly Number of Days HI between 80F-90F"
      number_of_days_max_heat_index_080_to_090_monthly_total@long_name     = "Monthly Number of Days HI between 80F-90F"
      number_of_days_max_heat_index_080_to_090_monthly_total@units         = "days"
      number_of_days_max_heat_index_080_to_090_monthly_total@cell_methods  = "doy: count(80F<HiMax<90F)"
      number_of_days_max_heat_index_080_to_090_monthly_total@heat_index_flag  = "Caution"



   number_of_days_max_heat_index_090_to_105_monthly_total =  heating_degree_days_monthly_total
      number_of_days_max_heat_index_090_to_105_monthly_total@description   = "Monthly Number of Days HI between 90F-105F"
      number_of_days_max_heat_index_090_to_105_monthly_total@long_name     = "Monthly Number of Days HI between 90F-105F"
      number_of_days_max_heat_index_090_to_105_monthly_total@units         = "days"
      number_of_days_max_heat_index_090_to_105_monthly_total@cell_methods  = "doy: count(90F<HiMax<105F)"
      number_of_days_max_heat_index_090_to_105_monthly_total@heat_index_flag  = "Extreme Caution"



   number_of_days_max_heat_index_105_to_130_monthly_total =  heating_degree_days_monthly_total
      number_of_days_max_heat_index_105_to_130_monthly_total@description   = "Monthly Number of Days HI between 105F-130F"
      number_of_days_max_heat_index_105_to_130_monthly_total@long_name     = "Monthly Number of Days HI between 105F-130F"
      number_of_days_max_heat_index_105_to_130_monthly_total@units         = "days"
      number_of_days_max_heat_index_105_to_130_monthly_total@cell_methods  = "doy: count(105F<HiMax<130F)"
      number_of_days_max_heat_index_105_to_130_monthly_total@heat_index_flag  = "Danger"



   number_of_days_max_heat_index_over_130_monthly_total =  heating_degree_days_monthly_total
      number_of_days_max_heat_index_over_130_monthly_total@description   = "Monthly Number of Days HI over 130F"
      number_of_days_max_heat_index_over_130_monthly_total@long_name     = "Monthly Number of Days HI over 130F"
      number_of_days_max_heat_index_over_130_monthly_total@units         = "days"
      number_of_days_max_heat_index_over_130_monthly_total@cell_methods  = "doy: count(HiMax>130F)"
      number_of_days_max_heat_index_over_130_monthly_total@heat_index_flag  = "Danger"




   ;
   ;  Annual Floating Point Values
   ;

   precip_daily_total_annual_total =  new((/ nyr, ny, nx /), float)
      precip_daily_total_annual_total!0             = "year"
      precip_daily_total_annual_total!1             = "lat"
      precip_daily_total_annual_total!2             = "lon"
      precip_daily_total_annual_total@standard_name = "precipitation_amount"
      precip_daily_total_annual_total@description   = "Annual Total Precip"
      precip_daily_total_annual_total@long_name     = "Annual Total Precip"
      precip_daily_total_annual_total@units         = "kg m-2"
      precip_daily_total_annual_total@cell_methods  = "hour: sum; doy: total"
      precip_daily_total_annual_total@crd_climate_change_index  = "PRCPTOT"

   t2m_daily_mean_annual_mean =  precip_daily_total_annual_total
      delete(t2m_daily_mean_annual_mean@standard_name)
      t2m_daily_mean_annual_mean@description   = "Annual Mean Daily 2-m Temperature"
      t2m_daily_mean_annual_mean@long_name     = "Annual Mean Daily 2-m Temperature"
      t2m_daily_mean_annual_mean@units         = "K"
      t2m_daily_mean_annual_mean@standard_name = "air_temperature"
      t2m_daily_mean_annual_mean@cell_methods  = "hour: mean; doy: mean"



   length_of_dry_periods_daily_sum_annual_mean = new((/ nyr, ny, nx /), float)
      length_of_dry_periods_daily_sum_annual_mean!0             = "year"
      length_of_dry_periods_daily_sum_annual_mean!1             = "lat"
      length_of_dry_periods_daily_sum_annual_mean!2             = "lon"
      length_of_dry_periods_daily_sum_annual_mean@description   = "Annual Mean Consecutive Dry Days"
      length_of_dry_periods_daily_sum_annual_mean@long_name    = "Annual Mean Consecutive Dry Days"
      length_of_dry_periods_daily_sum_annual_mean@units         = "days"
      length_of_dry_periods_daily_sum_annual_mean@cell_methods  = "hour: sum; doy: maxrun(P<1mm)"


   length_of_wet_periods_daily_sum_annual_mean = new((/ nyr, ny, nx /), float)
      length_of_wet_periods_daily_sum_annual_mean!0             = "year"
      length_of_wet_periods_daily_sum_annual_mean!1             = "lat"
      length_of_wet_periods_daily_sum_annual_mean!2             = "lon"
      length_of_wet_periods_daily_sum_annual_mean@description   = "Annual Mean Consecutive Wet Days"
      length_of_wet_periods_daily_sum_annual_mean@long_name     = "Annual Mean Consecutive Wet Days"
      length_of_wet_periods_daily_sum_annual_mean@units         = "days"
      length_of_wet_periods_daily_sum_annual_mean@cell_methods  = "hour: sum; doy: maxrun(P>1mm)"




   ;
   ;  Annual Integer Values
   ;



   length_of_dry_periods_daily_sum_annual_max = new((/ nyr, ny, nx /), integer)
      length_of_dry_periods_daily_sum_annual_max!0             = "year"
      length_of_dry_periods_daily_sum_annual_max!1             = "lat"
      length_of_dry_periods_daily_sum_annual_max!2             = "lon"
      length_of_dry_periods_daily_sum_annual_max@description   = "Annual Max Consecutive Dry Days"
      length_of_dry_periods_daily_sum_annual_max@long_name    = "Annual Max Consecutive Dry Days"
      length_of_dry_periods_daily_sum_annual_max@units         = "days"
      length_of_dry_periods_daily_sum_annual_max@cell_methods  = "hour: sum; doy: maxrun(P<1mm)"
      length_of_dry_periods_daily_sum_annual_max@crd_climate_change_index  = "CDD"


   length_of_wet_periods_daily_sum_annual_max = new((/ nyr, ny, nx /), integer)
      length_of_wet_periods_daily_sum_annual_max!0             = "year"
      length_of_wet_periods_daily_sum_annual_max!1             = "lat"
      length_of_wet_periods_daily_sum_annual_max!2             = "lon"
      length_of_wet_periods_daily_sum_annual_max@description   = "Annual Max Consecutive Wet Days"
      length_of_wet_periods_daily_sum_annual_max@long_name     = "Annual Max Consecutive Wet Days"
      length_of_wet_periods_daily_sum_annual_max@units         = "days"
      length_of_wet_periods_daily_sum_annual_max@cell_methods  = "hour: sum; doy: maxrun(P>1mm)"
      length_of_wet_periods_daily_sum_annual_max@crd_climate_change_index  = "CWD"




   number_of_frost_days_annual_total =  new((/ nyr, ny, nx /), integer)
      number_of_frost_days_annual_total!0             = "year"
      number_of_frost_days_annual_total!1             = "lat"
      number_of_frost_days_annual_total!2             = "lon"
      number_of_frost_days_annual_total@description   = "Annual Number of Frost Days"
      number_of_frost_days_annual_total@long_name     = "Annual Number of Frost Days"
      number_of_frost_days_annual_total@units         = "days"
      number_of_frost_days_annual_total@cell_methods  = "doy: count(TN<273.15K)"
      number_of_frost_days_annual_total@crd_climate_change_index  = "FD"

   number_of_summer_days_annual_total =  number_of_frost_days_annual_total
      number_of_summer_days_annual_total@description   = "Annual Number of Summer Days"
      number_of_summer_days_annual_total@long_name     = "Annual Number of Summer Days"
      number_of_summer_days_annual_total@units         = "days"
      number_of_summer_days_annual_total@cell_methods  = "doy: count(TX>298.15K)"
      number_of_summer_days_annual_total@crd_climate_change_index  = "SD"


   heating_degree_days_annual_total = number_of_frost_days_annual_total
      heating_degree_days_annual_total@description   = "Annual Heating Degree Days"
      heating_degree_days_annual_total@long_name     = "Annual Heating Degree Days"
      heating_degree_days_annual_total@units         = "HDD"
      heating_degree_days_annual_total@cell_methods  = "hour: mean; doy: accrued_excedence(avg_daily_t<291.48K)"

   cooling_degree_days_annual_total = number_of_frost_days_annual_total
      cooling_degree_days_annual_total@description   = "Annual Cooling Degree Days"
      cooling_degree_days_annual_total@long_name     = "Annual Cooling Degree Days"
      cooling_degree_days_annual_total@units         = "CDD"
      cooling_degree_days_annual_total@cell_methods  = "hour: mean; doy: accrued_excedence(avg_daily_t>291.48K)"



   growing_degree_days_annual_total = number_of_frost_days_annual_total
      growing_degree_days_annual_total@description   = "Annual Growing Degree Days"
      growing_degree_days_annual_total@long_name     = "Annual Growing Degree Days"
      growing_degree_days_annual_total@units         = "GDD"
      growing_degree_days_annual_total@cell_methods  = "hour: mean; doy: accrued_excedence(avg_daily_t>283.15K)"

   maize_growing_degree_days_annual_total = number_of_frost_days_annual_total
      maize_growing_degree_days_annual_total@description   = "Annual Maize Growing Degree Days"
      maize_growing_degree_days_annual_total@long_name     = "Annual Maize Growing Degree Days"
      maize_growing_degree_days_annual_total@units         = "MGDD"
      maize_growing_degree_days_annual_total@cell_methods  = "hour: mean; doy: accrued_excedence(283.15K<avg_daily_t<303.15K)"



   number_of_icing_days_annual_total =  number_of_frost_days_annual_total
      number_of_icing_days_annual_total@description   = "Annual Number of Icing Days"
      number_of_icing_days_annual_total@long_name     = "Annual Number of Icing Days"
      number_of_icing_days_annual_total@units         = "days"
      number_of_icing_days_annual_total@cell_methods  = "doy: count(TX<273.15K)"
      number_of_icing_days_annual_total@crd_climate_change_index  = "ID"

   number_of_tropical_nights_annual_total =  number_of_frost_days_annual_total
      number_of_tropical_nights_annual_total@description   = "Annual Number of Tropical Nights"
      number_of_tropical_nights_annual_total@long_name     = "Annual Number of Tropical Nights"
      number_of_tropical_nights_annual_total@units         = "days"
      number_of_tropical_nights_annual_total@cell_methods  = "doy: count(TN>293.15K)"
      number_of_tropical_nights_annual_total@crd_climate_change_index  = "TN"

   number_of_days_precip_ge_010mm_annual_total  =  number_of_frost_days_annual_total
      number_of_days_precip_ge_010mm_annual_total@description   = "Annual Number of Days with Precip over 10 mm"
      number_of_days_precip_ge_010mm_annual_total@long_name     = "Annual Number of Days with Precip over 10 mm"
      number_of_days_precip_ge_010mm_annual_total@units         = "days"
      number_of_days_precip_ge_010mm_annual_total@cell_methods  = "doy: count(P>010mm)"
      number_of_days_precip_ge_010mm_annual_total@crd_climate_change_index  = "R10mm"

   number_of_days_precip_ge_020mm_annual_total  =  number_of_frost_days_annual_total
      number_of_days_precip_ge_020mm_annual_total@description   = "Annual Number of Days with Precip over 20 mm"
      number_of_days_precip_ge_020mm_annual_total@long_name     = "Annual Number of Days with Precip over 20 mm"
      number_of_days_precip_ge_020mm_annual_total@units         = "days"
      number_of_days_precip_ge_020mm_annual_total@cell_methods  = "doy: count(P>20mm)"
      number_of_days_precip_ge_020mm_annual_total@crd_climate_change_index  = "R20mm"


   number_of_days_precip_ge_050mm_annual_total  =  number_of_frost_days_annual_total
      number_of_days_precip_ge_050mm_annual_total@description   = "Annual Number of Days with Precip over 50 mm"
      number_of_days_precip_ge_050mm_annual_total@long_name     = "Annual Number of Days with Precip over 50 mm"
      number_of_days_precip_ge_050mm_annual_total@units         = "days"
      number_of_days_precip_ge_050mm_annual_total@cell_methods  = "doy: count(P>50mm)"
      number_of_days_precip_ge_050mm_annual_total@crd_climate_change_index  = "R50mm"

 
   number_of_days_precip_ge_075mm_annual_total  =  number_of_frost_days_annual_total
      number_of_days_precip_ge_075mm_annual_total@description   = "Annual Number of Days with Precip over 75 mm"
      number_of_days_precip_ge_075mm_annual_total@long_name     = "Annual Number of Days with Precip over 75 mm"
      number_of_days_precip_ge_075mm_annual_total@units         = "days"
      number_of_days_precip_ge_075mm_annual_total@cell_methods  = "doy: count(P>75mm)"
      number_of_days_precip_ge_075mm_annual_total@crd_climate_change_index  = "R75mm"


   number_of_days_precip_ge_100mm_annual_total  =  number_of_frost_days_annual_total
      number_of_days_precip_ge_100mm_annual_total@description   = "Annual Number of Days with Precip over 100 mm"
      number_of_days_precip_ge_100mm_annual_total@long_name     = "Annual Number of Days with Precip over 100 mm"
      number_of_days_precip_ge_100mm_annual_total@units         = "days"
      number_of_days_precip_ge_100mm_annual_total@cell_methods  = "doy: count(P>100mm)"
      number_of_days_precip_ge_100mm_annual_total@crd_climate_change_index  = "R100mm"


   length_of_growing_season_annual_total =  number_of_frost_days_annual_total
      length_of_growing_season_annual_total@description   = "Length of Growing Season"
      length_of_growing_season_annual_total@long_name     = "Length of Growing Season"
      length_of_growing_season_annual_total@units         = "days"
      length_of_growing_season_annual_total@cell_methods  = "doy: maxrun(avg_daily_temp>278.15)>12"
      length_of_growing_season_annual_total@crd_climate_change_index  = "GSL"



   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  Annual and Daily Land Mask
   ;

   mask2d                  = tofloat(temp_mask(:,:))
         mask2d@_FillValue = default_fillvalue("float")
         mask2d            = where(mask2d .eq. 0,     \
                                         mask2d@_FillValue, \
                                         1)

   imask2d                  = temp_mask(:,:)
         imask2d@_FillValue = default_fillvalue("integer")
         imask2d            = where(imask2d .eq. 0,         \
                                    imask2d@_FillValue,   \
                                    1)

   mask3d                  = tofloat(conform_dims((/n_soil, ny, nx /), temp_mask(:,:), (/1,2/)))
         mask3d@_FillValue = default_fillvalue("float")
         mask3d            = where(mask3d .eq. 0,       \
                                   mask3d@_FillValue,   \
                                   1 )
   ;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  Begin Processing with Scenario Loop
   ;
   ;

   ;
   ; create dummy array to hold last days's final rain
   ;

   hour_wrf = ispan(0, 18, 6)

   n_wrfhours = dimsizes(hour_wrf)

   last_rain   = new( (/         24, ny,nx /), float)  ; final total rain from "yesterday"
   last_sfroff = new( (/             ny,nx /), float)
   last_udroff = new( (/             ny,nx /), float)

   ldom_calendar = toint((/ -1,  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /))


   do s = 0, 0; n_scenarios-1

      ;
      ; Counters and Buckets Reset for new simulation
      ;

      mm  = toint(0)
      yy  = toint(0)

      last_rain(:,:,:)    = 0
      last_sfroff(:,:)    = 0
      last_udroff(:,:)    = 0


      
      heating_degree_days_monthly_total(:,:,:)       = 0
      cooling_degree_days_monthly_total(:,:,:)       = 0

      growing_degree_days_monthly_total(:,:,:)       = 0
      maize_growing_degree_days_monthly_total(:,:,:) = 0

    
      heating_degree_days_annual_total(:,:,:)       = 0
      cooling_degree_days_annual_total(:,:,:)       = 0

      growing_degree_days_annual_total(:,:,:)       = 0
      maize_growing_degree_days_annual_total(:,:,:) = 0


      ;
      ; Time Control
      ;

      missing_aux_days = ""

      year_start = year_starts(s)
      year_end   = year_ends(s)

      start_time = cd_inv_calendar(year_start, 01, 01, 00, 00, 00, time_units, calendar)
      end_time   = cd_inv_calendar(  year_end, 12, 31, 00, 00, 00, time_units, calendar)

      time = ispan(toint(start_time),toint(end_time),1)
         time@units = time_units
         time@calendar = calendar@calendar

      nt   = dimsizes(time)

      month = new(nm, integer)
         month@description   = "month"
         month@standard_name = "time"
         delete(month@_FillValue)

      year = ispan(year_start,year_end,1)
         year@description   = "Year"
         year@long_name     = "Year"
         year@standard_name = "time"
         year!0             = "year"
         year@units         = year_units
         ;year@calendar     = calendar@calendar


      calendar_year_in_monthly_data = new(nm, integer)
         calendar_year_in_monthly_data@description = "Calendar Year In Monthly Records"
         calendar_year_in_monthly_data@long_name   = "Calendar Year In Monthly Records"
         calendar_year_in_monthly_data!0           = "month"

 
      calendar_month_in_monthly_data = new(nm, integer)
         calendar_year_in_monthly_data@description    = "Calendar Month In Monthly Records"
         calendar_month_in_monthly_data@long_name     = "Calendar Month In Monthly Records"
         calendar_month_in_monthly_data!0             = "month"



      nnm = toint(0)
      do yyyy = year_start, year_end
         do mm = 1, 12
            calendar_year_in_monthly_data(nnm) = yyyy
            calendar_month_in_monthly_data(nnm) = mm
            month(nnm) = toint(cd_inv_calendar(  yyyy, mm, 00, 00, 00, 00, time_units, calendar))
            nnm = nnm + 1
         end do
      end do
      month!0             = "month"
      month@units         = time_units
      month@calendar      = calendar@calendar

      ;
      ; Loop through Simulations (daily time step)
      ;


      mm = 0
      do t = 0, nt-1

          ;
          ; Calendar Control
          ;

          YMDHS = cd_calendar(time(t),-5)

          calendar_day    = YMDHS(0,2)
          calendar_month  = YMDHS(0,1)
          calendar_year   = YMDHS(0,0)
          doy             = day_of_year(calendar_year,calendar_month,calendar_day)
          last_day        = ldom_calendar(calendar_month)

          ndy_year = 365 ; day_of_year(calendar_year, 12, 31)


          ;
          ; Reset Monthly Aggregations (Hourly Intervals)
          ;

          if (calendar_day .eq. 1) then

             max_temp_daily = new( (/ last_day, ny, nx /), float )
             min_temp_daily = new( (/ last_day, ny, nx /), float )
             avg_temp_daily = new( (/ last_day, ny, nx /), float )
             sum_prec_daily = new( (/ last_day, ny, nx /), float )
             max_prec_daily = new( (/ last_day, ny, nx /), float )

             avg_uwind_daily = new( (/ last_day, ny, nx /), float )
             avg_vwind_daily = new( (/ last_day, ny, nx /), float )
             avg_mwind_daily = new( (/ last_day, ny, nx /), float )

             avg_vpd_daily   = new( (/ last_day, ny, nx /), float )
             avg_sphum_daily = new( (/ last_day, ny, nx /), float )

             sum_sfroff_daily = new( (/ last_day, ny, nx /), float )
             sum_udroff_daily = new( (/ last_day, ny, nx /), float )
             avg_swe_daily    = new( (/ last_day, ny, nx /), float )

             avg_hi_daily    = new( (/ last_day, ny, nx /), float )
             max_hi_daily    = new( (/ last_day, ny, nx /), float )
             min_hi_daily    = new( (/ last_day, ny, nx /), float )


             snowh_max        = new( (/ last_day, ny, nx /), float )
             avg_smois_daily  = new( (/ last_day, n_soil, ny, nx /), float )
             avg_vswc_daily   = new( (/ last_day, n_soil, ny, nx /), float )
             avg_stemp_daily  = new( (/ last_day, n_soil, ny, nx /), float )

          end if

          ;
          ; Reset Annual Aggregations
          ;

          if ((calendar_day .eq. 1).and.(calendar_month .eq. 1)) then

             max_temp_daily_ann = new( (/ ndy_year, ny, nx /), float )
             min_temp_daily_ann = new( (/ ndy_year, ny, nx /), float )
             avg_temp_daily_ann = new( (/ ndy_year, ny, nx /), float )
             sum_prec_daily_ann = new( (/ ndy_year, ny, nx /), float )




          end if

          ;
          ; File Control for Aux Files
          ;

          aux_product_time_string = sprinti("%0.4i-",         YMDHS(0,0)) + \
                                    sprinti("%0.2i-",         YMDHS(0,1)) + \
                                    sprinti("%0.2i_00:00:00", YMDHS(0,2))

          aux_model_file = root_dir             + \
                           scenario(s)          + \
                           "/aux/auxhist1_d01_" + \
                           aux_product_time_string


          print("===> " + calendar_year  + "-"   + \
                      calendar_month + "-"   + \
                      calendar_day   + "/"   + \
                      last_day       + " : " + \
                      mm             + "/"   + \
                      (nm-1)         + " : " + \
                      yy             + "/"   + \
                      (nyr-1)  )
              me =    systemfunc("ls -al "+aux_model_file)
              print(" "+me)
              delete(me)


          is_aux_here = isfilepresent(aux_model_file)

          if (is_aux_here) then

               af = addfile(aux_model_file, "r")

               Time_s := chartostring( af->Times )

               nt_auxes = dimsizes(Time_s)

               if (nt_auxes .lt. 24) then

                  print("Missing Time Steps in AUX File! nt="+nt_auxes)

                  missing_aux_days = missing_aux_days + " " +aux_model_file

               else

                  ;
                  ; Process Hourly Accumulated Buckets
                  ;

                  I_RAINC  = af->I_RAINC( :, lat_s:lat_n, lon_w:lon_e)
                  I_RAINNC = af->I_RAINNC(:, lat_s:lat_n, lon_w:lon_e)
                  RAINC    = af->RAINC(   :, lat_s:lat_n, lon_w:lon_e)
                  RAINNC   = af->RAINNC(  :, lat_s:lat_n, lon_w:lon_e)

                  ;
                  ; xxx2 are the original fields from the file -- the last timestep
                  ;     is extracted later down pipe to hold the last accumuated
                  ;     bucket for "yesturday"

                  rain2 = (/ I_RAINC  * bucket_mm + RAINC  /)  +  \
                          (/ I_RAINNC * bucket_mm + RAINNC /)

                  rain = (/ I_RAINC  * bucket_mm + RAINC  /)  +  \
                          (/ I_RAINNC * bucket_mm + RAINNC /)  -  \
                          (/ last_rain /)

                  rain(1:23,:,:) = (/ rain(1:23,:,:) - rain(0:22,:,:) /)

                  last_rain =  conform_dims((/24, ny, nx /), rain2(23,:,:), (/1,2/))

                  ;
                  ; Process Hourly Un-accumulated Fields
                  ;

                  t2m  = (/   af->T2(:, lat_s:lat_n, lon_w:lon_e)  /)  ; 2-m air temp
                  q2m  = (/   af->Q2(:, lat_s:lat_n, lon_w:lon_e)  /)  ; 2-m spec hum
                  psfc = (/ af->PSFC(:, lat_s:lat_n, lon_w:lon_e)  /)  ; 2-m pres

                  e2m  =  q2m * psfc / ( 0.622 + q2m  * 0.378)

                  es2m =  (/ 610.94 * exp( 17.625 * (t2m-273.15) / ( (t2m-273.15) + 243.04)   ) /)

                  rh2m = (/ (e2m/es2m) * 100.0 /)

                  edef2m = es2m - e2m

                  es2m = where(es2m .ge. 0, es2m, 0)

                  rh2m = where(rh2m.gt.100,100,rh2m)
                  rh2m = where(rh2m.lt.0, 0,rh2m)

                  u10m  = (/  af->U10(:, lat_s:lat_n, lon_w:lon_e) /)   ; 10-m u-wind component
                  v10m  = (/  af->V10(:, lat_s:lat_n, lon_w:lon_e) /)   ; 10-m v-wind component
                  m10m  = (/ sqrt( u10m * u10m +  v10m * v10m ) /) ; 10-m wind vector magnitude

                  delete(af)
               end if

          else

               print("Missing AUX File!")

          end if


          ;
          ; 3-hr File Inner Loop Goes Here When Ready
          ;    Design:  Create a loop from hh = 0, 21, 3
          ;                  File Control
          ;                  Input 3-hrly Field Here.
          ;

          sfroff = new( (/ n_wrfhours,         ny, nx /), float )
          udroff = new( (/ n_wrfhours,         ny, nx /), float )
          swe    = new( (/ n_wrfhours,         ny, nx /), float )
          smois  = new( (/ n_wrfhours, n_soil, ny, nx /), float )
          vswc   = new( (/ n_wrfhours, n_soil, ny, nx /), float )
          stemp  = new( (/ n_wrfhours, n_soil, ny, nx /), float )



          snowh  = new( (/ n_wrfhours,         ny, nx /), float )

          swe   = where(swe   .lt.0, swe@_FillValue,   swe)
          vswc  = where(vswc  .lt.0, vswc@_FillValue,  vswc)
          stemp = where(stemp .lt.0, stemp@_FillValue, stemp)
          smois = where(smois .lt.0, smois@_FillValue, smois)



          do hh = 0, n_wrfhours-1


             wrf_product_time_string = sprinti("%0.4i-",        YMDHS(0,0)) + \
                                       sprinti("%0.2i-",        YMDHS(0,1)) + \
                                       sprinti("%0.2i_",        YMDHS(0,2)) + \
                                       sprinti("%0.2i:00:00", hour_wrf(hh))

             wrf_model_file = root_dir                + \
                              scenario(s)             + \
                              "/plevs/wrfout_d01_"    + \
                              wrf_product_time_string + \
                              "_PLEV"

             me =    systemfunc("ls -al "+wrf_model_file)



             print(" "+me)
             delete(me)

             if (.not. isfilepresent(wrf_model_file+".nc")) then
                print("!!!!!!  Warning File is not supported  !!!!!! ")
             end if




             wf = addfile(wrf_model_file+".nc", "r")

             if (.not. isfilevar(wf,"SFROFF")) then
                print("!!!!!!  SFROFF is Missing in File  !!!!!! ")
             end if





             if ((nt_auxes .lt. 24) .or. (.not. is_aux_here)) then

                 rain  = new( (/ 24, ny, nx /), float )
                 rain2 = new( (/ 24, ny, nx /), float )

                 t2m  = new( (/ 24, ny, nx /), float ) ; 2-m air temp
                 q2m  = new( (/ 24, ny, nx /), float ) ; 2-m spec hum
                 psfc = new( (/ 24, ny, nx /), float ) ; 2-m pres
                 u10m = new( (/ 24, ny, nx /), float ) ;
                 v10m = new( (/ 24, ny, nx /), float ) ;

                 do tttt = 0, 6-1

                    print("-------> downscaling by time wrf to aux")
                    print("          hh,t=" + hh + "," + tttt)
                    print("          rain")
                    printVarSummary(rain((/ hh_aux_plevl(hh,tttt) /), :, :))
                    print("          wf->RAINC(0,lat_s:lat_n, lon_w:lon_e)")
                    printVarSummary(wf->RAINC(0,lat_s:lat_n, lon_w:lon_e))

                    if (isfilevar(wf, "RAINSC")) then
                       rain((/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->RAINC(0,lat_s:lat_n, lon_w:lon_e)  +  \
                                                                    wf->RAINSC(0,lat_s:lat_n, lon_w:lon_e) +  \
                                                                    wf->RAINNC(0,lat_s:lat_n, lon_w:lon_e)    /)
                    else
                       rain((/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->RAINC(0,lat_s:lat_n, lon_w:lon_e)  +  \
                                                                    wf->RAINNC(0,lat_s:lat_n, lon_w:lon_e)    /)
                    end if


                    rain((/ hh_aux_plevl(hh,tttt) /), :, :) = (/ rain((/ hh_aux_plevl(hh,tttt) /), :, :) / 6.0 /)


                    t2m( (/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->T2(0,lat_s:lat_n, lon_w:lon_e) /)

                    q2m( (/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->Q2(0,lat_s:lat_n, lon_w:lon_e) /)

                    psfc((/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->PSFC(0,lat_s:lat_n, lon_w:lon_e) /)


                    u10m( (/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->U10(0,lat_s:lat_n, lon_w:lon_e) /)

                    v10m( (/ hh_aux_plevl(hh,tttt) /), :, :) = (/ wf->V10(0,lat_s:lat_n, lon_w:lon_e) /)

                 end do

                 m10m  = (/ sqrt( u10m * u10m + v10m * v10m ) /) ; 10-m wind vector magnitude

                 e2m  =  q2m * psfc / ( 0.622 + q2m  * 0.378)

                 es2m =  (/ 610.94 * exp( 17.625 * (t2m-273.15) / ( (t2m-273.15) + 243.04)   ) /)

                 rh2m = (/ (e2m/es2m) * 100.0 /)

                 edef2m = es2m - e2m

                 es2m = where(es2m .ge. 0, es2m, 0)

                 rh2m = where(rh2m .gt. 100, 100, rh2m)
                 rh2m = where(rh2m .lt.   0,   0, rh2m)

             end if

             sfroff(hh,:,:)  = (/ wf->SFROFF(0, lat_s:lat_n, lon_w:lon_e) /) - last_sfroff
             udroff(hh,:,:)  = (/ wf->UDROFF(0, lat_s:lat_n, lon_w:lon_e) /) - last_udroff

             swe(hh,:,:)   = (/ wf->SNOW(   0, lat_s:lat_n, lon_w:lon_e) /)
             snowh(hh,:,:) = (/ wf->SNOWH(  0, lat_s:lat_n, lon_w:lon_e) /)

             smois(hh,:,:,:) = (/ wf->SMCREL( 0, :, lat_s:lat_n, lon_w:lon_e) /)
             vswc( hh,:,:,:) = (/ wf->SMOIS(  0, :, lat_s:lat_n, lon_w:lon_e) /)
             stemp(hh,:,:,:) = (/ wf->TSLB(   0, :, lat_s:lat_n, lon_w:lon_e) /)

             if (hh .eq. (n_wrfhours-1)) then
                last_sfroff = (/ wf->SFROFF(0, lat_s:lat_n, lon_w:lon_e) /)
                last_udroff = (/ wf->UDROFF(0, lat_s:lat_n, lon_w:lon_e) /)
             end if

          end do

          ;
          ; Aggregate Daily Sums, Means, Maxs and Mins
          ;


          heat_index  = (/ heat_index_nws(t2m, rh2m, (/ 1, 2 /), False) /)



          max_temp_daily(calendar_day-1, :, :) = (/ dim_max_n(t2m, 0) /)
          min_temp_daily(calendar_day-1, :, :) = (/ dim_min_n(t2m, 0) /)
          avg_temp_daily(calendar_day-1, :, :) = (/ dim_avg_n(t2m, 0) /)

          sum_prec_daily(calendar_day-1, :, :) = (/ dim_sum_n(rain, 0) /)
          max_prec_daily(calendar_day-1, :, :) = (/ dim_max_n(rain, 0) /)

          avg_uwind_daily(calendar_day-1, :, :) = (/ dim_avg_n(u10m, 0) /)
          avg_vwind_daily(calendar_day-1, :, :) = (/ dim_avg_n(v10m, 0) /)
          avg_mwind_daily(calendar_day-1, :, :) = (/ dim_avg_n(m10m, 0) /)

          avg_vpd_daily(  calendar_day-1, :, :) =  (/ dim_avg_n(edef2m, 0) /)
          avg_sphum_daily(calendar_day-1, :, :) =  (/ dim_avg_n(   q2m, 0) /)


          sum_sfroff_daily(calendar_day-1,    :, :)  = (/ dim_sum_n(sfroff, 0) /)
          sum_udroff_daily(calendar_day-1,    :, :)  = (/ dim_sum_n(udroff, 0) /)


          avg_swe_daily(   calendar_day-1,    :, :)  = (/ dim_avg_n(swe, 0)    /)
          avg_smois_daily( calendar_day-1, :, :, :)  = (/ dim_avg_n(smois, 0)  /)
          avg_vswc_daily(  calendar_day-1, :, :, :)  = (/ dim_avg_n(vswc, 0)   /)
          avg_stemp_daily( calendar_day-1, :, :, :)  = (/ dim_avg_n(stemp, 0)  /)

          snowh_max(calendar_day-1,    :, :)  = (/ dim_avg_n(snowh, 0) /)


          avg_hi_daily(calendar_day-1, :, :) =  (/ dim_avg_n(   heat_index, 0) /)
          max_hi_daily(calendar_day-1, :, :) =  (/ dim_max_n(   heat_index, 0) /)
          min_hi_daily(calendar_day-1, :, :) =  (/ dim_min_n(   heat_index, 0) /)


          max_temp_daily_ann(doy-1, :, :) = max_temp_daily(calendar_day-1, :, :)
          min_temp_daily_ann(doy-1, :, :) = min_temp_daily(calendar_day-1, :, :)
          avg_temp_daily_ann(doy-1, :, :) = avg_temp_daily(calendar_day-1, :, :)
          sum_prec_daily_ann(doy-1, :, :) = sum_prec_daily(calendar_day-1, :, :)

          ;
          ; Destroy Daily 1-hrly and 3-hrly Inputs
          ;

          delete(rain)
          delete(t2m)
          delete(rain2)
          delete(u10m)
          delete(v10m)
          delete(m10m)

          delete(q2m)
          delete(psfc)
          delete(e2m)
          delete(es2m)
          delete(edef2m)


          delete(sfroff)
          delete(udroff)
          delete(swe)
          delete(smois)
          delete(vswc)
          delete(stemp)
          delete(snowh)

          delete(rh2m)


          ;
          ; Aggregate Monthly Sums, Means, Maxs, Mins, Indices and Day-counts
          ;

          if (calendar_day .eq. last_day)  then

             print("aggregating month "+ calendar_year + ":" + calendar_month + ":" + calendar_day)

             t2m_daily_mean_monthly_mean(mm,:,:)      = (/ dim_avg_n(avg_temp_daily, 0) /)
             t2m_daily_min_monthly_mean(mm,:,:)       = (/ dim_avg_n(min_temp_daily, 0) /)
             t2m_daily_max_monthly_mean(mm,:,:)       = (/ dim_avg_n(max_temp_daily, 0) /)
             t2m_daily_min_monthly_min(mm,:,:)        = (/ dim_min_n(min_temp_daily, 0) /)
             t2m_daily_min_monthly_max(mm,:,:)        = (/ dim_max_n(min_temp_daily, 0) /)
             t2m_daily_max_monthly_min(mm,:,:)        = (/ dim_min_n(max_temp_daily, 0) /)
             t2m_daily_max_monthly_max(mm,:,:)        = (/ dim_max_n(max_temp_daily, 0) /)

             q2m_daily_mean_monthly_mean(mm,:,:)         = (/ dim_avg_n(avg_sphum_daily, 0) /)
             e2m_deficit_daily_mean_monthly_mean(mm,:,:) = (/ dim_avg_n(avg_vpd_daily, 0) /)


             sum_prec_daily = where((sum_prec_daily .ge. 1),sum_prec_daily,0)

             precip_daily_total_monthly_total(mm,:,:) = (/ dim_sum_n(sum_prec_daily, 0) /)
             precip_daily_total_monthly_max(mm,:,:)   = (/ dim_max_n(sum_prec_daily, 0) /)
             precip_hourly_max_monthly_max(mm,:,:)    = (/ dim_max_n(max_prec_daily, 0) /)


             delete_me = dim_num_n((sum_prec_daily .ge. 1.0), 0)  ; number of days over 1mm
                delete_me@_FillValue = default_fillvalue("integer")
             printVarSummary(delete_me)

             bad_zero_gridcells = delete_me

             number_of_zero_gridcells = num((delete_me .eq. 0))
             print(number_of_zero_gridcells)

             bad_zero_gridcells = where(delete_me .eq. 0, delete_me@_FillValue, delete_me)
             printVarSummary(bad_zero_gridcells)

             simple_precipitation_intensity_index_monthly_value(mm,:,:) =  (/ (dim_sum_n(sum_prec_daily, 0)  /  bad_zero_gridcells )  /)


             delete(delete_me)
             delete(number_of_zero_gridcells)
             delete(bad_zero_gridcells)


             non_zero_sum_prec_daily                         = (/ where(sum_prec_daily .lt. 1.0, sum_prec_daily@_FillValue, sum_prec_daily) /)
             non_zero_sum_prec_daily_mean                    = (/ dim_avg_n(non_zero_sum_prec_daily, 0) /)
             precip_nonzero_daily_total_monthly_mean(mm,:,:) = (/ where(ismissing(non_zero_sum_prec_daily_mean), 0, non_zero_sum_prec_daily_mean) /)
             delete(non_zero_sum_prec_daily)
             delete(non_zero_sum_prec_daily_mean)


             wet_days_monthly                                = toint( sum_prec_daily )
             wet_days_monthly(:,:,:)                         = (/ where(sum_prec_daily .ge. 1.0, 1, 0) /) 
             dry_days_monthly                                = toint( sum_prec_daily )
             dry_days_monthly(:,:,:)                         = (/ where(sum_prec_daily .lt. 1.0, 1, 0) /) 

             wet_run_monthly    = dim_numrun_n(wet_days_monthly, 1, 0)
             dry_run_monthly    = dim_numrun_n(dry_days_monthly, 1, 0)


             length_of_dry_periods_daily_sum_monthly_mean(mm,:,:) = (/  avgrunlen(dry_run_monthly)   /)
             length_of_wet_periods_daily_sum_monthly_mean(mm,:,:) = (/  avgrunlen(wet_run_monthly)   /)

             length_of_dry_periods_daily_sum_monthly_max(mm,:,:) = (/  maxrunlen(dry_run_monthly)   /)
             length_of_wet_periods_daily_sum_monthly_max(mm,:,:) =  (/  maxrunlen(wet_run_monthly)   /)

             delete(wet_days_monthly)
             delete(dry_days_monthly)
             delete(wet_run_monthly)
             delete(dry_run_monthly)

             t2m_daily_mean_daily_range(mm,:,:)       = (/ dim_avg_n( (/  max_temp_daily - min_temp_daily /), 0) /)

             u10m_wind_speed_daily_mean_monthly_mean(mm,:,:) = (/ dim_avg_n(avg_uwind_daily, 0) /)
             v10m_wind_speed_daily_mean_monthly_mean(mm,:,:) = (/ dim_avg_n(avg_vwind_daily, 0) /)
             m10m_wind_speed_daily_mean_monthly_mean(mm,:,:) = (/ dim_avg_n(avg_mwind_daily, 0) /)

             sfroff_daily_total_monthly_total(mm,:,:) = (/ dim_sum_n(sum_sfroff_daily, 0) *mask2d /)
             ugroff_daily_total_monthly_total(mm,:,:) = (/ dim_sum_n(sum_udroff_daily, 0)* mask2d /)

             swe_daily_mean_monthly_mean(mm,:,:)                    = (/ dim_avg_n(avg_swe_daily,   0)* mask2d /)

             soil_mois_daily_mean_monthly_mean(mm,:,:,:)            = (/ dim_avg_n(avg_smois_daily, 0)* mask3d /)
             vol_soil_wat_content_daily_mean_monthly_mean(mm,:,:,:) = (/ dim_avg_n(avg_vswc_daily,  0)* mask3d /)
             soil_temp_daily_mean_monthly_mean(mm,:,:,:)            = (/ dim_avg_n(avg_stemp_daily, 0)* mask3d /)

             number_of_days_with_snow_cover_monthly_total(mm, :, :) = (/ dim_num_n((snowh_max .gt. 0.01), 0) * imask2d /)



             cdd_deleteme = (/ round(avg_temp_daily - 291.48,3) /)
             cdd_deleteme = where(cdd_deleteme .lt. 0, 0, cdd_deleteme)

             hdd_deleteme = (/ round(291.48 - avg_temp_daily,3) /)
             hdd_deleteme = where(hdd_deleteme .lt. 0, 0, hdd_deleteme)

             gdd_deleteme = (/ round(avg_temp_daily - 283.15,3)  /)
             gdd_deleteme = where(gdd_deleteme .lt. 0, 0, gdd_deleteme)

             mdd_deleteme = (/ round(avg_temp_daily - 283.15,3)/)
             mdd_deleteme = where((avg_temp_daily .gt. 303.15), 0, mdd_deleteme)
             mdd_deleteme = where((mdd_deleteme       .lt.      0), 0, mdd_deleteme)


             heating_degree_days_monthly_total(mm,:,:)       = (/ dim_sum_n(hdd_deleteme, 0)  /)
             cooling_degree_days_monthly_total(mm,:,:)       = (/ dim_sum_n(cdd_deleteme, 0)  /)

             growing_degree_days_monthly_total(mm,:,:)       = (/ dim_sum_n(gdd_deleteme, 0) * imask2d /)
             maize_growing_degree_days_monthly_total(mm,:,:) = (/ dim_sum_n(mdd_deleteme, 0) * imask2d /)
 


             delete(hdd_deleteme)
             delete(cdd_deleteme)
             delete(gdd_deleteme)
             delete(mdd_deleteme)


             number_of_frost_days_monthly_total(mm,:,:) = (/ dim_num_n((min_temp_daily .lt. 273.15), 0) /)
             number_of_icing_days_monthly_total(mm,:,:) = (/ dim_num_n((max_temp_daily .lt. 273.15), 0) /)

             number_of_summer_days_monthly_total(mm,:,:)     = (/ dim_num_n((max_temp_daily .gt. 298.15), 0) /)
             number_of_tropical_nights_monthly_total(mm,:,:) = (/ dim_num_n((min_temp_daily .gt. 293.15), 0) /)

             number_of_days_precip_ge_010mm_monthly_total(mm,:,:) = (/ dim_num_n((sum_prec_daily .ge. 10.), 0)  /)
             number_of_days_precip_ge_020mm_monthly_total(mm,:,:) = (/ dim_num_n((sum_prec_daily .ge. 20.), 0)  /)
             number_of_days_precip_ge_050mm_monthly_total(mm,:,:) = (/ dim_num_n((sum_prec_daily .ge. 50.), 0)  /)
             number_of_days_precip_ge_075mm_monthly_total(mm,:,:) = (/ dim_num_n((sum_prec_daily .ge. 75.), 0)  /)
             number_of_days_precip_ge_100mm_monthly_total(mm,:,:) = (/ dim_num_n((sum_prec_daily .ge. 100.), 0)  /)





             nws_heat_index_daily_mean_monthly_mean(mm,:,:) = (/ dim_avg_n(avg_hi_daily, 0) /)
             nws_heat_index_daily_mean_monthly_max(mm,:,:)  = (/ dim_max_n(avg_hi_daily, 0) /)
             nws_heat_index_daily_max_monthly_max(mm,:,:)   = (/ dim_max_n(max_hi_daily, 0) /)

             number_of_days_max_heat_index_080_to_090_monthly_total(mm,:,:) = (/ dim_num_n((max_hi_daily .ge.  80.) .and. (max_hi_daily .lt.  90.), 0)  /)
             number_of_days_max_heat_index_090_to_105_monthly_total(mm,:,:) = (/ dim_num_n((max_hi_daily .ge.  90.) .and. (max_hi_daily .lt. 105.), 0)  /)
             number_of_days_max_heat_index_105_to_130_monthly_total(mm,:,:) = (/ dim_num_n((max_hi_daily .ge. 105.) .and. (max_hi_daily .lt. 130.), 0)  /)
             number_of_days_max_heat_index_over_130_monthly_total(mm,:,:)   = (/ dim_num_n((max_hi_daily .ge. 130.), 0)  /)


             ;
             ; Destroy Monthly Temporary Aggregates
             ;

             delete(max_temp_daily)
             delete(min_temp_daily)
             delete(avg_temp_daily)
             delete(sum_prec_daily)
             delete(max_prec_daily)
             delete(avg_uwind_daily)
             delete(avg_vwind_daily)
             delete(avg_mwind_daily)
             delete(sum_sfroff_daily)
             delete(sum_udroff_daily)
             delete(avg_swe_daily)
             delete(avg_smois_daily)
             delete(avg_vswc_daily)
             delete(avg_stemp_daily)
             delete(snowh_max)
             delete(avg_vpd_daily)
             delete(avg_sphum_daily)
             delete(avg_hi_daily)
             delete(min_hi_daily)
             delete(max_hi_daily)



             mm = mm + 1

          end if

          ;
          ; Aggregate Annual Sums, Means, Maxs, Mins, Indices and Day-counts
          ;

          if ((calendar_day .eq. last_day).and.(calendar_month .eq. 12))  then

             print("aggregating year "+ calendar_year + ":" + calendar_month + ":" + calendar_day)


             sum_prec_daily_ann = where((sum_prec_daily_ann .ge. 1),sum_prec_daily_ann,0)


             t2m_daily_mean_annual_mean(yy,:,:)         = (/ dim_avg_n(avg_temp_daily_ann, 0) /)
             precip_daily_total_annual_total(yy,:,:) = (/ dim_sum_n(sum_prec_daily_ann, 0) /)




             cdd_deleteme = (/ round(avg_temp_daily_ann - 291.48,3) /)
             cdd_deleteme = where(cdd_deleteme .lt. 0, 0, cdd_deleteme)

             hdd_deleteme = (/ round(291.48 - avg_temp_daily_ann,3) /)
             hdd_deleteme = where(hdd_deleteme .lt. 0, 0, hdd_deleteme)

             gdd_deleteme = (/ round(avg_temp_daily_ann - 283.15,3)  /)
             gdd_deleteme = where(gdd_deleteme .lt. 0, 0, gdd_deleteme)

             mdd_deleteme = (/ round(avg_temp_daily_ann - 283.15,3)/)
             mdd_deleteme = where((avg_temp_daily_ann .gt. 303.15), 0, mdd_deleteme)
             mdd_deleteme = where((mdd_deleteme       .lt.      0), 0, mdd_deleteme)




             heating_degree_days_annual_total(yy,:,:)   = (/ dim_sum_n(hdd_deleteme, 0) /)
             cooling_degree_days_annual_total(yy,:,:)   = (/ dim_sum_n(cdd_deleteme, 0) /)

             growing_degree_days_annual_total(yy,:,:)       = (/ dim_sum_n(gdd_deleteme, 0) * imask2d /)
             maize_growing_degree_days_annual_total(yy,:,:) = (/ dim_sum_n(mdd_deleteme, 0) * imask2d /)






             delete(hdd_deleteme)
             delete(cdd_deleteme)
             delete(gdd_deleteme)
             delete(mdd_deleteme)




             avg_gs_run_ann        = toint(avg_temp_daily_ann)
             avg_gs_run_ann(:,:,:) = 0
             avg_gs_run_ann        = where((avg_temp_daily_ann .gt. 278.15), 1, 0)
             gseason_run_ann       = dim_numrun_n(avg_gs_run_ann, 1, 0)



             wet_days_annual                                = toint( sum_prec_daily_ann )
             wet_days_annual(:,:,:)                         = (/ where(sum_prec_daily_ann .ge. 1.0, 1, 0) /) 
             dry_days_annual                                = toint( sum_prec_daily_ann )
             dry_days_annual(:,:,:)                         = (/ where(sum_prec_daily_ann .lt. 1.0, 1, 0) /) 

             wet_run_annual    = dim_numrun_n(wet_days_annual, 1, 0)
             dry_run_annual    = dim_numrun_n(dry_days_annual, 1, 0)


             length_of_dry_periods_daily_sum_annual_mean(yy,:,:) = (/  avgrunlen(dry_run_annual)   /)
             length_of_wet_periods_daily_sum_annual_mean(yy,:,:) = (/  avgrunlen(wet_run_annual)   /)

             length_of_dry_periods_daily_sum_annual_max(yy,:,:) = (/  maxrunlen(dry_run_annual)   /)
             length_of_wet_periods_daily_sum_annual_max(yy,:,:) = (/  maxrunlen(wet_run_annual)   /)

             delete(wet_days_annual)
             delete(dry_days_annual)
             delete(wet_run_annual)
             delete(dry_run_annual)


             length_of_growing_season_annual_total(yy,:,:) = (/  maxrunlen(gseason_run_ann) * imask2d  /)

             number_of_frost_days_annual_total(yy,:,:) = (/ dim_num_n((min_temp_daily_ann .lt. 273.15), 0) /)
             number_of_icing_days_annual_total(yy,:,:) = (/ dim_num_n((max_temp_daily_ann .lt. 273.15), 0) /)

             number_of_summer_days_annual_total(yy,:,:)     = (/ dim_num_n((max_temp_daily_ann .gt. 298.15), 0) /)
             number_of_tropical_nights_annual_total(yy,:,:) = (/ dim_num_n((min_temp_daily_ann .gt. 293.15), 0) /)

             number_of_days_precip_ge_010mm_annual_total(yy,:,:) = (/ dim_num_n((sum_prec_daily_ann .ge.  10.), 0)  /)
             number_of_days_precip_ge_020mm_annual_total(yy,:,:) = (/ dim_num_n((sum_prec_daily_ann .ge.  20.), 0)  /)
             number_of_days_precip_ge_050mm_annual_total(yy,:,:) = (/ dim_num_n((sum_prec_daily_ann .ge.  50.), 0)  /)
             number_of_days_precip_ge_075mm_annual_total(yy,:,:) = (/ dim_num_n((sum_prec_daily_ann .ge.  75.), 0)  /)
             number_of_days_precip_ge_100mm_annual_total(yy,:,:) = (/ dim_num_n((sum_prec_daily_ann .ge. 100.), 0)  /)

             ;
             ; Destroy Monthly Temporary Aggregates & Masks
             ;

             delete(max_temp_daily_ann)
             delete(min_temp_daily_ann)
             delete(avg_temp_daily_ann)
             delete(sum_prec_daily_ann)
             delete(avg_gs_run_ann)
             delete(gseason_run_ann)
             




             yy = yy + 1

          end if



      end do

      ;
      ; File Output
      ;

      length_of_growing_season_annual_total = where(length_of_growing_season_annual_total .le. 12, 0, length_of_growing_season_annual_total)

      length_of_growing_season_annual_total = where(length_of_growing_season_annual_total .lt. 0,     \
                                                    length_of_growing_season_annual_total@_FillValue, \
                                                    length_of_growing_season_annual_total)

      simple_precipitation_intensity_index_monthly_value = where(ismissing(simple_precipitation_intensity_index_monthly_value),0,simple_precipitation_intensity_index_monthly_value)



      soil_mois_daily_mean_monthly_mean            = where(soil_mois_daily_mean_monthly_mean .lt. 0, soil_mois_daily_mean_monthly_mean@_FillValue, soil_mois_daily_mean_monthly_mean)
      vol_soil_wat_content_daily_mean_monthly_mean = where(vol_soil_wat_content_daily_mean_monthly_mean .lt. 0, vol_soil_wat_content_daily_mean_monthly_mean@_FillValue, vol_soil_wat_content_daily_mean_monthly_mean)
      soil_temp_daily_mean_monthly_mean            = where(soil_temp_daily_mean_monthly_mean .lt. 0, soil_temp_daily_mean_monthly_mean@_FillValue, soil_temp_daily_mean_monthly_mean)



      length_of_dry_periods_daily_sum_annual_mean = where(ismissing(length_of_dry_periods_daily_sum_annual_mean),0,length_of_dry_periods_daily_sum_annual_mean)
      length_of_wet_periods_daily_sum_annual_mean = where(ismissing(length_of_wet_periods_daily_sum_annual_mean),0,length_of_wet_periods_daily_sum_annual_mean)

      length_of_dry_periods_daily_sum_annual_max = where(ismissing(length_of_dry_periods_daily_sum_annual_max),0,length_of_dry_periods_daily_sum_annual_max)
      length_of_wet_periods_daily_sum_annual_max = where(ismissing(length_of_wet_periods_daily_sum_annual_max),0,length_of_wet_periods_daily_sum_annual_max)



      length_of_dry_periods_daily_sum_monthly_mean = where(ismissing(length_of_dry_periods_daily_sum_monthly_mean),0,length_of_dry_periods_daily_sum_monthly_mean)
      length_of_wet_periods_daily_sum_monthly_mean = where(ismissing(length_of_wet_periods_daily_sum_monthly_mean),0,length_of_wet_periods_daily_sum_monthly_mean)

      length_of_dry_periods_daily_sum_monthly_max = where(ismissing(length_of_dry_periods_daily_sum_monthly_max),0,length_of_dry_periods_daily_sum_monthly_max)
      length_of_wet_periods_daily_sum_monthly_max = where(ismissing(length_of_wet_periods_daily_sum_monthly_max),0,length_of_wet_periods_daily_sum_monthly_max)



   lat2d_for_spei = conform_dims((/nm, ny, nx /), lat, (/ 1 /))
      lat2d_for_spei!0             = "month"
      lat2d_for_spei!1             = "lat"
      lat2d_for_spei!2             = "lon"
      lat2d_for_spei@description   = "latitude"
      lat2d_for_spei@long_name      = "latitude"
      lat2d_for_spei@standard_name = "latitude"
      lat2d_for_spei@units         = "degrees_north"


      if(1 .gt. 2) then

      standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value(:,:,:) = (/ speidx(precip_daily_total_monthly_total, \
                                                                                                        t2m_daily_mean_monthly_mean,      \
                                                                                                        lat2d_for_spei,                   \
                                                                                                        3,                                \
                                                                                                        False,                            \
                                                                                                        0 ) /)

      standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value(:,:,:) = (/ speidx(precip_daily_total_monthly_total, \
                                                                                                        t2m_daily_mean_monthly_mean,      \
                                                                                                        lat2d_for_spei,                   \
                                                                                                        6,                                \
                                                                                                        False,                            \
                                                                                                        0 ) /)

      standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value(:,:,:) = (/ speidx(precip_daily_total_monthly_total, \
                                                                                                        t2m_daily_mean_monthly_mean,      \
                                                                                                        lat2d_for_spei,                   \
                                                                                                        12,                                \
                                                                                                        False,                             \
                                                                                                        0 ) /)

      standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value(:,:,:) = (/ speidx(precip_daily_total_monthly_total, \
                                                                                                        t2m_daily_mean_monthly_mean,      \
                                                                                                        lat2d_for_spei,                   \
                                                                                                        24,                                \
                                                                                                        False,                             \
                                                                                                        0 ) /)

      standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value(:,:,:) = (/ speidx(precip_daily_total_monthly_total, \
                                                                                                        t2m_daily_mean_monthly_mean,      \
                                                                                                        lat2d_for_spei,                   \
                                                                                                        36,                                \
                                                                                                        False,                            \
                                                                                                        0 ) /)


end if

      delete(lat2d_for_spei)
      output_file = "./summary_" + scenario(s) + ".nc"

      system("rm -fv " + output_file)


      of = addfile(output_file, "c")

      print("Outputting Data " + output_file)
      print(of)

         of@title                 = "Climate Extreme Indicies for WRF Ensemble Run: "+scenario(s)
         of@institution           = "South Dakota School of Mines (aggregated dataset), National Center for Atmospheric Research (original WRF RCM runs)"
         of@source                = "Weather Reserach and Forecasting Model v3.5"
         of@radiation_scheme      = radiation(s)
         of@cumulus_scheme        = convection(s)
         of@microphysics_scheme   = microphysics(s)
         of@boundary_layer_scheme = boundary_layer(s)
         of@start_year            = year_starts(s)
         of@end_year              = year_ends(s)
         of@scenario_code         = scenario(s)
         of@missing_aux_files     = missing_aux_days

         of->month          = month
         of->year           = year
         of->lat            = lat
         of->lon            = lon
         of->soil_depth     = soil_depth
         of->soil_thickness = soil_thickness

         of->calendar_year_in_monthly_data = calendar_year_in_monthly_data
         of->calendar_month_in_monthly_data = calendar_month_in_monthly_data

         of->terrain_height                   = terrain_height
         of->terrain_height_variance          = terrain_height_variance
         of->soil_texture                     = soil_texture
         of->land_classification              = land_classification
         of->land_mask                        = land_mask

         of->t2m_daily_mean_monthly_mean      = t2m_daily_mean_monthly_mean
         of->t2m_daily_min_monthly_mean       = t2m_daily_min_monthly_mean
         of->t2m_daily_max_monthly_mean       = t2m_daily_max_monthly_mean
         of->t2m_daily_min_monthly_min        = t2m_daily_min_monthly_min
         of->t2m_daily_min_monthly_max        = t2m_daily_min_monthly_max
         of->t2m_daily_max_monthly_min        = t2m_daily_max_monthly_min
         of->t2m_daily_max_monthly_max        = t2m_daily_max_monthly_max
         of->t2m_daily_mean_daily_range       = t2m_daily_mean_daily_range


         of->precip_daily_total_monthly_total        = precip_daily_total_monthly_total
         of->precip_daily_total_monthly_max          = precip_daily_total_monthly_max
         of->precip_hourly_max_monthly_max           = precip_hourly_max_monthly_max
         of->precip_nonzero_daily_total_monthly_mean = precip_nonzero_daily_total_monthly_mean

         of->length_of_dry_periods_daily_sum_monthly_mean = length_of_dry_periods_daily_sum_monthly_mean
         of->length_of_wet_periods_daily_sum_monthly_mean = length_of_wet_periods_daily_sum_monthly_mean

         of->length_of_dry_periods_daily_sum_monthly_max = length_of_dry_periods_daily_sum_monthly_max
         of->length_of_wet_periods_daily_sum_monthly_max = length_of_wet_periods_daily_sum_monthly_max


       ;  of->standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value    = standardized_precipitation_evapotranspiration_index_03monthrange_monthly_value
       ;  of->standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value    = standardized_precipitation_evapotranspiration_index_06monthrange_monthly_value
       ;  of->standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value    = standardized_precipitation_evapotranspiration_index_12monthrange_monthly_value
       ;  of->standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value    = standardized_precipitation_evapotranspiration_index_24monthrange_monthly_value
       ;  of->standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value    = standardized_precipitation_evapotranspiration_index_36monthrange_monthly_value

   

         of->q2m_daily_mean_monthly_mean          = q2m_daily_mean_monthly_mean
         of->e2m_deficit_daily_mean_monthly_mean = e2m_deficit_daily_mean_monthly_mean

         of->u10m_wind_speed_daily_mean_monthly_mean = u10m_wind_speed_daily_mean_monthly_mean
         of->v10m_wind_speed_daily_mean_monthly_mean = v10m_wind_speed_daily_mean_monthly_mean
         of->m10m_wind_speed_daily_mean_monthly_mean = m10m_wind_speed_daily_mean_monthly_mean

         of->heating_degree_days_monthly_total       = heating_degree_days_monthly_total
         of->cooling_degree_days_monthly_total       = cooling_degree_days_monthly_total
         of->growing_degree_days_monthly_total       = growing_degree_days_monthly_total
         of->maize_growing_degree_days_monthly_total = maize_growing_degree_days_monthly_total

         of->number_of_frost_days_monthly_total   = number_of_frost_days_monthly_total
         of->number_of_icing_days_monthly_total   = number_of_icing_days_monthly_total

         of->number_of_summer_days_monthly_total     = number_of_summer_days_monthly_total
         of->number_of_tropical_nights_monthly_total = number_of_tropical_nights_monthly_total


         of->sfroff_daily_total_monthly_total = sfroff_daily_total_monthly_total
         of->ugroff_daily_total_monthly_total = ugroff_daily_total_monthly_total


         of->swe_daily_mean_monthly_mean                  = swe_daily_mean_monthly_mean
         of->soil_mois_daily_mean_monthly_mean            = soil_mois_daily_mean_monthly_mean
         of->vol_soil_wat_content_daily_mean_monthly_mean = vol_soil_wat_content_daily_mean_monthly_mean
         of->soil_temp_daily_mean_monthly_mean            = soil_temp_daily_mean_monthly_mean

         of->number_of_days_precip_ge_010mm_monthly_total = number_of_days_precip_ge_010mm_monthly_total
         of->number_of_days_precip_ge_020mm_monthly_total = number_of_days_precip_ge_020mm_monthly_total
         of->number_of_days_precip_ge_050mm_monthly_total = number_of_days_precip_ge_050mm_monthly_total
         of->number_of_days_precip_ge_075mm_monthly_total = number_of_days_precip_ge_075mm_monthly_total
         of->number_of_days_precip_ge_100mm_monthly_total = number_of_days_precip_ge_100mm_monthly_total


         of->simple_precipitation_intensity_index_monthly_value = simple_precipitation_intensity_index_monthly_value

         of->nws_heat_index_daily_mean_monthly_mean = nws_heat_index_daily_mean_monthly_mean
         of->nws_heat_index_daily_mean_monthly_max  = nws_heat_index_daily_mean_monthly_max
         of->nws_heat_index_daily_max_monthly_max   = nws_heat_index_daily_max_monthly_max

         of->number_of_days_max_heat_index_080_to_090_monthly_total = number_of_days_max_heat_index_080_to_090_monthly_total
         of->number_of_days_max_heat_index_090_to_105_monthly_total = number_of_days_max_heat_index_090_to_105_monthly_total
         of->number_of_days_max_heat_index_105_to_130_monthly_total = number_of_days_max_heat_index_105_to_130_monthly_total
         of->number_of_days_max_heat_index_over_130_monthly_total   = number_of_days_max_heat_index_over_130_monthly_total



         of->number_of_days_with_snow_cover_monthly_total = number_of_days_with_snow_cover_monthly_total

         of->t2m_daily_mean_annual_mean         =  t2m_daily_mean_annual_mean
         of->precip_daily_total_annual_total    =  precip_daily_total_annual_total

         of->heating_degree_days_annual_total         = heating_degree_days_annual_total
         of->cooling_degree_days_annual_total         = cooling_degree_days_annual_total
         of->growing_degree_days_annual_total         = growing_degree_days_annual_total
         of->maize_growing_degree_days_annual_total   = maize_growing_degree_days_annual_total
         of->length_of_growing_season_annual_total    = length_of_growing_season_annual_total

         of->number_of_frost_days_annual_total   = number_of_frost_days_annual_total
         of->number_of_icing_days_annual_total   = number_of_icing_days_annual_total

         of->number_of_summer_days_annual_total     = number_of_summer_days_annual_total
         of->number_of_tropical_nights_annual_total = number_of_tropical_nights_annual_total

         of->number_of_days_precip_ge_010mm_annual_total = number_of_days_precip_ge_010mm_annual_total
         of->number_of_days_precip_ge_020mm_annual_total = number_of_days_precip_ge_020mm_annual_total
         of->number_of_days_precip_ge_050mm_annual_total = number_of_days_precip_ge_050mm_annual_total
         of->number_of_days_precip_ge_075mm_annual_total = number_of_days_precip_ge_075mm_annual_total
         of->number_of_days_precip_ge_100mm_annual_total = number_of_days_precip_ge_100mm_annual_total


         of->length_of_dry_periods_daily_sum_annual_mean = length_of_dry_periods_daily_sum_annual_mean
         of->length_of_wet_periods_daily_sum_annual_mean = length_of_wet_periods_daily_sum_annual_mean

         of->length_of_dry_periods_daily_sum_annual_max = length_of_dry_periods_daily_sum_annual_max
         of->length_of_wet_periods_daily_sum_annual_max = length_of_wet_periods_daily_sum_annual_max




         print(of)

      delete(of)

      print("Done with Year " + calendar_year)

   end do



end
