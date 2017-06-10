
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/glade/u/home/capehart/heat_stress.ncl"

begin

    ;setfileoption("nc", "FileStructure", "Advanced")
    ;setfileoption("nc", "Format",  "NetCDF4")


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;
   ;  User Modification Area
   ;
   ;   User selects Repositories, Years and File Locations
   ;

   ;
   ; Repository, and Date Runs
   ;


   scenario = (/  "rktm_1990-2000", \
                  "rktm_2020-2030", \
                  "rktm_2030-2040", \
                  "rktm_2050-2060", \
                  "rktm_2080-2090" /)


   ;
   ; Physics Options
   ;

   radiation      = (/ "RRTM", \
                       "RRTM", \
                       "RRTM", \
                       "RRTM", \
                       "RRTM" /)

   convection     = (/ "Kain-Fritsch", \
                       "Kain-Fritsch", \
                       "Kain-Fritsch", \
                       "Kain-Fritsch", \
                       "Kain-Fritsch" /)

   microphysics   = (/ "Thompson", \
                       "Thompson", \
                       "Thompson", \
                       "Thompson", \
                       "Thompson" /)

   boundary_layer = (/ "MYJ", \
                       "MYJ", \
                       "MYJ", \
                       "MYJ", \
                       "MYJ" /)
   ;
   ; Start and End Years
   ;

   year_starts = (/ 1990, 2020, 2030, 2050, 2080 /)
   year_ends   = (/ 2000, 2030, 2040, 2060, 2090 /)

   nyr = 11   ; # of years to run... currently each has 11 yr

   nt = 11*365
   n_scenarios = dimsizes(scenario)

   ;
   ; Clipbox Latitude and Longitudes CONUS
   ;


   far_south =   23.00 ; degrees_north
   far_north =   55.00 ; degrees_north

   far_west  = -129.00 ; degrees_east
   far_east  =  -52.00 ; degrees_east


   ;
   ; Clipbox Latitude and Longitudes UNR/FSD/ABR CWAa
   ;


   far_south =   41.00 ; degrees_north
   far_north =   46.90 ; degrees_north

   far_west  = -106.3 ; degrees_east
   far_east  =  -94.4 ; degrees_east



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

   t2m_daily_mean = new((/ nt, ny, nx /), float)
      t2m_daily_mean!0             = "time"
      t2m_daily_mean!1             = "lat"
      t2m_daily_mean!2             = "lon"
      t2m_daily_mean@description   = "Avg Daily 2-m Temp"
      t2m_daily_mean@long_name     = "Avg Daily 2-m Temp"
      t2m_daily_mean@units         = "K"
      t2m_daily_mean@standard_name = "air_temperature"
      t2m_daily_mean@cell_methods  = "hour: mean"

   t2m_daily_max = new((/ nt, ny, nx /), float)
      t2m_daily_max!0             = "time"
      t2m_daily_max!1             = "lat"
      t2m_daily_max!2             = "lon"
      t2m_daily_max@description   = "Max Daily 2-m Temp"
      t2m_daily_max@long_name     = "Max Daily 2-m Temp"
      t2m_daily_max@units         = "K"
      t2m_daily_max@standard_name = "air_temperature"
      t2m_daily_max@cell_methods  = "hour: mean"


   t2m_daily_min = new((/ nt, ny, nx /), float)
      t2m_daily_min!0             = "time"
      t2m_daily_min!1             = "lat"
      t2m_daily_min!2             = "lon"
      t2m_daily_min@description   = "Min Daily 2-m Temp"
      t2m_daily_min@long_name     = "Min Daily 2-m Temp"
      t2m_daily_min@units         = "K"
      t2m_daily_min@standard_name = "air_temperature"
      t2m_daily_min@cell_methods  = "hour: mean"


   td2m_daily_mean = new((/ nt, ny, nx /), float)
      td2m_daily_mean!0             = "time"
      td2m_daily_mean!1             = "lat"
      td2m_daily_mean!2             = "lon"
      td2m_daily_mean@description   = "Avg Daily 2-m Dew Point Temp"
      td2m_daily_mean@long_name     = "Avg Daily 2-m Dew Point Temp"
      td2m_daily_mean@units         = "K"
      td2m_daily_mean@standard_name = "dew_point_temperature"
      td2m_daily_mean@cell_methods  = "hour: mean"

   td2m_daily_max = new((/ nt, ny, nx /), float)
      td2m_daily_max!0             = "time"
      td2m_daily_max!1             = "lat"
      td2m_daily_max!2             = "lon"
      td2m_daily_max@description   = "Max Daily 2-m Dew Point Temp"
      td2m_daily_max@long_name     = "Max Daily 2-m Dew Point Temp"
      td2m_daily_max@units         = "K"
      td2m_daily_max@standard_name = "dew_point_temperature"
      td2m_daily_max@cell_methods  = "hour: mean"


   td2m_daily_min = new((/ nt, ny, nx /), float)
      td2m_daily_min!0             = "time"
      td2m_daily_min!1             = "lat"
      td2m_daily_min!2             = "lon"
      td2m_daily_min@description   = "Min Daily 2-m Dew Point Temp"
      td2m_daily_min@long_name     = "Min Daily 2-m Dew Point Temp"
      td2m_daily_min@units         = "K"
      td2m_daily_min@standard_name = "dew_point_temperature"
      td2m_daily_min@cell_methods  = "hour: mean"



   precip_total_daily = new((/ nt, ny, nx /), float)
      precip_total_daily!0             = "time"
      precip_total_daily!1             = "lat"
      precip_total_daily!2             = "lon"
      precip_total_daily@description   = "Daily Total Precip"
      precip_total_daily@long_name     = "Daily Total Precip"
      precip_total_daily@units         = "kg m-2"
      precip_total_daily@cell_methods  = "hour: sum"
      precip_total_daily@standard_name = "precipitation_amount"

   precip_hourly_max = new((/ nt, ny, nx /), float)
      precip_hourly_max!0             = "time"
      precip_hourly_max!1             = "lat"
      precip_hourly_max!2             = "lon"
      precip_hourly_max@description   = "Daily Max Hrly Precip"
      precip_hourly_max@long_name     = "Daily Max Hrly Precip"
      precip_hourly_max@units         = "kg m-2"
      precip_hourly_max@cell_methods  = "hour: max"
      precip_hourly_max@standard_name = "precipitation_amount"

   m10m_wind_speed_daily_mean = new((/ nt, ny, nx /), float)
      m10m_wind_speed_daily_mean!0             = "time"
      m10m_wind_speed_daily_mean!1             = "lat"
      m10m_wind_speed_daily_mean!2             = "lon"
      m10m_wind_speed_daily_mean@description   = "Mean Daily 10-m Wind Speed"
      m10m_wind_speed_daily_mean@long_name     = "Mean Daily 10-m Wind Speed"
      m10m_wind_speed_daily_mean@units         = "m s-1"
      m10m_wind_speed_daily_mean@standard_name = "wind_speed"
      m10m_wind_speed_daily_mean@cell_methods  = "hour: mean"


   m10m_wind_speed_daily_max = new((/ nt, ny, nx /), float)
      m10m_wind_speed_daily_max!0             = "time"
      m10m_wind_speed_daily_max!1             = "lat"
      m10m_wind_speed_daily_max!2             = "lon"
      m10m_wind_speed_daily_max@description   = "Max Daily 10-m Wind Speed"
      m10m_wind_speed_daily_max@long_name     = "Max Daily 10-m Wind Speed"
      m10m_wind_speed_daily_max@units         = "m s-1"
      m10m_wind_speed_daily_max@standard_name = "wind_speed"
      m10m_wind_speed_daily_max@cell_methods  = "hour: max"


   u10m_wind_speed_daily_mean = new((/ nt, ny, nx /), float)
      u10m_wind_speed_daily_mean!0             = "time"
      u10m_wind_speed_daily_mean!1             = "lat"
      u10m_wind_speed_daily_mean!2             = "lon"
      u10m_wind_speed_daily_mean@description   = "Mean Daily 10-m Eastward Wind Speed"
      u10m_wind_speed_daily_mean@long_name     = "Mean Daily 10-m Eastward Wind Speed"
      u10m_wind_speed_daily_mean@units         = "m s-1"
      u10m_wind_speed_daily_mean@standard_name = "eastward_wind"
      u10m_wind_speed_daily_mean@cell_methods  = "hour: mean"


   v10m_wind_speed_daily_mean = new((/ nt, ny, nx /), float)
      v10m_wind_speed_daily_mean!0             = "time"
      v10m_wind_speed_daily_mean!1             = "lat"
      v10m_wind_speed_daily_mean!2             = "lon"
      v10m_wind_speed_daily_mean@description   = "Mean Daily 10-m Northward Wind Speed"
      v10m_wind_speed_daily_mean@long_name     = "Mean Daily 10-m Northward Wind Speed"
      v10m_wind_speed_daily_mean@units         = "m s-1"
      v10m_wind_speed_daily_mean@standard_name = "northward_wind"
      v10m_wind_speed_daily_mean@cell_methods  = "hour: mean"





   thi_daily_mean = new((/ nt, ny, nx /), float)
      thi_daily_mean!0             = "time"
      thi_daily_mean!1             = "lat"
      thi_daily_mean!2             = "lon"
      thi_daily_mean@description   = "Avg Daily Temperature Humidity Index"
      thi_daily_mean@long_name     = "Avg Daily Temperature Humidity Index"
      thi_daily_mean@cell_methods  = "hour: mean"

   thi_daily_max = new((/ nt, ny, nx /), float)
      thi_daily_max!0             = "time"
      thi_daily_max!1             = "lat"
      thi_daily_max!2             = "lon"
      thi_daily_max@description   = "Max Daily Temperature Humidity Index"
      thi_daily_max@long_name     = "Max Daily Temperature Humidity Index"
      thi_daily_max@cell_methods  = "hour: mean"


   thi_daily_min = new((/ nt, ny, nx /), float)
      thi_daily_min!0             = "time"
      thi_daily_min!1             = "lat"
      thi_daily_min!2             = "lon"
      thi_daily_min@description   = "Min Daily Temperature Humidity Index"
      thi_daily_min@long_name     = "Min Daily Temperature Humidity Index"
      thi_daily_min@cell_methods  = "hour: mean"



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


   ldom_calendar = toint((/ -1,  31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /))


   do s = 0,  n_scenarios-1

      ;
      ; Counters and Buckets Reset for new simulation
      ;

      mm  = toint(0)
      yy  = toint(0)

      last_rain(:,:,:)    = 0






      ;
      ; Time Control
      ;

      missing_aux_days = ""

      year_start = year_starts(s)
      year_end   = year_ends(s)

      start_time = cd_inv_calendar(year_start, 01, 01, 00, 00, 00, time_units, calendar)
      end_time   = cd_inv_calendar(  year_end, 12, 31, 00, 00, 00, time_units, calendar)

      time = ispan(toint(start_time),toint(end_time),1)
         time@long_name   = "time"
         time@description = "time"
         time@units       = time_units
         time@calendar    = calendar@calendar
         time!0           = "time"

      nt   = dimsizes(time)

      day = time

      aux_file_status = ispan(toint(start_time),toint(end_time),1)
         aux_file_status(:)          = 1
         aux_file_status!0           = "time"
         aux_file_status@long_name   = "Aux File Presence String"
         aux_file_status@description = "Aux File Presence String"
         aux_file_status@comment1    = "1 = aux file present; 0 = aux file missing"
         aux_file_status@comment2    = "When Aux File is Missing we distribute the 6-hrly data over the 24-hr daily period"



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


      calendar_year_in_monthly_data = new(nt, integer)
         calendar_year_in_monthly_data@description = "Calendar Year In Monthly Records"
         calendar_year_in_monthly_data@long_name   = "Calendar Year In Monthly Records"
         calendar_year_in_monthly_data!0           = "day"

 
      calendar_month_in_monthly_data = new(nt, integer)
         calendar_year_in_monthly_data@description    = "Calendar Month In Monthly Records"
         calendar_month_in_monthly_data@long_name     = "Calendar Month In Monthly Records"
         calendar_month_in_monthly_data!0             = "day"



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
          calendar_year_in_monthly_data(t)  = calendar_year
          calendar_month_in_monthly_data(t) = calendar_month



          ndy_year = 365 ; day_of_year(calendar_year, 12, 31)

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

                  ; e2m goes up wiht temp

                  ; max value of e2m =

                  es2m =  (/ 610.94 * exp( 17.625 * (t2m-273.15) / ( (t2m-273.15) + 243.04)   ) /)

                  rh2m = (/ (e2m/es2m) * 100.0 /)
                  es2m = where(es2m .ge. 0, es2m, 0)

                  rh2m = where(rh2m.gt.100,100,rh2m)
                  rh2m = where(rh2m.lt.  0,  0,rh2m)

                  e2m  = where(rh2m.ge.100, es2m, e2m)
                  e2m  = where(rh2m.le.  0,    0, e2m)

                  td2m = dewtemp_trh(t2m,rh2m)

                  thi  = 41.2 + (t2m  - 273.15) + 0.36 * (td2m  - 273.15)


                  u10m  = (/  af->U10(:, lat_s:lat_n, lon_w:lon_e) /)   ; 10-m u-wind component
                  v10m  = (/  af->V10(:, lat_s:lat_n, lon_w:lon_e) /)   ; 10-m v-wind component
                  m10m  = (/ sqrt( u10m * u10m +  v10m * v10m ) /) ; 10-m wind vector magnitude

                  delete(af)


                  t2m_daily_mean(t,:,:) = (/ dim_avg_n(t2m,0) /)
                  t2m_daily_min(t,:,:)  = (/ dim_min_n(t2m,0) /)
                  t2m_daily_max(t,:,:)  = (/ dim_max_n(t2m,0) /)


                  thi_daily_mean(t,:,:) = (/ dim_avg_n(thi,0) /)
                  thi_daily_min(t,:,:)  = (/ dim_min_n(thi,0) /)
                  thi_daily_max(t,:,:)  = (/ dim_max_n(thi,0) /)


                 td2m_daily_mean(t,:,:) = (/ dim_avg_n(td2m,0) /)
                 td2m_daily_min(t,:,:)  = (/ dim_min_n(td2m,0) /)
                 td2m_daily_max(t,:,:)  = (/ dim_max_n(td2m,0) /)



                  precip_total_daily(t,:,:) = (/ dim_sum_n(rain,0) /)
                  precip_hourly_max(t,:,:)  = (/ dim_max_n(rain,0) /)

                  m10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(m10m,0) /)
                  m10m_wind_speed_daily_max(t,:,:)  =  (/ dim_max_n(m10m,0) /)
                  u10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(u10m,0) /)
                  v10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(v10m,0) /)

               end if

          else

               print("Missing AUX File!")
               missing_aux_days = missing_aux_days + " " +aux_model_file
               aux_file_status(t) = 0

          end if









          if ((nt_auxes .lt. 24) .or. (.not. is_aux_here)) then
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


                 rain  = new( (/ 24, ny, nx /), float )
                 rain2 = new( (/ 24, ny, nx /), float )

                 t2m   = new( (/ 24, ny, nx /), float ) ; 2-m air temp
                 psfc  = new( (/ 24, ny, nx /), float ) ; 2-m air temp
                 q2m   = new( (/ 24, ny, nx /), float ) ; 2-m air temp



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

                  es2m =  (/ 610.94 * exp( 17.625 * (t2m-273.15) / ( (t2m-273.15) + 243.04)   ) /)

                  rh2m = (/ (e2m/es2m) * 100.0 /)
                  es2m = where(es2m .ge. 0, es2m, 0)

                  rh2m = where(rh2m.gt.100,100,rh2m)
                  rh2m = where(rh2m.lt.  0,  0,rh2m)

                  e2m  = where(rh2m.ge.100, es2m, e2m)
                  e2m  = where(rh2m.le.  0,    0, e2m)

                  td2m = dewtemp_trh(t2m,rh2m)

                  thi  = 41.2 + (t2m  - 273.15) + 0.36 * (td2m  - 273.15)

                 thi_daily_mean(t,:,:) = (/ dim_avg_n(thi,0) /)
                 thi_daily_min(t,:,:)  = (/ dim_min_n(thi,0) /)
                 thi_daily_max(t,:,:)  = (/ dim_max_n(thi,0) /)

                 td2m_daily_mean(t,:,:) = (/ dim_avg_n(td2m,0) /)
                 td2m_daily_min(t,:,:)  = (/ dim_min_n(td2m,0) /)
                 td2m_daily_max(t,:,:)  = (/ dim_max_n(td2m,0) /)

                 t2m_daily_mean(t,:,:) = (/ dim_avg_n(t2m,0) /)
                 t2m_daily_min(t,:,:)  = (/ dim_min_n(t2m,0) /)
                 t2m_daily_max(t,:,:)  = (/ dim_max_n(t2m,0) /)

                 precip_total_daily(t,:,:) = (/ dim_sum_n(rain,0) /)
                 precip_hourly_max(t,:,:)  = (/ dim_max_n(rain,0) /)

                 m10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(m10m,0) /)
                 m10m_wind_speed_daily_max(t,:,:)  =  (/ dim_max_n(m10m,0) /)
                 u10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(u10m,0) /)
                 v10m_wind_speed_daily_mean(t,:,:) =  (/ dim_avg_n(v10m,0) /)
                                aux_file_status(t) = 0

                delete(wf)
             end do
          end if




          ;
          ; Destroy Daily 1-hrly and 3-hrly Inputs
          ;

          delete(rain)
          delete(t2m)
          delete(rain2)
          delete(u10m)
          delete(v10m)
          delete(m10m)


      end do

      ;
      ; File Output
      ;

      output_file = "./daily_summary_" + scenario(s) + ".nc"

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

         of->time           = time
         of->lat            = lat
         of->lon            = lon

         of->aux_file_status = aux_file_status

         of->calendar_year_in_monthly_data  = calendar_year_in_monthly_data
         of->calendar_month_in_monthly_data = calendar_month_in_monthly_data

         of->terrain_height                   = terrain_height
         of->terrain_height_variance          = terrain_height_variance
         of->soil_texture                     = soil_texture
         of->land_classification              = land_classification
         of->land_mask                        = land_mask


         of->t2m_daily_mean      = t2m_daily_mean
         of->t2m_daily_max       = t2m_daily_max
         of->t2m_daily_min       = t2m_daily_min

         of->td2m_daily_mean      = td2m_daily_mean
         of->td2m_daily_max       = td2m_daily_max
         of->td2m_daily_min       = td2m_daily_min

         of->thi_daily_mean      = thi_daily_mean
         of->thi_daily_max       = thi_daily_max
         of->thi_daily_min       = thi_daily_min



         of->precip_total_daily  = precip_total_daily
         of->precip_hourly_max   = precip_hourly_max

         of->m10m_wind_speed_daily_mean      = m10m_wind_speed_daily_mean
         of->u10m_wind_speed_daily_mean      = u10m_wind_speed_daily_mean
         of->v10m_wind_speed_daily_mean      = v10m_wind_speed_daily_mean
         of->m10m_wind_speed_daily_max       = m10m_wind_speed_daily_max


         print(of)

      delete(of)

      print("Done with Year " + calendar_year)

   end do



end
