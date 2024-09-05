module atm_import_export

  use shr_kind_mod  , only: r8 => shr_kind_r8, cl=>shr_kind_cl
!   use camsrfexch,     only: cam_in_t, cam_out_t
  use cam_logfile,       only : iulog ! <-- for debugging
  

  implicit none
!   type(cam_out_t), pointer :: cam_out(:) ! WT line

contains

  subroutine atm_import( x2a, cam_in, cam_out, restart_init )

    !-----------------------------------------------------------------------
    use cam_cpl_indices
    use camsrfexch,     only: cam_in_t, cam_out_t
    use phys_grid ,     only: get_ncols_p, get_rlat_p, get_rlon_p
    use ppgrid    ,     only: begchunk, endchunk       
    use shr_const_mod,  only: shr_const_stebol
    use seq_drydep_mod, only: n_drydep
    use co2_cycle     , only: c_i, co2_readFlux_ocn, co2_readFlux_fuel
    use co2_cycle     , only: co2_transport, co2_time_interp_ocn, co2_time_interp_fuel
    use co2_cycle     , only: data_flux_ocn, data_flux_fuel
    use physconst     , only: mwco2
    use time_manager  , only: is_first_step


   !Water isotopes:

   use water_tracer_vars, only: wtrc_nsrfvap, wtrc_iasrfvap, wtrc_indices, wtrc_species
   use water_tracers    , only: wtrc_ratio



    !
    ! Arguments
    !
    real(r8)      , intent(in)    :: x2a(:,:)
    type(cam_in_t), intent(inout) :: cam_in(begchunk:endchunk)
   type(cam_out_t), intent(in) :: cam_out(begchunk:endchunk)
    logical, optional, intent(in) :: restart_init
    !
    ! Local variables
    !		
    integer            :: i,lat,n,c,ig,j  ! indices
    integer            :: ncols         ! number of columns
    logical, save      :: first_time = .true.
    integer, parameter :: ndst = 2
    integer, target    :: spc_ndx(ndst)
    integer, pointer   :: dst_a5_ndx, dst_a7_ndx
    integer, pointer   :: dst_a1_ndx, dst_a3_ndx
    logical :: overwrite_flds
    !water tracers:
    real(r8) :: R  !water tracer ratio

    real(r8)           :: wtlat
    real(r8)           :: wtlon
    real(r8), parameter:: radtodeg = 180.0_r8/SHR_CONST_PI

    !-----------------------------------------------------------------------
    overwrite_flds = .true.
    ! don't overwrite fields if invoked during the initialization phase 
    ! of a 'continue' or 'branch' run type with data from .rs file
    if (present(restart_init)) overwrite_flds = .not. restart_init

    ! ccsm sign convention is that fluxes are positive downward

    ig=1
    do c=begchunk,endchunk
       ncols = get_ncols_p(c) 

       ! initialize constituent surface fluxes to zero
       ! NOTE:overwrite_flds is .FALSE. for the first restart
       ! time step making cflx(:,1)=0.0 for the first restart time step.
       ! cflx(:,1) should not be zeroed out, start the second index of cflx from 2.
       cam_in(c)%cflx(:,2:) = 0._r8 
                                               
       do i =1,ncols                                                               
          if (overwrite_flds) then
             ! Prior to this change, "overwrite_flds" was always .true. therefore wsx and wsy were always updated.
             ! Now, overwrite_flds is .false. for the first time step of the restart run. Move wsx and wsy out of 
             ! this if-condition so that they are still updated everytime irrespective of the value of overwrite_flds.

             ! Move lhf to this if-block so that it is not overwritten to ensure BFB restarts when qneg4 correction 
             ! occurs at the restart time step
             ! Modified by Wuyin Lin
             cam_in(c)%shf(i)    = -x2a(index_x2a_Faxx_sen, ig)     
             cam_in(c)%cflx(i,1) = -x2a(index_x2a_Faxx_evap,ig)                
             cam_in(c)%lhf(i)    = -x2a(index_x2a_Faxx_lat, ig)     
          endif
         ! WT-block begins 
         !Need to define lat/lon for water tracers:
         wtlat = get_rlat_p(c,i)*radtodeg
         wtlon = get_rlon_p(c,i)*radtodeg
         ! WT-bock ends
         !Need to set this before doing water tracers:
         ! cam_in(c)%landfrac(i)  =  x2a(index_x2a_Sf_lfrac, ig) ! land fraction
         ! write(iulog,*) 'qbot shape', SHAPE(cam_out(c)%qbot(:,:))
         ! write(iulog,*) ''

         do j = 1, wtrc_nsrfvap
         if(j .eq. 1) then !Normal water vapour (total)
            cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = - x2a(index_x2a_Faxx_evap, ig)
         else !water tag
            if( -x2a(index_x2a_Faxx_evap,ig) .lt. 0._r8) then !dew/frost?
               ! calculate surface vapour ratio
               R = wtrc_ratio(j,cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(j))),&
              cam_out(c)%qbot(i,wtrc_indices(wtrc_iasrfvap(1))))
               cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = R*-x2a(index_x2a_Faxx_evap,ig)
            else !Sources of water vapour tags

              if(j .eq. 2) then
                  cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)*0.5_r8

              else if(j .eq. 3) then
                  !Latitude band from 90S to 80S, LAT85S
                  if((wtlat >= -90._r8) .and. (wtlat <= -80._r8)) then
                    cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig) 
                  else
                    cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 4) then
                !Latitude band from 80S to 70S, LAT75S
                 if((wtlat > -80._r8) .and. (wtlat <= -70._r8)) then
                    cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                 else
                    cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                 end if

               else if(j .eq. 5) then
                  !Latitude band from 70S to 60S, LAT65S
                  if((wtlat > -70._r8) .and. (wtlat <= -60._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 6) then
                  !Latitude band from 60S to 50S, LAT55S
                  if((wtlat > -60._r8) .and. (wtlat <= -50._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 7) then
                  !Latitude band from 50S to 40S, LAT45S
                  if((wtlat > -50._r8) .and. (wtlat <= -40._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 8) then
                  !Latitude band from 40S to 30S, LAT35S
                  if((wtlat > -40._r8) .and. (wtlat <= -30._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 9) then
                  !Latitude band from 30S to 20S, LAT25S
                  if((wtlat > -30._r8) .and. (wtlat <= -20._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 10) then
                  !Latitude band from 20S to 10S, LAT15S
                  if((wtlat > -20._r8) .and. (wtlat <= -10._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 11) then
                  !Latitude band from 10S to Eq, LAT05S
                  if((wtlat > -10._r8) .and. (wtlat <= 0._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 12) then
                  !Latitude band from Eq to 10N, LAT05N
                  if((wtlat > 0._r8) .and. (wtlat <= 10._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 13) then
                  !Latitude band from 10N to 20N, LAT15N
                  if((wtlat > 10._r8) .and. (wtlat <= 20._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 14) then
                  !Latitude band from 20N to 30N, LAT25N
                  if((wtlat > 20._r8) .and. (wtlat <= 30._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               
               else if(j .eq. 15) then
                  !Latitude band from 30N to 40N, LAT35N
                  if((wtlat > 30._r8) .and. (wtlat <= 40._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 16) then
                  !Latitude band from 40N to 50N, LAT45N
                  if((wtlat > 40._r8) .and. (wtlat <= 50._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 17) then
                  !Latitude band from 50N to 60N, LAT55N
                  if((wtlat > 50._r8) .and. (wtlat <= 60._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 18) then
                  !Latitude band from 60N to 70N, LAT65N
                  if((wtlat > 60._r8) .and. (wtlat <= 70._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 19) then
                  !Latitude band from 70N to 80N, LAT75N
                  if((wtlat > 70._r8) .and. (wtlat <= 80._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 20) then
                  !Latitude band from 80N to 90N, LAT85N
                  if((wtlat > 80._r8) .and. (wtlat <= 90._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 21) then
                  !Longitude band from 0E to 10E, LON05E
                  if((wtlon >= 0._r8) .and. (wtlon <= 10._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 22) then
                  !Longitude band from 10E to 20E, LON15E
                  if((wtlon > 10._r8) .and. (wtlon <= 20._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 23) then
                  !Longitude band from 20E to 30E, LON25E
                  if((wtlon > 20._r8) .and. (wtlon <= 30._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               
               else if(j .eq. 24) then
                  !Longitude band from 30E to 40E, LON35E
                  if((wtlon > 30._r8) .and. (wtlon <= 40._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 25) then
                  !Longitude band from 40E to 50E, LON45E
                  if((wtlon > 40._r8) .and. (wtlon <= 50._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 26) then
                  !Longitude band from 50E to 60E, LON55E
                  if((wtlon > 50._r8) .and. (wtlon <= 60._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 27) then
                  !Longitude band from 60E to 70E, LON65E
                  if((wtlon > 60._r8) .and. (wtlon <= 70._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 28) then
                  !Longitude band from 70E to 80E, LON75E
                  if((wtlon > 70._r8) .and. (wtlon <= 80._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 29) then
                  !Longitude band from 80E to 90E, LON85E
                  if((wtlon > 80._r8) .and. (wtlon <= 90._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 30) then
                  !Longitude band from 90E to 100E, LON95E
                  if((wtlon > 90._r8) .and. (wtlon <= 100._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 31) then
                  !Longitude band from 100E to 110E, LON105E
                  if((wtlon > 100._r8) .and. (wtlon <= 110._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 32) then
                  !Longitude band from 110E to 120E, LON115E
                  if((wtlon > 110._r8) .and. (wtlon <= 120._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 33) then
                  !Longitude band from 120E to 130E, LON125E
                  if((wtlon > 120._r8) .and. (wtlon <= 130._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               
               else if(j .eq. 34) then
                  !Longitude band from 130E to 140E, LON135E
                  if((wtlon > 130._r8) .and. (wtlon <= 140._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 35) then
                  !Longitude band from 140E to 150E, LON145E
                  if((wtlon > 140._r8) .and. (wtlon <= 150._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 36) then
                  !Longitude band from 150E to 160E, LON155E
                  if((wtlon > 150._r8) .and. (wtlon <= 160._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 37) then
                  !Longitude band from 160E to 170E, LON165E
                  if((wtlon > 160._r8) .and. (wtlon <= 170._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 38) then
                  !Longitude band from 170E to 180E, LON175E
                  if((wtlon > 170._r8) .and. (wtlon <= 180._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 39) then
                  !Longitude band from 180E to 190E, LON185E
                  if((wtlon > 180._r8) .and. (wtlon <= 190._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 40) then
                  !Longitude band from 190E to 200E, LON195E
                  if((wtlon > 190._r8) .and. (wtlon <= 200._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               else if(j .eq. 41) then
                  !Longitude band from 200E to 210E, LON205E
                  if((wtlon > 200._r8) .and. (wtlon <= 210._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 42) then
                  !Longitude band from 210E to 220E, LON215E
                  if((wtlon > 210._r8) .and. (wtlon <= 220._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 43) then
                  !Longitude band from 220E to 230E, LON225E
                  if((wtlon > 220._r8) .and. (wtlon <= 230._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               
               else if(j .eq. 44) then
                  !Longitude band from 230E to 240E, LON235E
                  if((wtlon > 230._r8) .and. (wtlon <= 240._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 45) then
                  !Longitude band from 240E to 250E, LON245E
                  if((wtlon > 240._r8) .and. (wtlon <= 250._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 46) then
                  !Longitude band from 250E to 260E, LON255E
                  if((wtlon > 250._r8) .and. (wtlon <= 260._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 47) then
                  !Longitude band from 260E to 270E, LON265E
                  if((wtlon > 260._r8) .and. (wtlon <= 270._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 48) then
                  !Longitude band from 270E to 280E, LON275E
                  if((wtlon > 270._r8) .and. (wtlon <= 280._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 49) then
                  !Longitude band from 280E to 290E, LON285E
                  if((wtlon > 280._r8) .and. (wtlon <= 290._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 50) then
                  !Longitude band from 290E to 300E, LON295E
                  if((wtlon > 290._r8) .and. (wtlon <= 300._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 51) then
                  !Longitude band from 300E to 310E, LON305E
                  if((wtlon > 300._r8) .and. (wtlon <= 310._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 52) then
                  !Longitude band from 310E to 320E, LON315E
                  if((wtlon > 310._r8) .and. (wtlon <= 320._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 53) then
                  !Longitude band from 320E to 330E, LON325E
                  if((wtlon > 320._r8) .and. (wtlon <= 330._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if
               
               else if(j .eq. 54) then
                  !Longitude band from 330E to 340E, LON335E
                  if((wtlon > 330._r8) .and. (wtlon <= 340._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 55) then
                  !Longitude band from 340E to 350E, LON345E
                  if((wtlon > 340._r8) .and. (wtlon <= 350._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               else if(j .eq. 56) then
                  !Longitude band from 350E to 360E, LON355E
                  if((wtlon > 350._r8) .and. (wtlon <= 360._r8)) then
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = -x2a(index_x2a_Faxx_evap,ig)
                  else
                     cam_in(c)%cflx(i,wtrc_indices(wtrc_iasrfvap(j))) = 0._r8
                  end if

               end if ! water tracer index j
            end if ! dew/frost or evap
         end if ! Normal water vapour H20
      end do




         

          if (index_x2a_Faoo_h2otemp /= 0) then
             cam_in(c)%h2otemp(i) = -x2a(index_x2a_Faoo_h2otemp,ig)
          end if
           
          cam_in(c)%wsx(i)    = -x2a(index_x2a_Faxx_taux,ig)     
          cam_in(c)%wsy(i)    = -x2a(index_x2a_Faxx_tauy,ig)     
          cam_in(c)%lwup(i)      = -x2a(index_x2a_Faxx_lwup,ig)    
          cam_in(c)%asdir(i)     =  x2a(index_x2a_Sx_avsdr, ig)  
          cam_in(c)%aldir(i)     =  x2a(index_x2a_Sx_anidr, ig)  
          cam_in(c)%asdif(i)     =  x2a(index_x2a_Sx_avsdf, ig)  
          cam_in(c)%aldif(i)     =  x2a(index_x2a_Sx_anidf, ig)
          cam_in(c)%ts(i)        =  x2a(index_x2a_Sx_t,     ig)  
          cam_in(c)%sst(i)       =  x2a(index_x2a_So_t,     ig)             
          cam_in(c)%snowhland(i) =  x2a(index_x2a_Sl_snowh, ig)  
          cam_in(c)%snowhice(i)  =  x2a(index_x2a_Si_snowh, ig)  
          cam_in(c)%tref(i)      =  x2a(index_x2a_Sx_tref,  ig)  
          cam_in(c)%qref(i)      =  x2a(index_x2a_Sx_qref,  ig)
          cam_in(c)%u10(i)       =  x2a(index_x2a_Sx_u10,   ig)
          cam_in(c)%icefrac(i)   =  x2a(index_x2a_Sf_ifrac, ig)  
          cam_in(c)%ocnfrac(i)   =  x2a(index_x2a_Sf_ofrac, ig)
          cam_in(c)%landfrac(i)  =  x2a(index_x2a_Sf_lfrac, ig)
          if ( associated(cam_in(c)%ram1) ) &
               cam_in(c)%ram1(i) =  x2a(index_x2a_Sl_ram1 , ig)
          if ( associated(cam_in(c)%fv) ) &
               cam_in(c)%fv(i)   =  x2a(index_x2a_Sl_fv   , ig)
          if ( associated(cam_in(c)%soilw) ) &
               cam_in(c)%soilw(i) =  x2a(index_x2a_Sl_soilw, ig)
          if ( associated(cam_in(c)%dstflx) ) then
             cam_in(c)%dstflx(i,1) = x2a(index_x2a_Fall_flxdst1, ig)
             cam_in(c)%dstflx(i,2) = x2a(index_x2a_Fall_flxdst2, ig)
             cam_in(c)%dstflx(i,3) = x2a(index_x2a_Fall_flxdst3, ig)
             cam_in(c)%dstflx(i,4) = x2a(index_x2a_Fall_flxdst4, ig)
          endif
          if ( associated(cam_in(c)%meganflx) ) then
             cam_in(c)%meganflx(i,1:shr_megan_mechcomps_n) = &
                  x2a(index_x2a_Fall_flxvoc:index_x2a_Fall_flxvoc+shr_megan_mechcomps_n-1, ig)
          endif

          ! dry dep velocities
          if ( index_x2a_Sl_ddvel/=0 .and. n_drydep>0 ) then
             cam_in(c)%depvel(i,:n_drydep) = &
                  x2a(index_x2a_Sl_ddvel:index_x2a_Sl_ddvel+n_drydep-1, ig)
          endif
          !
          ! fields needed to calculate water isotopes to ocean evaporation processes
          !
          cam_in(c)%ustar(i) = x2a(index_x2a_So_ustar,ig)
          cam_in(c)%re(i)    = x2a(index_x2a_So_re   ,ig)
          cam_in(c)%ssq(i)   = x2a(index_x2a_So_ssq  ,ig)
          !
          ! bgc scenarios
          !
          if (index_x2a_Fall_fco2_lnd /= 0) then
             cam_in(c)%fco2_lnd(i) = -x2a(index_x2a_Fall_fco2_lnd,ig)
          end if
          if (index_x2a_Faoo_fco2_ocn /= 0) then
             cam_in(c)%fco2_ocn(i) = -x2a(index_x2a_Faoo_fco2_ocn,ig)
          end if
          if (index_x2a_Faoo_fdms_ocn /= 0) then
             cam_in(c)%fdms(i)     = -x2a(index_x2a_Faoo_fdms_ocn,ig)
          end if

          ig=ig+1

       end do
    end do

    ! Get total co2 flux from components,
    ! Note - co2_transport determines if cam_in(c)%cflx(i,c_i(1:4)) is allocated

    if (co2_transport().and.overwrite_flds) then

       ! Interpolate in time for flux data read in
       if (co2_readFlux_ocn) then
          call co2_time_interp_ocn
       end if
       if (co2_readFlux_fuel) then
          call co2_time_interp_fuel
       end if
       
       ! from ocn : data read in or from coupler or zero
       ! from fuel: data read in or zero
       ! from lnd : through coupler or zero
       do c=begchunk,endchunk
          ncols = get_ncols_p(c)                                                 
          do i=1,ncols                                                               
             
             ! all co2 fluxes in unit kgCO2/m2/s ! co2 flux from ocn 
             if (index_x2a_Faoo_fco2_ocn /= 0) then
                cam_in(c)%cflx(i,c_i(1)) = cam_in(c)%fco2_ocn(i)
             else if (co2_readFlux_ocn) then 
                ! convert from molesCO2/m2/s to kgCO2/m2/s
! The below section involves a temporary workaround for fluxes from data (read in from a file)
! There is an issue with infld that does not allow time-varying 2D files to be read correctly.
! The work around involves adding a singleton 3rd dimension offline and reading the files as 
! 3D fields.  Once this issue is corrected, the old implementation can be reinstated.
! This is the case for both data_flux_ocn and data_flux_fuel
!++BEH  vvv old implementation vvv
!                cam_in(c)%cflx(i,c_i(1)) = &
!                     -data_flux_ocn%co2flx(i,c)*(1._r8- cam_in(c)%landfrac(i)) &
!                     *mwco2*1.0e-3_r8
!       ^^^ old implementation ^^^   ///    vvv new implementation vvv
                cam_in(c)%cflx(i,c_i(1)) = &
                     -data_flux_ocn%co2flx(i,1,c)*(1._r8- cam_in(c)%landfrac(i)) &
                     *mwco2*1.0e-3_r8
!--BEH  ^^^ new implementation ^^^
             else
                cam_in(c)%cflx(i,c_i(1)) = 0._r8
             end if
             
             ! co2 flux from fossil fuel
             if (co2_readFlux_fuel) then
!++BEH  vvv old implementation vvv
!                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,c)
!       ^^^ old implementation ^^^   ///    vvv new implementation vvv
                cam_in(c)%cflx(i,c_i(2)) = data_flux_fuel%co2flx(i,1,c)
!--BEH  ^^^ new implementation ^^^
             else
                cam_in(c)%cflx(i,c_i(2)) = 0._r8
             end if
             
             ! co2 flux from land (cpl already multiplies flux by land fraction)
             if (index_x2a_Fall_fco2_lnd /= 0) then
                cam_in(c)%cflx(i,c_i(3)) = cam_in(c)%fco2_lnd(i)
             else
                cam_in(c)%cflx(i,c_i(3)) = 0._r8
             end if
             
             ! merged co2 flux
             cam_in(c)%cflx(i,c_i(4)) = cam_in(c)%cflx(i,c_i(1)) + &
                                        cam_in(c)%cflx(i,c_i(2)) + &
                                        cam_in(c)%cflx(i,c_i(3))
          end do
       end do
    end if
    !
    ! if first step, determine longwave up flux from the surface temperature 
    !
    if (first_time) then
       if (is_first_step()) then
          do c=begchunk, endchunk
             ncols = get_ncols_p(c)
             do i=1,ncols
                cam_in(c)%lwup(i) = shr_const_stebol*(cam_in(c)%ts(i)**4)
             end do
          end do
       end if
       first_time = .false.
    end if

  end subroutine atm_import

  !===============================================================================

  subroutine atm_export( cam_out, a2x )

    !-------------------------------------------------------------------
    use camsrfexch, only: cam_out_t
    use phys_grid , only: get_ncols_p
    use ppgrid    , only: begchunk, endchunk       
    use cam_cpl_indices
    use phys_control, only: phys_getopts
    use lnd_infodata, only: precip_downscaling_method
    !
    ! Arguments
    !
    type(cam_out_t), intent(in)    :: cam_out(begchunk:endchunk) 
    real(r8)       , intent(inout) :: a2x(:,:)
    !
    ! Local variables
    !
    integer :: avsize, avnat
    integer :: i,m,c,n,ig       ! indices
    integer :: ncols            ! Number of columns
    logical :: linearize_pbl_winds, export_gustiness
    !-----------------------------------------------------------------------

    call phys_getopts(linearize_pbl_winds_out=linearize_pbl_winds, &
                      export_gustiness_out=export_gustiness)

    ! Copy from component arrays into chunk array data structure
    ! Rearrange data from chunk structure into lat-lon buffer and subsequently
    ! create attribute vector

    ig=1
    do c=begchunk, endchunk
       ncols = get_ncols_p(c)
       do i=1,ncols
          a2x(index_a2x_Sa_pslv   ,ig) = cam_out(c)%psl(i)
          a2x(index_a2x_Sa_z      ,ig) = cam_out(c)%zbot(i)   
          a2x(index_a2x_Sa_u      ,ig) = cam_out(c)%ubot(i)   
          a2x(index_a2x_Sa_v      ,ig) = cam_out(c)%vbot(i)
          if (linearize_pbl_winds) then
             a2x(index_a2x_Sa_wsresp ,ig) = cam_out(c)%wsresp(i)
             a2x(index_a2x_Sa_tau_est,ig) = cam_out(c)%tau_est(i)
          end if
          if (export_gustiness) then
             a2x(index_a2x_Sa_ugust  ,ig) = cam_out(c)%ugust(i)
          end if
          a2x(index_a2x_Sa_tbot   ,ig) = cam_out(c)%tbot(i)   
          a2x(index_a2x_Sa_ptem   ,ig) = cam_out(c)%thbot(i)  
          a2x(index_a2x_Sa_pbot   ,ig) = cam_out(c)%pbot(i)   
          a2x(index_a2x_Sa_shum   ,ig) = cam_out(c)%qbot(i,1) 
	  a2x(index_a2x_Sa_dens   ,ig) = cam_out(c)%rho(i)

          if (trim(adjustl(precip_downscaling_method)) == "FNM") then
             !if the land model's precip downscaling method is FNM, export uovern to the coupler
             a2x(index_a2x_Sa_uovern ,ig) = cam_out(c)%uovern(i)
          end if
          a2x(index_a2x_Faxa_swnet,ig) = cam_out(c)%netsw(i)      
          a2x(index_a2x_Faxa_lwdn ,ig) = cam_out(c)%flwds(i)  
          a2x(index_a2x_Faxa_rainc,ig) = (cam_out(c)%precc(i)-cam_out(c)%precsc(i))*1000._r8
          a2x(index_a2x_Faxa_rainl,ig) = (cam_out(c)%precl(i)-cam_out(c)%precsl(i))*1000._r8
          a2x(index_a2x_Faxa_snowc,ig) = cam_out(c)%precsc(i)*1000._r8
          a2x(index_a2x_Faxa_snowl,ig) = cam_out(c)%precsl(i)*1000._r8
          a2x(index_a2x_Faxa_swndr,ig) = cam_out(c)%soll(i)   
          a2x(index_a2x_Faxa_swvdr,ig) = cam_out(c)%sols(i)   
          a2x(index_a2x_Faxa_swndf,ig) = cam_out(c)%solld(i)  
          a2x(index_a2x_Faxa_swvdf,ig) = cam_out(c)%solsd(i)  

          ! aerosol deposition fluxes
          a2x(index_a2x_Faxa_bcphidry,ig) = cam_out(c)%bcphidry(i)
          a2x(index_a2x_Faxa_bcphodry,ig) = cam_out(c)%bcphodry(i)
          a2x(index_a2x_Faxa_bcphiwet,ig) = cam_out(c)%bcphiwet(i)
          a2x(index_a2x_Faxa_ocphidry,ig) = cam_out(c)%ocphidry(i)
          a2x(index_a2x_Faxa_ocphodry,ig) = cam_out(c)%ocphodry(i)
          a2x(index_a2x_Faxa_ocphiwet,ig) = cam_out(c)%ocphiwet(i)
          a2x(index_a2x_Faxa_dstwet1,ig)  = cam_out(c)%dstwet1(i)
          a2x(index_a2x_Faxa_dstdry1,ig)  = cam_out(c)%dstdry1(i)
          a2x(index_a2x_Faxa_dstwet2,ig)  = cam_out(c)%dstwet2(i)
          a2x(index_a2x_Faxa_dstdry2,ig)  = cam_out(c)%dstdry2(i)
          a2x(index_a2x_Faxa_dstwet3,ig)  = cam_out(c)%dstwet3(i)
          a2x(index_a2x_Faxa_dstdry3,ig)  = cam_out(c)%dstdry3(i)
          a2x(index_a2x_Faxa_dstwet4,ig)  = cam_out(c)%dstwet4(i)
          a2x(index_a2x_Faxa_dstdry4,ig)  = cam_out(c)%dstdry4(i)

          if (index_a2x_Sa_co2prog /= 0) then
             a2x(index_a2x_Sa_co2prog,ig) = cam_out(c)%co2prog(i) ! atm prognostic co2
          end if
          if (index_a2x_Sa_co2diag /= 0) then
             a2x(index_a2x_Sa_co2diag,ig) = cam_out(c)%co2diag(i) ! atm diagnostic co2
          end if

          ig=ig+1
       end do
    end do
    
  end subroutine atm_export 

end module atm_import_export
