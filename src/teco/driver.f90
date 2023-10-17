module driver

    use datatypes
    use vegetation
    use soil
    use transfer
    use mcmc_mod
    use io_mod

    implicit none
    integer :: unit_h, unit_d, unit_m, unit_y
    logical :: do_out_csv

    contains
    subroutine teco_simu(vegn, do_out_csv_0)
        implicit none
        type(vegn_tile_type), intent(inout) :: vegn
        logical, intent(in) :: do_out_csv_0
        integer year0, first_year                               ! year0: record the current year to judge whether a new year
        real    Difference                                      ! GPP-Rauto-NPP. Jian: no sure whether for balance? 
        ! real    RaLeaf,RaStem,RaRoot                            ! for summary the automatic respiration of leaf, stem, root
        ! integer dlayer                                          ! to run cycle of each layer.
        ! real    Q_soil                                          ! total soil carbon
        ! real    RECOh                                           ! ecosytem respiration
        ! real    ETh, Th, Eh                                     ! record hourly ET, transp, evap. Jian: why not use original variable?
        ! real    INTh,ROh,DRAINh,LEh,SHh                         ! 
        ! real    VPDh, LWH
        real    :: esat1, eairP
        integer :: iclim, iyear, iday, ihour
        integer :: iTotHourly, iTotDaily, iTotMonthly, iTotYearly
        integer :: daysOfyear, daysOfmonth(12), hoursOfYear, hoursOfmonth
        integer :: ipft, iostat
        real    :: snow_depth_e, RH, radsol
        character(len=:), allocatable :: csv_fileName
        character(2000) :: header_csv

        ! Jian: start the cycle of the forcing data
        first_year = forcing(1)%year
        ! initilize the output 
        do_out_csv = do_out_csv_0
        if(do_out_csv) then
            allocate(character(len=200+len(outDir_csv)) :: csv_fileName)
            call def_header(header_csv)
            if(do_out_hr) then
                unit_h = 3987
                call def_csv_fileName(outDir_csv, "Hourly", csv_fileName)
                open(newunit=unit_h, file=csv_fileName, status='replace', action='write', iostat=iostat)
                write(unit_h, *) adjustl(trim(header_csv))
            endif
            if(do_out_day)then
                unit_d = 3988
                call def_csv_fileName(outDir_csv, "Daily", csv_fileName)
                open(newunit=unit_d, file=csv_fileName, status='replace', action='write', iostat=iostat)
                write(unit_d, *) adjustl(trim(header_csv))
            endif
            if(do_out_mon) then
                unit_m = 3989
                call def_csv_fileName(outDir_csv, "Monthly", csv_fileName)
                open(newunit=unit_m, file=csv_fileName, status='replace', action='write', iostat=iostat)
                write(unit_m, *) adjustl(trim(header_csv))
            endif
            if(do_out_yr)then
                unit_y = 3990
                call def_csv_fileName(outDir_csv, "Yearly", csv_fileName)
                open(newunit=unit_y, file=csv_fileName, status='replace', action='write', iostat=iostat)
                write(unit_y, *) adjustl(trim(header_csv))
            endif
            deallocate(csv_fileName)
        endif
        
        do iclim = 1, nforcing 
            ! if(iclim > 10) stop
            if (iclim .eq. 1) then
                year0       = first_year             ! Jian: record whether it is a new year.
                iTotHourly  = 1
                iTotDaily   = 1
                iTotMonthly = 1
                iTotYearly  = 1
            endif
            iyear = forcing(iclim)%year                      ! force%year
            iday  = forcing(iclim)%doy                    
            ihour = forcing(iclim)%hour
            ! if it is a new year
            if ((iday .eq. 1) .and. (ihour .eq. 0)) call init_yearly(vegn)
            if (do_simu .and. (iday .eq. 1) .and. (ihour .eq. 0)) write(*,*)iyear
            if (do_spruce) then
                if ((iyear .eq. 1974) .and. (iday .eq. 1) .and. (ihour .eq. 0))then
                    ! 1974 remove 99% of tree biomass
                    do ipft = 1, vegn%npft
                        vegn%allSp(ipft)%QC(1)    = 0.1 * vegn%allSp(ipft)%QC(1)
                        vegn%allSp(ipft)%QC(2)    = 0.1 * vegn%allSp(ipft)%QC(2)
                        vegn%allSp(ipft)%QC(3)    = 0.1 * vegn%allSp(ipft)%QC(3)
                        vegn%allSp(ipft)%QN(1)    = 0.1 * vegn%allSp(ipft)%QN(1)
                        vegn%allSp(ipft)%QN(2)    = 0.1 * vegn%allSp(ipft)%QN(2)
                        vegn%allSp(ipft)%QN(3)    = 0.1 * vegn%allSp(ipft)%QN(3)
                        vegn%allSp(ipft)%bmleaf   = 0.1 * vegn%allSp(ipft)%bmleaf
                        vegn%allSp(ipft)%bmstem   = 0.1 * vegn%allSp(ipft)%bmstem
                        vegn%allSp(ipft)%bmroot   = 0.1 * vegn%allSp(ipft)%bmroot
                        vegn%allSp(ipft)%nsc      = 0.1 * vegn%allSp(ipft)%nsc
                        vegn%allSp(ipft)%nsn      = 0.1 * vegn%allSp(ipft)%nsn
                        vegn%allSp(ipft)%storage  = 0.1 * vegn%allSp(ipft)%storage
                        vegn%allSp(ipft)%lai      = vegn%allSp(ipft)%LAIMIN!0.1 * lai
                        vegn%allSp(ipft)%stor_use = 0.1 * vegn%allSp(ipft)%stor_use
                    enddo
                    st%QC(1)      = 0.1 * st%QC(1)
                    st%QC(2)      = 0.1 * st%QC(2)
                    st%QC(3)      = 0.1 * st%QC(3)
                    st%QN(1)      = 0.1 * st%QN(1)
                    st%QN(2)      = 0.1 * st%QN(2)
                    st%QN(3)      = 0.1 * st%QN(3)
                    vegn%bmleaf   = 0.1 * vegn%bmleaf
                    vegn%bmstem   = 0.1 * vegn%bmstem
                    vegn%bmroot   = 0.1 * vegn%bmroot
                    vegn%nsc      = 0.1 * vegn%nsc
                    vegn%nsn      = 0.1 * vegn%nsn
                    vegn%storage  = 0.1 * vegn%storage
                    vegn%lai      = vegn%LAIMIN!0.1 * lai
                    vegn%stor_use = 0.1 * vegn%stor_use
                endif
            endif

            ! leap year
            if ((iday .eq. 1) .and. (ihour .eq. 0))then
                if (do_leap) then
                    if(iday .eq. 1) call isLeap_update_daysOfyear(iyear, daysOfyear)
                else
                    daysOfyear = 365
                endif
                ! for update the results of monthly and yearly
                call update_hoursOfYear_daysOfmonth_initMonthly(iday, ihour, &
                        daysOfyear, daysOfmonth, hoursOfYear, hoursOfmonth, iTotMonthly)
            endif

            call update_summary_monthly(iday, ihour, daysOfmonth, iTotMonthly)

            ! initialize the daily variables to run hourly simulaiton.
            if (ihour .eq. 0) then
                ! a new day simulation.
                if (do_snow) then 
                    if (iyear .eq. first_year .and. iday .eq. 1.) then
                        st%ta     = -12.85                      ! since changed the ta criteria (0. to 1.e-10)) in calculating melt
                        st%rain_d = 0.                          ! dbmemo
                    endif
                    call snow_d(iday)                               ! Jian: update snow_dsim snow_d(rain_d,lat,days,ta,snow_dsim,fa,fsub,rho_snow,melt,dcount,decay_m)                            
                    snow_depth_e = st%snow_dsim
                endif
                do ipft = 1, vegn%npft
                    vegn%allSp(ipft)%StemSap = AMIN1(vegn%allSp(ipft)%Stemmax,vegn%allSp(ipft)%SapS*vegn%allSp(ipft)%bmStem)            ! Stemmax and SapS were input from parameter file, what are they? Unit? Maximum stem biomass? -JJJJJJJJJJJJJJJJJJJJJJ 
                    vegn%allSp(ipft)%RootSap = AMIN1(vegn%allSp(ipft)%Rootmax,vegn%allSp(ipft)%SapR*vegn%allSp(ipft)%bmRoot)
                    vegn%allSp(ipft)%NSCmax  = 0.05*(vegn%allSp(ipft)%StemSap+vegn%allSp(ipft)%RootSap+vegn%allSp(ipft)%QC(1))          ! Jian: update the NSCmax each step? and fixed NSCmin  = 5.? 
                enddo
                if(st%Ta.gt.5.0) st%GDD5 = st%GDD5+st%Ta
                call init_daily()                                 ! Jian: initilize the daily data.
            endif

            ! forcing data --------------------------------------------------------------------------------
            ! Tair  = forcing(iclim)%Tair                      ! Tair
            ! Tsoil = forcing(iclim)%Tsoil                     ! SLT
            co2ca = forcing(iclim)%CO2*1.0E-6                       ! CO2 concentration,ppm-->1.0E-6
            if (co2ca .lt. 0) co2ca = 380.0*1.0E-6                  ! Jian: if no CO2 (-9999), then use the default value 
            forcing(iclim)%Tair  = forcing(iclim)%Tair  + Ttreat    ! Jian: whether it has the treatment
            forcing(iclim)%Tsoil = forcing(iclim)%Tsoil + Ttreat
            if (CO2treat .ne. 0.) co2ca = CO2treat*1.0E-6 
            ! ----------------------------------------------------------                
            RH       = forcing(iclim)%RH
            st%Dair  = forcing(iclim)%VPD                      ! air water vapour defficit? Unit Pa
            radsol   = forcing(iclim)%PAR                      ! unit ? PAR actually  Jian: or incoming shortwave/longwave radiation?
            st%dpatm = forcing(iclim)%PBOT
            if (do_ndep) st%N_deposit = forcing(iclim)%Ndep*3600
            ! Ajust some unreasonable values 
            RH      = AMAX1(0.01,AMIN1(99.99,RH))                ! relative humidity
            esat1   = 610.78*exp(17.27*forcing(iclim)%Tair/(forcing(iclim)%Tair + 237.3))      ! intermediate parameter
            eairP   = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure
            st%Dair = esat1-eairP                                ! Jian: confused that SPRUCE has the VPD data, why calculate it again?
            radsol = AMAX1(radsol,0.01)

            ! intially added for soil thermal/ soil water
            if (do_snow) then
                st%snow_depth = snow_depth_e
            else
                st%snow_depth = snow_in(iclim)                  ! read from input file
            endif
            if (st%snow_depth .lt. 0.0) st%snow_depth = 0.0   
            st%snow_depth = st%snow_depth*100.                        ! change from m to cm  
            
            ! Jian: G and Esoil?
            if (do_soilphy) then 
                GOTO 160
            endif
            if(radsol.gt.10.0) then
                st%G = -25.0
            else
                st%G = 20.5
            endif
            if (isnan(st%G)) then
                print *, "st%G is nan", st%G
                stop
            endif
            do ipft = 1, vegn%npft
                vegn%allSp(ipft)%Esoil = (0.05/vegn%npft)* vegn%allSp(ipft)%Esoil
                if (ipft .eq. 1) then
                    st%Esoil = vegn%allSp(ipft)%Esoil
                else
                    st%Esoil = st%Esoil + vegn%allSp(ipft)%esoil
                endif
            enddo
            ! st%Esoil = 0.05*radsol
            if(radsol.LE.10.0) then
                do ipft = 1, vegn%npft
                    vegn%allSp(ipft)%Esoil = (0.5/vegn%npft)* st%G
                    if (ipft .eq. 1) then
                        st%Esoil = vegn%allSp(ipft)%Esoil
                    else
                        st%Esoil = st%Esoil + vegn%allSp(ipft)%esoil
                    endif
                enddo
                ! st%Esoil = 0.5*st%G
            endif
160 continue        
            ! for daily mean conditions 
            st%ta     = st%ta + forcing(iclim)%Tair/24.0                             ! sum of a day, for calculating daily mean temperature, snow_d and soilwater
            st%rain_d = st%rain_d + forcing(iclim)%rain                                
            ! calculating scaling factor of NSC
            do ipft = 1, vegn%npft
                if(vegn%allSp(ipft)%NSC.le.vegn%allSp(ipft)%NSCmin) vegn%allSp(ipft)%fnsc=0.0
                if(vegn%allSp(ipft)%NSC.ge.vegn%allSp(ipft)%NSCmax) vegn%allSp(ipft)%fnsc=1.0
                if((vegn%allSp(ipft)%NSC.lt.vegn%allSp(ipft)%NSCmax).and.(vegn%allSp(ipft)%NSC.gt.vegn%allSp(ipft)%NSCmin))then 
                    vegn%allSp(ipft)%fnsc = (vegn%allSp(ipft)%NSC-vegn%allSp(ipft)%NSCmin) / &
                                            (vegn%allSp(ipft)%NSCmax-vegn%allSp(ipft)%NSCmin)
                endif
                ! update vcmx0 and eJmx0 according to C/N of leaves
                ! print*,"update snvxmax:", vegn%allSp(ipft)%Vcmax0,vegn%allSp(ipft)%SNvcmax*1.0e-6
                vegn%allSp(ipft)%Vcmx0 = vegn%allSp(ipft)%Vcmax0*vegn%allSp(ipft)%SNvcmax*1.0e-6
                vegn%allSp(ipft)%eJmx0 = 1.67*vegn%allSp(ipft)%Vcmx0 ! Weng 02/21/2011 Medlyn et al. 2002 
                vegn%allSp(ipft)%eJmx0 = vegn%allSp(ipft)%JV*vegn%allSp(ipft)%Vcmx0   ! added for acclimation study,replace 1.67 with JV Feb 19 2019 Shuang  
            enddo
            ! print*,"before LAI: ", vegn%allSp(1)%LAI
            call vegn_canopy(vegn, forcing(iclim))      ! run canopy module
            ! run soil water processes
            call soilwater(vegn, forcing(iclim))                      
            st%ET = st%evap + vegn%transp
            
            ! ! Jian: to update module
            call respiration(vegn, forcing(iclim))
            ! THE Third Part: update LAI
            call vegn_plantgrowth(vegn, forcing(iclim))

            ! THE Fourth PART: simulating C influx allocation in pools
            ! print*, "before TCS_CN: ", vegn%allSp(1)%npp
            call TCS_CN(vegn, forcing(iclim)) 
            ! print*, "after TCS_CN: ", vegn%allSp(1)%npp  
            ! if (do_matrix) call matrix_struct() 
            call methane(vegn, forcing(iclim))        !update single value of Rh_pools,Tsoil,zwt,wsc 
            ! update NSC
            do ipft = 1, vegn%npft
                vegn%allSp(ipft)%Rauto = vegn%allSp(ipft)%Rmain + vegn%allSp(ipft)%Rgrowth + vegn%allSp(ipft)%Rnitrogen
                vegn%allSp(ipft)%NSC   = vegn%allSp(ipft)%NSC   + vegn%allSp(ipft)%GPP     - vegn%allSp(ipft)%Rauto - &
                                          (vegn%allSp(ipft)%NPP - vegn%allSp(ipft)%add)    - vegn%allSp(ipft)%store
                Difference = vegn%allSp(ipft)%GPP - vegn%allSp(ipft)%Rauto - vegn%allSp(ipft)%NPP
                if(vegn%allSp(ipft)%NSC < 0)then
                    vegn%allSp(ipft)%bmstem = vegn%allSp(ipft)%bmstem + vegn%allSp(ipft)%NSC/0.48
                    vegn%allSp(ipft)%NPP    = Amax1(vegn%allSp(ipft)%NPP + vegn%allSp(ipft)%NSC, 0.) 
                    vegn%allSp(ipft)%NSN    = vegn%allSp(ipft)%NSN    - vegn%allSp(ipft)%NSC/vegn%allSp(ipft)%CN(2)
                    vegn%allSp(ipft)%NSC    = 0.
                endif
                ! vegn%allSp(ipft)%NSN     = 0.35
                vegn%allSp(ipft)%bmleaf  = vegn%allSp(ipft)%QC(1)/0.48
                vegn%allSp(ipft)%bmstem  = vegn%allSp(ipft)%QC(2)/0.48
                vegn%allSp(ipft)%bmroot  = vegn%allSp(ipft)%QC(3)/0.48
                vegn%allSp(ipft)%bmplant = vegn%allSp(ipft)%bmleaf + vegn%allSp(ipft)%bmroot + vegn%allSp(ipft)%bmstem
                ! print *, "check LAI: ",vegn%allSp(ipft)%bmleaf, vegn%allSp(ipft)%SLA
                vegn%allSp(ipft)%LAI     = vegn%allSp(ipft)%bmleaf*vegn%allSp(ipft)%SLA
                ! summary
                ! if (vegn%allSp(ipft)%gpp>0) print*, "gpp > 0",  vegn%allSp(ipft)%gpp
                if (ipft .eq. 1) then
                    vegn%Rauto   = vegn%allSp(ipft)%Rauto
                    vegn%gpp     = vegn%allSp(ipft)%gpp
                    vegn%npp     = vegn%allSp(ipft)%npp
                    vegn%NPP_L   = vegn%allSp(ipft)%NPP_L
                    vegn%NPP_W   = vegn%allSp(ipft)%NPP_W
                    vegn%NPP_R   = vegn%allSp(ipft)%NPP_R
                    vegn%Rgrowth = vegn%allSp(ipft)%Rgrowth
                    st%Rnitrogen = vegn%allSp(ipft)%Rnitrogen
                    vegn%NSC     = vegn%allSp(ipft)%NSC
                    vegn%NSN     = vegn%allSp(ipft)%NSN
                else
                    vegn%Rauto = vegn%Rauto + vegn%allSp(ipft)%Rauto
                    vegn%gpp   = vegn%gpp   + vegn%allSp(ipft)%gpp
                    vegn%npp   = vegn%npp   + vegn%allSp(ipft)%npp
                    vegn%NPP_L = vegn%NPP_L + vegn%allSp(ipft)%NPP_L
                    vegn%NPP_W = vegn%NPP_W + vegn%allSp(ipft)%NPP_W
                    vegn%NPP_R = vegn%NPP_R + vegn%allSp(ipft)%NPP_R
                    vegn%Rgrowth = vegn%Rgrowth + vegn%allSp(ipft)%Rgrowth
                    st%Rnitrogen = st%Rnitrogen + vegn%allSp(ipft)%Rnitrogen
                    vegn%NSC     = vegn%NSC + vegn%allSp(ipft)%NSC
                    vegn%NSN     = vegn%NSN + vegn%allSp(ipft)%NSN
                endif
            enddo 
            ! print*,"end LAI: ", vegn%allSp(1)%LAI, vegn%allSp(1)%bmleaf, vegn%allSp(1)%SLA

            ! Rhetero=Rh_f + Rh_c + Rh_Micr + Rh_Slow + Rh_Pass
            st%Rhetero = st%Rh_pools(1) + st%Rh_pools(2) + st%Rh_pools(3) &
                &      + st%Rh_pools(4) + st%Rh_pools(5)

            st%NEE     = vegn%Rauto+st%Rhetero - vegn%GPP
            

            ! call updateHourly(vegn, iclim, iyear, iday, ihour)    ! hourly simulation
            call init_hourly()
            call updateOutVars(vegn, outvars_h, 1, iyear, iday, ihour)

            ! call updateDaily(vegn, iTotDaily, iyear, iday, ihour)
            call updateOutVars(vegn, outvars_d, 24, iyear, iday, ihour)
            ! call updateMonthly(vegn, iTotMonthly, hoursOfmonth, iyear, iday, ihour)
            call updateOutVars(vegn, outvars_m, hoursOfmonth, iyear, iday, ihour)
            ! call updateYearly(vegn, iTotYearly, hoursOfYear, iyear, iday, ihour)
            call updateOutVars(vegn, outvars_y, hoursOfYear, iyear, iday, ihour)

            if (do_mcmc) call GetSimuData(iyear, iday, ihour, vegn, iclim, iTotDaily, iTotMonthly, iTotYearly)

            if(do_out_csv) then
                if(do_out_hr) call write_data_csv(unit_h, outVars_h)
            endif
            
            if (ihour .eq. 23) then
                if(do_out_csv) then 
                    if(do_out_day) call write_data_csv(unit_d, outVars_d)
                endif
                iTotDaily = iTotDaily + 1 
            endif
                 
            if (iclim < nforcing)then
                if (forcing(iclim+1)%year>iyear) then            
                    year0        = iyear                      ! update the record of year (year0)
                    if(do_out_csv) then 
                        if(do_out_yr) call write_data_csv(unit_y, outVars_y)
                    endif
                    iTotYearly   = iTotYearly + 1
                    do ipft = 1, vegn%npft
                        vegn%allSp(ipft)%storage      = vegn%allSp(ipft)%accumulation
                        vegn%allSp(ipft)%stor_use     = vegn%allSp(ipft)%Storage/times_storage_use
                        vegn%allSp(ipft)%accumulation = 0.0
                        vegn%allSp(ipft)%onset        = 0
                    enddo
                endif
            else
                year0        = iyear                          ! update the record of year (year0)
                if(do_out_csv) then
                    if(do_out_yr) call write_data_csv(unit_y, outVars_y)
                endif
                do ipft = 1, vegn%npft
                    vegn%allSp(ipft)%storage      = vegn%allSp(ipft)%accumulation
                    vegn%allSp(ipft)%stor_use     = vegn%allSp(ipft)%Storage/times_storage_use
                    vegn%allSp(ipft)%accumulation = 0.0
                    vegn%allSp(ipft)%onset        = 0
                enddo
            endif
        enddo
        if(do_out_csv)then
            if(do_out_hr)  close(unit_h)
            if(do_out_day) close(unit_d)
            if(do_out_mon) close(unit_m)
            if(do_out_yr)  close(unit_y)
        endif
    end subroutine teco_simu

    subroutine isLeap_update_daysOfyear(iyear, daysOfyear)    
        implicit none
        integer, intent(in) :: iyear
        integer, intent(inout) :: daysOfyear
        if (MOD(iyear, 4) .eq. 0)then
            if (MOD(iyear, 100) .eq. 0)then
                if (MOD(iyear, 400) .eq. 0)then
                    daysOfyear = 366
                else
                    daysOfyear = 365
                endif
            else
                daysOfyear = 366
            endif
        else
            daysOfyear = 365
        endif
        return
    end subroutine isLeap_update_daysOfyear

    subroutine update_hoursOfYear_daysOfmonth_initMonthly(iday, ihour, &
        daysOfyear, daysOfmonth, hoursOfYear, hoursOfmonth, iTotMonthly)
        implicit none
        integer, intent(in)    :: iday, ihour
        integer, intent(inout) :: daysOfyear, daysOfmonth(12)
        integer, intent(inout) :: hoursOfYear, hoursOfmonth, iTotMonthly
        if (daysOfyear .eq. 365) then ! common year
            hoursOfYear = 365*24
            daysOfmonth = (/31,59,90,120,151,181,212,243,273,304,334,365/)
        else
            hoursOfYear = 366*24
            daysOfmonth = (/31,60,91,121,152,182,213,244,274,305,335,366/)
        endif
        ! hours of month
        ! January:
        if (iday .eq. 1)then 
            hoursOfmonth = (daysOfmonth(1)-0)*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Feburay:
        if (iday .eq. daysOfmonth(1)+1)then
            hoursOfmonth = (daysOfmonth(2)-daysOfmonth(1))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! March
        if (iday .eq. daysOfmonth(2)+1)then
            hoursOfmonth = (daysOfmonth(3)-daysOfmonth(2))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! April
        if (iday .eq. daysOfmonth(3)+1)then
            hoursOfmonth = (daysOfmonth(4)-daysOfmonth(3))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! May
        if (iday .eq. daysOfmonth(4)+1)then
            hoursOfmonth = (daysOfmonth(5)-daysOfmonth(4))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! June
        if (iday .eq. daysOfmonth(5)+1)then
            hoursOfmonth = (daysOfmonth(6)-daysOfmonth(5))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! July
        if(iday .eq. daysOfmonth(6)+1)then
            hoursOfmonth = (daysOfmonth(7)-daysOfmonth(6))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Auguest
        if(iday .eq. daysOfmonth(7)+1)then
            hoursOfmonth = (daysOfmonth(8)-daysOfmonth(7))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! Septemble
        if(iday .eq. daysOfmonth(8)+1)then
            hoursOfmonth = (daysOfmonth(9)-daysOfmonth(8))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! October
        if(iday .eq. daysOfmonth(9)+1)then
            hoursOfmonth = (daysOfmonth(10)-daysOfmonth(9))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! November
        if(iday .eq. daysOfmonth(10)+1)then
            hoursOfmonth = (daysOfmonth(11)-daysOfmonth(10))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        ! December
        if(iday .eq. daysOfmonth(11)+1)then
            hoursOfmonth = (daysOfmonth(12)-daysOfmonth(11))*24
            if (ihour .eq. 0) call init_monthly()
        endif
        return
    end subroutine update_hoursOfYear_daysOfmonth_initMonthly

    subroutine update_summary_monthly(iday, ihour, daysOfmonth, iTotMonthly)
        implicit none
        integer, intent(in) :: iday, ihour, daysOfmonth(12)
        integer, intent(inout) :: iTotMonthly
        ! January
        if ((iday .eq. daysOfmonth(1))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(2))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(3))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(4))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(5))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(6))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(7))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(8))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(9))  .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(10)) .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(11)) .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
        if ((iday .eq. daysOfmonth(12)) .and. (ihour .eq. 23)) then
            if(do_out_csv) then
                if(do_out_mon) call write_data_csv(unit_m, outVars_m)
            endif
            iTotMonthly = iTotMonthly + 1
        endif
    end subroutine update_summary_monthly
    
end module driver