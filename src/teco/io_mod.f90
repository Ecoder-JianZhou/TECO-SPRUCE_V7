!#define USE_NETCDF

module io_mod
    use datatypes
#ifdef USE_NETCDF
        use netcdf
#endif
    implicit none
    CHARACTER(len=4) :: str_startyr, str_endyr
    real convert_g2kg, convert_h2s
    
    contains

    subroutine updateOutVars(vegn, outvars, ntime, iyear, iday, ihour)
        implicit none
        integer, intent(in) :: ntime, iyear, iday, ihour
        type(outvars_data_type), intent(inout) :: outVars
        type(vegn_tile_type), intent(in) :: vegn
        integer :: ipft, npft
        ! integer iTotHourly
        outVars%year = iyear
        outVars%doy  = iday
        outVars%hour = ihour
        ! stop
        convert_g2kg = 1 !0.001
        convert_h2s  = 1!1/3600.
        if (allocated(outvars%allSpec)) then
            npft = size(outvars%allSpec)
            do ipft = 1, npft
                ! carbon fluxes (Kg C m-2 s-1)
                outvars%allSpec(ipft)%gpp      = outvars%allSpec(ipft)%gpp     + &
                                                    vegn%allSp(ipft)%gpp*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%npp      = outvars%allSpec(ipft)%npp     + &
                                                    vegn%allSp(ipft)%npp*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%nppLeaf  = outvars%allSpec(ipft)%nppLeaf + &
                                                    vegn%allSp(ipft)%NPP_L*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%nppWood  = outvars%allSpec(ipft)%nppWood + &
                                                    vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%nppStem  = outvars%allSpec(ipft)%nppStem + &
                                                    vegn%allSp(ipft)%NPP_W*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%nppRoot  = outvars%allSpec(ipft)%nppRoot + &
                                                    vegn%allSp(ipft)%NPP_R*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%nppOther = outvars%allSpec(ipft)%nppOther + &
                                                    vegn%allSp(ipft)%NSC*convert_g2kg*convert_h2s/ntime    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
                outvars%allSpec(ipft)%ra       = outvars%allSpec(ipft)%ra      + &
                                                    vegn%allSp(ipft)%Rauto*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%raLeaf   = outvars%allSpec(ipft)%raLeaf  + &
                                                    vegn%allSp(ipft)%RmLeaf*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%raStem   = outvars%allSpec(ipft)%raStem  + &
                                                    vegn%allSp(ipft)%RmStem*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%raRoot   = outvars%allSpec(ipft)%raRoot  + &
                                                    vegn%allSp(ipft)%RmRoot*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%rMaint   = outvars%allSpec(ipft)%rMaint  + &
                                                    vegn%allSp(ipft)%Rmain*convert_g2kg*convert_h2s/ntime
                outvars%allSpec(ipft)%rGrowth  = outvars%allSpec(ipft)%rGrowth + &
                                                    vegn%allSp(ipft)%Rgrowth*convert_g2kg*convert_h2s/ntime
                ! Carbon Pools  (KgC m-2)
                outvars%allSpec(ipft)%cLeaf    = outvars%allSpec(ipft)%cLeaf  + &
                                                    vegn%allSp(ipft)%QC(1)*convert_g2kg/ntime
                outvars%allSpec(ipft)%cStem    = outvars%allSpec(ipft)%cStem  + &
                                                    vegn%allSp(ipft)%QC(2)*convert_g2kg/ntime
                outvars%allSpec(ipft)%cRoot    = outvars%allSpec(ipft)%cRoot  + &
                                                    vegn%allSp(ipft)%QC(3)*convert_g2kg/ntime
                ! Nitrogen pools (kgN m-2)
                outvars%allSpec(ipft)%nLeaf    = outvars%allSpec(ipft)%nLeaf  + &
                                                    vegn%allSp(ipft)%QN(1)*convert_g2kg/ntime
                outvars%allSpec(ipft)%nStem    = outvars%allSpec(ipft)%nStem  + &
                                                    vegn%allSp(ipft)%QN(2)*convert_g2kg/ntime
                outvars%allSpec(ipft)%nRoot    = outvars%allSpec(ipft)%nRoot  + &
                                                    vegn%allSp(ipft)%QN(3)*convert_g2kg/ntime
                ! water fluxes (kg m-2 s-1)
                outvars%allSpec(ipft)%tran     = outvars%allSpec(ipft)%tran   + &
                                                    vegn%allSp(ipft)%transp*convert_g2kg*convert_h2s/ntime
                ! other
                outvars%allSpec(ipft)%lai      = outvars%allSpec(ipft)%lai    + &
                                                    vegn%allSp(ipft)%lai/ntime
            enddo
        endif
        ! carbon fluxes (KgC m-2 s-1) Jian: TECO unit is gC m-2 h-1
        outvars%gpp             = outvars%gpp      + vegn%gpp*convert_g2kg*convert_h2s/ntime
        outvars%npp             = outvars%npp      + vegn%npp*convert_g2kg*convert_h2s/ntime
        outvars%nppLeaf         = outvars%nppLeaf  + vegn%NPP_L*convert_g2kg*convert_h2s/ntime
        outvars%nppWood         = outvars%nppWood  + vegn%NPP_W*convert_g2kg*convert_h2s/ntime  
        outvars%nppStem         = outvars%nppStem  + vegn%NPP_W*convert_g2kg*convert_h2s/ntime 
        outvars%nppRoot         = outvars%nppRoot  + vegn%NPP_R*convert_g2kg*convert_h2s/ntime
        outvars%nppOther        = outvars%nppOther + vegn%NSC*convert_g2kg*convert_h2s/ntime 
        outvars%ra              = outvars%ra       + vegn%Rauto*convert_g2kg*convert_h2s/ntime
        outvars%raLeaf          = outvars%raLeaf   + vegn%Rmleaf*convert_g2kg*convert_h2s/ntime
        outvars%raStem          = outvars%raStem   + vegn%Rmstem*convert_g2kg*convert_h2s/ntime
        outvars%raRoot          = outvars%raRoot   + vegn%Rmroot*convert_g2kg*convert_h2s/ntime
        outvars%raOther         = outvars%raOther  + st%Rnitrogen *convert_g2kg*convert_h2s/ntime
        outvars%rMaint          = outvars%rMaint   + vegn%Rmain *convert_g2kg*convert_h2s/ntime
        outvars%rGrowth         = outvars%rGrowth  + vegn%Rgrowth *convert_g2kg*convert_h2s/ntime 
        outvars%rh              = outvars%rh       + st%Rhetero *convert_g2kg*convert_h2s/ntime 
        outvars%nbp             = outvars%nbp      + &
                                    (vegn%gpp - st%Rhetero - vegn%Rauto) *convert_g2kg*convert_h2s/ntime   
        outvars%wetlandCH4      = outvars%wetlandCH4 + st%simuCH4 *convert_g2kg*convert_h2s/ntime   
        outvars%wetlandCH4prod  = outvars%wetlandCH4prod + st%Pro_sum *convert_g2kg*convert_h2s/ntime 
        outvars%wetlandCH4cons  = outvars%wetlandCH4cons + st%Oxi_sum *convert_g2kg*convert_h2s/ntime 
        ! Carbon Pools  (KgC m-2)
        outvars%cLeaf           = outvars%cLeaf + st%QC(1)*convert_g2kg/ntime
        outvars%cStem           = outvars%cStem + st%QC(2)*convert_g2kg/ntime
        outvars%cRoot           = outvars%cRoot + st%QC(3)*convert_g2kg/ntime
        outvars%cOther          = outvars%cOther+ vegn%NSC*convert_g2kg/ntime
        outvars%cLitter         = outvars%cLitter + st%QC(4)*convert_g2kg/ntime
        outvars%cLitterCwd      = outvars%cLitterCwd + st%QC(5)*convert_g2kg/ntime
        outvars%cSoil           = outvars%cSoil + (st%QC(6) + st%QC(7) + st%QC(8))*convert_g2kg/ntime
        outvars%cSoilLevels(:)  = outvars%cSoilLevels(:) + (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)/ntime
        outvars%cSoilFast       = outvars%cSoilFast + st%QC(6)*convert_g2kg/ntime 
        outvars%cSoilSlow       = outvars%cSoilSlow + st%QC(7)*convert_g2kg/ntime 
        outvars%cSoilPassive    = outvars%cSoilPassive + st%QC(8)*convert_g2kg/ntime 
        outvars%CH4(:)          = outvars%CH4(:) + st%CH4*convert_g2kg/ntime 
        ! Nitrogen fluxes (kgN m-2 s-1)
        outvars%fBNF            = outvars%fBNF + st%N_fixation*convert_g2kg*convert_h2s/ntime 
        outvars%fN2O            = outvars%fN2O + &
                                    (st%N_transfer+st%N_uptake+st%N_fixation)*convert_g2kg*convert_h2s/ntime
        outvars%fNloss          = outvars%fNloss + &
                                    (vegn%N_leaf+vegn%N_wood+vegn%N_root)*convert_g2kg*convert_h2s/ntime
        outvars%fNnetmin        = outvars%fNnetmin + st%fNnetmin*convert_g2kg*convert_h2s/ntime 
        outvars%fNdep           = outvars%fNdep + st%N_deposit*convert_g2kg*convert_h2s/ntime 
        ! Nitrogen pools (kgN m-2)
        outvars%nLeaf           = outvars%nLeaf      + st%QN(1)*convert_g2kg/ntime
        outvars%nStem           = outvars%nStem      + st%QN(2)*convert_g2kg/ntime
        outvars%nRoot           = outvars%nRoot      + st%QN(3)*convert_g2kg/ntime
        outvars%nOther          = outvars%nOther     + vegn%NSN*convert_g2kg/ntime
        outvars%nLitter         = outvars%nLitter    + st%QN(4)*convert_g2kg/ntime
        outvars%nLitterCwd      = outvars%nLitterCwd + st%QN(5)*convert_g2kg/ntime
        outvars%nSoil           = outvars%nSoil      + (st%QN(6)+st%QN(7)+st%QN(8))*convert_g2kg/ntime
        outvars%nMineral        = outvars%nMineral   + st%QNminer*convert_g2kg/ntime 
        ! energy fluxes (W m-2)
        outvars%hfls            = outvars%hfls + st%Hsoil/ntime ! Sensible heat flux;
        outvars%hfss            = outvars%hfss + st%Esoil/ntime ! Latent heat flux;
        outvars%SWnet           = outvars%SWnet + 0/ntime       ! Net shortwave radiation;
        outvars%LWnet           = outvars%LWnet + 0/ntime       ! Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        outvars%ec              = outvars%ec    + 0/ntime       !evap*convert_g2kg*convert_h2s/ntime        ! Canopy evaporation;
        outvars%tran            = outvars%tran  + vegn%transp*convert_g2kg*convert_h2s/ntime      ! Canopy transpiration;
        outvars%es              = outvars%es    + st%evap*convert_g2kg*convert_h2s/ntime ! Soil evaporation
        outvars%hfsbl           = outvars%hfsbl + st%sublim*convert_g2kg*convert_h2s/ntime ! Snow sublimation
        outvars%mrro            = outvars%mrro  + st%runoff*convert_g2kg*convert_h2s/ntime
        ! outvars%mrros         = forcing(iforcing)%Rain    
        outvars%mrrob           = outvars%mrrob + 0/ntime ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        outvars%mrso(:)         = outvars%mrso(:) + st%liq_water*1000/ntime  ! Kg m-2, soil moisture in each soil layer
        outvars%tsl(:)          = outvars%tsl(:) + (st%tsoil_layer(1:10)+273.15)/ntime                            ! K, soil temperature in each soil layer Jian: not sure the tsoil_layer is correct or not
        ! outvars%tsland        = forcing(iforcing)%Tair+273.15                                   ! K, surface temperature
        outvars%wtd             = outvars%wtd +  (st%zwt/1000)/ntime                                       ! m, Water table depth
        outvars%snd             = outvars%snd +  (st%snow_depth/100)/ntime                               ! m, Total snow depth, Jian: change from m to cm in code, and now change from cm to m
        outvars%lai             = outvars%lai +  (vegn%LAI)/ntime                                           ! m2 m-2, Leaf area index
       
    end subroutine updateOutVars

! ========================================================================================================
!  This part is used to write the outputs to the csv-format file
!       open_file_csv, write_data_csv, close_data_csv
!  ---------------------------------------------------------------
    subroutine def_header(header_csv)
        implicit none
        character(*), intent(inout) :: header_csv
        integer :: ipft

        ! Write header line
        header_csv = "year,doy,hour,"
        do ipft = 1, count_pft
            header_csv = adjustl(trim(header_csv))//"gpp_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nee_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"npp_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nppLeaf_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nppWood_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nppStem_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nppRoot_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nppOther_"//adjustl(trim(spec_names(ipft)))//","    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
            header_csv = adjustl(trim(header_csv))//"ra_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"raLeaf_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"raStem_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"raRoot_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"raOther_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"rMaint_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"rGrowth_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nbp_"//adjustl(trim(spec_names(ipft)))//","
            ! Carbon Pools  (KgC m-2)
            header_csv = adjustl(trim(header_csv))//"cLeaf_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"cStem_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"cRoot_"//adjustl(trim(spec_names(ipft)))//","
            ! Nitrogen pools (kgN m-2)
            header_csv = adjustl(trim(header_csv))//"nLeaf_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nStem_"//adjustl(trim(spec_names(ipft)))//","
            header_csv = adjustl(trim(header_csv))//"nRoot_"//adjustl(trim(spec_names(ipft)))//","
            ! water fluxes (kg m-2 s-1)
            header_csv = adjustl(trim(header_csv))//"tran_"//adjustl(trim(spec_names(ipft)))//","
            ! other
            header_csv = adjustl(trim(header_csv))//"lai_"//adjustl(trim(spec_names(ipft)))//"," 
        enddo
        header_csv = adjustl(trim(header_csv))//"gpp,nee,npp,nppLeaf,nppWood,nppStem,nppRoot,nppOther,ra,&
            raLeaf,raStem,raRoot,raOther,rMaint,rGrowth,rh,nbp,wetlandCH4,wetlandCH4prod,&
            wetlandCH4cons,cLeaf,cStem,cRoot,cOther,cLitter,cLitterCwd,cSoil,&
            cSoilLevels_1,cSoilLevels_2,cSoilLevels_3,cSoilLevels_4,cSoilLevels_5,&
            cSoilLevels_6,cSoilLevels_7,cSoilLevels_8,cSoilLevels_9,cSoilLevels_10,&
            cSoilFast,cSoilSlow,cSoilPassive,CH4_1,CH4_2,CH4_3,CH4_4,CH4_5,CH4_6,CH4_7,&
            CH4_8,CH4_9,CH4_10,fBNF,fN2O,fNloss,fNnetmin,fNdep,nLeaf,nStem,nRoot,nOther,&
            nLitter,nLitterCwd,nSoil,nMineral,hfls,hfss,SWnet,LWnet,ec,tran,es,hfsbl,&
            mrro,mrros,mrrob,mrso_1,mrso_2,mrso_3,mrso_4,mrso_5,mrso_6,mrso_7,mrso_8,&
            mrso_9,mrso_10,tsl_1,tsl_2,tsl_3,tsl_4,tsl_5,tsl_6,tsl_7,tsl_8,tsl_9,tsl_10,&
            &tsland,wtd,snd,lai"
        return
    end subroutine def_header

    subroutine def_csv_fileName(out_path, str_freq, csv_fileName)
        implicit none
        character(*), intent(in) :: out_path, str_freq
        character(*), intent(inout) :: csv_fileName

        csv_fileName = adjustl(trim(out_path))//"/TECO-SPRUCE_"//adjustl(trim(case_name))//"_"//str_freq//".csv"
! windows
#if defined(WIN32) || defined(WIN64) || defined(_WIN32) || defined(_WIN64)
        csv_fileName = adjustl(trim(out_path))//"\TECO-SPRUCE_"//adjustl(trim(case_name))//"_"//str_freq//".csv"
#endif
        return
    end subroutine def_csv_fileName

    subroutine write_data_csv(unit, outVars)
        implicit none
        integer, intent(in) :: unit
        type(outvars_data_type), intent(in) :: outVars
        integer :: ipft, ilayer, nformat
        character(len=2500) :: format_string
        
        ! write the date
        nformat       = 24*count_pft+98 
        format_string = '((i4,",")(i3,",")(i2,",")' // repeat('(f15.4, ",")', nformat-1) // 'f15.4)'
        write(unit, adjustl(trim(format_string)))outVars%year, outVars%doy, outVars%hour, &
            (outVars%allSpec(ipft)%gpp,     outVars%allSpec(ipft)%nee,      outVars%allSpec(ipft)%npp,     &     
            outVars%allSpec(ipft)%nppLeaf,  outVars%allSpec(ipft)%nppWood,  outVars%allSpec(ipft)%nppStem, &
            outVars%allSpec(ipft)%nppRoot,  outVars%allSpec(ipft)%nppOther, outVars%allSpec(ipft)%ra,      & 
            outVars%allSpec(ipft)%raLeaf,   outVars%allSpec(ipft)%raStem,   outVars%allSpec(ipft)%raRoot,  & 
            outVars%allSpec(ipft)%raOther,  outVars%allSpec(ipft)%rMaint,   outVars%allSpec(ipft)%rGrowth, &
            outVars%allSpec(ipft)%nbp,      outVars%allSpec(ipft)%cLeaf,    outVars%allSpec(ipft)%cStem,   & 
            outVars%allSpec(ipft)%cRoot,    outVars%allSpec(ipft)%nLeaf,    outVars%allSpec(ipft)%nStem,   &  
            outVars%allSpec(ipft)%nRoot,    outVars%allSpec(ipft)%tran,     outVars%allSpec(ipft)%lai,     &    
            ipft = 1, count_pft),& 
            outVars%gpp,     outVars%nee,        outVars%npp,            outVars%nppLeaf,  &        
            outVars%nppWood, outVars%nppStem,    outVars%nppRoot,        outVars%nppOther, &   
            outVars%ra,      outVars%raLeaf,     outVars%raStem,         outVars%raRoot,   &  
            outVars%raOther, outVars%rMaint,     outVars%rGrowth,        outVars%rh,       &
            outVars%nbp,     outVars%wetlandCH4, outVars%wetlandCH4prod, outVars%wetlandCH4cons,   &
            outVars%cLeaf,   outVars%cStem,      outVars%cRoot,          outVars%cOther,   &
            outVars%cLitter, outVars%cLitterCwd, outVars%cSoil,           &
            (outVars%cSoilLevels(ilayer), ilayer = 1, nlayers),           &
            outVars%cSoilFast,  outVars%cSoilSlow, outVars%cSoilPassive,  &
            (outVars%CH4(ilayer), ilayer = 1, nlayers),                   &
            outVars%fBNF,       outVars%fN2O,   outVars%fNloss,   outVars%fNnetmin,   outVars%fNdep,   &
            outVars%nLeaf,      outVars%nStem,  outVars%nRoot,    outVars%nOther,     outVars%nLitter, &
            outVars%nLitterCwd, outVars%nSoil,  outVars%nMineral, outVars%hfls,       outVars%hfss,    &
            outVars%SWnet,      outVars%LWnet,  outVars%ec,       outVars%tran,       outVars%es,      &
            outVars%hfsbl,      outVars%mrro,   outVars%mrros,    outVars%mrrob,                       &     
            (outVars%mrso(ilayer), ilayer = 1, nlayers), &
            (outVars%tsl(ilayer),  ilayer = 1, nlayers), &
            outVars%tsland,     outVars%wtd,    outVars%snd,      outVars%lai
        return
    end subroutine write_data_csv

! #ifdef USE_NETCDF
!     subroutine write_outputs_nc(out_path, outVars, nSimuLen, str_freq)
!         ! Daily and monthly
!         ! carbon flux (KgC m-2 s-1): gpp, npp, nppLeaf, nppWood, nppRoot, nppOther,
!         !              ra, raLeaf, raStem, raRoot, raOther, rMaint, rGrowth, rh
!         !              nbp (=gpp - Rh - Ra - other losses)
!         !              wetlandCH4, wetlandCH4prod, wetlandCH4cons
!         ! carbon pools (KgC m-2): cLeaf, cStem, cRoot, cOther, cLitter (excluding coarse wood debris), cLitterCwd
!         !              cSoil, cSoilLevels, cSoilPools (soil organic carbon for each pool), CH4 (Methane concentration)
!         ! Nitrogen flux (KgN m-2 s-1) : fBNF(biological nitrogen fixation), fN2O, fNloss, fNnetmin, fNdep
!         ! Nitrogen pools (KgN m-2): nleaf, nStem, nRoot, nOther, nLitter, nLitterCwd, nSoil, nMineral
!         ! Energy Fluxes (W m-2): hfls(sensible heat flux), hfss(Latent heat flux), SWnet (Net Shortwave radiation), LWnet(Net Longwave radiation)
!         ! Water Fluxes  (Kg m-2 s-1): ec(canopy evaporation), tran(canopy transpiration), es(soil evaporation), hfsbl (snow sublimation), mrro(total runoff),
!         !                 mrros (surface runoff), mrrob(subsurface runoff)
!         ! other         : mrso (soil moisture in each soil layer, Kg m-2), tsl(soil temperature in each soil layer, K), tsland(surface temperature, K),
!         !                 wtd (Water table depth, m), snd (total snow depth, m), lai(m2 m-2) 
!         ! ===================================================================================================================================================
!         ! carbon fluxes variables
!         ! ----------:-----------:----------:-----------------------
!         implicit none
!         character(*), intent(in) :: out_path, str_freq
!         type(outvars_data_type), allocatable, intent(in) :: outVars
!         integer, intent(in) :: nSimuLen
!         integer :: ipft

!         write(str_startyr,"(I4)")forcing(1)%year
!         write(str_endyr,"(I4)")forcing(nforcing)%year

!         if (allocated(outVars%allSpec)) then
!             do ipft = 1, count_pft
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%gpp,      "gpp_"//adjustl(trim(spec_names(ipft))),     &
!                     "kgC m-2 s-1", "gross primary productivity",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nee,      "nee_"//adjustl(trim(spec_names(ipft))),      &
!                     "kgC m-2 s-1", "net ecosystem exchange",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%npp,      "npp_"//adjustl(trim(spec_names(ipft))),      &
!                     "kgC m-2 s-1", "net primary productivity",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nppLeaf,  "nppLeaf_"//adjustl(trim(spec_names(ipft))),  &
!                     "kgC m-2 s-1", "NPP allocated to leaf tissues",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nppWood,  "nppWood_"//adjustl(trim(spec_names(ipft))),  &
!                     "kgC m-2 s-1", "NPP allocated to above ground woody tissues",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nppStem,  "nppStem_"//adjustl(trim(spec_names(ipft))),  &
!                     "kgC m-2 s-1", "NPP allocated to stem tissues",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nppRoot,  "nppRoot_"//adjustl(trim(spec_names(ipft))),  &
!                     "kgC m-2 s-1", "NPP allocated to root tissues",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nppOther, "nppOther_"//adjustl(trim(spec_names(ipft))), &
!                     "kgC m-2 s-1", "NPP allocated to other plant organs (reserves, fruits, exudates)",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%ra, "ra_"//adjustl(trim(spec_names(ipft))), &
!                     "kgC m-2 s-1", "Plant Autotrophic Respiration",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%raLeaf,   "raLeaf_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2 s-1", "Ra from leaves",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%raStem,   "raStem_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2 s-1", "Ra from above ground woody tissues",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%raRoot,   "raRoot_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2 s-1", "Ra from fine roots",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%raOther,   "raOther_"//adjustl(trim(spec_names(ipft))), &
!                     "kgC m-2 s-1", "Ra from other plant organs (reserves, fruits, exudates)",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%rMaint,    "rMaint_"//adjustl(trim(spec_names(ipft))),  &
!                     "kgC m-2 s-1", "Maintenance respiration",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%rGrowth,   "rGrowth_"//adjustl(trim(spec_names(ipft))), &
!                     "kgC m-2 s-1", "Growth respiration",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nbp,       "nbp_"//adjustl(trim(spec_names(ipft))),     &
!                     "kgC m-2 s-1", "Net Biome productivity (NBP = GPP - Rh - Ra - other losses)",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%cLeaf,     "cLeaf_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2", "Carbon biomass in leaves",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%cStem,     "cStem_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2", "Carbon above ground woody biomass",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%cRoot,     "cRoot_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgC m-2", "Carbon biomass in roots",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nLeaf,     "nLeaf_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgN m-2", "Nitrogen in leaves",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nStem,     "nStem_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgN m-2", "Nitrogen in stems",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%nRoot,     "nRoot_"//adjustl(trim(spec_names(ipft))),   &
!                     "kgN m-2", "Nirogen in roots",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%tran,      "tran_"//adjustl(trim(spec_names(ipft))),    &
!                     "kg m-2 s-1", "Canopy transpiration",str_freq,1)
!                 call write_nc(out_path, nSimuLen, outVars%allSpec(ipft)%lai,       "lai_"//adjustl(trim(spec_names(ipft))),     &
!                     "m2 m-2", "Leaf area index",str_freq,1)
!             enddo
!         endif

!         ! outputs
!         call write_nc(out_path, nSimuLen, outVars%gpp,"gpp","kgC m-2 s-1", "gross primary productivity",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%npp,"npp","kgC m-2 s-1", "Total net primary productivity",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nppLeaf,"nppLeaf","kgC m-2 s-1", "NPP allocated to leaf tissues",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nppWood,"nppWood","kgC m-2 s-1", &
!             & "NPP allocated to above ground woody tissues",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nppStem,"nppStem","kgC m-2 s-1", "NPP allocated to stem tissues",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nppRoot,"nppRoot","kgC m-2 s-1", "NPP allocated to root tissues",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nppOther,"nppOther","kgC m-2 s-1", &
!             & "NPP allocated to other plant organs (reserves, fruits, exudates)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%ra,"ra","kgC m-2 s-1", "Plant Autotrophic Respiration",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%raLeaf,"raLeaf","kgC m-2 s-1", "Ra from leaves",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%raStem,"raStem","kgC m-2 s-1", "Ra from above ground woody tissues",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%raRoot,"raRoot","kgC m-2 s-1", "Ra from fine roots",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%raOther,"raOther","kgC m-2 s-1", &
!             & "Ra from other plant organs (reserves, fruits, exudates)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%rMaint,"rMaint","kgC m-2 s-1", "Maintenance respiration",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%rGrowth,"rGrowth","kgC m-2 s-1", "Growth respiration",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%rh,"rh","kgC m-2 s-1", "Heterotrophic respiration rate",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nbp,"nbp","kgC m-2 s-1", &
!             &"Net Biome productivity (NBP = GPP - Rh - Ra - other losses)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%wetlandCH4,"wetlandCH4","kgC m-2 s-1", "Net fluxes of CH4",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%wetlandCH4prod,"wetlandCH4prod","kgC m-2 s-1", "CH4 production",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%wetlandCH4cons,"wetlandCH4cons","kgC m-2 s-1", "CH4 consumption",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cLeaf,"cLeaf","kgC m-2", "Carbon biomass in leaves",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cStem,"cStem","kgC m-2", "Carbon above ground woody biomass",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cRoot,"cRoot","kgC m-2", "Carbon biomass in roots",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cOther,"cOther","kgC m-2", &
!             & "Carbon biomass in other plant organs (reserves, fruits)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cLitter,"cLitter","kgC m-2", &
!             & "Carbon in litter (excluding coarse woody debris)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cLitterCwd,"cLitterCwd","kgC m-2", "Carbon in coarse woody debris",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cSoil,"cSoil","kgC m-2", "Total soil organic carbon",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cSoilLevels,"cSoilLevels","kgC m-2", &
!             & "Depth-specific soil organic carbon",str_freq,nlayers)
!         call write_nc(out_path, nSimuLen, outVars%cSoilFast,"cSoilFast","kgC m-2", "Fast soil organic carbon",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cSoilSlow,"cSoilSlow","kgC m-2 s-1", "Slow soil organic carbon",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%cSoilPassive,"cSoilPassive","kgC m-2 s-1", &
!             & "Passive soil organic carbon",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%CH4,"CH4","kgC m-2 s-1", "Methane concentration",str_freq,nlayers)
!         call write_nc(out_path, nSimuLen, outVars%fBNF,"fBNF","kgN m-2 s-1", "biological nitrogen fixation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%fN2O,"fN2O","kgN m-2 s-1", &
!             & "loss of nitrogen through emission of N2O",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%fNloss,"fNloss","kgN m-2 s-1", &
!             & "Total loss of nitrogen to the atmosphere and from leaching",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%fNnetmin,"fNnetmin","kgN m-2 s-1", "net mineralization of N",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%fNdep,"fNdep","kgN m-2 s-1", "Nitrogen deposition",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nLeaf,"nLeaf","kgN m-2", "Nitrogen in leaves",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nStem,"nStem","kgN m-2", "Nitrogen in stems",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nRoot,"nRoot","kgN m-2", "Nirogen in roots",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nOther,     "nOther","kgN m-2", &
!             & "nitrogen in other plant organs (reserves, fruits)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nLitter,    "nLitter","kgN m-2", &
!             & "Nitrogen in litter (excluding coarse woody debris)",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nLitterCwd, "nLitterCwd","kgN m-2", &
!             & "Nitrogen in coarse woody debris",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nSoil,"nSoil","kgN m-2", "Nitrogen in soil organic matter",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%nMineral,"nMineral","kgN m-2", "Mineral nitrogen pool",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%hfls,"hfls","W m-2", "Sensible heat flux",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%hfss,"hfss","W m-2", "Latent heat flux",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%SWnet,"SWnet","W m-2", "Net shortwave radiation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%LWnet,"LWnet","W m-2", "Net longwave radiation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%ec,"ec","kg m-2 s-1", "Canopy evaporation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%tran,"tran","kg m-2 s-1", "Canopy transpiration",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%es,"es","kg m-2 s-1", "Soil evaporation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%hfsbl,"hfsbl","kg m-2 s-1", "Snow sublimation",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%mrro,"mrro","kg m-2 s-1", "Total runoff",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%mrros,"mrros","kg m-2 s-1", "Surface runoff",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%mrrob,"mrrob","kg m-2 s-1", "Subsurface runoff",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%mrso,"mrso","kg m-2", "soil moisture in each soil layer",str_freq,nlayers)
!         call write_nc(out_path, nSimuLen, outVars%tsl,"tsl","K", "soil temperature in each soil layer",str_freq,nlayers)
!         call write_nc(out_path, nSimuLen, outVars%tsland,"tsland","K", "surface temperature",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%wtd,"wtd","m", "Water table depth",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%snd,"snd","m", "Total snow depth",str_freq,1)
!         call write_nc(out_path, nSimuLen, outVars%lai,"lai","m2 m-2", "Leaf area index",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_gdd5_h,"GDD5","m2 m-2", "GDD5",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_onset_h,"onset","m2 m-2", "onset",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_storage_h,"storage","m2 m-2", "onset",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_add_h,"add","m2 m-2", "onset",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_accumulation_h,"accumulation","m2 m-2", "accumulation",str_freq,1)
!         ! call write_nc(out_path, nSimuLen, all_test_h,"test_gpp","m2 m-2", "test_gpp",str_freq,9)
        
!     end subroutine write_outputs_nc

!     subroutine write_nc(outfile, nSimuLen, data, varName, unit, description, str_freq, nSoilLayer)
!         IMPLICIT NONE
!         real(kind=4), Dimension(nSimuLen, nSoilLayer), intent(in) :: data
!         integer(kind=4) :: nSoilLayer
!         integer(KIND=4) :: ncid, timid, dp_dimid, timvarid
!         integer(kind=4) :: varid
!         integer(kind=4), intent(in) :: nSimuLen
!         CHARACTER(LEN=*), INTENT(IN) :: outfile, str_freq
!         CHARACTER(len=*), intent(in) :: varName, unit, description
!         character(len=:), allocatable :: nc_fileName
!         character(len=100) :: timeUnit
!         integer itime
!         real, dimension(nSimuLen) :: time_values 
!         integer :: start(1), count(1)
        
!         allocate(character(len=200+len(outfile)) :: nc_fileName)
!         nc_fileName = adjustl(trim(outfile))//"/"//adjustl(trim(varName))//"_"//str_freq//"_TECO-SPRUCE_"//&
!             & adjustl(trim(case_name))//"_"//adjustl(trim(str_startyr))//"-"//adjustl(trim(str_endyr))//".nc"   
        
!         !Create the netCDF file.
!         CALL check(nf90_create(nc_fileName, NF90_CLOBBER, ncid))

!         !Define the dimensions.
!         ! CALL check(nf90_def_dim(ncid, "nSimu", nfreq,    simuid))
!         CALL check(nf90_def_dim(ncid, "time",  nSimuLen, timid))
    
!         if (nSoilLayer>1)then
!             call check(nf90_def_dim(ncid, "depth", nSoilLayer, dp_dimid))
!             CALL check(nf90_def_var(ncid = ncid, name = varName,  xtype = NF90_FLOAT, &
!                 & dimids = (/timid, dp_dimid/),  varID =varid))
!         else
!             CALL check(nf90_def_var(ncid = ncid, name = varName,  xtype = NF90_FLOAT, &
!                 & dimids = (/timid/),  varID =varid))
!         endif

!         call check(nf90_def_var(ncid, "time",  NF90_DOUBLE, timid,  timvarid))
!         !Define data variable
        
!         !Add attributes
!         if (str_freq .eq. "hourly") then
!             timeUnit = "hours since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
!         else if (str_freq .eq. "daily") then
!             timeUnit = "days since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
!         else if (str_freq .eq. "monthly") then
!             timeUnit = "months since "//adjustl(trim(str_startyr))//"-01-01 00:00:00"
!         end if
        
!         call check(nf90_put_att(ncid,timvarid,"units",adjustl(trim(timeUnit))))
!         CALL check(nf90_put_att(ncid,varid,"units",unit))
!         CALL check(nf90_put_att(ncid,varid,"description",description))
!         CALL check(nf90_enddef(ncid)) 
!         !End Definitions

!         !Write Data
!         ! if (nSoilLayer>1)then
!         !     do i = 1, nSoilLayer
!         !         CALL check(nf90_put_var(ncid, varid, data, start=[1,i], count=[nSimuLen,1]))
!         !     enddo
!         ! else

!         do itime = 1, nSimuLen
!             time_values(itime) = itime-1
!         enddo
!         start = 1
!         count = nSimuLen

!         CALL check(nf90_put_var(ncid, timvarid, time_values,start,count))
!         CALL check(nf90_put_var(ncid, varid, data))
        
!         CALL check(nf90_close(ncid))
!     end subroutine write_nc

!     ! check (ever so slightly modified from www.unidata.ucar.edu)
!     subroutine check(istatus)
!         ! use netcdf
!         implicit none
!         integer, intent(in) :: istatus
!         if(istatus /= nf90_noerr) then
!             write(*,*) trim(adjustl(nf90_strerror(istatus)))
!         end if
!     end subroutine check
! #endif
end module io_mod