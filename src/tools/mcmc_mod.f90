module mcmc_mod
    ! some functions that both driver module and MCMC module.
    use datatypes
    implicit none

    ! parameters and observation files

    integer npar, nDAsimu, ncov, nRand, nSpecParams,upgraded,iDAsimu
    real search_scale
    logical :: do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr

    real mc_lat, mc_Longitude, mc_wsmax, mc_wsmin
    real mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax
    real mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n
    real mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha
    real mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F
    real mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive
    real mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0
    real mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre
    real mc_Tveg, mc_Tpro_me, mc_Toxi
    real mc_f, mc_bubprob, mc_Vmaxfraction
    real mc_Q10rh, mc_JV, mc_Entrpy
    real mc_etaL, mc_etaW, mc_etaR
    real mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S
    real mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M

    ! type params_mcmc
    !     real lat, Longitude, wsmax, wsmin
    !     real LAIMAX, LAIMIN, rdepth, Rootmax, Stemmax
    !     real SapR, SapS, SLA, GLmax, GRmax, Gsmax, stom_n
    !     real a1, Ds0, Vcmx0, extkU, xfang, alpha
    !     real Tau_Leaf, Tau_Wood, Tau_Root, Tau_F
    !     real Tau_C,  Tau_Micro, Tau_SlowSOM, Tau_Passive
    !     real gddonset, Q10, Rl0, Rs0, Rr0
    !     real r_me, Q10pro, kCH4, Omax, CH4_thre
    !     real Tveg, Tpro_me, Toxi
    !     real f, bubprob, Vmaxfraction
    !     real Q10rh, JV, Entrpy
    !     real etaL, etaW, etaR
    !     real f_F2M, f_C2M, f_C2S, f_M2S
    !     real f_M2P, f_S2P, f_S2M, f_P2M
    ! end type params_mcmc

    type params_mcmc
        real, allocatable :: parval(:)
        real, allocatable :: parmin(:)
        real, allocatable :: parmax(:)
    end type params_mcmc
    type(params_mcmc), allocatable :: mc_parvals(:)

    type params_DApar
        real, allocatable :: DAparmin(:)
        real, allocatable :: DAparmax(:)
        real, allocatable :: DApar(:)
        real, allocatable :: DApar_old(:)
        integer, allocatable :: DAparidx(:)
        real, allocatable :: gamma(:,:)
        real, allocatable :: gamnew(:,:)
        real, allocatable :: coefhistory(:,:)
        real, allocatable :: coefnorm(:)
        real, allocatable :: coefac(:)
    end type params_DApar
    type(params_DApar), allocatable :: mc_DApar(:)          ! all variables for DATA ASSIMILATION

    ! type(nml_params_data_type) :: in_params, in_parval, in_parval_min, in_parval_max
    ! real, allocatable :: parval(:), parmin(:), parmax(:)
    character(20), allocatable :: parnames(:)

    ! observational file path
    character(500) :: obsfile_ANPP_Shrub_y
    character(500) :: obsfile_ANPP_Tree_y
    character(500) :: obsfile_NPP_sphag_y
    character(500) :: obsfile_BNPP_y        ! tree + shrub
    character(500) :: obsfile_er_d          ! shrub + sphag.
    character(500) :: obsfile_er_h          ! shrub + sphag.
    character(500) :: obsfile_gpp_d         ! Shrub + sphag.
    character(500) :: obsfile_nee_d         ! Shrub + sphag.
    character(500) :: obsfile_nee_h         ! shrub + sphag.
    character(500) :: obsfile_LAI_d         ! tree  + Shrub
    !
    character(500) :: obsfile_leaf_mass_shrub_y
    character(500) :: obsfile_stem_mass_shrub_y
    character(500) :: obsfile_leaf_resp_shrub_d 
    character(500) :: obsfile_leaf_resp_tree_d 
    ! methane
    character(500) :: obsfile_ch4_d 
    character(500) :: obsfile_ch4_h 
    ! 
    character(500) :: obsfile_CN_shag_d 
    character(500) :: obsfile_photo_shrub_d 
    character(500) :: obsfile_photo_tree_d  

    ! variables for calculating the cost in MCMC processes
    type interCostVariable
        character(300) :: filepath
        logical :: existOrNot
        real, allocatable :: obsData(:,:)
        real, allocatable :: mdData(:,:)
        integer :: mc_itime
    end type interCostVariable

    type allCostVariables
    ! default variables, you can add the variable names here. (year, doy, hour, value, std.)
        ! carbon flux 
        type(interCostVariable) :: ANPP_Shrub_y
        type(interCostVariable) :: ANPP_Tree_y
        type(interCostVariable) :: NPP_sphag_y
        type(interCostVariable) :: BNPP_y        ! tree + shrub
        type(interCostVariable) :: er_d          ! shrub + sphag.
        type(interCostVariable) :: er_h          ! shrub + sphag.
        type(interCostVariable) :: gpp_d         ! Shrub + sphag.
        type(interCostVariable) :: nee_d         ! Shrub + sphag.
        type(interCostVariable) :: nee_h         ! shrub + sphag.
        type(interCostVariable) :: LAI_d         ! tree  + Shrub
        !
        type(interCostVariable) :: leaf_mass_shrub_y
        type(interCostVariable) :: stem_mass_shrub_y
        type(interCostVariable) :: leaf_resp_shrub_d 
        type(interCostVariable) :: leaf_resp_tree_d 
        ! methane
        type(interCostVariable) :: ch4_d 
        type(interCostVariable) :: ch4_h 
        ! 
        type(interCostVariable) :: CN_shag_d 
        type(interCostVariable) :: photo_shrub_d 
        type(interCostVariable) :: photo_tree_d 
    end type allCostVariables

    type(allCostVariables) :: vars4MCMC      ! define a allCostVariables first

    ! variables for marking the cycle number
    integer mc_itime_gpp_d, mc_itime_nee_d, mc_itime_reco_d
    integer mc_itime_gpp_h, mc_itime_nee_h, mc_itime_reco_h
    integer mc_itime_ch4_h, mc_itime_cleaf, mc_itime_cwood
    integer mc_itime_anpp_y, mc_itime_bnpp_y, mc_itime_lai_h
    integer mc_itime_npp_y, mc_itime_reco_y
    integer mc_iyear,  mc_iday, mc_ihour

    type mcmc_spec_outvars_type
        ! carbon fluxes (Kg C m-2 s-1)
        real, allocatable :: gpp(:, :)
        real, allocatable :: nee(:, :)
        real, allocatable :: npp(:, :)
        real, allocatable :: nppLeaf(:, :)
        real, allocatable :: nppWood(:, :)
        real, allocatable :: nppStem(:, :)
        real, allocatable :: nppRoot(:, :)
        real, allocatable :: nppOther(:, :)    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real, allocatable :: ra(:, :)
        real, allocatable :: raLeaf(:, :)
        real, allocatable :: raStem(:, :)
        real, allocatable :: raRoot(:, :)
        real, allocatable :: raOther(:, :)
        real, allocatable :: rMaint(:, :)
        real, allocatable :: rGrowth(:, :)
        real, allocatable :: nbp(:, :)
        ! Carbon Pools  (KgC m-2)
        real, allocatable :: cLeaf(:, :)
        real, allocatable :: cStem(:, :)
        real, allocatable :: cRoot(:, :)
        ! Nitrogen pools (kgN m-2)
        real, allocatable :: nLeaf(:, :)
        real, allocatable :: nStem(:, :)
        real, allocatable :: nRoot(:, :)
        ! real, allocatable :: nOther(:)
        ! water fluxes (kg m-2 s-1)
        real, allocatable :: tran(:, :)
        ! other
        real, allocatable :: lai(:, :)                     ! m2 m-2, Leaf area index
    end type mcmc_spec_outvars_type

    type mcmc_outVars_type
        type(mcmc_spec_outvars_type), allocatable :: allSpec(:)
        ! carbon fluxes (Kg C m-2 s-1)
        real, allocatable :: gpp(:, :)
        real, allocatable :: nee(:, :)
        real, allocatable :: npp(:, :)
        real, allocatable :: nppLeaf(:, :)
        real, allocatable :: nppWood(:, :)
        real, allocatable :: nppStem(:, :)
        real, allocatable :: nppRoot(:, :)
        real, allocatable :: nppOther(:, :)           ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real, allocatable :: ra(:, :)
        real, allocatable :: raLeaf(:, :)
        real, allocatable :: raStem(:, :)
        real, allocatable :: raRoot(:, :)
        real, allocatable :: raOther(:, :)
        real, allocatable :: rMaint(:, :)
        real, allocatable :: rGrowth(:, :)            ! maintenance respiration and growth respiration
        real, allocatable :: rh(:, :)
        real, allocatable :: nbp(:, :)                ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        real, allocatable :: wetlandCH4(:, :)
        real, allocatable :: wetlandCH4prod(:, :)
        real, allocatable :: wetlandCH4cons(:, :)     ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        real, allocatable :: cLeaf(:, :)
        real, allocatable :: cStem(:, :)
        real, allocatable :: cRoot(:, :)
        real, allocatable :: cOther(:, :)              ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        real, allocatable :: cLitter(:, :)
        real, allocatable :: cLitterCwd(:, :)          ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        real, allocatable :: cSoil(:, :)
        real, allocatable :: cSoilLevels(:, :, :)
        real, allocatable :: cSoilFast(:, :)
        real, allocatable :: cSoilSlow(:, :)
        real, allocatable :: cSoilPassive(:, :)           ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        real, allocatable :: CH4(:, :, :)          ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        real, allocatable :: fBNF(:, :)
        real, allocatable :: fN2O(:, :)
        real, allocatable :: fNloss(:, :)
        real, allocatable :: fNnetmin(:, :)
        real, allocatable :: fNdep(:, :)                   ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        real, allocatable :: nLeaf(:, :)
        real, allocatable :: nStem(:, :)
        real, allocatable :: nRoot(:, :)
        real, allocatable :: nOther(:, :)
        real, allocatable :: nLitter(:, :)
        real, allocatable :: nLitterCwd(:, :)
        real, allocatable :: nSoil(:, :)
        real, allocatable :: nMineral(:, :)                ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        real, allocatable :: hfls(:, :)
        real, allocatable :: hfss(:, :)
        real, allocatable :: SWnet(:, :)
        real, allocatable :: LWnet(:, :)                   ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        real, allocatable :: ec(:, :)
        real, allocatable :: tran(:, :)
        real, allocatable :: es(:, :)                      ! Canopy evaporation; Canopy transpiration; Soil evaporation
        real, allocatable :: hfsbl(:, :)                   ! Snow sublimation
        real, allocatable :: mrro(:, :)
        real, allocatable :: mrros(:, :)
        real, allocatable :: mrrob(:, :)                   ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        real, allocatable :: mrso(:, :, :)           ! Kg m-2, soil moisture in each soil layer
        real, allocatable :: tsl(:, :, :)            ! K, soil temperature in each soil layer
        real, allocatable :: tsland(:, :)                  ! K, surface temperature
        real, allocatable :: wtd(:, :)                     ! m, Water table depth
        real, allocatable :: snd(:, :)                     ! m, Total snow depth
        real, allocatable :: lai(:, :)                     ! m2 m-2, Leaf area index            
    end type mcmc_outVars_type

    type(mcmc_outVars_type) sel_paramsets_outs_h
    type(mcmc_outVars_type) sel_paramsets_outs_d
    type(mcmc_outVars_type) sel_paramsets_outs_m
    ! total simulation outputs
    type(mcmc_outVars_type) tot_paramsets_outs_h
    type(mcmc_outVars_type) tot_paramsets_outs_d
    type(mcmc_outVars_type) tot_paramsets_outs_m

    contains

    subroutine mcmc_functions_init()
        implicit none
        vars4MCMC%ANPP_Shrub_y%mc_itime = 1
        vars4MCMC%ANPP_Tree_y%mc_itime  = 1
        vars4MCMC%NPP_sphag_y%mc_itime  = 1
        vars4MCMC%BNPP_y%mc_itime       = 1
        vars4MCMC%er_d%mc_itime         = 1
        vars4MCMC%er_h%mc_itime         = 1
        vars4MCMC%gpp_d%mc_itime        = 1
        vars4MCMC%nee_d%mc_itime        = 1
        vars4MCMC%nee_h%mc_itime        = 1
        vars4MCMC%LAI_d%mc_itime        = 1
        
        vars4MCMC%leaf_mass_shrub_y%mc_itime = 1
        vars4MCMC%stem_mass_shrub_y%mc_itime = 1
        vars4MCMC%leaf_resp_shrub_d%mc_itime = 1
        vars4MCMC%leaf_resp_tree_d%mc_itime  = 1 

        vars4MCMC%ch4_d%mc_itime = 1
        vars4MCMC%ch4_h%mc_itime = 1

        vars4MCMC%CN_shag_d%mc_itime     = 1
        vars4MCMC%photo_shrub_d%mc_itime = 1
        vars4MCMC%photo_tree_d%mc_itime  = 1

        mc_iyear = 1
        mc_iday  = 1
        mc_ihour = 1
    end subroutine mcmc_functions_init

    subroutine readConfsNml()
        integer ipft, npft
        integer io
        character(20) :: parnames_1, parnames_2, parnames_3, parnames_4, parnames_5 
        character(20) :: parnames_6, parnames_7, parnames_8, parnames_9, parnames_10

        character(20) :: parnames_11, parnames_12, parnames_13, parnames_14, parnames_15 
        character(20) :: parnames_16, parnames_17, parnames_18, parnames_19, parnames_20

        character(20) :: parnames_21, parnames_22, parnames_23, parnames_24, parnames_25 
        character(20) :: parnames_26, parnames_27, parnames_28, parnames_29, parnames_30

        character(20) :: parnames_31, parnames_32, parnames_33, parnames_34, parnames_35 
        character(20) :: parnames_36, parnames_37, parnames_38, parnames_39, parnames_40

        character(20) :: parnames_41, parnames_42, parnames_43, parnames_44, parnames_45 
        character(20) :: parnames_46, parnames_47, parnames_48, parnames_49, parnames_50

        character(20) :: parnames_51, parnames_52, parnames_53, parnames_54, parnames_55 
        character(20) :: parnames_56, parnames_57, parnames_58, parnames_59, parnames_60 

        namelist /nml_obsfiles/ obsfile_ANPP_Shrub_y, obsfile_ANPP_Tree_y, obsfile_NPP_sphag_y, &
            obsfile_BNPP_y, obsfile_er_d, obsfile_er_h, obsfile_gpp_d, obsfile_nee_d, &
            obsfile_nee_h, obsfile_LAI_d, obsfile_leaf_mass_shrub_y, obsfile_stem_mass_shrub_y, &
            obsfile_leaf_resp_shrub_d, obsfile_leaf_resp_tree_d, obsfile_ch4_d, obsfile_ch4_h,  & 
            obsfile_CN_shag_d, obsfile_photo_shrub_d, obsfile_photo_tree_d

        namelist /nml_param_names/parnames_1, parnames_2, parnames_3, parnames_4, parnames_5, & 
                parnames_6, parnames_7, parnames_8, parnames_9, parnames_10, &
                parnames_11, parnames_12, parnames_13, parnames_14, parnames_15, & 
                parnames_16, parnames_17, parnames_18, parnames_19, parnames_20, &
                parnames_21, parnames_22, parnames_23, parnames_24, parnames_25, &
                parnames_26, parnames_27, parnames_28, parnames_29, parnames_30, &
                parnames_31, parnames_32, parnames_33, parnames_34, parnames_35, &
                parnames_36, parnames_37, parnames_38, parnames_39, parnames_40, &
                parnames_41, parnames_42, parnames_43, parnames_44, parnames_45, & 
                parnames_46, parnames_47, parnames_48, parnames_49, parnames_50, &
                parnames_51, parnames_52, parnames_53, parnames_54, parnames_55, &
                parnames_56, parnames_57, parnames_58, parnames_59, parnames_60

        namelist /nml_mcmc_settings/ nDAsimu, search_scale, ncov, nRand, &
                do_mc_out_hr, do_mc_out_day, do_mc_out_mon, do_mc_out_yr, nSpecParams

        allocate(parnames(60))
        print *, "mcmc_configfile: ", adjustl(trim("configs/"))//adjustl(trim(mcmc_configfile))
        open(145, file=adjustl(trim("configs/"))//adjustl(trim(mcmc_configfile)))
        read(145, nml=nml_mcmc_settings, iostat=io)
        read(145, nml=nml_obsfiles,      iostat=io)
        read(145, nml=nml_param_names,   iostat=io)
        close(145)

        parnames(1)  = parnames_1
        parnames(2)  = parnames_2
        parnames(3)  = parnames_3
        parnames(4)  = parnames_4
        parnames(5)  = parnames_5 
        parnames(6)  = parnames_6
        parnames(7)  = parnames_7
        parnames(8)  = parnames_8
        parnames(9)  = parnames_9
        parnames(10) = parnames_10
        parnames(11)  = parnames_11
        parnames(12)  = parnames_12
        parnames(13)  = parnames_13
        parnames(14)  = parnames_14
        parnames(15)  = parnames_15 
        parnames(16)  = parnames_16
        parnames(17)  = parnames_17
        parnames(18)  = parnames_18
        parnames(19)  = parnames_19
        parnames(20)  = parnames_20
        parnames(21)  = parnames_21
        parnames(22)  = parnames_22
        parnames(23)  = parnames_23
        parnames(24)  = parnames_24
        parnames(25)  = parnames_25 
        parnames(26)  = parnames_26
        parnames(27)  = parnames_27
        parnames(28)  = parnames_28
        parnames(29)  = parnames_29
        parnames(30)  = parnames_30
        parnames(31)  = parnames_31
        parnames(32)  = parnames_32
        parnames(33)  = parnames_33
        parnames(34)  = parnames_34
        parnames(35)  = parnames_35 
        parnames(36)  = parnames_36
        parnames(37)  = parnames_37
        parnames(38)  = parnames_38
        parnames(39)  = parnames_39
        parnames(40)  = parnames_40
        parnames(41)  = parnames_41
        parnames(42)  = parnames_42
        parnames(43)  = parnames_43
        parnames(44)  = parnames_44
        parnames(45)  = parnames_45 
        parnames(46)  = parnames_46
        parnames(47)  = parnames_47
        parnames(48)  = parnames_48
        parnames(49)  = parnames_49
        parnames(50)  = parnames_50        
        parnames(51)  = parnames_51
        parnames(52)  = parnames_52
        parnames(53)  = parnames_53
        parnames(54)  = parnames_54
        parnames(55)  = parnames_55 
        parnames(56)  = parnames_56
        parnames(57)  = parnames_57
        parnames(58)  = parnames_58
        parnames(59)  = parnames_59
        parnames(60)  = parnames_60

        ! give the filepath to each variable
        vars4MCMC%ANPP_Shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ANPP_Shrub_y))
        vars4MCMC%ANPP_Tree_y%filepath  = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ANPP_Tree_y))
        vars4MCMC%NPP_sphag_y%filepath  = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_NPP_sphag_y))
        vars4MCMC%BNPP_y%filepath       = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_BNPP_y))        ! tree + shrub
        vars4MCMC%er_d%filepath         = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_er_d))          ! shrub + sphag.
        vars4MCMC%er_h%filepath         = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_er_h))          ! shrub + sphag.
        vars4MCMC%gpp_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_gpp_d))         ! Shrub + sphag.
        vars4MCMC%nee_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_nee_d))         ! Shrub + sphag.
        vars4MCMC%nee_h%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_nee_h))         ! shrub + sphag.
        vars4MCMC%LAI_d%filepath        = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_LAI_d))         ! tree  + Shrub
        !
        vars4MCMC%leaf_mass_shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_mass_shrub_y))
        vars4MCMC%stem_mass_shrub_y%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_stem_mass_shrub_y))
        vars4MCMC%leaf_resp_shrub_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_resp_shrub_d))
        vars4MCMC%leaf_resp_tree_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_leaf_resp_tree_d)) 
        ! methane
        vars4MCMC%ch4_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ch4_d))
        vars4MCMC%ch4_h%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_ch4_h))
        ! 
        vars4MCMC%CN_shag_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_CN_shag_d))
        vars4MCMC%photo_shrub_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_photo_shrub_d)) 
        vars4MCMC%photo_tree_d%filepath = adjustl(trim(inDir))//"/"//adjustl(trim(obsfile_photo_tree_d))
        
        return
    end subroutine readConfsNml

    subroutine readParamNml(param_nml_file, in_params, init_params, arr_parval, arr_parmin, arr_parmax)
    ! default nml file name of "TECO_MCMC_configs.nml"
        implicit none
        character(*), intent(in) :: param_nml_file
        type(nml_params_data_type), intent(inout)    :: in_params
        type(nml_initValue_data_type), intent(inout) :: init_params
        real, intent(inout) :: arr_parval(:), arr_parmin(:), arr_parmax(:)
        integer io

        ! ! site special variables that are read from namelist file
        ! type spec_data_type
        real :: lat, lon
        real :: wsmax, wsmin
        real :: LAImin, LAImax
        real :: SLAx, rdepth
        real :: Rootmax, Stemmax
        real :: SapR, SapS
        real :: GLmax, GRmax, Gsmax
        real :: stom_n, a1, Ds0                     ! a1 recalculated in vegn model
        real :: Vcmax0                              ! Jian: Vcmax0 and Vcmx0 is same? Vcmax0 is Vcmx0 in consts
        real :: extkU, xfang, alpha               
        real :: Tau_Leaf, Tau_Wood, Tau_Root        ! turnover rate of plant carbon pools : leaf, wood, root  
        real :: Tau_F, Tau_C                        ! turnover rate of litter carbon pools: fine, coarse 
        real :: Tau_Micro, Tau_slowSOM, Tau_Passive ! turnover rate of soil carbon pools  : fast, slow, passive 
        real :: gddonset
        real :: Q10, Q10rh                          ! Q10rh modified from Ma et al.,2023 for aclimate study, change in transfer module of Q10h
        real :: Rl0, Rs0, Rr0
        ! added for parameters in methane module   
        real :: r_me, Q10pro
        real :: kCH4, Omax
        real :: CH4_thre
        real :: Tveg, Toxi 
        real :: Tpro_me
        ! add based on Ma et al., 2022
        real :: f, bubprob, Vmaxfraction  
        ! add based on Ma et al., 2023
        real :: JV, Entrpy              ! J/mol/K (entropy term, for Jmax & Vcmax)
        real :: etaL, etaW, etaR        ! etaL and etaR are not used. ! the percentage of fine litter of the litters from plant parts
        real :: f_F2M, f_C2M, f_C2S
        real :: f_M2S, f_M2P, f_S2P
        real :: f_S2M, f_P2M
        real :: hmax, hl0, LAIMAX0, la0
        ! ----------------------------------------------------
        namelist/nml_params/ lat, lon, wsmax, wsmin, LAImax, LAImin, rdepth,    &
            Rootmax, Stemmax, SapR, SapS, SLAx, GLmax, GRmax, Gsmax, stom_n,    &
            a1, Ds0, Vcmax0, extkU, xfang, alpha, Tau_Leaf, Tau_Wood, Tau_Root, &
            Tau_F, Tau_C, Tau_Micro, Tau_SlowSOM, Tau_Passive, gddonset, Q10,   &
            Q10rh, Rl0, Rs0, Rr0, r_me, Q10pro, kCH4, Omax, CH4_thre, Tveg,     &
            Tpro_me, Toxi, f, bubprob, Vmaxfraction, JV, Entrpy, etaL, etaW,    &
            etaR, f_F2M, f_C2M, f_C2S, f_M2S, f_M2P, f_S2P, f_S2M, f_P2M, hmax, hl0, LAIMAX0, la0 
        ! ------------------------------------------------------------------------
        ! parameters that are needed to be initilized.
        real :: QC(8), CN0(8)               ! leaf,wood,root,fine lit.,coarse lit.,Micr,Slow,Pass
        real :: NSCmin, Storage, nsc
        real :: accumulation, SNvcmax 
        real :: N_deposit, alphaN, NSN
        real :: QNminer, N_deficit
        real :: thksl(10)                   ! thickness of every soil layer
        real :: FRLEN(10)                   ! ratio of roots in every layer, Oak Ridge FACE: Shuang
        real :: liq_water(10)               ! unit m
        real :: fwsoil, topfws, omega       ! update in soilwater module
        real :: zwt, infilt
        real :: sftmp, Tsnow, Twater
        real :: Tice, G, snow_dsim
        real :: dcount, dcount_soil
        real :: ice_tw   
        real :: Tsoill(10), ice(10)
        real :: shcap_snow                  ! tuneice worker better
        real :: condu_snow, condu_b         ! yuanyuan soil thermal version value  ... int: this par is not sensitive to CWE
        real :: depth_ex, diff_s, diff_snow ! int diffusivity of snow not sensitive for ice
        real :: albedo_snow, resht     
        real :: thd_snow_depth, b_bound     ! tuneice  not sensitive for ice, b_bound no used
        real :: infilt_rate, fa, fsub
        real :: rho_snow, decay_m           ! aging factor on snow melting
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        real :: CH4_V(10), CH4(10), Vp(10)  ! assume in the very beginning no bubbles exist in the first three layers (30cm)
        real :: bubble_methane_tot, Nbub
        real :: depth_1                     ! calculate soil depth unit cm
        ! -----------------------------------------------------------------
        namelist/nml_initial_values/ QC, CN0, NSCmin, Storage, nsc, accumulation, SNvcmax, &
            N_deposit, alphaN, NSN, QNminer, N_deficit, thksl, FRLEN, liq_water, fwsoil,   &
            topfws, omega, zwt, infilt, sftmp, Tsnow, Twater, Tice, G, snow_dsim, dcount,  &
            dcount_soil, ice_tw, Tsoill, ice, shcap_snow, condu_snow, condu_b, depth_ex,   & 
            diff_s, diff_snow, albedo_snow, resht, thd_snow_depth, b_bound, infilt_rate,   &
            fa, fsub, rho_snow, decay_m, CH4_V, CH4, Vp, bubble_methane_tot, Nbub, depth_1
        ! -----------------------------------------------------------------------------------

        ! -------------------------------------------------------------

        namelist /nml_parval/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin,            &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax,           &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha,             &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M
        namelist /nml_parmin/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin, &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax, &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha, &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M
        namelist /nml_parmax/ mc_lat, mc_Longitude, mc_wsmax, mc_wsmin, &                                                    
                mc_LAIMAX, mc_LAIMIN, mc_rdepth, mc_Rootmax, mc_Stemmax, &                        
                mc_SapR, mc_SapS, mc_SLA, mc_GLmax, mc_GRmax, mc_Gsmax, mc_stom_n, &            
                mc_a1, mc_Ds0, mc_Vcmx0, mc_extkU, mc_xfang, mc_alpha, &                    
                mc_Tau_Leaf, mc_Tau_Wood, mc_Tau_Root, mc_Tau_F, &
                mc_Tau_C,  mc_Tau_Micro, mc_Tau_SlowSOM, mc_Tau_Passive, &            
                mc_gddonset, mc_Q10, mc_Rl0, mc_Rs0, mc_Rr0, &        
                mc_r_me, mc_Q10pro, mc_kCH4, mc_Omax, mc_CH4_thre, &
                mc_Tveg, mc_Tpro_me, mc_Toxi, &
                mc_f, mc_bubprob, mc_Vmaxfraction, &                                    
                mc_Q10rh, mc_JV, mc_Entrpy, &                                
                mc_etaL, mc_etaW, mc_etaR, &
                mc_f_F2M, mc_f_C2M, mc_f_C2S, mc_f_M2S, &
                mc_f_M2P, mc_f_S2P, mc_f_S2M, mc_f_P2M

        
        print *, "# read parameters nml file: ", param_nml_file
        open(343, file = param_nml_file)
        read(343, nml  = nml_params,         iostat=io)
        read(343, nml  = nml_initial_values, iostat=io)
        read(343, nml  = nml_parval,         iostat=io)
        call giveValues2par(arr_parval)
        read(343, nml=nml_parmin,         iostat=io)
        call giveValues2par(arr_parmin)
        read(343, nml=nml_parmax,         iostat=io)
        call giveValues2par(arr_parmax)
        close(343)
        ! update the parameters in in_params and init_params
        in_params%lat         = lat
        in_params%lon         = lon
        in_params%wsmax       = wsmax
        in_params%wsmin       = wsmin
        in_params%LAImin      = LAImin
        in_params%LAImax      = LAImax
        in_params%SLAx        = SLAx
        in_params%rdepth      = rdepth
        in_params%Rootmax     = Rootmax
        in_params%Stemmax     = Stemmax
        in_params%SapR        = SapR
        in_params%SapS        = SapS
        in_params%GLmax       = GLmax
        in_params%GRmax       = GRmax
        in_params%Gsmax       = Gsmax
        in_params%stom_n      = stom_n
        in_params%a1          = a1
        in_params%Ds0         = Ds0
        in_params%Vcmax0      = Vcmax0
        in_params%extkU       = extkU
        in_params%xfang       = xfang
        in_params%alpha       = alpha            
        in_params%Tau_Leaf    = Tau_Leaf
        in_params%Tau_Wood    = Tau_Wood
        in_params%Tau_Root    = Tau_Root 
        in_params%Tau_F       = Tau_F
        in_params%Tau_C       = Tau_C  
        in_params%Tau_Micro   = Tau_Micro
        in_params%Tau_slowSOM = Tau_slowSOM
        in_params%Tau_Passive = Tau_Passive
        in_params%gddonset    = gddonset
        in_params%Q10         = Q10
        in_params%Q10rh       = Q10rh    
        in_params%Rl0         = Rl0
        in_params%Rs0         = Rs0
        in_params%Rr0         = Rr0
        ! added for parameters in methane module   
        in_params%r_me        = r_me
        in_params%Q10pro      = Q10pro
        in_params%kCH4        = kCH4
        in_params%Omax        = Omax
        in_params%CH4_thre    = CH4_thre
        in_params%Tveg        = Tveg
        in_params%Toxi        = Toxi
        in_params%Tpro_me     = Tpro_me
        ! add based on Ma et al., 2022
        in_params%f            = f
        in_params%bubprob      = bubprob
        in_params%Vmaxfraction = Vmaxfraction
        ! add based on Ma et al., 2023
        in_params%JV     = JV
        in_params%Entrpy = Entrpy 
        in_params%etaL   = etaL
        in_params%etaW   = etaW
        in_params%etaR   = etaR
        in_params%f_F2M  = f_F2M
        in_params%f_C2M  = f_C2M
        in_params%f_C2S  = f_C2S
        in_params%f_M2S  = f_M2S
        in_params%f_M2P  = f_M2P
        in_params%f_S2P  = f_S2P
        in_params%f_S2M  = f_S2M
        in_params%f_P2M  = f_P2M
        in_params%hmax   = hmax 
        in_params%hl0    = hl0
        in_params%LAIMAX0 = LAIMAX0
        in_params%la0     = la0
        ! =====================================
        init_params%QC             = QC
        init_params%CN0            = CN0            
        init_params%NSCmin         = NSCmin
        init_params%Storage        = Storage
        init_params%nsc            = nsc
        init_params%accumulation   = accumulation
        init_params%SNvcmax        = SNvcmax
        init_params%N_deposit      = N_deposit
        init_params%alphaN         = alphaN
        init_params%NSN            = NSN
        init_params%QNminer        = QNminer
        init_params%N_deficit      = N_deficit
        init_params%thksl          = thksl 
        init_params%FRLEN          = FRLEN
        init_params%liq_water      = liq_water
        init_params%fwsoil         = fwsoil
        init_params%topfws         = topfws
        init_params%omega          = omega
        init_params%zwt            = zwt
        init_params%infilt         = infilt
        init_params%sftmp          = sftmp
        init_params%Tsnow          = Tsnow
        init_params%Twater         = Twater
        init_params%Tice           = Tice
        init_params%G              = G
        init_params%snow_dsim      = snow_dsim
        init_params%dcount         = dcount
        init_params%dcount_soil    = dcount_soil
        init_params%ice_tw         = ice_tw
        init_params%Tsoill         = Tsoill
        init_params%ice            = ice
        init_params%shcap_snow     = shcap_snow
        init_params%condu_snow     = condu_snow
        init_params%condu_b        = condu_b
        init_params%depth_ex       = depth_ex
        init_params%diff_s         = diff_s
        init_params%diff_snow      = diff_snow
        init_params%albedo_snow    = albedo_snow
        init_params%resht          = resht
        init_params%thd_snow_depth = thd_snow_depth
        init_params%b_bound        = b_bound
        init_params%infilt_rate    = infilt_rate
        init_params%fa             = fa
        init_params%fsub           = fsub
        init_params%rho_snow       = rho_snow
        init_params%decay_m        = decay_m
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        init_params%CH4_V              = CH4_V
        init_params%CH4                = CH4
        init_params%Vp                 = Vp
        init_params%bubble_methane_tot = bubble_methane_tot
        init_params%Nbub               = Nbub
        init_params%depth_1            = depth_1
    end subroutine readParamNml

    subroutine readObsData()
        implicit none
        call readObsData_var(vars4MCMC%ANPP_Shrub_y)
        call readObsData_var(vars4MCMC%ANPP_Tree_y)
        call readObsData_var(vars4MCMC%NPP_sphag_y)
        call readObsData_var(vars4MCMC%BNPP_y)        ! tree + shrub
        call readObsData_var(vars4MCMC%er_d)          ! shrub + sphag.
        call readObsData_var(vars4MCMC%er_h)          ! shrub + sphag.
        call readObsData_var(vars4MCMC%gpp_d)         ! Shrub + sphag.
        call readObsData_var(vars4MCMC%nee_d)         ! Shrub + sphag.
        call readObsData_var(vars4MCMC%nee_h)         ! shrub + sphag.
        call readObsData_var(vars4MCMC%LAI_d)         ! tree  + Shrub
        !
        call readObsData_var(vars4MCMC%leaf_mass_shrub_y)
        call readObsData_var(vars4MCMC%stem_mass_shrub_y)
        call readObsData_var(vars4MCMC%leaf_resp_shrub_d) 
        call readObsData_var(vars4MCMC%leaf_resp_tree_d) 
        ! methane
        call readObsData_var(vars4MCMC%ch4_d) 
        call readObsData_var(vars4MCMC%ch4_h) 
        ! 
        call readObsData_var(vars4MCMC%CN_shag_d) 
        call readObsData_var(vars4MCMC%photo_shrub_d) 
        call readObsData_var(vars4MCMC%photo_tree_d) 
    end subroutine readObsData

    subroutine readObsData_var(var_obsData)
        implicit none
        type(interCostVariable), intent(inout) :: var_obsData
        logical toExistOrNot
        integer toCountLines

        INQUIRE(FILE=var_obsData%filepath, EXIST=toExistOrNot)
        var_obsData%existOrNot = toExistOrNot
        if (var_obsData%existOrNot) then
            call ReadLineNumFromFile(var_obsData%filepath, toCountLines)
            allocate(var_obsData%obsData(toCountLines, 5))
            call ReadObsDataFromFile(var_obsData%filepath, toCountLines, var_obsData%obsData)
            allocate(var_obsData%mdData(toCountLines, 4))
        endif
        return
    end subroutine readObsData_var

    subroutine renewMDpars(parval, re_in_params)
        implicit none
        real, intent(in) :: parval(:)
        type(nml_params_data_type), intent(inout) :: re_in_params

        re_in_params%lat         = parval(1)
        re_in_params%lon         = parval(2)
        re_in_params%wsmax       = parval(3)
        re_in_params%wsmin       = parval(4)                                            
        re_in_params%LAIMAX      = parval(5)
        re_in_params%LAIMIN      = parval(6)
        re_in_params%rdepth      = parval(7)
        re_in_params%Rootmax     = parval(8)
        re_in_params%Stemmax     = parval(9)                                    
        re_in_params%SapR        = parval(10)
        re_in_params%SapS        = parval(11)
        re_in_params%SLAx        = parval(12)
        re_in_params%GLmax       = parval(13)
        re_in_params%GRmax       = parval(14)
        re_in_params%Gsmax       = parval(15)
        re_in_params%stom_n      = parval(16)         
        re_in_params%a1          = parval(17)
        re_in_params%Ds0         = parval(18)
        re_in_params%Vcmax0      = parval(19)
        re_in_params%extkU       = parval(20)
        re_in_params%xfang       = parval(21)
        re_in_params%alpha       = parval(22)    
        re_in_params%Tau_Leaf    = parval(23)
        re_in_params%Tau_Wood    = parval(24)
        re_in_params%Tau_Root    = parval(25)
        re_in_params%Tau_F       = parval(26)
        re_in_params%Tau_C       = parval(27)
        re_in_params%Tau_Micro   = parval(28)
        re_in_params%Tau_SlowSOM = parval(29)
        re_in_params%Tau_Passive = parval(30)    
        re_in_params%gddonset    = parval(31)
        re_in_params%Q10         = parval(32)
        re_in_params%Rl0         = parval(33)     
        re_in_params%Rs0         = parval(34)    
        re_in_params%Rr0         = parval(35)                    
        re_in_params%r_me        = parval(36)
        re_in_params%Q10pro      = parval(37)
        re_in_params%kCH4        = parval(38)
        re_in_params%Omax         = parval(39)
        re_in_params%CH4_thre     = parval(40)
        re_in_params%Tveg         = parval(41)
        re_in_params%Tpro_me      = parval(42)
        re_in_params%Toxi         = parval(43)        
        re_in_params%f            = parval(44)
        re_in_params%bubprob      = parval(45)
        re_in_params%Vmaxfraction = parval(46)                                    
        re_in_params%Q10rh        = parval(47)
        re_in_params%JV           = parval(48)
        re_in_params%Entrpy       = parval(49)                    
        re_in_params%etaL         = parval(50)
        re_in_params%etaW         = parval(51)
        re_in_params%etaR         = parval(52)
        re_in_params%f_F2M        = parval(53)
        re_in_params%f_C2M        = parval(54)
        re_in_params%f_C2S        = parval(55)
        re_in_params%f_M2S        = parval(56)
        re_in_params%f_M2P        = parval(57)
        re_in_params%f_S2P        = parval(58)
        re_in_params%f_S2M        = parval(59)
        re_in_params%f_P2M        = parval(60)
        return
    end subroutine renewMDpars


    ! subroutine giveValues2var(filepath, existOrNot, data)
    !     implicit none
    !     character(500) filepath
    !     logical existOrNot
    !     real, allocatable :: data(:, :)
    !     integer count_lines

    !     INQUIRE(FILE=filepath, EXIST=existOrNot)
    !     if(existOrNot)then
    !         call ReadLineNumFromFile(filepath, count_lines)
    !         allocate(data(count_lines, 5))
    !         call ReadObsDataFromFile(filepath, count_lines, data)
    !     end if
    !     return
    ! end subroutine giveValues2var

    subroutine giveValues2par(arr_par)
        implicit none
        real, intent(inout) :: arr_par(:)

        arr_par(1)  = mc_lat
        arr_par(2)  = mc_Longitude 
        arr_par(3)  = mc_wsmax 
        arr_par(4)  = mc_wsmin                                                      
        arr_par(5)  = mc_LAIMAX
        arr_par(6)  = mc_LAIMIN    
        arr_par(7)  = mc_rdepth    
        arr_par(8)  = mc_Rootmax    
        arr_par(9)  = mc_Stemmax                                            
        arr_par(10) = mc_SapR    
        arr_par(11) = mc_SapS     
        arr_par(12) = mc_SLA        
        arr_par(13) = mc_GLmax    
        arr_par(14) = mc_GRmax    
        arr_par(15) = mc_Gsmax    
        arr_par(16) = mc_stom_n                                            
        arr_par(17) = mc_a1       
        arr_par(18) = mc_Ds0        
        arr_par(19) = mc_Vcmx0    
        arr_par(20) = mc_extkU    
        arr_par(21) = mc_xfang    
        arr_par(22) = mc_alpha                         
        arr_par(23) = mc_Tau_Leaf   
        arr_par(24) = mc_Tau_Wood   
        arr_par(25) = mc_Tau_Root   
        arr_par(26) = mc_Tau_F       
        arr_par(27) = mc_Tau_C       
        arr_par(28) = mc_Tau_Micro   
        arr_par(29) = mc_Tau_SlowSOM 
        arr_par(30) = mc_Tau_Passive                             
        arr_par(31) = mc_gddonset    
        arr_par(32) = mc_Q10         
        arr_par(33) = mc_Rl0        
        arr_par(34) = mc_Rs0        
        arr_par(35) = mc_Rr0                            
        arr_par(36) = mc_r_me   
        arr_par(37) = mc_Q10pro   
        arr_par(38) = mc_kCH4    
        arr_par(39) = mc_Omax   
        arr_par(40) = mc_CH4_thre 
        arr_par(41) = mc_Tveg  
        arr_par(42) = mc_Tpro_me 
        arr_par(43) = mc_Toxi               
        arr_par(44) = mc_f    
        arr_par(45) = mc_bubprob  
        arr_par(46) = mc_Vmaxfraction                                        
        arr_par(47) = mc_Q10rh  
        arr_par(48) = mc_JV   
        arr_par(49) = mc_Entrpy                                
        arr_par(50) = mc_etaL   
        arr_par(51) = mc_etaW  
        arr_par(52) = mc_etaR   
        arr_par(53) = mc_f_F2M   
        arr_par(54) = mc_f_C2M  
        arr_par(55) = mc_f_C2S 
        arr_par(56) = mc_f_M2S  
        arr_par(57) = mc_f_M2P 
        arr_par(58) = mc_f_S2P  
        arr_par(59) = mc_f_S2M  
        arr_par(60) = mc_f_P2M 

    end subroutine giveValues2par

    subroutine GetSimuData(get_iyear, get_iday, get_ihour, in_vegn, nHr, nDay, nMon, nYr)
        implicit none
        type(vegn_tile_type), intent(in) :: in_vegn
        integer, intent(in) :: nHr, nDay, nMon, nYr

        integer get_iyear, get_iday, get_ihour
        integer i, ipft, npft
        ! vars4MCMC%
        mc_iyear = get_iyear
        mc_iday  = get_iday
        mc_ihour = get_ihour + 1

        call GetSimuData_var(vars4MCMC%ANPP_Shrub_y,(outVars_y%allSpec(2)%nppLeaf + outVars_y%allSpec(2)%nppStem)*24*365)
        call GetSimuData_var(vars4MCMC%ANPP_Tree_y, (outVars_y%allSpec(1)%nppLeaf + outVars_y%allSpec(1)%nppStem)*24*365)
        call GetSimuData_var(vars4MCMC%NPP_sphag_y,  outVars_y%allSpec(3)%npp*24*365)
        call GetSimuData_var(vars4MCMC%BNPP_y,       (outVars_y%allSpec(1)%nppRoot + outVars_y%allSpec(2)%nppRoot)*24*365)        ! tree + shrub
        call GetSimuData_var(vars4MCMC%er_d,         (outVars_d%allSpec(2)%ra + outVars_d%allSpec(3)%ra + outVars_d%rh)*24)          ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%er_h,         (outVars_h%allSpec(2)%ra + outVars_h%allSpec(3)%ra + outVars_h%rh))          ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%gpp_d,        (outVars_d%allSpec(2)%gpp + outVars_d%allSpec(3)%gpp)*24)         ! Shrub + sphag.
        call GetSimuData_var(vars4MCMC%nee_d,        (outVars_d%allSpec(2)%gpp + outVars_d%allSpec(3)%gpp - &
                                                    (outVars_d%allSpec(2)%ra  + outVars_d%allSpec(3)%ra  + outVars_d%rh))*24)         ! Shrub + sphag.
        call GetSimuData_var(vars4MCMC%nee_h,        (outVars_h%allSpec(2)%gpp + outVars_h%allSpec(3)%gpp - &
                                                    (outVars_h%allSpec(2)%ra  + outVars_h%allSpec(3)%ra  + outVars_h%rh))*24)         ! shrub + sphag.
        call GetSimuData_var(vars4MCMC%LAI_d,       (outVars_d%allSpec(1)%LAI + outVars_d%allSpec(2)%LAI)/2)         ! tree  + Shrub
        !
        call GetSimuData_var(vars4MCMC%leaf_mass_shrub_y, outVars_y%allSpec(2)%cleaf*0.48)
        call GetSimuData_var(vars4MCMC%stem_mass_shrub_y, outVars_y%allSpec(2)%cStem*0.48)
        call GetSimuData_var(vars4MCMC%leaf_resp_shrub_d, outVars_d%allSpec(2)%ra*24) 
        call GetSimuData_var(vars4MCMC%leaf_resp_tree_d,  outVars_d%allSpec(1)%ra*24) 
        ! methane
        call GetSimuData_var(vars4MCMC%ch4_d, sum(outVars_d%ch4)*24) 
        call GetSimuData_var(vars4MCMC%ch4_h, sum(outVars_h%ch4)) 
        ! 
        call GetSimuData_var(vars4MCMC%CN_shag_d,    ((outVars_d%allSpec(3)%cleaf + outVars_d%allSpec(3)%cStem + &
                                                      outVars_d%allSpec(3)%cRoot)/(outVars_d%allSpec(3)%nLeaf + &
                                                      outVars_d%allSpec(3)%nStem + outVars_d%allSpec(3)%nRoot))) 
        call GetSimuData_var(vars4MCMC%photo_shrub_d, outVars_d%allSpec(2)%gpp*24) 
        call GetSimuData_var(vars4MCMC%photo_tree_d,  outVars_d%allSpec(1)%gpp*24) 
        ! ----------------------------------------------------------------------------
        call update_mcmc_tot_outputs(tot_paramsets_outs_d, outVars_d, nDay)
    end subroutine GetSimuData

    subroutine GetSimuData_var(var_obsData, var_mdData)
        implicit none
        type(interCostVariable), intent(inout) :: var_obsData
        real, intent(in) :: var_mdData

        if(var_obsData%existOrNot)then  ! if the observation file is existed
            if(var_obsData%mc_itime <= size(var_obsData%obsData, dim=1))then  ! still have observation data not being matched
                do while(var_obsData%obsData(var_obsData%mc_itime, 1) .lt. forcing(1)%year) ! some observation is beyond the range of simulation
                    var_obsData%mdData(var_obsData%mc_itime, 4) = -9999
                    var_obsData%mc_itime = var_obsData%mc_itime + 1
                enddo

                if(var_obsData%obsData(var_obsData%mc_itime, 1) .eq. mc_iyear .and. &
                   var_obsData%obsData(var_obsData%mc_itime, 2) .eq. mc_iday  .and. &
                   var_obsData%obsData(var_obsData%mc_itime, 3) .eq. mc_ihour) then
                        var_obsData%mdData(var_obsData%mc_itime, 1) = mc_iyear
                        var_obsData%mdData(var_obsData%mc_itime, 2) = mc_iday
                        var_obsData%mdData(var_obsData%mc_itime, 3) = mc_ihour
                        var_obsData%mdData(var_obsData%mc_itime, 4) = var_mdData
                        var_obsData%mc_itime = var_obsData%mc_itime + 1
                endif
            endif
        endif
    end subroutine GetSimuData_var

    subroutine update_mcmc_tot_outputs(tot_mcmc_outputs, simu_outputs, inTime)
        implicit none
        type(mcmc_outVars_type), intent(inout) :: tot_mcmc_outputs 
        type(outvars_data_type), intent(in)    :: simu_outputs
        integer, intent(in) :: inTime
        integer :: ipft

        do ipft = 1, count_pft
            tot_mcmc_outputs%allSpec(ipft)%gpp(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%gpp
            tot_mcmc_outputs%allSpec(ipft)%nee(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nee
            tot_mcmc_outputs%allSpec(ipft)%npp(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%npp
            tot_mcmc_outputs%allSpec(ipft)%nppLeaf(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nppLeaf
            tot_mcmc_outputs%allSpec(ipft)%nppWood(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nppWood
            tot_mcmc_outputs%allSpec(ipft)%nppStem(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nppStem
            tot_mcmc_outputs%allSpec(ipft)%nppRoot(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nppRoot
            tot_mcmc_outputs%allSpec(ipft)%nppOther(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nppOther    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
            tot_mcmc_outputs%allSpec(ipft)%ra(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%ra
            tot_mcmc_outputs%allSpec(ipft)%raLeaf(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%raLeaf
            tot_mcmc_outputs%allSpec(ipft)%raStem(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%raStem
            tot_mcmc_outputs%allSpec(ipft)%raRoot(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%raRoot
            tot_mcmc_outputs%allSpec(ipft)%raOther(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%raOther
            tot_mcmc_outputs%allSpec(ipft)%rMaint(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%rMaint
            tot_mcmc_outputs%allSpec(ipft)%rGrowth(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%rGrowth
            tot_mcmc_outputs%allSpec(ipft)%nbp(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nbp
            ! Carbon Pools  (KgC m-2)
            tot_mcmc_outputs%allSpec(ipft)%cLeaf(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%cLeaf
            tot_mcmc_outputs%allSpec(ipft)%cStem(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%cStem
            tot_mcmc_outputs%allSpec(ipft)%cRoot(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%cRoot
            ! Nitrogen pools (kgN m-2)
            tot_mcmc_outputs%allSpec(ipft)%nLeaf(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nLeaf
            tot_mcmc_outputs%allSpec(ipft)%nStem(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nStem
            tot_mcmc_outputs%allSpec(ipft)%nRoot(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%nRoot
            ! tot_mcmc_outputs%allSpec(ipft)%nOther(:)
            ! water fluxes (kg m-2 s-1)
            tot_mcmc_outputs%allSpec(ipft)%tran(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%tran
            ! other
            tot_mcmc_outputs%allSpec(ipft)%lai(iDAsimu, inTime)   = simu_outputs%allSpec(ipft)%lai
        enddo

        tot_mcmc_outputs%gpp(iDAsimu, inTime)            = simu_outputs%gpp    
        tot_mcmc_outputs%nee(iDAsimu, inTime)            = simu_outputs%nee
        tot_mcmc_outputs%npp(iDAsimu, inTime)            = simu_outputs%npp
        tot_mcmc_outputs%nppLeaf(iDAsimu, inTime)        = simu_outputs%nppLeaf
        tot_mcmc_outputs%nppWood(iDAsimu, inTime)        = simu_outputs%nppWood
        tot_mcmc_outputs%nppStem(iDAsimu, inTime)        = simu_outputs%nppStem
        tot_mcmc_outputs%nppRoot(iDAsimu, inTime)        = simu_outputs%nppRoot
        tot_mcmc_outputs%nppOther(iDAsimu, inTime)       = simu_outputs%nppOther
        tot_mcmc_outputs%ra(iDAsimu, inTime)             = simu_outputs%ra
        tot_mcmc_outputs%raLeaf(iDAsimu, inTime)         = simu_outputs%raLeaf
        tot_mcmc_outputs%raStem(iDAsimu, inTime)         = simu_outputs%raStem
        tot_mcmc_outputs%raRoot(iDAsimu, inTime)         = simu_outputs%raRoot
        tot_mcmc_outputs%raOther(iDAsimu, inTime)        = simu_outputs%raOther
        tot_mcmc_outputs%rMaint(iDAsimu, inTime)         = simu_outputs%rMaint
        tot_mcmc_outputs%rGrowth(iDAsimu, inTime)        = simu_outputs%rGrowth
        tot_mcmc_outputs%rh(iDAsimu, inTime)             = simu_outputs%rh
        tot_mcmc_outputs%nbp(iDAsimu, inTime)            = simu_outputs%nbp
        tot_mcmc_outputs%wetlandCH4(iDAsimu, inTime)     = simu_outputs%wetlandCH4
        tot_mcmc_outputs%wetlandCH4prod(iDAsimu, inTime) = simu_outputs%wetlandCH4prod
        tot_mcmc_outputs%wetlandCH4cons(iDAsimu, inTime) = simu_outputs%wetlandCH4cons
        ! Carbon Pools  (KgC m-2)
        tot_mcmc_outputs%cLeaf(iDAsimu, inTime)          = simu_outputs%cLeaf
        tot_mcmc_outputs%cStem(iDAsimu, inTime)          = simu_outputs%cStem
        tot_mcmc_outputs%cRoot(iDAsimu, inTime)          = simu_outputs%cRoot
        tot_mcmc_outputs%cOther(iDAsimu, inTime)         = simu_outputs%cOther
        tot_mcmc_outputs%cLitter(iDAsimu, inTime)        = simu_outputs%cLitter
        tot_mcmc_outputs%cLitterCwd(iDAsimu, inTime)     = simu_outputs%cLitterCwd
        tot_mcmc_outputs%cSoil(iDAsimu, inTime)          = simu_outputs%cSoil
        tot_mcmc_outputs%cSoilLevels(iDAsimu, inTime, :) = simu_outputs%cSoilLevels
        tot_mcmc_outputs%cSoilFast(iDAsimu, inTime)      = simu_outputs%cSoilFast
        tot_mcmc_outputs%cSoilSlow(iDAsimu, inTime)      = simu_outputs%cSoilSlow
        tot_mcmc_outputs%cSoilPassive(iDAsimu, inTime)   = simu_outputs%cSoilPassive
        tot_mcmc_outputs%CH4(iDAsimu, inTime, :)         = simu_outputs%CH4
        ! Nitrogen fluxes (kgN m-2 s-1)
        tot_mcmc_outputs%fBNF(iDAsimu, inTime)           = simu_outputs%fBNF
        tot_mcmc_outputs%fN2O(iDAsimu, inTime)           = simu_outputs%fN2O
        tot_mcmc_outputs%fNloss(iDAsimu, inTime)         = simu_outputs%fNloss
        tot_mcmc_outputs%fNnetmin(iDAsimu, inTime)       = simu_outputs%fNnetmin
        tot_mcmc_outputs%fNdep(iDAsimu, inTime)          = simu_outputs% fNdep
        ! Nitrogen pools (kgN m-2)
        tot_mcmc_outputs%nLeaf(iDAsimu, inTime)          = simu_outputs%nLeaf
        tot_mcmc_outputs%nStem(iDAsimu, inTime)          = simu_outputs%nStem
        tot_mcmc_outputs%nRoot(iDAsimu, inTime)          = simu_outputs%nRoot
        tot_mcmc_outputs%nOther(iDAsimu, inTime)         = simu_outputs%nOther
        tot_mcmc_outputs%nLitter(iDAsimu, inTime)        = simu_outputs%nLitter
        tot_mcmc_outputs%nLitterCwd(iDAsimu, inTime)     = simu_outputs%nLitterCwd
        tot_mcmc_outputs%nSoil(iDAsimu, inTime)          = simu_outputs%nSoil
        tot_mcmc_outputs%nMineral(iDAsimu, inTime)       = simu_outputs%nMineral
        ! energy fluxes (W m-2)
        tot_mcmc_outputs%hfls(iDAsimu, inTime)           = simu_outputs%hfls
        tot_mcmc_outputs%hfss(iDAsimu, inTime)           = simu_outputs%hfss
        tot_mcmc_outputs%SWnet(iDAsimu, inTime)          = simu_outputs%SWnet
        tot_mcmc_outputs%LWnet(iDAsimu, inTime)          = simu_outputs%LWnet
        ! water fluxes (kg m-2 s-1)
        tot_mcmc_outputs%ec(iDAsimu, inTime)             = simu_outputs%ec
        tot_mcmc_outputs%tran(iDAsimu, inTime)           = simu_outputs%tran
        tot_mcmc_outputs%es(iDAsimu, inTime)             = simu_outputs%es
        tot_mcmc_outputs%hfsbl(iDAsimu, inTime)          = simu_outputs%hfsbl
        tot_mcmc_outputs%mrro(iDAsimu, inTime)           = simu_outputs%mrro
        tot_mcmc_outputs%mrros(iDAsimu, inTime)          = simu_outputs%mrros
        tot_mcmc_outputs%mrrob(iDAsimu, inTime)          = simu_outputs%mrrob
        ! other
        tot_mcmc_outputs%mrso(iDAsimu, inTime, :)        = simu_outputs%mrso
        tot_mcmc_outputs%tsl(iDAsimu, inTime, :)         = simu_outputs%tsl
        tot_mcmc_outputs%tsland(iDAsimu, inTime)         = simu_outputs%tsland
        tot_mcmc_outputs%wtd(iDAsimu, inTime)            = simu_outputs%wtd
        tot_mcmc_outputs%snd(iDAsimu, inTime)            = simu_outputs%snd
        tot_mcmc_outputs%lai(iDAsimu, inTime)            = simu_outputs%lai
        return
    end subroutine update_mcmc_tot_outputs


    ! subroutine ReadLineNumFromFile(filepath, count_lines)
    !     implicit none
    !     character(len=*), intent(in) :: filepath
    !     character(len=100) header, line
    !     integer STAT, count_lines

    !     open(38, file=trim(filepath), status="old", action="read", iostat=STAT) ! open file
    !     read(38, '(a100)') header           ! read the header of the file
    !     count_lines = 0                     ! initilize the count_lines
    !     do while(.TRUE.)
    !         read(38, *, iostat=STAT) line   ! read each line
    !         if(STAT .ne. 0) exit            ! until the end of the file
    !         count_lines = count_lines + 1   ! recording the count of the lines
    !     enddo
    !     return
    ! end subroutine ReadLineNumFromFile

    subroutine ReadObsDataFromFile(filepath, count_lines, resData)
        ! Jian: note that this subroutine is used to read the observational data. 
        ! The observational file must be .txt format, and with 5 columns: year, doy, hour, value, std.
        implicit none
        character(len=*), intent(in) :: filepath
        character(len=100) header
        integer STAT, count_lines, iline, n
        real resData(count_lines, 5), readData(5) ! 5 colunms: year, doy, hour, value, std.

        OPEN(34, FILE=trim(filepath), status='old', ACTION='read', IOSTAT=STAT) ! open file
        read(34, '(a100)') header
        iline = 1
        do
            read(34,*,iostat=STAT, end=567) (readData(n), n = 1, 5)
            if(STAT .ne. 0) exit
            resData(iline, :) = readData
            iline = iline + 1
        end do
567     continue
        close(34)
        return
    end subroutine ReadObsDataFromFile

end module mcmc_mod