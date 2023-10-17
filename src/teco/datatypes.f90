module datatypes
    implicit none
    ! settings
    character(100) :: teco_configfile  ! the TECO config file
    character(50)  :: case_name        ! define the case name
    logical :: do_simu                 ! simulation mode
    logical :: do_mcmc                 ! MCMC for data assimilation mode
    logical :: do_spinup               ! spinup mode
    logical :: do_matrix               ! whether run matrix or not
    logical :: do_restart              ! whether read restart file or not
    ! simulation selections
    logical :: do_snow                 ! do soil snow process or not
    logical :: do_soilphy              ! do soil physics or not
    logical :: do_EBG                  ! run EBG or not based on Ma et al., 2022
    logical :: do_ndep                 ! N deposit
    logical :: do_leap                 ! judge leap year or not
    ! output selections
    logical :: do_out_hr
    logical :: do_out_day
    logical :: do_out_mon
    logical :: do_out_yr
    ! 
    logical :: do_spruce = .False.    ! set the spruce site handle in 1974 for cutting the vegetation
    integer :: dtimes                 ! 24: hourly simulation
    ! set the input and output path and files
    character(200) :: inDir
    character(200) :: outDir
    ! input files
    character(300) :: climfile
    character(300) :: watertablefile
    character(300) :: snowdepthfile
    character(300) :: in_restartfile
    character(300) :: mcmc_configfile
    character(300) :: spinup_configfile
    ! above settings from the nml file --------------------------------------------------
    
    ! output path
    character(250) :: outDir_nc     = "results_nc_format"
    character(250) :: outDir_csv    = "results_csv_format"
    character(250) :: outDir_mcmc   = "results_mcmc"
    character(250) :: outDir_spinup = "results_spinup"

    ! other paths
    character(300) :: outdir_case, outDir_h, outDir_d, outDir_m, outDir_y
    character(300) :: outDir_mcmc_h, outDir_mcmc_d, outDir_mcmc_m, outfile_restart, restartfile

    integer :: count_pft
    character(300) :: file_site_params
    character(300), allocatable :: spec_names(:), files_pft_params(:)

    ! experiment settings
    real :: Ttreat   = 0.        ! Temperature treatment, warming in air and soil temperature
    real :: CO2treat = 0.        ! CO2 treatmant, up to CO2treat, not add to Ca. CO2
    real :: N_fert   = 0.        ! 5.6 ! (11.2 gN m-2 yr-1, in spring, Duke Forest FACE)
    ! ---------------------------------------------------------------------------------------
    integer :: nHours, nDays, nMonths, nYears         ! how many iterations hourly, daily, monthly and yearly
    ! --------------------------
    
    integer, parameter :: nlayers = 10                ! how many layers
    real,    parameter :: pi      = 3.1415926
    ! physical constants
    real,    parameter :: tauL(3) = (/0.1, 0.425, 0.00/)  ! leaf transmittance for vis, for NIR, for thermal
    real,    parameter :: rhoL(3) = (/0.1, 0.425, 0.00/)  ! leaf reflectance for vis, for NIR, for thermal
    real,    parameter :: emleaf  = 0.96
    real,    parameter :: emsoil  = 0.94
    real,    parameter :: Rconst  = 8.314                 ! universal gas constant (J/mol)
    real,    parameter :: sigma   = 5.67e-8               ! Steffan Boltzman constant (W/m2/K4)
    real,    parameter :: cpair   = 1010.                 ! heat capapcity of air (J/kg/K)
    real,    parameter :: Patm    = 101325. !1.e5         ! atmospheric pressure  (Pa)
    real,    parameter :: Trefk   = 293.2                 ! reference temp K for Kc, Ko, Rd
    real,    parameter :: H2OLv0  = 2.501e6               ! latent heat H2O (J/kg)
    real,    parameter :: AirMa   = 29.e-3                ! mol mass air (kg/mol)
    real,    parameter :: H2OMw   = 18.e-3                ! mol mass H2O (kg/mol)
    real,    parameter :: chi     = 0.93                  ! gbH/gbw
    real,    parameter :: Dheat   = 21.5e-6               ! molecular diffusivity for heat
    ! plant parameters
    real,    parameter :: gsw0    = 1.0e-2                ! g0 for H2O in BWB model
    real,    parameter :: theta   = 0.9
    real,    parameter :: wleaf   = 0.01                  ! leaf width (m)
    ! thermodynamic parameters for Kc and Ko (Leuning 1990)
    real,    parameter :: conKc0  = 302.e-6               ! mol mol^-1
    real,    parameter :: conKo0  = 256.e-3               ! mol mol^-1
    real,    parameter :: Ekc     = 59430.                ! J mol^-1
    real,    parameter :: Eko     = 36000.                ! J mol^-1
    ! Erd = 53000.                                        ! J mol^-1
    real,    parameter :: o2ci    = 210.e-3               ! mol mol^-1
    ! thermodynamic parameters for Vcmax & Jmax (Eq 9, Harley et al, 1992; #1392)
    real,    parameter :: Eavm    = 116300.               ! J/mol  (activation energy)
    real,    parameter :: Edvm    = 202900.               ! J/mol  (deactivation energy)
    real,    parameter :: Eajm    = 79500.                ! J/mol  (activation energy) 
    real,    parameter :: Edjm    = 201000.               ! J/mol  (deactivation energy)
    ! parameters for temperature dependence of gamma* (revised from von Caemmerer et al 1993)
    real,    parameter :: gam0    = 28.0e-6               ! mol mol^-1 @ 20C = 36.9 @ 25C
    real,    parameter :: gam1    = .0509
    real,    parameter :: gam2    = .0010
    real,    parameter :: times_storage_use=3*720.        ! 720 hours, 30 days
    real :: rhoS(3) = (/0.1, 0.3,   0.00/)                ! soil reflectance for vis, for NIR, for thermal, update in vegetation?
    ! end of consts parameters -------------------------------------------------------------------

    ! climate data type
    type forcing_data_type
        integer :: year
        integer :: doy
        integer :: hour
        real    :: Tair
        real    :: Tsoil
        real    :: RH                   ! Jian: RH seems confused in forcing and soil respiration
        real    :: VPD
        real    :: Rain
        real    :: WS
        real    :: PAR
        real    :: CO2
        real    :: PBOT                 ! unit patm Pa dynamic atmosphere pressure
        real    :: Ndep
    end type forcing_data_type
    type(forcing_data_type), allocatable, save :: forcing(:)
    integer :: nforcing
    real :: co2ca
    ! alternative input variable
    real,DIMENSION(:), ALLOCATABLE :: snow_in       ! if not run snow process, then read from the input file

    type site_data_type   ! data use in this model, but the common parameters in this site
        real :: lat, lon
        real :: wsmax
        real :: wsmin
        real :: extkU
        real :: rdepth
        real :: Q10rh
        real :: Q10pro
        real :: r_me
        real :: Toxi
        real :: Omax
        real :: kCH4
        real :: CH4_thre
        real :: Tveg
        real :: f
        real :: bubprob
        real :: Vmaxfraction  
        real :: etaW
        real :: tauC(8)
        real :: Tpro_me
        real :: f_F2M, f_C2M, f_C2S
        real :: f_M2S, f_M2P, f_S2P, f_S2M, f_P2M
        ! parameters need input initilized values
        real :: G
        real :: CN(8), CN0(8), QC(8), QN(8)
        real :: N_deposit, alphaN, QNminer, N_deficit
        real :: THKSL(nlayers), FRLEN(10)
        real :: liq_water(10), fwsoil, topfws, omega, zwt
        real :: infilt, sftmp, Tsnow, Twater, Tice, snow_dsim
        real :: dcount, dcount_soil
        real :: ice_tw, Tsoill(10), ice(10)
        real :: shcap_snow, condu_snow, condu_b
        real :: depth_ex, diff_s, diff_snow, albedo_snow
        real :: resht, thd_snow_depth, fa, fsub
        real :: rho_snow, decay_m
        real :: CH4_V(nlayers), CH4(nlayers), Vp(nlayers) 
        real :: bubble_methane_tot, Nbub
        ! parameters used in the cycle
        real :: GDD5
        ! some environmental variables
        real :: Dair, raero, ta, rain_d
        real :: tsoil_layer(11)
        real :: sublim
        real :: dpatm
        ! state variables
        real :: depth(nlayers)
        real :: wcl(10), wsc(10)
        real ::  scalW
        real :: Rsoilab1, Rsoilab2, Rsoilab3, Rsoilabs
        real :: melt, runoff
        real :: evap, ET
        real :: water_tw,  snow_depth
        ! flux 
        real :: NPP_L, NPP_W, NPP_R, NEE
        ! energy
        real :: Esoil, Hsoil
        ! soil flux
        real :: Rhetero, Rh_pools(5), OutC(8), OutN(8)
        ! methane
        real :: simuCH4, Pro_sum, Oxi_sum
        real :: pwater(nlayers), presP(nlayers)
        real :: methanebP(nlayers), methaneP(nlayers)
        ! N cycle
        real :: Rnitrogen, N_miner, N_transfer, N_uptake, N_fixation
        real :: N_leach, N_vol, N_loss, N_immob
        real :: fNnetmin
    end type site_data_type
    type(site_data_type) :: st

    ! pft data type for different species
    type spec_data_type
        real :: LAI
        ! special paramaters for different species
        real :: LAImin
        real :: LAImax
        real :: SLA
        real :: xfang
        real :: Vcmx0 
        real :: Vcmax0    ! fixed value from namelist
        real :: eJmx0
        real :: Entrpy
        real :: gddonset
        real :: stom_n
        real :: alpha
        real :: Ds0
        real :: Rl0
        real :: Rs0
        real :: Rr0
        real :: Q10
        real :: SapS, SapR
        real :: hmax            ! in plant growth hmax = 24.19   ! m
        real :: hl0             ! in plant growth hl0  = 0.00019  ! m2/kg C
        real :: LAIMAX0         ! in plant growth LAIMAX0 = 8.    ! maybe the LAImax
        real :: la0             ! in plant growht la0     = 0.2
        real :: GLmax, GSmax, GRmax
        ! flux 
        real :: gpp, npp
        real :: transp, plantup(10)
        real :: evap
        real :: RmLeaf, RmStem, RmRoot, Rmain, Rauto
        real :: Rgrowth
        real :: NPP_L, NPP_W, NPP_R
        ! states
        real :: QC(8), QN(8), CN0(8), CN(8), tauC(8), OutC(8)            ! leaf, stem, root, Litm, Lits
        real :: bmleaf, bmstem, bmroot, bmplant 
        real :: StemSap, RootSap, NSC, NSCmax, NSCmin, add
        real :: Stemmax, Rootmax
        real :: storage, stor_use, accumulation, store
        real :: L_fall
        ! special cycle 
        real :: Rh_pools(2)
        ! real :: f_F2M, f_C2M, f_C2S
        ! N scalar
        real :: SNvcmax, SNgrowth, SNRauto, JV
        real :: fnsc, NSN, N_miner, N_transfer, N_uptake, N_fixation
        real :: OutN(8), alphaN, N_deficit, Rnitrogen
        real :: N_leaf, N_wood, N_root, N_demand
        ! 
        integer :: onset
        real :: alpha_L, alpha_W, alpha_R
        ! energy
        real :: QLleaf, Esoil, Hsoil
    end type spec_data_type

    type vegn_tile_type
        integer :: npft
        type(spec_data_type), allocatable :: allSp(:)
        real :: LAI
        real :: LAImin
        real :: LAImax
        ! total flux
        real :: gpp, npp
        real :: NPP_L, NPP_W, NPP_R
        real :: transp, plantup(10)
        real :: evap
        real :: RmLeaf, RmStem, RmRoot, Rmain
        real :: Rauto, Rgrowth
        ! states
        real :: bmleaf, bmstem, bmroot, bmplant
        real :: NSC, NSN, storage, stor_use
        real :: N_leaf, N_wood, N_root
        ! energy
        real :: QLleaf, Rsoilab1, Rsoilab2, Esoil, Hsoil
    end type vegn_tile_type

    ! some parameters, may be in species, may be in cycle

    ! just one time output
    ! outputs
    type spec_outvars_type
        ! carbon fluxes (Kg C m-2 s-1)
        real :: gpp
        real :: nee
        real :: npp
        real :: nppLeaf
        real :: nppWood
        real :: nppStem
        real :: nppRoot
        real :: nppOther    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real :: ra
        real :: raLeaf
        real :: raStem
        real :: raRoot
        real :: raOther
        real :: rMaint
        real :: rGrowth
        real :: nbp
        ! Carbon Pools  (KgC m-2)
        real :: cLeaf
        real :: cStem
        real :: cRoot
        ! Nitrogen pools (kgN m-2)
        real :: nLeaf
        real :: nStem
        real :: nRoot
        ! real :: nOther(:)
        ! water fluxes (kg m-2 s-1)
        real :: tran
        ! other
        real :: lai                     ! m2 m-2, Leaf area index
    end type spec_outvars_type

    ! total outputs
    type outvars_data_type
        integer :: year
        integer :: doy
        integer :: hour
        type(spec_outvars_type), allocatable :: allSpec(:)
        ! carbon fluxes (Kg C m-2 s-1)
        real :: gpp
        real :: nee
        real :: npp
        real :: nppLeaf
        real :: nppWood
        real :: nppStem
        real :: nppRoot
        real :: nppOther           ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
        real :: ra
        real :: raLeaf
        real :: raStem
        real :: raRoot
        real :: raOther
        real :: rMaint
        real :: rGrowth            ! maintenance respiration and growth respiration
        real :: rh
        real :: nbp                ! heterotrophic respiration. NBP(net biome productivity) = GPP - Rh - Ra - other losses  
        real :: wetlandCH4
        real :: wetlandCH4prod
        real :: wetlandCH4cons     ! wetland net fluxes of CH4, CH4 production, CH4 consumption
        ! Carbon Pools  (KgC m-2)
        real :: cLeaf
        real :: cStem
        real :: cRoot
        real :: cOther              ! cOther: carbon biomass in other plant organs(reserves, fruits), Jian: maybe NSC storage in TECO?
        real :: cLitter
        real :: cLitterCwd          ! litter (excluding coarse woody debris), Jian: fine litter in TECO?, cLitterCwd: carbon in coarse woody debris
        real :: cSoil
        real :: cSoilLevels(nlayers)
        real :: cSoilFast
        real :: cSoilSlow
        real :: cSoilPassive        ! cSoil: soil organic carbon (Jian: total soil carbon); cSoilLevels(depth-specific soil organic carbon, Jian: depth?); cSoilPools (different pools without depth)
        real :: CH4(nlayers)              ! methane concentration
        ! Nitrogen fluxes (kgN m-2 s-1)
        real :: fBNF
        real :: fN2O
        real :: fNloss
        real :: fNnetmin
        real :: fNdep               ! fBNF: biological nitrogen fixation; fN2O: loss of nitrogen through emission of N2O; fNloss:Total loss of nitrogen to the atmosphere and from leaching; net mineralizaiton and deposition of N
        ! Nitrogen pools (kgN m-2)
        real :: nLeaf
        real :: nStem
        real :: nRoot
        real :: nOther
        real :: nLitter
        real :: nLitterCwd
        real :: nSoil
        real :: nMineral                ! nMineral: Mineral nitrogen pool
        ! energy fluxes (W m-2)
        real :: hfls
        real :: hfss
        real :: SWnet
        real :: LWnet                   ! Sensible heat flux; Latent heat flux; Net shortwave radiation; Net longwave radiation
        ! water fluxes (kg m-2 s-1)
        real :: ec
        real :: tran
        real :: es                      ! Canopy evaporation; Canopy transpiration; Soil evaporation
        real :: hfsbl                   ! Snow sublimation
        real :: mrro
        real :: mrros
        real :: mrrob                   ! Total runoff; Surface runoff; Subsurface runoff
        ! other
        real :: mrso(nlayers)           ! Kg m-2, soil moisture in each soil layer
        real :: tsl(nlayers)            ! K, soil temperature in each soil layer
        real :: tsland                  ! K, surface temperature
        real :: wtd                     ! m, Water table depth
        real :: snd                     ! m, Total snow depth
        real :: lai                     ! m2 m-2, Leaf area index 
    end type outvars_data_type
    type(outvars_data_type) :: outVars_h, outVars_d, outVars_m, outVars_y
    type(outvars_data_type), allocatable :: tot_outVars_h(:) 
    type(outvars_data_type), allocatable :: tot_outVars_d(:)
    type(outvars_data_type), allocatable :: tot_outVars_m(:)
    type(outvars_data_type), allocatable :: tot_outVars_y(:)

    ! parameters from the namelist file
    type nml_params_data_type
        real :: lat, lon
        real :: wsmax, wsmin
        real :: LAImin, LAImax
        real :: SLAx, rdepth
        real :: Rootmax, Stemmax
        real :: SapR, SapS
        real :: GLmax, GRmax, Gsmax
        real :: stom_n, a1, Ds0
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
    end type nml_params_data_type

    type nml_initValue_data_type
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
        real :: thd_snow_depth, b_bound     ! tuneice  not sensitive for ice
        real :: infilt_rate, fa, fsub       ! infilt_rate no used
        real :: rho_snow, decay_m           ! aging factor on snow melting
        ! methane module. update: Shuang methane bog species even more shallowly rooted than the tundra. add initials for methane module Shuang version
        real :: CH4_V(10), CH4(10), Vp(10)  ! assume in the very beginning no bubbles exist in the first three layers (30cm)
        real :: bubble_methane_tot, Nbub
        real :: depth_1                     ! calculate soil depth unit cm
    end type nml_initValue_data_type

contains

    subroutine read_teco_configs()
    ! read the config file of TECO model, which is used for setting the simulaiton
        implicit none
        integer io
        character(300) :: spec_names_0(10), files_pft_params_0(10)
        namelist /nml_teco_settings/ case_name, do_simu, do_mcmc, do_spinup, do_matrix, &
            do_restart, do_snow, do_soilphy, do_EBG, do_ndep, do_leap, do_out_hr,       &
            do_out_day, do_out_mon, do_out_yr, dtimes, inDir, outDir, climfile,         &
            watertablefile, snowdepthfile, in_restartfile, mcmc_configfile,             &
            spinup_configfile, file_site_params, count_pft, spec_names_0, files_pft_params_0
        namelist /nml_exps/ Ttreat, CO2treat, N_fert
        print *, "# read TECO config nml file ...", teco_configfile
        open(388, file = teco_configfile)
        read(388, nml  = nml_teco_settings,  iostat=io)
        read(388, nml  = nml_exps,           iostat=io)
        close(388)
        allocate(spec_names(count_pft))
        allocate(files_pft_params(count_pft))
        spec_names = spec_names_0(1:count_pft)
        files_pft_params = files_pft_params_0(1:count_pft)
    end subroutine read_teco_configs

    subroutine read_parameters_nml(param_nml_file, in_params, init_params)
    ! read the parameter file and the initilized parameters, return in_params and init_params
        implicit none
        character(*), intent(in) :: param_nml_file
        type(nml_params_data_type), intent(inout)    :: in_params
        type(nml_initValue_data_type), intent(inout) :: init_params
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
        print *, "# read parameters nml file: ", param_nml_file
        open(343, file = param_nml_file)
        read(343, nml  = nml_params,         iostat=io)
        read(343, nml  = nml_initial_values, iostat=io)
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
    end subroutine read_parameters_nml

    subroutine initilize(site_params, vegn_params, vegn)
        implicit none
        character(*), intent(in) :: site_params, vegn_params(:)
        type(nml_params_data_type)    :: in_params
        type(nml_initValue_data_type) :: init_params
        type(vegn_tile_type), intent(inout) :: vegn
        call read_parameters_nml(adjustl(trim("configs/"//adjustl(trim(site_params)))), &
                in_params, init_params)
        call initilize_site(in_params, init_params)
        call initilize_vegn(vegn, vegn_params)
        return
    end subroutine initilize

    subroutine initilize_site(in_params, init_params)
        implicit none
        type(nml_params_data_type), intent(inout)    :: in_params
        type(nml_initValue_data_type), intent(inout) :: init_params
        integer :: i
        ! site based parameters
        st%lat          = in_params%lat
        st%lon          = in_params%lon
        st%wsmax        = in_params%wsmax
        st%wsmin        = in_params%wsmin
        st%extkU        = in_params%extkU
        st%rdepth       = in_params%rdepth
        st%Q10rh        = in_params%Q10rh
        st%Q10pro       = in_params%Q10pro
        st%r_me         = in_params%r_me
        st%Toxi         = in_params%Toxi
        st%Omax         = in_params%Omax
        st%kCH4         = in_params%kCH4
        st%CH4_thre     = in_params%CH4_thre
        st%Tveg         = in_params%Tveg
        st%bubprob      = in_params%bubprob
        st%Vmaxfraction = in_params%Vmaxfraction
        st%etaW         = in_params%etaW
        st%f            = in_params%f
        st%tauC(1)      = in_params%Tau_Leaf*8760.
        st%tauC(2)      = in_params%Tau_Wood*8760.
        st%tauC(3)      = in_params%Tau_Root*8760.
        st%tauC(4)      = in_params%Tau_F*8760.
        st%tauC(5)      = in_params%Tau_C*8760.
        st%tauC(6)      = in_params%Tau_Micro*8760.
        st%tauC(7)      = in_params%Tau_slowSOM*8760.
        st%tauC(8)      = in_params%Tau_Passive*8760.
        st%Tpro_me      = in_params%Tpro_me
        st%f_F2M        = in_params%f_F2M
        st%f_C2M        = in_params%f_C2M
        st%f_C2S        = in_params%f_C2S
        st%f_M2S        = in_params%f_M2S
        st%f_M2P        = in_params%f_M2P
        st%f_S2P        = in_params%f_S2P
        st%f_S2M        = in_params%f_S2M
        st%f_P2M        = in_params%f_P2M
        ! initilize data
        st%G            = init_params%G
        st%QC           = init_params%QC
        st%CN0          = init_params%CN0
        st%CN           = init_params%CN0
        st%QN           = init_params%QC/init_params%CN0
        st%N_deposit    = init_params%N_deposit/8760.   ! Nitrogen input (gN/h/m2, )
        st%alphaN       = init_params%alphaN
        st%QNminer      = init_params%QNminer
        st%N_deficit    = init_params%N_deficit
        st%THKSL        = init_params%thksl
        st%FRLEN        = init_params%FRLEN
        st%liq_water    = init_params%liq_water
        st%fwsoil       = init_params%fwsoil
        st%topfws       = init_params%topfws
        st%omega        = init_params%omega
        st%zwt          = init_params%zwt
        st%infilt       = init_params%infilt
        st%sftmp        = init_params%sftmp
        st%Tsnow        = init_params%Tsnow
        st%Twater       = init_params%Twater
        st%Tice         = init_params%Tice
        st%snow_dsim    = init_params%snow_dsim
        st%dcount       = init_params%dcount
        st%dcount_soil  = init_params%dcount_soil
        st%ice_tw       = init_params%ice_tw
        st%Tsoill       = init_params%Tsoill
        st%ice          = init_params%ice
        st%shcap_snow   = init_params%shcap_snow
        st%condu_snow   = init_params%condu_snow
        st%condu_b      = init_params%condu_b
        st%depth_ex     = init_params%depth_ex
        st%diff_s       = init_params%diff_s
        st%diff_snow    = init_params%diff_snow
        st%albedo_snow  = init_params%albedo_snow
        st%resht        = init_params%resht
        st%thd_snow_depth = init_params%thd_snow_depth
        st%fa           = init_params%fa
        st%fsub         = init_params%fsub
        st%rho_snow     = init_params%rho_snow
        st%decay_m      = init_params%decay_m
        st%CH4_V        = init_params%CH4_V
        st%CH4          = init_params%CH4
        st%Vp           = init_params%Vp
        st%bubble_methane_tot = init_params%bubble_methane_tot
        st%Nbub         = st%Nbub
        st%depth(1)     = init_params%depth_1
        do i = 1, 10
            st%wcl(i)   = st%wsmax/100.
        enddo
        do i = 2, nlayers
            st%depth(i)=st%depth(i-1)+st%THKSL(i)
        enddo
        do i=1,nlayers
            if (st%depth(i) .le. (-st%zwt)*0.1) then
                st%pwater(i) = 1000*9.81*(st%depth(i)*0.01-(-st%zwt)*0.001)
            else
                st%pwater(i) = 0.
            endif
            st%presP(i) = 101325 + st%pwater(i)  ! unit Pa
            st%methanebP(i) = st%f * st%presP(i) * st%Vp(i)/(8.3144621 * (st%Tsoill(i)+273.15))  !unit mol/layer
            st%methaneP(i) = st%CH4(i)/12
            ! gC/layer  /12   unit molC/layer
        enddo
    end subroutine initilize_site

    subroutine initilize_vegn(vegn, files_vegn_params)
        implicit none
        type(vegn_tile_type), intent(inout) :: vegn
        character(*), intent(in) :: files_vegn_params(:)
        type(nml_params_data_type)    :: in_params
        type(nml_initValue_data_type) :: init_params
        integer :: ipft, npft

        npft = size(files_vegn_params)
        allocate(vegn%allSp(npft))
        vegn%npft = npft
        if(allocated(vegn%allSp)) then
            do ipft = 1, npft
                call read_parameters_nml(adjustl(trim("configs/"//adjustl(trim(files_vegn_params(ipft))))), &
                        in_params, init_params)
                call initilize_spec(vegn%allSp(ipft), in_params, init_params)
                if (ipft .eq. 1) then
                    vegn%LAImax = vegn%allSp(ipft)%LAImax
                    vegn%LAImin = vegn%allSp(ipft)%LAImin
                else
                    vegn%LAImax = AMAX1(vegn%LAImax, vegn%allSp(ipft)%LAImax)
                    vegn%LAImin = AMAX1(vegn%LAImin, vegn%allSp(ipft)%LAImin)
                endif
            enddo
        endif
        
    end subroutine initilize_vegn

    subroutine initilize_spec(spec, in_params, init_params)
        implicit none
        type(spec_data_type), intent(inout)          :: spec
        type(nml_params_data_type), intent(inout)    :: in_params
        type(nml_initValue_data_type), intent(inout) :: init_params

        spec%LAImin = in_params%LAImin
        spec%LAImax = in_params%LAImax
        spec%SLA    = in_params%SLAx/10000.          ! Convert unit from cm2/g to m2/g
        spec%Rootmax = in_params%Rootmax
        spec%Stemmax = in_params%Stemmax
        spec%SapR    = in_params%SapR
        spec%SapS    = in_params%SapS
        spec%GLmax   = in_params%GLmax/8760.          ! growth rates of plant. Jian: per year to per hour ?
        spec%GRmax   = in_params%GRmax/8760.
        spec%Gsmax   = in_params%Gsmax/8760.
        spec%stom_n  = in_params%stom_n
        spec%Ds0     = in_params%Ds0
        spec%Vcmax0  = in_params%Vcmax0
        spec%xfang   = in_params%xfang
        spec%alpha   = in_params%alpha
        spec%gddonset = in_params%gddonset
        spec%Q10     = in_params%Q10
        spec%Rl0     = in_params%Rl0
        spec%Rs0     = in_params%Rs0
        spec%Rr0     = in_params%Rr0
        spec%JV      = in_params%JV
        spec%Entrpy  = in_params%Entrpy
        spec%alphaN  = init_params%alphaN
        ! 
        spec%tauC(1)      = in_params%Tau_Leaf*8760.
        spec%tauC(2)      = in_params%Tau_Wood*8760.
        spec%tauC(3)      = in_params%Tau_Root*8760.
        ! st%tauC(4)      = in_params%Tau_F
        ! st%tauC(5)      = in_params%Tau_C
        ! st%tauC(6)      = in_params%Tau_Micro
        ! st%tauC(7)      = in_params%Tau_slowSOM
        ! st%tauC(8)      = in_params%Tau_Passive
        spec%QC      = init_params%QC
        spec%CN0     = init_params%CN0
        spec%CN      = init_params%CN0
        spec%QN      = init_params%QC/init_params%CN0
        spec%NSCmin  = init_params%NSCmin
        spec%storage = init_params%Storage
        spec%stor_use= spec%Storage/times_storage_use
        spec%nsc     = init_params%nsc
        spec%accumulation = init_params%accumulation
        spec%SNvcmax = init_params%SNvcmax
        spec%NSN     = init_params%NSN
        spec%N_deficit = init_params%N_deficit
        spec%LAI     = spec%LAImin
        spec%bmleaf  = spec%QC(1)/0.48
        spec%bmstem  = spec%QC(2)/0.48
        spec%bmroot  = spec%QC(3)/0.48
        spec%hmax    = in_params%hmax        ! in plant growth hmax = 24.19   ! m
        spec%hl0     = in_params%hl0             ! in plant growth hl0  = 0.00019  ! m2/kg C
        spec%LAIMAX0 = in_params%LAIMAX0        ! in plant growth LAIMAX0 = 8.    ! maybe the LAImax
        spec%la0     = in_params%la0        ! in plant growht la0     = 0.2
        return
    end subroutine initilize_spec

    ! subroutine assign_outVars(outVars, ntime, nspec)
    !     type(outvars_data_type), intent(inout) :: outVars
    !     integer, intent(in) :: ntime, nspec
    !     integer :: itime
    !     allocate(outVars(ntime))
    !     do itime = 1, ntime
    !         allocate(outVars(itime)%allSpec(nspec))
    !     enddo
    !     return
    ! end subroutine assign_outVars

    subroutine init_hourly()
        implicit none
        if (do_out_hr) call init_outVars(outVars_h)
    end subroutine init_hourly

    subroutine init_daily()
        implicit none
        if (do_out_day) call init_outVars(outVars_d)
    end subroutine init_daily

    subroutine init_monthly()
        implicit none
        if(do_out_mon) call init_outVars(outVars_m)
    end subroutine init_monthly

    subroutine init_yearly(vegn)
        implicit none
        type(vegn_tile_type), intent(inout) :: vegn
        integer ipft
        st%GDD5 = 0.
        do ipft = 1, vegn%npft
            vegn%allSp(ipft)%onset = 0
        enddo 
        if (do_out_yr) call init_outVars(outVars_y)
    end subroutine init_yearly

    subroutine init_outVars(outVars)
        implicit none
        type(outvars_data_type), intent(inout) :: outVars
        integer :: npft, ipft
        outVars%year = 0
        outVars%doy  = 0
        outVars%hour = 0
        if (allocated(outVars%allSpec))then
            npft = size(outVars%allSpec)
            do ipft = 1, npft
                ! carbon fluxes (Kg C m-2 s-1)
                outVars%allSpec(ipft)%gpp      = 0.
                outVars%allSpec(ipft)%nee      = 0.
                outVars%allSpec(ipft)%npp      = 0.
                outVars%allSpec(ipft)%nppLeaf  = 0.
                outVars%allSpec(ipft)%nppWood  = 0.
                outVars%allSpec(ipft)%nppStem  = 0.
                outVars%allSpec(ipft)%nppRoot  = 0.
                outVars%allSpec(ipft)%nppOther = 0.    ! According to SPRUCE-MIP, stem means above ground woody tissues which is different from wood tissues.
                outVars%allSpec(ipft)%ra       = 0.
                outVars%allSpec(ipft)%raLeaf   = 0.
                outVars%allSpec(ipft)%raStem   = 0.
                outVars%allSpec(ipft)%raRoot   = 0.
                outVars%allSpec(ipft)%raOther  = 0.
                outVars%allSpec(ipft)%rMaint   = 0.
                outVars%allSpec(ipft)%rGrowth  = 0.
                outVars%allSpec(ipft)%nbp      = 0.
                ! Carbon Pools  (KgC m-2)
                outVars%allSpec(ipft)%cLeaf    = 0.
                outVars%allSpec(ipft)%cStem    = 0.
                outVars%allSpec(ipft)%cRoot    = 0.
                ! Nitrogen pools (kgN m-2)
                outVars%allSpec(ipft)%nLeaf    = 0.
                outVars%allSpec(ipft)%nStem    = 0.
                outVars%allSpec(ipft)%nRoot    = 0.
                ! water fluxes (kg m-2 s-1)
                outVars%allSpec(ipft)%tran     = 0.
                ! other
                outVars%allSpec(ipft)%lai      = 0. 
            enddo
        endif
        outVars%gpp              = 0.
        outVars%nee              = 0.
        outVars%npp              = 0.
        outVars%nppLeaf          = 0.
        outVars%nppWood          = 0.
        outVars%nppStem          = 0.
        outVars%nppRoot          = 0.
        outVars%nppOther         = 0.  
        outVars%ra               = 0.
        outVars%raLeaf           = 0.
        outVars%raStem           = 0.
        outVars%raRoot           = 0.
        outVars%raOther          = 0.
        outVars%rMaint           = 0.
        outVars%rGrowth          = 0.
        outVars%rh               = 0.
        outVars%nbp              = 0.
        outVars%wetlandCH4       = 0.
        outVars%wetlandCH4prod   = 0.
        outVars%wetlandCH4cons   = 0. 
        ! Carbon Pools  (KgC m-2)
        outVars%cLeaf            = 0.
        outVars%cStem            = 0.
        outVars%cRoot            = 0.
        outVars%cOther           = 0.
        outVars%cLitter          = 0.
        outVars%cLitterCwd       = 0.  
        outVars%cSoil            = 0.
        outVars%cSoilLevels(:)   = 0.
        outVars%cSoilFast        = 0.
        outVars%cSoilSlow        = 0.
        outVars%cSoilPassive     = 0. 
        outVars%CH4(:)           = 0.
        ! Nitrogen fluxes (kgN m-2 s-1)
        outVars%fBNF             = 0.
        outVars%fN2O             = 0.
        outVars%fNloss           = 0.
        outVars%fNnetmin         = 0.
        outVars%fNdep            = 0.  
        ! Nitrogen pools (kgN m-2)
        outVars%nLeaf            = 0.
        outVars%nStem            = 0.
        outVars%nRoot            = 0.
        outVars%nOther           = 0.
        outVars%nLitter          = 0.
        outVars%nLitterCwd       = 0.
        outVars%nSoil            = 0.
        outVars%nMineral         = 0. 
        ! energy fluxes (W m-2)
        outVars%hfls             = 0.
        outVars%hfss             = 0.
        outVars%SWnet            = 0.
        outVars%LWnet            = 0.
        ! water fluxes (kg m-2 s-1)
        outVars%ec               = 0.
        outVars%tran             = 0.
        outVars%es               = 0.   
        outVars%hfsbl            = 0.  
        outVars%mrro             = 0.
        outVars%mrros            = 0.
        outVars%mrrob            = 0.   
        ! other
        outVars%mrso(:)          = 0.  
        outVars%tsl(:)           = 0.
        outVars%tsland           = 0.                 
        outVars%wtd              = 0.           
        outVars%snd              = 0.           
        outVars%lai              = 0.
    end subroutine init_outVars

    subroutine get_forcingdata()
        implicit none
        integer STAT, COUNT
        character(150) commts
        ! define variable for each line
        integer :: tmp_yr, tmp_doy, tmp_h
        real    :: tmp_Ta, tmp_Ts,  tmp_rh, tmp_vpd, tmp_rain, tmp_ws 
        real    :: tmp_par, tmp_co2, tmp_pbot, tmp_ndep

        call ReadLineNumFromFile(climfile, nforcing)  ! get the line number

        allocate(forcing(nforcing))                   ! allocate the array

        COUNT = 0
        OPEN(1,FILE=climfile,status='old',ACTION='read',IOSTAT=STAT)
        read(1,'(a160)') commts
        DO WHILE (.TRUE.)
            COUNT=COUNT+1
            READ(1,*,IOSTAT=STAT, end=993) tmp_yr, tmp_doy, tmp_h,   &
                tmp_Ta,  tmp_Ts,  tmp_rh, tmp_vpd, tmp_rain, tmp_ws, & 
                tmp_par, tmp_co2, tmp_pbot, tmp_ndep
            IF(STAT .NE. 0) EXIT
            forcing(COUNT)%year  = tmp_yr
            forcing(COUNT)%doy   = tmp_doy
            forcing(COUNT)%hour  = tmp_h
            forcing(COUNT)%Tair  = tmp_Ta
            forcing(COUNT)%Tsoil = tmp_Ts
            forcing(COUNT)%RH    = tmp_rh
            forcing(COUNT)%VPD   = tmp_vpd
            forcing(COUNT)%Rain  = tmp_rain
            forcing(COUNT)%WS    = tmp_ws
            forcing(COUNT)%PAR   = tmp_par
            forcing(COUNT)%CO2   = tmp_co2
            forcing(COUNT)%PBOT  = tmp_pbot
            forcing(COUNT)%Ndep  = tmp_ndep
        ENDDO
993     continue
        CLOSE(1)
    end subroutine get_forcingdata

    subroutine get_snowdepth()
        implicit none
        ! real temp_snow_depth(max_nlines)
        integer STAT, COUNT, nrow
        character(50) commts

        ! integer m,n,istat1,lines,yr_length
        real snow_depth_read
        integer tmp_yr, tmp_doy, tmp_hr

        call ReadLineNumFromFile(snowdepthfile, nrow)  ! get the line number
        allocate(snow_in(nrow))

        open(11,file = snowdepthfile, status ='old',ACTION='read', IOSTAT=STAT)
        read(11,'(a160)') commts ! skip 2 lines of input met data file
        COUNT = 0
        do
            COUNT = COUNT + 1
            read (11,*,IOSTAT=STAT, end=1018) tmp_yr,tmp_doy,tmp_hr,snow_depth_read
            IF(STAT .NE. 0) EXIT
            snow_in(COUNT)=snow_depth_read     
        enddo
1018    continue
        close(11)    ! close snow depth file
        return
    end subroutine get_snowdepth

    subroutine ReadLineNumFromFile(filepath, count_lines)
        implicit none
        character(len=*), intent(in) :: filepath
        character(len=100) header, line
        integer STAT, count_lines
        print*, "file path: ", trim(filepath)
        open(38, file=trim(filepath), status="old", action="read", iostat=STAT) ! open file
        read(38, '(a100)') header           ! read the header of the file
        count_lines = 0                     ! initilize the count_lines
        do while(.TRUE.)
            read(38, *, iostat=STAT) line   ! read each line
            if(STAT .ne. 0) exit            ! until the end of the file
            count_lines = count_lines + 1   ! recording the count of the lines
        enddo
        close(38)
        return
    end subroutine ReadLineNumFromFile

    subroutine deallocate_results(outVars, npft)
        implicit none
        type(outvars_data_type), intent(inout) :: outVars
        integer, intent(in) :: npft
        integer :: ipft

        do ipft = 1, npft
            if(allocated(outVars%allSpec)) deallocate(outVars%allSpec)
        enddo
    end subroutine deallocate_results

    subroutine deallocate_date_type()
        if (allocated(forcing)) deallocate(forcing)
        if (allocated(snow_in)) deallocate(snow_in)
        if (allocated(spec_names)) deallocate(spec_names)
        if (allocated(files_pft_params)) deallocate(files_pft_params)
    end subroutine deallocate_date_type
end module datatypes