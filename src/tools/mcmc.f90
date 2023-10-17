module mcmc
    use driver
    use datatypes
    use mcmc_mod
    use MCMC_outputs

    implicit none

    integer npar4DA, ipar, covexist
    
    real fact_rejet
    real J_last, J_new, accept_rate
    integer new, reject
    logical do_cov2createNewPars, do_cov

    type(nml_params_data_type),allocatable     :: mc_in_params(:)     !   in_params for MCMC
    type(nml_initValue_data_type), allocatable :: mc_init_params(:)   ! init_params for MCMC

    contains
    subroutine init_mcmc(files_vegn_params, vegn)
    ! *** files_vegn_params : array of files for vegetation parameters
    ! ***              vegn : return the vegn data
    ! *** put values to the variables for MCMC
        implicit none
        real, allocatable :: temp_parmin(:), temp_parmax(:), temp_parval(:)
        integer, allocatable :: temp_paridx(:)
        ! different PFTs
        type(vegn_tile_type), intent(inout) :: vegn
        character(*), intent(in) :: files_vegn_params(:)

        integer :: ipft, npft

        ! read the nml file of MCMC configs (eg. TECO_MCMC_configs.nml)
        call readConfsNml()     
        npar = nSpecParams
        npft = size(files_vegn_params)
        ! initilize the parameters and initial values in TECO model
        allocate(vegn%allSp(npft))
        vegn%npft = npft

        allocate(mc_parvals(npft))      ! parval, parmin, parmax
        allocate(mc_in_params(npft))
        allocate(mc_init_params(npft))

        if(allocated(vegn%allSp)) then
            do ipft = 1, count_pft
                allocate(mc_parvals(ipft)%parval(npar), mc_parvals(ipft)%parmin(npar), mc_parvals(ipft)%parmax(npar))
                call readParamNml(adjustl(trim("configs/"//adjustl(trim(files_vegn_params(ipft))))), &
                    mc_in_params(ipft), mc_init_params(ipft), &
                    mc_parvals(ipft)%parval, mc_parvals(ipft)%parmin, mc_parvals(ipft)%parmax)
                call initilize_site(mc_in_params(ipft), mc_init_params(ipft))  ! Jian: this version not separate the site parameters and pft parameters
                call initilize_spec(vegn%allSp(ipft), mc_in_params(ipft), mc_init_params(ipft))
                if (ipft .eq. 1) then
                    vegn%LAImax = vegn%allSp(ipft)%LAImax
                    vegn%LAImin = vegn%allSp(ipft)%LAImin
                else
                    vegn%LAImax = AMAX1(vegn%LAImax, vegn%allSp(ipft)%LAImax)
                    vegn%LAImin = AMAX1(vegn%LAImin, vegn%allSp(ipft)%LAImin)
                endif
            enddo
        endif ! finish the initilize parameters and the mc_parvals

        ! read the observational data
        call readObsData() ! return a type array of vars4MCMC
        ! handle the parameters for MCMC
        allocate(mc_DApar(npft)) 
        allocate(temp_parmin(npar), temp_parmax(npar))  ! allocate the temporary parmin value
        allocate(temp_paridx(npar), temp_parval(npar))  ! mark the index of parameters for MCMC
        ! allocate(MDparval(npar))                        ! record the parameters set for model simulation
        ! MDparval = parval                               ! parameters for running model
        do ipft = 1, npft
            npar4DA  = 0 ! record the number of parameters for data assimilation
            do ipar = 1, npar
                if (mc_parvals(ipft)%parmin(ipar) .ne. mc_parvals(ipft)%parmax(ipar)) then
                    npar4DA              = npar4DA + 1
                    temp_paridx(npar4DA) = ipar
                    temp_parmin(npar4DA) = mc_parvals(ipft)%parmin(ipar)
                    temp_parmax(npar4DA) = mc_parvals(ipft)%parmax(ipar)
                    temp_parval(npar4DA) = mc_parvals(ipft)%parval(ipar)
                endif
            enddo
            allocate(mc_DApar(ipft)%DAparmin(npar4DA),  mc_DApar(ipft)%DAparmax(npar4DA), mc_DApar(ipft)%DApar(npar4DA), &
                     mc_DApar(ipft)%DApar_old(npar4DA), mc_DApar(ipft)%DAparidx(npar4DA))
            ! allocate(DAparmin(npar4DA), DAparmax(npar4DA), DAparidx(npar4DA))
            ! allocate(DApar(npar4DA),    DApar_old(npar4DA))

            mc_DApar(ipft)%DAparmin  = temp_parmin(:npar4DA)
            mc_DApar(ipft)%DAparmax  = temp_parmax(:npar4DA)
            mc_DApar(ipft)%DAparidx  = temp_paridx(:npar4DA)
            mc_DApar(ipft)%DApar     = temp_parval(:npar4DA)
            mc_DApar(ipft)%DApar_old = mc_DApar(ipft)%DApar   ! mark as old parameters
        enddo

        deallocate(temp_parmin, temp_parmax, temp_parval, temp_paridx)

        ! give some values to the parameters for MCMC
        covexist      = 0
        fact_rejet    = 2.4/sqrt(real(npar4DA))

        ! record
        do ipft = 1, npft
            allocate(mc_DApar(ipft)%coefhistory(ncov, npar4DA))
            ! create the coefnorm for generating the new parameters
            allocate(mc_DApar(ipft)%coefnorm(npar4DA)) 
            allocate(mc_DApar(ipft)%coefac(npar4DA))
            do ipar = 1, npar4DA
                mc_DApar(ipft)%coefnorm(ipar) = 0.5
                mc_DApar(ipft)%coefac(ipar)   = mc_DApar(ipft)%coefnorm(ipar)
            enddo
            allocate(mc_DApar(ipft)%gamnew(npar4DA, npar4DA))
        enddo
        J_last = 9000000.0
        ! init the outputs
        call init_mcmc_outputs(nDAsimu, npar4DA)
        do_cov = .False.
        do_cov2createNewPars = .False.
    end subroutine init_mcmc

    subroutine run_mcmc(vegn)
        implicit none
        type(vegn_tile_type), intent(inout) :: vegn
        integer temp_upgraded, ipft, npft
        real rand
        
        print *, "# Start to run mcmc ..."
        npft = count_pft
        ! call generate_newPar()
        call generate_rand_newpar()
        upgraded = 0
        new = 0
        do iDAsimu = 1, nDAsimu
            write(*,*) iDAsimu, "/", nDAsimu, J_last, J_new, upgraded, accept_rate
            call mcmc_functions_init()  ! initialize the mc_itime ... variables
                
            ! generate parameters 
            call generate_newPar()
            do ipft = 1, npft   
                ! update the parameters
                do ipar = 1, npar4DA
                    mc_parvals(ipft)%parval(mc_DApar(ipft)%DAparidx(ipar)) = mc_DApar(ipft)%DApar(ipar)
                enddo
                call renewMDpars(mc_parvals(ipft)%parval, mc_in_params(ipft))          ! call update parameters in TECO model
            enddo
            
            ! call initialize()           ! initialize the TECO model 
            if(allocated(vegn%allSp)) then
                do ipft = 1, count_pft
                    call initilize_site(mc_in_params(ipft), mc_init_params(ipft))  ! Jian: this version not separate the site parameters and pft parameters
                    call initilize_spec(vegn%allSp(ipft), mc_in_params(ipft), mc_init_params(ipft))
                    if (ipft .eq. 1) then
                        vegn%LAImax = vegn%allSp(ipft)%LAImax
                        vegn%LAImin = vegn%allSp(ipft)%LAImin
                    else
                        vegn%LAImax = AMAX1(vegn%LAImax, vegn%allSp(ipft)%LAImax)
                        vegn%LAImin = AMAX1(vegn%LAImin, vegn%allSp(ipft)%LAImin)
                    endif
                enddo
            endif ! finish ! initialize the TECO model

            call teco_simu(vegn, .False.)            ! run the model
            if (iDAsimu .eq. nDAsimu) call teco_simu(vegn, .True.)
            
            temp_upgraded = upgraded
            call costFuncObs()          ! calculate the cost between observations and simulations

            ! if upgraded is updated in costFuncObs
            if (upgraded .gt. temp_upgraded) then
                new =  new + 1  ! new is for what?
                if (covexist .eq. 1)then
                    do ipft = 1, npft
                        mc_DApar(ipft)%coefac = mc_DApar(ipft)%coefnorm                           ! coefac is old parameter sets? coef is the new one; coefnorm 
                        mc_DApar(ipft)%coefhistory(new, :) = mc_DApar(ipft)%coefnorm              ! coefhistory used to create new coef matrix
                    enddo
                else
                    do ipft = 1, npft
                        do ipar = 1, npar4DA
                            mc_DApar(ipft)%coefnorm(ipar) = (mc_DApar(ipft)%DApar(ipar)-mc_DApar(ipft)%DAparmin(ipar))&
                                                           /(mc_DApar(ipft)%DAparmax(ipar)-mc_DApar(ipft)%DAparmin(ipar))
                        enddo
                    enddo
                endif
                do ipft = 1, npft
                    mc_DApar(ipft)%coefhistory(new, :) = mc_DApar(ipft)%coefnorm 
                enddo
                if(new .ge. ncov)new=0
                ! update the parameters sets
                do ipft = 1, npft
                    arr_params_set(ipft)%tot_paramsets(upgraded,:) = mc_DApar(ipft)%DApar
                enddo
                ! if (do_mc_out_hr) then
                !     call mcmc_update_outputs(upgraded, tot_paramsets_outs_h, outVars_h)
                ! endif
                ! if (do_mc_out_day)then
                !     call mcmc_update_outputs(upgraded, tot_paramsets_outs_d, outVars_d)
                ! endif
                ! if (do_mc_out_mon) then
                !     call mcmc_update_outputs(upgraded, tot_paramsets_outs_m, outVars_m)
                ! endif
            else
                reject = reject + 1
            endif

            ! updates of the covariance matrix
            if(do_cov)then
                do ipft = 1, npft
                    
                    if (.not. do_cov2createNewPars .and. mod(upgraded, ncov).eq.0 .and. upgraded .ne. 0)then
                        do_cov2createNewPars = .True.
                        mc_DApar(ipft)%coefac   = mc_DApar(ipft)%coefnorm          ! coefnorm: normized values between min and max values
                        call varcov(mc_DApar(ipft)%coefhistory, mc_DApar(ipft)%gamnew, npar4DA, ncov) !
                        if(.not.(all(mc_DApar(ipft)%gamnew==0.)))then
                            mc_DApar(ipft)%gamma = mc_DApar(ipft)%gamnew
                            call racine_mat(mc_DApar(ipft)%gamma, mc_DApar(ipft)%gamnew, npar4DA)
                            mc_DApar(ipft)%gamma = mc_DApar(ipft)%gamnew
                        endif
                    endif
                    

                    if(mod(upgraded, ncov).eq.0 .and. covexist.eq.1 .and. upgraded .ne. 0)then
                        call varcov(mc_DApar(ipft)%coefhistory, mc_DApar(ipft)%gamnew, npar4DA, ncov)
                        if(.not.(all(mc_DApar(ipft)%gamnew==0.)))then
                            mc_DApar(ipft)%gamma = mc_DApar(ipft)%gamnew
                            call racine_mat(mc_DApar(ipft)%gamma, mc_DApar(ipft)%gamnew, npar4DA)
                            mc_DApar(ipft)%gamma = mc_DApar(ipft)%gamnew
                        endif
                    endif
                enddo
            endif
        enddo

        ! summary
        call mcmc_param_outputs(upgraded, npar4DA, parnames)!, mc_DApar)
    end subroutine run_mcmc

    subroutine generate_newPar()
        ! This subroutine is used to generate the new parameters to run MCMC
        ! Based on the Shuang's code, it need to use the coef to generate the new parameters.
        implicit none
        ! real, intent(in) :: par_old(:), par_min(:), par_max(:)
        ! real, intent(inout) :: par_new(:) 
        integer igenPar, parflag, ipft, npft
        real rand_harvest, rand

        call random_seed()

        ! DApar_old = DApar                   ! mark as old parameters 
        npft = count_pft
        if (do_cov2createNewPars) then
            do ipft = 1, npft
                parflag = 1                 ! mark
                do while(parflag .gt. 0)    ! create the new coefnorm
                    ! create the new coefnorm based on old one of coefac
                    call gengaussvect(fact_rejet*mc_DApar(ipft)%gamma, mc_DApar(ipft)%coefac, &
                         mc_DApar(ipft)%coefnorm, npar4DA)          ! generate the new cov parameters
                    parflag = 0
                    do igenPar = 1, npar4DA                                                 ! check the cov 
                        if(mc_DApar(ipft)%coefnorm(igenPar).lt.0. .or. mc_DApar(ipft)%coefnorm(igenPar).gt.1.)then
                            parflag=parflag+1
                            ! write(*,*)'out of range',parflag
                        endif
                    enddo
                enddo
                ! create the new parameters from 
                do ipar = 1, npar4DA
                    mc_DApar(ipft)%DApar(ipar) = mc_DApar(ipft)%DAparmin(ipar) + &
                        mc_DApar(ipft)%coefnorm(ipar) * (mc_DApar(ipft)%DAparmax(ipar)-mc_DApar(ipft)%DAparmin(ipar))
                enddo
            enddo
        else ! do not run cov to create new parameters, just random selections
            do ipft = 1, npft
                do igenPar = 1, npar4DA     ! for each parameters
    999             continue
                    call random_number(rand_harvest)    
                    rand = rand_harvest - 0.5           ! create a random number in [-0.5, 0.5]
                    mc_DApar(ipft)%DApar(igenPar) = mc_DApar(ipft)%DApar_old(igenPar) + &
                        rand*(mc_DApar(ipft)%DAparmax(igenPar) - mc_DApar(ipft)%DAparmin(igenPar)) * search_scale   ! create new parameter
                    if((mc_DApar(ipft)%DApar(igenPar) .gt. mc_DApar(ipft)%DAparmax(igenPar)) &
                        &   .or. (mc_DApar(ipft)%DApar(igenPar) .lt. mc_DApar(ipft)%DAparmin(igenPar))) then 
                        goto 999                  ! judge the range of new parameter
                    endif
                enddo
            enddo
        endif
        return   ! mainly return the DApar, meanwhile update the coefnorm
    end subroutine generate_newPar

    subroutine generate_rand_newpar()
        integer igenPar, parflag, ipft
        real rand_harvest, rand

        call random_seed()
        do ipft = 1, count_pft
            do igenPar = 1, npar4DA     ! for each parameters
1999             continue
                call random_number(rand_harvest)    
                rand = rand_harvest - 0.5           ! create a random number in [-0.5, 0.5]
                mc_DApar(ipft)%DApar(igenPar) = mc_DApar(ipft)%DAparmin(igenPar) + &
                    rand_harvest*(mc_DApar(ipft)%DAparmax(igenPar) - mc_DApar(ipft)%DAparmin(igenPar))   ! create new parameter
                if((mc_DApar(ipft)%DApar(igenPar) .gt. mc_DApar(ipft)%DAparmax(igenPar)) &
                    &   .or. (mc_DApar(ipft)%DApar(igenPar) .lt. mc_DApar(ipft)%DAparmin(igenPar))) then 
                    goto 1999                  ! judge the range of new parameter
                endif
            enddo
        enddo
        return
    end subroutine generate_rand_newpar

    subroutine costFuncObs()
        implicit none
        real J_cost, delta_J, cs_rand
        integer :: ipft, npft
        
        J_new = 0

        ! ANPP_Shrub_y
        if(vars4MCMC%ANPP_Shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
                 vars4MCMC%ANPP_Shrub_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/200
        endif
        ! print*, "J_new1: ", J_new
        ! print*, "test1:",vars4MCMC%ANPP_Shrub_y%mdData(:,4), vars4MCMC%ANPP_Shrub_y%obsData(:,4),&
        !          vars4MCMC%ANPP_Shrub_y%obsData(:,5)

        ! ANPP_Tree_y
        if(vars4MCMC%ANPP_Tree_y%existOrNot)then
            call CalculateCost(vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%obsData(:,4),&
                 vars4MCMC%ANPP_Tree_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/800
        endif
        ! print*, "J_new2: ", J_new
        ! print*, "test1:",vars4MCMC%ANPP_Tree_y%mdData(:,4), vars4MCMC%ANPP_Tree_y%mdData(:,1),&
        !          vars4MCMC%ANPP_Tree_y%mdData(:,2),vars4MCMC%ANPP_Tree_y%mdData(:,3)

        ! NPP_sphag_y
        if(vars4MCMC%NPP_sphag_y%existOrNot)then
            call CalculateCost(vars4MCMC%NPP_sphag_y%mdData(:,4), vars4MCMC%NPP_sphag_y%obsData(:,4),&
                 vars4MCMC%NPP_sphag_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/2000
        endif
        ! print*, "J_new3: ", J_new

        ! BNPP_y        ! tree + shrub
        if(vars4MCMC%BNPP_y%existOrNot)then
            call CalculateCost(vars4MCMC%BNPP_y%mdData(:,4), vars4MCMC%BNPP_y%obsData(:,4),&
                 vars4MCMC%BNPP_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/400
        endif
        ! print*, "J_new4: ", J_new
        ! er_d          ! shrub + sphag.
        if(vars4MCMC%er_d%existOrNot)then
            call CalculateCost(vars4MCMC%er_d%mdData(:,4), vars4MCMC%er_d%obsData(:,4),&
                 vars4MCMC%er_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new5: ", J_new
        ! er_h          ! shrub + sphag.

        if(vars4MCMC%er_h%existOrNot)then
            call CalculateCost(vars4MCMC%er_h%mdData(:,4), vars4MCMC%er_h%obsData(:,4),&
                 vars4MCMC%er_h%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new6: ", J_new
        ! gpp_d         ! Shrub + sphag.
        if(vars4MCMC%gpp_d%existOrNot)then
            call CalculateCost(vars4MCMC%gpp_d%mdData(:,4), vars4MCMC%gpp_d%obsData(:,4),&
                 vars4MCMC%gpp_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new7: ", J_new
        ! nee_d         ! Shrub + sphag.
        if(vars4MCMC%nee_d%existOrNot)then
            call CalculateCost(vars4MCMC%nee_d%mdData(:,4), vars4MCMC%nee_d%obsData(:,4),&
                 vars4MCMC%nee_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost/30
        endif
        ! print*, "J_new8: ", J_new
        ! nee_h         ! shrub + sphag.
        if(vars4MCMC%nee_h%existOrNot)then
            call CalculateCost(vars4MCMC%nee_h%mdData(:,4), vars4MCMC%nee_h%obsData(:,4),&
                 vars4MCMC%nee_h%obsData(:,5), J_cost)
            J_new = J_new + J_cost/50
        endif
        ! print*, "J_new9: ", J_new
        ! LAI_d         ! tree  + Shrub
        if(vars4MCMC%LAI_d%existOrNot)then
            call CalculateCost(vars4MCMC%LAI_d%mdData(:,4), vars4MCMC%LAI_d%obsData(:,4),&
                 vars4MCMC%LAI_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new10: ", J_new

        ! leaf_mass_shrub_y
        if(vars4MCMC%leaf_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_mass_shrub_y%mdData(:,4), vars4MCMC%leaf_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%leaf_mass_shrub_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/8000
        endif
        ! print*, "J_new11: ", J_new

        ! stem_mass_shrub_y
        if(vars4MCMC%stem_mass_shrub_y%existOrNot)then
            call CalculateCost(vars4MCMC%stem_mass_shrub_y%mdData(:,4), vars4MCMC%stem_mass_shrub_y%obsData(:,4),&
                 vars4MCMC%stem_mass_shrub_y%obsData(:,5), J_cost)
            J_new = J_new + J_cost/8000
        endif
        ! print*, "J_new12: ", J_new

        ! leaf_resp_shrub_d 
        if(vars4MCMC%leaf_resp_shrub_d%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_resp_shrub_d%mdData(:,4), vars4MCMC%leaf_resp_shrub_d%obsData(:,4),&
                 vars4MCMC%leaf_resp_shrub_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new13: ", J_new

        ! leaf_resp_tree_d 
        if(vars4MCMC%leaf_resp_tree_d%existOrNot)then
            call CalculateCost(vars4MCMC%leaf_resp_tree_d%mdData(:,4), vars4MCMC%leaf_resp_tree_d%obsData(:,4),&
                 vars4MCMC%leaf_resp_tree_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new14: ", J_new
        ! ch4_d 
        if(vars4MCMC%ch4_d%existOrNot)then
            call CalculateCost(vars4MCMC%ch4_d%mdData(:,4), vars4MCMC%ch4_d%obsData(:,4),&
                 vars4MCMC%ch4_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost/80000
        endif
        ! print*, "J_new15: ", J_new
        ! ch4_h
        if(vars4MCMC%ch4_h%existOrNot)then
            call CalculateCost(vars4MCMC%ch4_h%mdData(:,4), vars4MCMC%ch4_h%obsData(:,4),&
                 vars4MCMC%ch4_h%obsData(:,5), J_cost)
            J_new = J_new + J_cost/400
        endif
        ! print*, "J_new16: ", J_new
        
        ! CN_shag_d 
        if(vars4MCMC%CN_shag_d%existOrNot)then
            call CalculateCost(vars4MCMC%CN_shag_d%mdData(:,4), vars4MCMC%CN_shag_d%obsData(:,4),&
                 vars4MCMC%CN_shag_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost/10
        endif
        ! print*, "J_new17: ", J_new

        ! photo_shrub_d 
        if(vars4MCMC%photo_shrub_d%existOrNot)then
            call CalculateCost(vars4MCMC%photo_shrub_d%mdData(:,4), vars4MCMC%photo_shrub_d%obsData(:,4),&
                 vars4MCMC%photo_shrub_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost/10
        endif
        ! print*, "J_new18: ", J_new

        ! photo_tree_d 
        if(vars4MCMC%photo_tree_d%existOrNot)then
            call CalculateCost(vars4MCMC%photo_tree_d%mdData(:,4), vars4MCMC%photo_tree_d%obsData(:,4),&
                 vars4MCMC%photo_tree_d%obsData(:,5), J_cost)
            J_new = J_new + J_cost
        endif
        ! print*, "J_new19: ", J_new
        ! ------------------------------------------------------------------------------------
        ! write(*,*) "here2",J_new
        if(J_new .eq. 0) then ! no data is available
            delta_J = -0.1
        else
            delta_J = J_new - J_last
        endif

        delta_J = delta_J/10

        call random_number(cs_rand)
        if(AMIN1(1.0, exp(-delta_J)) .gt. cs_rand)then
            upgraded = upgraded + 1
            J_last = J_new
        endif
        accept_rate = real(upgraded)/real(iDAsimu)
    end subroutine costFuncObs

    subroutine CalculateCost(datMod4MCMC, datObs4MCMC, stdObs4MCMC, JCost)
        ! calculate the cost of the observation and simulation, and update the number of updated
        implicit none
        real, intent(in) :: datMod4MCMC(:), datObs4MCMC(:), stdObs4MCMC(:)
        integer nLine, iLine, nCost
        real JCost, dObsSimu, std4cal

        nLine = size(datObs4MCMC)
        nCost = 0
        JCost = 0.

        do iLine = 1, nLine
            if(datObs4MCMC(iLine) .gt. -999 .and. datMod4MCMC(iLine) .gt. -999)then
                nCost    = nCost + 1   
                dObsSimu = datMod4MCMC(iLine) - datObs4MCMC(iLine) 
                if (stdObs4MCMC(iLine) < 0.001) then 
                    std4cal = 0.5
                else
                    std4cal = stdObs4MCMC(iLine)
                endif
                JCost    = JCost + (dObsSimu*dObsSimu)/(2*std4cal)
            endif
        enddo
        if(nCost .gt. 0) JCost=JCost/real(nCost)
        return ! JCost
    end subroutine CalculateCost

    subroutine racine_mat(M, Mrac,npara)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Square root of a matrix							  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,i, nrot
        real M(npara,npara),Mrac(npara,npara)
        real valpr(npara),vectpr(npara,npara)
        Mrac=0.
        call jacobi(M,npara,npara,valpr,vectpr,nrot)
        do i=1,npara
            if(valpr(i).ge.0.) then
                Mrac(i,i)=sqrt(valpr(i))
            else
                print*, 'WARNING!!! Square root of the matrix is undefined.'
                print*, ' A negative eigenvalue has been set to zero - results may be wrong'
                Mrac=M
                return
            endif
        enddo
        Mrac=matmul(matmul(vectpr, Mrac),transpose(vectpr))

    end subroutine racine_mat

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Extraction of the eigenvalues and the eigenvectors !!
    !! of a matrix (Numerical Recipes)					  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE jacobi(a,n,np,d,v,nrot)
        INTEGER :: n,np,nrot
        REAL :: a(np,np),d(np),v(np,np)
        INTEGER, PARAMETER :: NMAX=500
        INTEGER :: i,ip,iq,j
        REAL :: c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
        
        do ip=1,n
            do iq=1,n
                v(ip,iq)=0.
            end do
            v(ip,ip)=1.
        end do
        
        do ip=1,n
            b(ip)=a(ip,ip)
            d(ip)=b(ip)
            z(ip)=0.
        end do
        
        nrot=0
        do i=1,50
            sm=0.
            do ip=1,n-1
                do iq=ip+1,n
                    sm=sm+abs(a(ip,iq))
                end do
            end do
            if(sm.eq.0.)return
            if(i.lt.4)then
                tresh=0.2*sm/n**2
            else
                tresh=0.
            endif
            do ip=1,n-1
                do iq=ip+1,n
                    g=100.*abs(a(ip,iq))
                    if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
                        a(ip,iq)=0.
                    else if(abs(a(ip,iq)).gt.tresh)then
                        h=d(iq)-d(ip)
                        if(abs(h)+g.eq.abs(h))then
                            t=a(ip,iq)/h
                        else
                            theta=0.5*h/a(ip,iq)
                            t=1./(abs(theta)+sqrt(1.+theta**2))
                            if(theta.lt.0.) then
                                t=-t
                            endif
                        endif
                        c=1./sqrt(1+t**2)
                        s=t*c
                        tau=s/(1.+c)
                        h=t*a(ip,iq)
                        z(ip)=z(ip)-h
                        z(iq)=z(iq)+h
                        d(ip)=d(ip)-h
                        d(iq)=d(iq)+h
                        a(ip,iq)=0.
                        do j=1,ip-1
                            g=a(j,ip)
                            h=a(j,iq)
                            a(j,ip)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        end do
                        do j=ip+1,iq-1
                            g=a(ip,j)
                            h=a(j,iq)
                            a(ip,j)=g-s*(h+g*tau)
                            a(j,iq)=h+s*(g-h*tau)
                        end do
                        do j=iq+1,n
                            g=a(ip,j)
                            h=a(iq,j)
                            a(ip,j)=g-s*(h+g*tau)
                            a(iq,j)=h+s*(g-h*tau)
                        end do
                        do j=1,n
                            g=v(j,ip)
                            h=v(j,iq)
                            v(j,ip)=g-s*(h+g*tau)
                            v(j,iq)=h+s*(g-h*tau)
                        end do
                        nrot=nrot+1
                    endif
                end do
            end do
            do ip=1,n
                b(ip)=b(ip)+z(ip)
                d(ip)=b(ip)
                z(ip)=0.
            end do
        end do
        print*, 'too many iterations in jacobi'
        return
    END subroutine jacobi

    subroutine gengaussvect(gamma_racine,xold,xnew,npara)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Generation of a random vector from a multivariate  !!
    !! normal distribution with mean zero and covariance  !!
    !! matrix gamma.									  !!
    !! Beware!!! In order to improve the speed of the	  !!
    !! algorithms, the subroutine use the Square root	  !!
    !! matrix of gamma									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara, i
        real gamma_racine(npara,npara)
        real x(npara),xold(npara),xnew(npara)
        
        do i=1,npara
            x(i)=rangauss(25)
        enddo
        
        x = matmul(gamma_racine, x)
        xnew = xold + x
    end subroutine gengaussvect

   real function rangauss(idum)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Generation of a random number from a standard	  !!
    !! normal distribution. (Numerical Recipes)           !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer idum
        real v1, v2, r, fac, gset
        real r_num
        integer :: iset
        
        ! data iset/0/
        iset = 0
        if(iset==0) then
1	        CALL random_number(r_num)
            v1=2.*r_num-1
            CALL random_number(r_num)
            v2=2.*r_num-1
            r=(v1)**2+(v2)**2
            if(r>=1) go to 1
            fac=sqrt(-2.*log(r)/r)
            gset=v1*fac
            rangauss=v2*fac
            iset=1
        else
            rangauss=gset
            iset=0
        end if
        return
    end function

    subroutine varcov(tab,varcovar,npara,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! variance matrix of a matrix of data				  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,ncov
        real tab(ncov,npara),tab2(ncov,npara)
        real varcovar(npara,npara)
        
        call centre(tab,tab2,npara,ncov)
        
        varcovar = matmul(transpose(tab2), tab2)*(1./real(ncov))
        
    end subroutine varcov



    subroutine centre(mat,mat_out,npara,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Compute the centered matrix, ie. the matrix minus  !!
    !! the column means									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer npara,i,ncov
        real mat(ncov,npara),mat_out(ncov,npara)
        ! real mean

        do i=1,npara
            mat_out(:,i) = mat(:,i) - mean(mat(:,i),ncov)
        enddo

    end subroutine centre

    real function mean(tab,ncov)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! mean of a vector									  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer ncov, incov
        real tab(ncov)
        real mean_tt
        mean_tt=0.
        do incov=1,ncov
        mean_tt=mean_tt+tab(incov)/real(ncov)
        enddo
        mean=mean_tt
    End Function mean

    ! real function mean(tab, ncov)
    !     integer ncov
    !     real tab(ncov)
    ! end function mean
    
    subroutine deallocate_mcmc()
    ! deallocate some variables and summary the information of MCMC
        implicit none
        integer :: ipft, npft
        
        

        if(allocated(mc_parvals)) then
            npft = size(mc_parvals)
            do ipft = 1, npft
                if(allocated(mc_parvals(ipft)%parval)) deallocate(mc_parvals(ipft)%parval)
                if(allocated(mc_parvals(ipft)%parmin)) deallocate(mc_parvals(ipft)%parmin)
                if(allocated(mc_parvals(ipft)%parmax)) deallocate(mc_parvals(ipft)%parmax)
            enddo
            deallocate(mc_parvals)
        endif

        ! if(allocated(MDparval)) deallocate(MDparval)
        if(allocated(mc_DApar)) then
            npft = size(mc_DApar)
            do ipft = 1, npft
                if(allocated(mc_DApar(ipft)%DAparmin))  deallocate(mc_DApar(ipft)%DAparmin)
                if(allocated(mc_DApar(ipft)%DAparmax))  deallocate(mc_DApar(ipft)%DAparmax)
                if(allocated(mc_DApar(ipft)%DAparidx))  deallocate(mc_DApar(ipft)%DAparidx)
                if(allocated(mc_DApar(ipft)%DApar))     deallocate(mc_DApar(ipft)%DApar)
                if(allocated(mc_DApar(ipft)%DApar_old)) deallocate(mc_DApar(ipft)%DApar_old)

                if(allocated(mc_DApar(ipft)%coefhistory)) deallocate(mc_DApar(ipft)%coefhistory)
                if(allocated(mc_DApar(ipft)%coefnorm))    deallocate(mc_DApar(ipft)%coefnorm)
                if(allocated(mc_DApar(ipft)%coefac))      deallocate(mc_DApar(ipft)%coefac)
        
                if(allocated(mc_DApar(ipft)%gamnew))   deallocate(mc_DApar(ipft)%gamnew)
            enddo
            deallocate(mc_DApar)
        endif
        

        if(allocated(vars4MCMC%ANPP_Shrub_y%obsData))  deallocate(vars4MCMC%ANPP_Shrub_y%obsData)
        if(allocated(vars4MCMC%ANPP_Tree_y%obsData))  deallocate(vars4MCMC%ANPP_Tree_y%obsData)
        if(allocated(vars4MCMC%NPP_sphag_y%obsData))  deallocate(vars4MCMC%NPP_sphag_y%obsData)
        if(allocated(vars4MCMC%BNPP_y%obsData))  deallocate(vars4MCMC%BNPP_y%obsData)        ! tree + shrub
        if(allocated(vars4MCMC%er_d%obsData))  deallocate(vars4MCMC%er_d%obsData)          ! shrub + sphag.
        if(allocated(vars4MCMC%er_h%obsData))  deallocate(vars4MCMC%er_h%obsData)          ! shrub + sphag.
        if(allocated(vars4MCMC%gpp_d%obsData))  deallocate(vars4MCMC%gpp_d%obsData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_d%obsData))  deallocate(vars4MCMC%nee_d%obsData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_h%obsData))  deallocate(vars4MCMC%nee_h%obsData)         ! shrub + sphag.
        if(allocated(vars4MCMC%LAI_d%obsData))  deallocate(vars4MCMC%LAI_d%obsData)         ! tree  + Shrub
        !
        if(allocated(vars4MCMC%leaf_mass_shrub_y%obsData))  deallocate(vars4MCMC%leaf_mass_shrub_y%obsData)
        if(allocated(vars4MCMC%stem_mass_shrub_y%obsData))  deallocate(vars4MCMC%stem_mass_shrub_y%obsData)
        if(allocated(vars4MCMC%leaf_resp_shrub_d%obsData))  deallocate(vars4MCMC%leaf_resp_shrub_d%obsData) 
        if(allocated(vars4MCMC%leaf_resp_tree_d%obsData))  deallocate(vars4MCMC%leaf_resp_tree_d%obsData) 
        ! methane
        if(allocated(vars4MCMC%ch4_d%obsData))  deallocate(vars4MCMC%ch4_d%obsData) 
        if(allocated(vars4MCMC%ch4_h%obsData))  deallocate(vars4MCMC%ch4_h%obsData) 
        ! 
        if(allocated(vars4MCMC%CN_shag_d%obsData))  deallocate(vars4MCMC%CN_shag_d%obsData) 
        if(allocated(vars4MCMC%photo_shrub_d%obsData))  deallocate(vars4MCMC%photo_shrub_d%obsData) 
        if(allocated(vars4MCMC%photo_tree_d%obsData))  deallocate(vars4MCMC%photo_tree_d%obsData) 
        ! ---------------------------------------------------------------------------------------------

        if(allocated(vars4MCMC%ANPP_Shrub_y%mdData))  deallocate(vars4MCMC%ANPP_Shrub_y%mdData)
        if(allocated(vars4MCMC%ANPP_Tree_y%mdData))  deallocate(vars4MCMC%ANPP_Tree_y%mdData)
        if(allocated(vars4MCMC%NPP_sphag_y%mdData))  deallocate(vars4MCMC%NPP_sphag_y%mdData)
        if(allocated(vars4MCMC%BNPP_y%mdData))  deallocate(vars4MCMC%BNPP_y%mdData)        ! tree + shrub
        if(allocated(vars4MCMC%er_d%mdData))  deallocate(vars4MCMC%er_d%mdData)          ! shrub + sphag.
        if(allocated(vars4MCMC%er_h%mdData))  deallocate(vars4MCMC%er_h%mdData)          ! shrub + sphag.
        if(allocated(vars4MCMC%gpp_d%mdData))  deallocate(vars4MCMC%gpp_d%mdData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_d%mdData))  deallocate(vars4MCMC%nee_d%mdData)         ! Shrub + sphag.
        if(allocated(vars4MCMC%nee_h%mdData))  deallocate(vars4MCMC%nee_h%mdData)         ! shrub + sphag.
        if(allocated(vars4MCMC%LAI_d%mdData))  deallocate(vars4MCMC%LAI_d%mdData)         ! tree  + Shrub
        !
        if(allocated(vars4MCMC%leaf_mass_shrub_y%mdData))  deallocate(vars4MCMC%leaf_mass_shrub_y%mdData)
        if(allocated(vars4MCMC%stem_mass_shrub_y%mdData))  deallocate(vars4MCMC%stem_mass_shrub_y%mdData)
        if(allocated(vars4MCMC%leaf_resp_shrub_d%mdData))  deallocate(vars4MCMC%leaf_resp_shrub_d%mdData) 
        if(allocated(vars4MCMC%leaf_resp_tree_d%mdData))  deallocate(vars4MCMC%leaf_resp_tree_d%mdData) 
        ! methane
        if(allocated(vars4MCMC%ch4_d%mdData))  deallocate(vars4MCMC%ch4_d%mdData) 
        if(allocated(vars4MCMC%ch4_h%mdData))  deallocate(vars4MCMC%ch4_h%mdData) 
        ! 
        if(allocated(vars4MCMC%CN_shag_d%mdData))  deallocate(vars4MCMC%CN_shag_d%mdData) 
        if(allocated(vars4MCMC%photo_shrub_d%mdData))  deallocate(vars4MCMC%photo_shrub_d%mdData) 
        if(allocated(vars4MCMC%photo_tree_d%mdData))  deallocate(vars4MCMC%photo_tree_d%mdData) 
        ! ---------------------------------------------------------------------------------------------

        if(allocated(parnames)) deallocate(parnames)

        ! in MCMC_outputs module
        do ipft = 1, npft
            if(allocated(arr_params_set(ipft)%tot_paramsets)) deallocate(arr_params_set(ipft)%tot_paramsets)
            if(allocated(arr_params_set(ipft)%sel_paramsets)) deallocate(arr_params_set(ipft)%sel_paramsets)
        enddo
        if(allocated(arr_params_set)) deallocate(arr_params_set)

        ! if (do_mc_out_hr)then
        !     call deallocate_mcmc_outs_type(sel_paramsets_outs_h)
        !     call deallocate_mcmc_outs_type(tot_paramsets_outs_h)
        ! endif
        if (do_mc_out_day)then
            call deallocate_mcmc_outs_type(sel_paramsets_outs_d)
            call deallocate_mcmc_outs_type(tot_paramsets_outs_d)
        endif
        ! if (do_mc_out_mon)then
        !     call deallocate_mcmc_outs_type(sel_paramsets_outs_m)
        !     call deallocate_mcmc_outs_type(tot_paramsets_outs_m)
        ! endif
    end subroutine deallocate_mcmc
end module mcmc
