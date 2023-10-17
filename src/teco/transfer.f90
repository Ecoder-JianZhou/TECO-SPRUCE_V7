module transfer
    use datatypes
    implicit none
    real S_omega    !  average values of the moisture scaling functions
    real S_t(5)     !  average values of temperature scaling functions

    contains
    subroutine TCS_CN(vegn, iforcing)    
        ! carbon transfer according to Xu et al. 2007 
        implicit none
        type(vegn_tile_type), intent(inout) :: vegn
        type(forcing_data_type), intent(in) :: iforcing

        ! real Q_plant, OutC(8), OutN(8)  
        real Q10h(5)                              
        real CNmin,CNmax,NSNmax,NSNmin              !,NSN
        real CN_foliage
        real N_immob,N_imm(5),Nfix0        ! ,N_deficit,N_fixation
        ! real N_transfer!,N_loss                      ! ,N_miner,N_uptake,N_deposit,N_leach,N_vol
        real Qroot0,Cfix0                           ! alphaN,
        real Scalar_N_flow,Scalar_N_T
        real ksye                                   ! NSC,fnsc,
        ! real SNvcmax,SNgrowth,SNRauto,SNrs
        real kappaVcmax
        real SNfine,SNcoarse,SNmicr,SNslow,SNpass
        real costCuptake,costCfix,costCreuse         ! Rnitrogen,
        real Creuse0,Nup0,N_deN0,LDON0
        real S_w_min    !  minimum decomposition rate at zero available water
        ! For test
        real ScNloss
        ! added for soil thermal
        real frac_soc(10)
        integer i, j, ipft
        ! real :: temp2
        
        frac_soc    = (/0.75,0.2,0.02,0.015,0.005,0.0,0.0,0.0,0.0,0.0/)
        ! temperature sensitivity of Rh_pools
        Q10h        = (/st%Q10rh,st%Q10rh,st%Q10rh,st%Q10rh,st%Q10rh/) 

        Qroot0      = 500.
        Nfix0       = 1./60.   ! maximum N fix ratio, N/C
        Nup0        = 0.02     ! nitrogen uptake rate
        Cfix0       = 12.      ! C cost per N for fixation
        ksye        = 0.05      ! C cost per N for uptake
        Creuse0     = 2.     ! C cost per N for resorption
        ScNloss     = 1.
        N_deN0      = 1.E-3*ScNloss   ! 1.E-3, 5.E-3, 10.E-3, 20.E-3
        LDON0       = 1.E-3*ScNloss
        st%Rnitrogen   = 0.
        ! for N scalars
        CNmin       = 40.0
        CNmax       = 200.0
        ! calculating soil scaling factors, S_omega and S_tmperature
        S_w_min = 0.08 !minimum decomposition rate at zero soil moisture
        S_omega = S_w_min + (1.-S_w_min) * Amin1(1.0, 0.3*st%omega)
        if (do_soilphy) then 
            S_t=(/0.0,0.0,0.0,0.0,0.0/)
            do i=1,5
                if(i.lt.3) then    ! couarse and fine litter use surface layer soil temperature
                    ! S_t(i)=Q10h(i)**((Tsoil-5.)/10.)  ! Oak
                    S_t(i)=Q10h(i)**((st%tsoil_layer(2)-10.)/10.)  ! Duke
                else 
                    do j=1,10       ! fast,slow and passive pool use weighed soil temperature in layers according to soc distribution
                        S_t(i)=S_t(i)+frac_soc(j)*Q10h(i)**((st%tsoil_layer(j+1)-10.)/10.)  ! Duke
                    enddo
                endif
            enddo
        else
            do i=1,5
                ! S_t(i)=Q10h(i)**((Tsoil-5.)/10.)  ! Oak
                S_t(i)=Q10h(i)**((iforcing%Tsoil-10.)/10.)  ! Duke
            enddo  
        endif 

        do ipft = 1, vegn%npft
            ! Max and min NSN pool
            NSNmax = vegn%allSp(ipft)%QN(1) + 0.2*vegn%allSp(ipft)%QN(2) + vegn%allSp(ipft)%QN(3)  ! 15.0
            NSNmin = 0.01

            ! Calculating NPP allocation and changes of each C pool
            ! print*, "NPP_L: ", vegn%allSp(ipft)%alpha_L, vegn%allSp(ipft)%NPP
            vegn%allSp(ipft)%NPP_L = vegn%allSp(ipft)%alpha_L * vegn%allSp(ipft)%NPP           ! NPP allocation
            vegn%allSp(ipft)%NPP_W = vegn%allSp(ipft)%alpha_W * vegn%allSp(ipft)%NPP
            vegn%allSp(ipft)%NPP_R = vegn%allSp(ipft)%alpha_R * vegn%allSp(ipft)%NPP
            ! N scalars on decomposition
            SNfine   = exp(-(vegn%allSp(ipft)%CN0(4)-vegn%allSp(ipft)%CN(4))/vegn%allSp(ipft)%CN0(4)) 
            SNcoarse = exp(-(vegn%allSp(ipft)%CN0(5)-vegn%allSp(ipft)%CN(5))/vegn%allSp(ipft)%CN0(5)) 
            ! the carbon leaving the pools
            ! print *, "outc:", vegn%allSp(ipft)%L_fall, vegn%allSp(ipft)%QC(2), vegn%allSp(ipft)%tauC(2)*S_omega
            vegn%allSp(ipft)%OutC(1)  = vegn%allSp(ipft)%L_fall
            vegn%allSp(ipft)%OutC(2)  = vegn%allSp(ipft)%QC(2)/vegn%allSp(ipft)%tauC(2)*S_omega !*exp(CN(2)/CN0(2)-1.) 
            vegn%allSp(ipft)%OutC(3)  = vegn%allSp(ipft)%QC(3)/vegn%allSp(ipft)%tauC(3)*S_omega
            ! vegn%allSp(ipft)%OutC(4)  = vegn%allSp(ipft)%QC(4)/vegn%allSp(ipft)%tauC(4)*S_omega* &
            !             S_T(1)*vegn%allSp(ipft)%CN(4)/vegn%allSp(ipft)%CN0(4)!*SNfine
            ! vegn%allSp(ipft)%OutC(5)  = vegn%allSp(ipft)%QC(5)/vegn%allSp(ipft)%tauC(5)*S_omega* & 
            !             S_T(2)*vegn%allSp(ipft)%CN(5)/vegn%allSp(ipft)%CN0(5)!*SNcoarse
            ! heterotrophic respiration from litter
            ! vegn%allSp(ipft)%Rh_pools(1) = OutC(4)* (1. - vegn%allSp(ipft)%f_F2M)
            ! vegn%allSp(ipft)%Rh_pools(2) = OutC(5)* (1. - vegn%allSp(ipft)%f_C2M - vegn%allSp(ipft)%f_C2S)

            do i=1,3
                vegn%allSp(ipft)%OutN(i) = vegn%allSp(ipft)%OutC(i)/vegn%allSp(ipft)%CN(i)
            enddo 
            vegn%allSp(ipft)%N_demand   = vegn%allSp(ipft)%NPP_L/vegn%allSp(ipft)%CN0(1) + &
                             vegn%allSp(ipft)%NPP_W/vegn%allSp(ipft)%CN0(2) + &
                             vegn%allSp(ipft)%NPP_R/vegn%allSp(ipft)%CN0(3) 
            ! summary
            if (ipft .eq. 1) then
                do i=1,3
                    st%OutC(i) = vegn%allSp(ipft)%OutC(i)
                    st%OutN(i) = vegn%allSp(ipft)%OutN(i)
                enddo
                st%NPP_L   = vegn%allSp(ipft)%NPP_L
                st%NPP_W   = vegn%allSp(ipft)%NPP_W
                st%NPP_R   = vegn%allSp(ipft)%NPP_R
                ! st%Rh_pools(1) = vegn%allSp(ipft)%Rh_pools(1)
                ! st%Rh_pools(2) = vegn%allSp(ipft)%Rh_pools(2)
                ! vegn%allSp(ipft)%N_demand   = vegn%allSp(ipft)%NPP_L/vegn%allSp(ipft)%CN0(1) + &
                !              vegn%allSp(ipft)%NPP_W/vegn%allSp(ipft)%CN0(2) + &
                !              vegn%allSp(ipft)%NPP_R/vegn%allSp(ipft)%CN0(3) 
            else
                do i = 1,3
                    st%OutC(i) = st%OutC(i) + vegn%allSp(ipft)%OutC(i)
                    st%OutN(i) = st%OutN(i) + vegn%allSp(ipft)%OutN(i)
                enddo
                st%NPP_L   = st%NPP_L + vegn%allSp(ipft)%NPP_L
                st%NPP_W   = st%NPP_L + vegn%allSp(ipft)%NPP_W
                st%NPP_R   = st%NPP_L + vegn%allSp(ipft)%NPP_R
                ! st%Rh_pools(1) = st%Rh_pools(1) + vegn%allSp(ipft)%Rh_pools(1)
                ! st%Rh_pools(2) = st%Rh_pools(2) + vegn%allSp(ipft)%Rh_pools(2)
                ! vegn%allSp(ipft)%N_demand   = vegn%allSp(ipft)%N_demand +                                       &
                !              vegn%allSp(ipft)%NPP_L/vegn%allSp(ipft)%CN0(1) + &
                !              vegn%allSp(ipft)%NPP_W/vegn%allSp(ipft)%CN0(2) + &
                !              vegn%allSp(ipft)%NPP_R/vegn%allSp(ipft)%CN0(3) 
            endif
        enddo

        ! N scalars on decomposition    
        SNmicr   = exp(-(st%CN0(6)-st%CN(6))/st%CN0(6)) 
        SNslow   = 1. !exp(-(CN0(7)-CNC(7))/CN0(7)) 
        SNpass   = exp(-(st%CN0(8)-st%CN(8))/st%CN0(8)) 
        ! the carbon leaving the pools
        ! OutC(1)  = L_fall
        ! OutC(2)  = QC(2)/tauC(2)*S_omega !*exp(CN(2)/CN0(2)-1.) 
        ! OutC(3)  = QC(3)/tauC(3)*S_omega
        st%OutC(4)  = st%QC(4)/st%tauC(4)*S_omega* S_T(1)*st%CN(4)/st%CN0(4)!*SNfine
        st%OutC(5)  = st%QC(5)/st%tauC(5)*S_omega* S_T(2)*st%CN(5)/st%CN0(5)!*SNcoarse
        st%OutC(6)  = st%QC(6)/st%tauC(6)*S_omega* S_T(3)!*SNmicr
        st%OutC(7)  = st%QC(7)/st%tauC(7)*S_omega* S_T(4)!*SNslow
        st%OutC(8)  = st%QC(8)/st%tauC(8)*S_omega* S_T(5)!*SNpass
    
        ! ! heterotrophic respiration from each pool
        st%Rh_pools(1) = st%OutC(4)* (1. - st%f_F2M)
        st%Rh_pools(2) = st%OutC(5)* (1. - st%f_C2M - st%f_C2S)
        st%Rh_pools(3) = st%OutC(6)* (1. - st%f_M2S - st%f_M2P)
        st%Rh_pools(4) = st%OutC(7)* (1. - st%f_S2P - st%f_S2M)
        st%Rh_pools(5) = st%OutC(8)* (1. - st%f_P2M)
        !========================================================================
        ! Nitrogen part
        ! nitrogen leaving the pools and resorption
        st%OutN(4) = st%OutC(4)/st%CN(4)
        st%OutN(5) = st%OutC(5)/st%CN(5) 
        st%OutN(6) = st%OutC(6)/st%CN(6)
        st%OutN(7) = st%OutC(7)/st%CN(7)
        st%OutN(8) = st%OutC(8)/st%CN(8)
        ! nitrogen mineralization
        st%N_miner = st%OutN(4)* (1. - st%f_F2M)  &
            &      + st%OutN(5)* (1. - st%f_C2M - st%f_C2S) &
            &      + st%OutN(6)* (1. - st%f_M2S - st%f_M2P) &
            &      + st%OutN(7)* (1. - st%f_S2P - st%f_S2M) &
            &      + st%OutN(8)* (1. - st%f_P2M)

        ! Nitrogen immobilization
        N_imm   = 0.
        N_immob = 0.
        if(st%QNminer>0)then
            do i=4,8
                if(st%CN(i)>st%CN0(i))then
                    N_imm(i-3) = Amin1(st%QC(i)/st%CN0(i)-st%QC(i)/st%CN(i),0.1*st%QNminer)           
                    N_immob    = N_immob+N_imm(i-3)
                endif
            enddo
        endif

        ! Let plant itself choose the strategy between using C to uptake
        ! or fix N2 by comparing C invest.
        ! ! N demand
        ! N_demand    = st%NPP_L/st%CN0(1)+st%NPP_W/st%CN0(2)+st%NPP_R/st%CN0(3) !+N_deficit
        ! Nitrogen input:
        st%N_transfer  = 0.
        st%N_uptake    = 0.
        st%N_fixation  = 0.
        do ipft = 1, vegn%npft
            vegn%allSp(ipft)%N_transfer = 0.
            vegn%allSp(ipft)%N_uptake   = 0.
            vegn%allSp(ipft)%N_fixation = 0.
            costCuptake    = 0.
            costCfix       = 0.
            costCreuse     = 0.
            ! vegn%allSp(ipft)%N_demand       = vegn%allSp(ipft)%NPP_L/vegn%allSp(ipft)%CN0(1)
            
            ! 1. Nitrogen resorption
            vegn%allSp(ipft)%N_transfer  = Amax1((vegn%allSp(ipft)%OutN(1) + vegn%allSp(ipft)%OutN(2) + &
                                                  vegn%allSp(ipft)%OutN(3))*vegn%allSp(ipft)%alphaN, 0.)
            costCreuse     = Amax1(Creuse0*vegn%allSp(ipft)%N_transfer, 0.)
            vegn%allSp(ipft)%N_demand       = vegn%allSp(ipft)%N_demand-vegn%allSp(ipft)%N_transfer
            If(vegn%allSp(ipft)%N_demand>0.0)then
                ! 2.  N uptake
                if(ksye/st%QNminer<Cfix0)then
                    vegn%allSp(ipft)%N_uptake = Amax1(AMIN1(vegn%allSp(ipft)%N_demand + vegn%allSp(ipft)%N_deficit,    &
                            &     st%QNminer*vegn%allSp(ipft)%QC(3)/(vegn%allSp(ipft)%QC(3)+Qroot0),  &
                            &     Nup0*vegn%allSp(ipft)%NSC/(ksye/st%QNminer)), 0.) 
                    costCuptake = Amax1(vegn%allSp(ipft)%N_uptake*(ksye/st%QNminer),0.)
                    vegn%allSp(ipft)%N_demand    = vegn%allSp(ipft)%N_demand-vegn%allSp(ipft)%N_uptake
                elseif(vegn%allSp(ipft)%NSN<24.*30.*vegn%allSp(ipft)%N_demand)then
                ! 3.  Nitrogen fixation
                    vegn%allSp(ipft)%N_fixation = Amax1(Amin1(vegn%allSp(ipft)%N_demand, &
                        vegn%allSp(ipft)%fnsc*Nfix0*vegn%allSp(ipft)%NSC), 0.)
                    costCfix      = Amax1(Cfix0*vegn%allSp(ipft)%N_fixation,0.)
                    vegn%allSp(ipft)%N_demand      = vegn%allSp(ipft)%N_demand-vegn%allSp(ipft)%N_fixation
                endif
            endif
            ! print *,"test_N_demand: ", vegn%allSp(ipft)%N_demand, vegn%allSp(ipft)%NPP_L,vegn%allSp(ipft)%CN0(1)
            vegn%allSp(ipft)%N_deficit = vegn%allSp(ipft)%N_deficit + vegn%allSp(ipft)%N_demand
            ! update NSN
            vegn%allSp(ipft)%NSN = vegn%allSp(ipft)%NSN      + vegn%allSp(ipft)%N_transfer + &
                                   vegn%allSp(ipft)%N_uptake + vegn%allSp(ipft)%N_fixation
            ! Total C cost for nitrogen
            vegn%allSp(ipft)%Rnitrogen = costCuptake + costCfix + costCreuse

            ! Nitrogen using, non-structural nitrogen pool, NSN
            vegn%allSp(ipft)%N_leaf = AMAX1(AMIN1(vegn%allSp(ipft)%NPP * vegn%allSp(ipft)%alpha_L/vegn%allSp(ipft)%CN(1) + &
                                            vegn%allSp(ipft)%QC(1)/vegn%allSp(ipft)%CN0(1) -                         &
                                            vegn%allSp(ipft)%QC(1)/vegn%allSp(ipft)%CN(1),  0.2*vegn%allSp(ipft)%NSN),0.)
            vegn%allSp(ipft)%N_wood = AMAX1(AMIN1(vegn%allSp(ipft)%NPP*vegn%allSp(ipft)%alpha_W/vegn%allSp(ipft)%CN(2),    &
                                            0.1*vegn%allSp(ipft)%NSN), 0.)
            vegn%allSp(ipft)%N_root = AMAX1(AMIN1(vegn%allSp(ipft)%NPP*vegn%allSp(ipft)%alpha_R/vegn%allSp(ipft)%CN(3) +   &
                                            vegn%allSp(ipft)%QC(3)/vegn%allSp(ipft)%CN0(3) -                         &
                                            vegn%allSp(ipft)%QC(3)/vegn%allSp(ipft)%CN(3),0.2*vegn%allSp(ipft)%NSN),0.)
            vegn%allSp(ipft)%NSN    = vegn%allSp(ipft)%NSN- &
                                      (vegn%allSp(ipft)%N_leaf+vegn%allSp(ipft)%N_wood+vegn%allSp(ipft)%N_root)

            if(vegn%allSp(ipft)%NSN < 0.) then
                print*, "vegn%allSp(ipft)%NSN < 0.", vegn%allSp(ipft)%NSN
                stop
            endif
            ! N_LF   = OutN(1)*(1.-alphaN)
            ! N_WF   = OutN(2)*(1.-alphaN)
            ! N_RF   = OutN(3)*(1.-alphaN)
            
            ! summary
            if (ipft .eq. 1) then
                st%N_transfer  = vegn%allSp(ipft)%N_transfer
                st%N_uptake    = vegn%allSp(ipft)%N_uptake
                st%N_fixation  = vegn%allSp(ipft)%N_fixation
            else
                st%N_transfer  = st%N_transfer + vegn%allSp(ipft)%N_transfer
                st%N_uptake    = st%N_uptake   + vegn%allSp(ipft)%N_uptake
                st%N_fixation  = st%N_fixation + vegn%allSp(ipft)%N_fixation
            endif
        enddo

        ! update QNminer
        st%N_immob = N_immob
        st%QNminer = st%QNminer + st%N_miner + st%N_deposit - (st%N_uptake + N_immob)
        ! Loss of mineralized N and dissolved organic N
        Scalar_N_flow = 0.5*st%runoff/st%rdepth
        ! Scalar_N_T=0.005*(Tsoil+273.)/(Tsoil+273+333.)

        ! commented line for soil thermal       
        ! Scalar_N_T=N_deN0*exp((Tsoil-25.)/10.)
        ! added lines for soil thermal
        if (do_soilphy) then 
            Scalar_N_T = 0.0 
            do j=1,10
                Scalar_N_T = Scalar_N_T + frac_soc(j)*N_deN0*exp((st%tsoil_layer(j+1)-25.)/10.)  
            enddo
        else
            Scalar_N_T=N_deN0*exp((iforcing%Tsoil-25.)/10.)
        endif  
        ! -------------------------------------------------------
        st%N_leach = Scalar_N_flow*st%QNminer+Scalar_N_flow*st%QN(6)*LDON0
        st%N_vol   = Scalar_N_T*st%QNminer
        st%N_loss  = st%N_leach + st%N_vol
        
        ! update QNminer
        st%QNminer  = st%QNminer - st%N_loss
        st%fNnetmin = st%N_miner + st%N_deposit - (st%N_uptake+N_immob)-st%N_loss
        ! update plant carbon pools, ! daily change of each pool size
        ! call matrix_struct()
        do ipft = 1, vegn%npft
            ! print *, "QC: ", vegn%allSp(ipft)%QC(1), vegn%allSp(ipft)%OutC(1),  vegn%allSp(ipft)%NPP_L
            vegn%allSp(ipft)%QC(1) = vegn%allSp(ipft)%QC(1) - vegn%allSp(ipft)%OutC(1) + vegn%allSp(ipft)%NPP_L
            vegn%allSp(ipft)%QC(2) = vegn%allSp(ipft)%QC(2) - vegn%allSp(ipft)%OutC(2) + vegn%allSp(ipft)%NPP_W
            vegn%allSp(ipft)%QC(3) = vegn%allSp(ipft)%QC(3) - vegn%allSp(ipft)%OutC(3) + vegn%allSp(ipft)%NPP_R
            ! update nitrogen pools
            vegn%allSp(ipft)%QN(1) = vegn%allSp(ipft)%QN(1) - vegn%allSp(ipft)%OutN(1) + vegn%allSp(ipft)%N_leaf
            vegn%allSp(ipft)%QN(2) = vegn%allSp(ipft)%QN(2) - vegn%allSp(ipft)%OutN(2) + vegn%allSp(ipft)%N_wood
            vegn%allSp(ipft)%QN(3) = vegn%allSp(ipft)%QN(3) - vegn%allSp(ipft)%OutN(3) + vegn%allSp(ipft)%N_root
            ! CN
            vegn%allSp(ipft)%CN(1) = vegn%allSp(ipft)%QC(1)/vegn%allSp(ipft)%QN(1)
            vegn%allSp(ipft)%CN(2) = vegn%allSp(ipft)%QC(2)/vegn%allSp(ipft)%QN(2)
            vegn%allSp(ipft)%CN(3) = vegn%allSp(ipft)%QC(3)/vegn%allSp(ipft)%QN(3)
            ! summary
            if (ipft .eq. 1) then
                st%QC(1) = vegn%allSp(ipft)%QC(1)
                st%QC(2) = vegn%allSp(ipft)%QC(2)
                st%QC(3) = vegn%allSp(ipft)%QC(3)
                ! Nitrogen
                st%QN(1) = vegn%allSp(ipft)%QN(1)
                st%QN(2) = vegn%allSp(ipft)%QN(2)
                st%QN(3) = vegn%allSp(ipft)%QN(3)
                vegn%N_leaf = vegn%allSp(ipft)%N_leaf
                vegn%N_wood = vegn%allSp(ipft)%N_wood
                vegn%N_root = vegn%allSp(ipft)%N_root
            else
                st%QC(1) = st%QC(1) + vegn%allSp(ipft)%QC(1)
                st%QC(2) = st%QC(2) + vegn%allSp(ipft)%QC(2)
                st%QC(3) = st%QC(3) + vegn%allSp(ipft)%QC(3)
                ! Nitrogen
                st%QN(1) = st%QN(1) + vegn%allSp(ipft)%QN(1)
                st%QN(2) = st%QN(2) + vegn%allSp(ipft)%QN(2)
                st%QN(3) = st%QN(3) + vegn%allSp(ipft)%QN(3)
                vegn%N_leaf = vegn%N_leaf + vegn%allSp(ipft)%N_leaf
                vegn%N_wood = vegn%N_wood + vegn%allSp(ipft)%N_wood
                vegn%N_root = vegn%N_root + vegn%allSp(ipft)%N_root
            endif
        enddo
        st%QC(4) = st%QC(4) - st%OutC(4) + st%OutC(1) + st%etaW*st%OutC(2)+st%OutC(3)
        st%QC(5) = st%QC(5) - st%OutC(5) + (1.-st%etaW)*st%OutC(2)
        st%QC(6) = st%QC(6) - st%OutC(6) + st%f_F2M*st%OutC(4)+st%f_C2M*st%OutC(5)     &         
            &   + st%f_S2M*st%OutC(7)+st%f_P2M * st%OutC(8)
        st%QC(7) = st%QC(7) - st%OutC(7)+st%f_C2S*st%OutC(5)+st%f_M2S*st%OutC(6)
        st%QC(8) = st%QC(8) - st%OutC(8)+st%f_M2P*st%OutC(6)+st%f_S2P*st%OutC(7)

        st%QN(4) = st%QN(4) - st%OutN(4) + N_imm(1) + (st%OutN(1) + st%etaW*st%OutN(2) + st%OutN(3))*(1.-st%alphaN)
        st%QN(5) = st%QN(5) - st%OutN(5) + N_imm(2) + (1.-st%etaW)*st%OutN(2)*(1.-st%alphaN)
        st%QN(6) = st%QN(6) - st%OutN(6) + N_imm(3) - Scalar_N_flow*st%QN(6)*LDON0  &
            &    + st%f_F2M*st%OutN(4)+st%f_C2M*st%OutN(5) + st%f_S2M*st%OutN(7)+st%f_P2M*st%OutN(8)
        st%QN(7) = st%QN(7) - st%OutN(7) + N_imm(4) + st%f_C2S*st%OutN(5) + st%f_M2S*st%OutN(6)
        st%QN(8) = st%QN(8) - st%OutN(8) + N_imm(5) + st%f_M2P*st%OutN(6) + st%f_S2P*st%OutN(7)
        ! st%QNplant = st%QN(1) + st%QN(2)+ st%QN(3)
        ! update C/N ratio
        st%CN      = st%QC/st%QN
        CN_foliage = (st%QC(1)+st%QC(3))/(st%QN(1)+st%QN(3))

        ! calculate N related scalars for Duke FACE
        do ipft = 1, vegn%npft
            kappaVcmax                = vegn%allSp(ipft)%CN0(1)/1.
            vegn%allSp(ipft)%SNvcmax  = exp(-kappaVcmax*(vegn%allSp(ipft)%CN(1)-vegn%allSp(ipft)%CN0(1))/vegn%allSp(ipft)%CN0(1)) ! /CN0(1) ! Duke
            ! print*, "snvcmax:",vegn%allSp(ipft)%SNvcmax, kappaVcmax, vegn%allSp(ipft)%CN(1), vegn%allSp(ipft)%CN0(1), &
            ! vegn%allSp(ipft)%QC(1), vegn%allSp(ipft)%QN(1)
            vegn%allSp(ipft)%SNvcmax  = AMAX1(AMIN1(vegn%allSp(ipft)%SNvcmax,1.),0.) 
            ! print*, "snvcmax2:", vegn%allSp(ipft)%SNvcmax
            vegn%allSp(ipft)%SNgrowth = exp(-(vegn%allSp(ipft)%CN(1)-vegn%allSp(ipft)%CN0(1))/vegn%allSp(ipft)%CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.25
            vegn%allSp(ipft)%SNRauto  = exp(-(vegn%allSp(ipft)%CN(1)-vegn%allSp(ipft)%CN0(1))/vegn%allSp(ipft)%CN0(1)) !  AMAX1((CNmax-CN_foliage)/(CNmax-CNmin),0.0)+0.5
            ! SNrs       = 1.
        enddo
        return
    end subroutine TCS_CN

end module transfer