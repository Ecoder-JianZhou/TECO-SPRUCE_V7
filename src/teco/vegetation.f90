module vegetation
   !=================================================================================
   !  main subroutines  :  canopy,  respiration, plantgrowth 
   !             canopy => yrday,   xlayers,     Tsoil_simu
   !            xlayers => Radiso,  goudriaan,   agsean_day,  agsean_ngt
   !         agsean_day => photosyn
   !           photosyn => ciandA
   ! functions:
   !           sinbet, esat, VJtemp, fJQres, EnzK
   !=================================================================================
   use datatypes
   use soil
   implicit none
   
   contains
   subroutine vegn_canopy(vegn, iforcing)
      implicit none
      type(vegn_tile_type), intent(inout) :: vegn
      type(forcing_data_type), intent(in) :: iforcing
      ! local variables
      integer :: doy, hour, ipft
      real    :: lat, radsol, fbeam, coszen, Radabv(2)
      real    :: Aleaf(2), Eleaf(2), Hleaf(2), gbleaf(2), gsleaf(2) ! from xlayers
      real    :: Acan1, Acan2, Ecan1, Ecan2, Rsoilab1, Rsoilab2 !, Esoil !from xlayers
      ! SUMMARY
      real    :: Acanop, Ecanop

      doy    = iforcing%doy
      hour   = iforcing%hour + 1
      radsol = iforcing%PAR
      radsol = AMAX1(radsol,0.01)
      call yrday(doy, hour, lat, radsol, fbeam)  ! calculate beam fraction in incoming solar radiation
      coszen = sinbet(doy, hour, lat)
      if(st%topfws.gt.0.5) then
         rhoS(2)=0.18
      else
         rhoS(2)=0.52-0.68*st%topfws
      endif
      Radabv(1) = 0.5 * radsol
      Radabv(2) = 0.5 * radsol

      do ipft = 1, vegn%npft
         call xlayers(vegn%allSp(ipft), iforcing, fbeam, coszen, Radabv, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf, &
            Acan1, Acan2, Ecan1, Ecan2, Rsoilab1, Rsoilab2)!, Esoil)
         Acanop = Acan1 + Acan2
         Ecanop = Ecan1 + Ecan2
         vegn%allSp(ipft)%gpp    = Acanop*3600.0*12.0                           ! every hour, g C m-2 h-1
         vegn%allSp(ipft)%transp = AMAX1(Ecanop*3600.0/(1.0e6*(2.501 - 0.00236*iforcing%Tair)), 0.) ! mm H2O /hour
         ! print *, "test_transp: ", Ecanop*3600.0/(1.0e6*(2.501 - 0.00236*iforcing%Tair)), Ecanop, Ecan1, Ecan2
         if (ipft .eq. 1) then
            vegn%gpp    = vegn%allSp(ipft)%gpp
            vegn%transp = vegn%allSp(ipft)%transp
            vegn%LAI    = (1./vegn%npft)*vegn%allSp(ipft)%LAI
            st%Esoil    = (1./vegn%npft)*vegn%allSp(ipft)%Esoil
            st%Rsoilab1 = (1./vegn%npft)*Rsoilab1
            st%Rsoilab2 = (1./vegn%npft)*Rsoilab2
            ! st%Rsoilab3 = (1/npft)*Rsoilab3
         else
            vegn%gpp    = vegn%gpp    + vegn%allSp(ipft)%gpp
            vegn%transp = vegn%transp + vegn%allSp(ipft)%transp
            vegn%LAI    = vegn%LAI    + (1./vegn%npft)*vegn%allSp(ipft)%LAI
            st%Esoil    = st%Esoil    + (1./vegn%npft)*vegn%allSp(ipft)%Esoil
            st%Rsoilab1 = st%Rsoilab1 + (1./vegn%npft)*Rsoilab1
            st%Rsoilab2 = st%Rsoilab2 + (1./vegn%npft)*Rsoilab2
         endif
      enddo 
      if (do_soilphy) call Tsoil_simu(vegn, iforcing) ! Jian: need modify
      st%evap   = AMAX1(st%Esoil *3600.0/(1.0e6*(2.501 - 0.00236*iforcing%Tair)), 0.)  ! Jian: need modify
      return
   end subroutine vegn_canopy

   ! -------- autotrophic respiration -----------------------------------------------------------------
   subroutine respiration(vegn, iforcing)
      ! calculate plant and soil respiration by the following equation:
      ! RD=BM*Rd*Q10**((T-25)/10) (Sun et al. 2005. Acta Ecologica Sinica)
      implicit none
      type(vegn_tile_type), intent(inout) :: vegn
      type(forcing_data_type), intent(in) :: iforcing

      integer :: ipft
      real conv 
      conv = 3600.*12./1000000.                  ! converter from "umol C /m2/s" to "gC/m2/hour"
      ! SNRauto =exp(-(CN(1)-CN0(1))/CN0(1))

      do ipft = 1, vegn%npft
         if (vegn%allSp(ipft)%LAI .gt. vegn%allSp(ipft)%LAImin) then
            
            vegn%allSp(ipft)%RmLeaf = vegn%allSp(ipft)%Rl0 * vegn%allSp(ipft)%SNRauto * vegn%allSp(ipft)%bmleaf  * 0.48 * &
               vegn%allSp(ipft)%SLA * 0.1 * vegn%allSp(ipft)%Q10**((iforcing%Tair - 10.)/10.)*vegn%allSp(ipft)%fnsc*conv
            vegn%allSp(ipft)%RmStem = vegn%allSp(ipft)%Rs0 * vegn%allSp(ipft)%SNRauto * vegn%allSp(ipft)%StemSap * 0.001* &
               vegn%allSp(ipft)%Q10**((iforcing%Tair - 25.)/10.)*vegn%allSp(ipft)%fnsc*conv
            vegn%allSp(ipft)%RmRoot = vegn%allSp(ipft)%Rr0 * vegn%allSp(ipft)%SNRauto * vegn%allSp(ipft)%RootSap * 0.001* & 
               vegn%allSp(ipft)%Q10**((iforcing%Tair - 25.)/10.)*vegn%allSp(ipft)%fnsc*conv
         else
            vegn%allSp(ipft)%RmLeaf = 0.3*vegn%allSp(ipft)%gpp
            vegn%allSp(ipft)%RmStem = 0.3*vegn%allSp(ipft)%gpp
            vegn%allSp(ipft)%RmRoot = 0.4*vegn%allSp(ipft)%gpp
         end if
         vegn%allSp(ipft)%Rmain = vegn%allSp(ipft)%Rmleaf + vegn%allSp(ipft)%Rmstem + vegn%allSp(ipft)%Rmroot
         if (vegn%allSp(ipft)%Rmain > 0.0015*vegn%allSp(ipft)%NSC) then             ! If Total autotropic respiration greater than 0.15% of Nonstructure Carbon, rescale.
            vegn%allSp(ipft)%Rmleaf = vegn%allSp(ipft)%Rmleaf/vegn%allSp(ipft)%Rmain*0.0015*vegn%allSp(ipft)%NSC
            vegn%allSp(ipft)%Rmstem = vegn%allSp(ipft)%Rmstem/vegn%allSp(ipft)%Rmain*0.0015*vegn%allSp(ipft)%NSC
            vegn%allSp(ipft)%Rmroot = vegn%allSp(ipft)%Rmstem/vegn%allSp(ipft)%Rmain*0.0015*vegn%allSp(ipft)%NSC
            vegn%allSp(ipft)%Rmain  = vegn%allSp(ipft)%Rmleaf + vegn%allSp(ipft)%Rmstem + vegn%allSp(ipft)%Rmroot
         end if
         ! summary the 
         if (ipft .eq. 1) then
            vegn%RmLeaf = vegn%allSp(ipft)%Rmleaf
            vegn%RmStem = vegn%allSp(ipft)%RmStem
            vegn%RmRoot = vegn%allSp(ipft)%RmRoot
            vegn%Rmain  = vegn%allSp(ipft)%Rmain
         else
            vegn%RmLeaf = vegn%RmLeaf +  vegn%allSp(ipft)%Rmleaf
            vegn%RmStem = vegn%RmStem +  vegn%allSp(ipft)%RmStem
            vegn%RmRoot = vegn%RmRoot +  vegn%allSp(ipft)%RmRoot
            vegn%Rmain  = vegn%Rmain  +  vegn%allSp(ipft)%Rmain
         endif
      enddo
      return
   end subroutine respiration

   ! ------------- plant growth model --------------------------------------------------------------------
   subroutine vegn_plantgrowth(vegn, iforcing)
      implicit none
      type(vegn_tile_type), intent(inout) :: vegn
      type(forcing_data_type), intent(in) :: iforcing

      integer ipft
      
      do ipft = 1, vegn%npft
         call plantgrowth(vegn%allSp(ipft), iforcing)
         ! if (ipft .eq. 1)then   
      enddo
      return
   end subroutine vegn_plantgrowth

   subroutine plantgrowth(spec, iforcing)
      implicit none  ! Jian: RS0 is different from Rs0 from data type
      type(spec_data_type), intent(inout) :: spec
      type(forcing_data_type), intent(in) :: iforcing

      ! real store
      real GrowthP, GrowthL, GrowthR, GrowthS
      real CNP0, GPmax, acP
      real bmL, bmR, bmP, bmS !, StemSap, RootSap
      real Rgroot, Rgleaf, Rgstem !, Rgrowth
      real, save :: addaccu = 0, GrowthLaccu = 0, GrowthSaccu = 0, GrowthRaccu = 0
      real scalT, SL_rs, SR_rs, Slai, phiN   !, scalS, scalN
      real RS, RS_0 ! , RSw
      real gamma_W, gamma_Wmax, gamma_T, gamma_Tmax, gamma_N
      real beta_T, Tcold, Twarm, Topt
      real bW, bT, W
      ! real L_add, NL_fall, NL_add
      ! real alpha_St
      real ht

      Twarm = 35.0
      Tcold = 5.0         ! Tcold=0.0       ! For SPRUCE
      Topt  = 30.
      phiN  = 0.33

      bmL = spec%bmleaf*0.48   ! Carbon
      bmR = spec%bmRoot*0.48
      bmS = spec%bmStem*0.48

      if (bmL .lt. spec%NSC/0.333) bmL = spec%NSC/0.333
      if (bmR .lt. spec%NSC/0.333) bmR = spec%NSC/0.333
      if (bmS .lt. spec%NSC/0.334) bmS = spec%NSC/0.334

      spec%StemSap = spec%SapS*bmS  ! Weng 12/05/2008
      spec%RootSap = spec%SapR*bmR
      if (spec%StemSap .lt. 0.001) spec%StemSap = 0.001
      if (spec%RootSap .lt. 0.001) spec%RootSap = 0.001

      bmP  = bmL + bmR + bmS                                        ! Plant C biomass
      acP  = bmL + spec%StemSap + bmS                               ! Plant available sapwood C
      CNp0 = bmP/(bmL/spec%CN0(1) + bmR/spec%CN0(3) + bmS/spec%CN0(2))            ! Plant CN ratio

      ! hmax = 24.19   ! m
      ! hl0  = 0.00019  ! m2/kg C
      ! LAIMAX0 = 6.
      ! LAIMAX0 = 8.
      ! la0     = 0.2
      ht          = spec%hmax*(1.-exp(-spec%hl0*bmP))                              ! Scaling plant C biomass to height
      spec%LAIMAX = AMAX1(spec%LAIMAX0*(1.-exp(-spec%la0*ht)), spec%LAIMIN + 0.1)  ! Scaling plant height to maximum LAI

      ! Phenology
      if ((st%GDD5 .gt. spec%gddonset) .and. spec%onset .eq. 0 .and. spec%storage .gt. spec%stor_use) then
         spec%onset = 1
      end if
      if ((spec%onset .eq. 1) .and. (spec%storage .gt. spec%stor_use)) then
         if (spec%LAI .lt. spec%LAIMAX) spec%add = spec%stor_use !Amin1(stor_use, 0.004*NSN*CNp0)
         ! if(LAI.lt.LAIMAX)add=stor_use/20.0
         spec%storage = spec%storage - spec%add
         ! add = 0
      else
         spec%add   = 0.0
         spec%onset = 0
      end if
      if (spec%accumulation .lt. (spec%NSCmax + 0.001*spec%RootSap)) then
         spec%store = AMAX1(0., 0.0001*spec%NSC)                        ! 0.5% of nonstructure carbon is stored
      else
         spec%store = 0.0
      end if
      spec%accumulation = spec%accumulation + spec%store

      ! Scalars for plant growth
      ! Sps=Amin1(1.0,3.33*AMAX1(0.0,1.0 - fnsc))
      ! Sps   = Sps*(1.-exp(-phiN*NSN))                                                        ! Sps is not assigned previous, something is wrong. -JJJJJJJJJJJJJJJJJJJJJ
      ! scalS = AMIN1(1.0, 2.*spec%fnsc) ! Ss
      RS_0  = 1.0
      RS    = bmR/bmL
      W     = AMIN1(1.0, 3.333 * st%omega)
      SL_rs = RS/(RS + RS_0*(2.-W))
      SR_rs = (RS_0*(2.-W))/(RS + RS_0*(2.-W))
      Slai  = amin1(1.0, 2.333*(spec%LAIMAX - spec%LAI)/(spec%LAIMAX - spec%LAIMIN))
      scalT = AMAX1(0.0, 1.0 - exp(-(iforcing%Tair - spec%gddonset/10.)/5.0))  !0.5 !
      st%scalW = AMIN1(0.5, AMAX1(0.333, 0.333 + st%omega))
      
      ! Plant growth and allocation, based on LM3V
      GPmax   = (spec%GLmax*bmL + spec%GSmax*spec%StemSap + spec%GRmax*bmR)                        !/acP
      ! Jian: limit the Growth > 0
      GrowthP = Amax1(AMIN1(GPmax*spec%fnsc*scalT*(1.-exp(-spec%NSN)),0.004*spec%NSC,0.004*spec%NSN*CNp0), 0.) !    ! Jian: GrowthP is used to calculate NPP, which means the minimum of group plant needs and 0.004*NSC, 0.004*NSN*CNp0
      GrowthL = MAX(0.0, GrowthP*0.5)                                          ! updated when QC leaf and wood changed due to the change of plot area for tree biomass
      GrowthR = Amax1(MIN(GrowthP*0.4, MAX(0.0, 0.75/st%scalW*bmL - bmR)), 0.)                  ! *c1/(1.+c1+c2)
      GrowthS = MAX(0.0, GrowthP - (GrowthL + GrowthR))                        ! *c2/(1.+c1+c2)

      spec%npp    = GrowthL + GrowthR + GrowthS + spec%add       ! Modified by Jiang Jiang 2015/10/13
      ! write(*,*) "test: NPP1 ", spec%npp, GrowthL, GrowthR, GrowthS, spec%add, spec%stor_use, spec%accumulation
      ! write(*,*) "test: NPP2 ", GPmax*spec%fnsc*scalT*(1.-exp(-spec%NSN)),0.004*spec%NSC,0.004*spec%NSN*CNp0
      ! write(*,*) "test: NPP3 ", GPmax,spec%fnsc,scalT,spec%NSN,spec%NSC,CNp0
      addaccu     = addaccu + spec%add
      GrowthLaccu = GrowthLaccu + GrowthL
      GrowthRaccu = GrowthRaccu + GrowthR
      GrowthSaccu = GrowthSaccu + GrowthS
      
      if (spec%npp .lt. 1E-18) then
         spec%alpha_L = 0.333
         spec%alpha_W = 0.333
         spec%alpha_R = 0.333
      else
         spec%alpha_L = (GrowthL + spec%add)/spec%NPP
         spec%alpha_W = GrowthS/spec%NPP
         spec%alpha_R = GrowthR/spec%NPP
      end if

      ! Carbon cost for growth
      ! Rgrowth,Rgroot,Rgleaf,Rgstem, 0.5 is from IBIS and Amthor, 1984
      Rgleaf  = 0.5*GrowthL
      Rgstem  = 0.5*GrowthS
      Rgroot  = 0.5*GrowthR
      spec%Rgrowth = Rgleaf + Rgstem + Rgroot
      ! Leaf litter
      gamma_Wmax = 0.12/24. ! maxmum leaf fall rate per hour
      gamma_Tmax = 0.12/24.
      bW = 4.0
      bT = 2.0
      if (iforcing%Tair .gt. (Tcold + 10.)) then
         beta_T = 1.
      else
         if (iforcing%Tair .gt. Tcold) beta_T = (iforcing%Tair - Tcold)/10.
         if (iforcing%Tair .LE. Tcold) beta_T = 0.0
      end if

      ! if (spec%tauC(1) < 8760.) then
      !    gamma_W = (1.-W)**bW*gamma_Wmax
      !    gamma_T = (1.-beta_T)**bT*gamma_Tmax
      ! else
      !    gamma_W = 0.
      !    gamma_T = 0.
      ! end if

      gamma_N = 1.0/spec%tauC(1)*st%scalW      ! Modify by Jiang Jiang 2015/10/20
      if (spec%LAI < spec%LAIMIN) then
         gamma_W = 0.
         gamma_T = 0.
         gamma_N = 0.
      end if
      ! print *, "L_fall: ", spec%bmleaf, 0.48*gamma_N, spec%tauC(1), st%scalW
      spec%L_fall = spec%bmleaf*0.48*gamma_N     ! L_fall=bmleaf*0.48*AMIN1((gamma_T+gamma_N),0.99)
      return
   end subroutine plantgrowth

   subroutine yrday(doy, hour, lat, radsol, fbeam)  
      ! Jian: This subroutine is used to calculate the fbeam.
      integer, intent(in) :: doy, hour
      real, intent(in)    :: lat, radsol
      real, intent(out)   :: fbeam
      ! parameters in this subroutine 
      real pidiv, slatx, sindec, cosdec
      real a, b, sinbet0, solext, tmprat, tmpR, tmpK, fdiff

      pidiv   = pi/180.0                ! divide the pi value
      slatx   = lat*pidiv               
      sindec  = -sin(23.4*pidiv)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec  = sqrt(1.-sindec*sindec)
      a       = sin(slatx)*sindec
      b       = cos(slatx)*cosdec
      sinbet0 = a + b*cos(2*pi*(hour - 12.)/24.)
      solext  = 1370.0*(1.0 + 0.033*cos(2.0*pi*(doy - 10.)/365.0))*sinbet0
      tmprat  = radsol/solext                                              ! radsol: par in forcing
      tmpR    = 0.847 - 1.61*sinbet0 + 1.04*sinbet0*sinbet0
      tmpK    = (1.47 - tmpR)/1.66
      if (tmprat .le. 0.22) fdiff = 1.0
      if (tmprat .gt. 0.22 .and. tmprat .le. 0.35) then
         fdiff = 1.0 - 6.4*(tmprat - 0.22)*(tmprat - 0.22)
      end if
      if (tmprat .gt. 0.35 .and. tmprat .le. tmpK) then
         fdiff = 1.47 - 1.66*tmprat
      end if
      if (tmprat .ge. tmpK) then
         fdiff = tmpR
      end if
      fbeam = 1.0 - fdiff
      if (fbeam .lt. 0.0) fbeam = 0.0
      return
   end subroutine yrday

   subroutine xlayers(spec, iforcing, fbeam, coszen, radabv, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf, &
      Acan1, Acan2, Ecan1, Ecan2, Rsoilab1, Rsoilab2)!, Esoil) 
      ! the multi-layered canopy model developed by
      ! Ray Leuning with the new radiative transfer scheme
      ! implemented by Y.P. Wang (from Sellers 1986)
      ! 12/Sept/96 (YPW) correction for mean surface temperature of sunlit
      ! and shaded leaves
      ! Tleaf,i=sum{Tleaf,i(n)*fslt*Gaussw(n)}/sum{fslt*Gaussw(n)}
      ! ----------------------------------------------------------------------
      type(spec_data_type), intent(inout) :: spec
      type(forcing_data_type), intent(in) :: iforcing
      real, intent(in) :: fbeam, coszen, Radabv(2)
      real, intent(inout) :: Aleaf(2), Eleaf(2), Hleaf(2), gbleaf(2), gsleaf(2)
      real :: Acan1, Acan2, Ecan1, Ecan2, Rsoilab1, Rsoilab2!, Esoil    ! outputs
      
      ! local vars
      real    :: Rnst1 , Rnst2 , Qcan1 , Qcan2  !net rad, sunlit, vis rad
      real    :: Rcan1 , Rcan2 , Hcan1 , Hcan2  !NIR rad !Sens heat
      real    :: Gbwc1 , Gbwc2 , Gswc1 , Gswc2  !Boundary layer conductance; Canopy conductance
      real    :: Tleaf1 , Tleaf2                !Leaf Temp
      real    :: xphi1, xphi2, funG 
      real    :: pi180, cozen15, cozen45, cozen75, xK15, xK45, xK75
      real    :: transd, extkn
      real    :: scatt(2)
      real    :: rhoc(3, 2)       !Goudriaan
      real    :: rhoch, rhoc15, rhoc45, rhoc75
      real    :: scalex, fslt, fshd
      real    :: RnStL(5), QcanL(5), RcanL(5), AcanL(5), EcanL(5), HcanL(5)
      real    :: layer1(5), layer2(5)
      real    :: FLAIT1, flait
      real    :: Gaussx(5), Gaussw(5), Gaussw_cum(5) 
      real    :: wind, raero, WILTPT, FILDCP, Dair, RH, esat1, eairP
      real    :: TairK, Tlk1, Tlk2, Rsoilab3, Rsoilabs 
      real    :: rhocp, H2OLv, slope, psyc, Cmolar, fw1, Rsoil, rLAI 
      ! real    :: Hsoil
      integer :: nw, ng
      ! transfer to other subroutine
      real    :: extkd, extkb, flai, kpr(3,2), reff(3,2)
      real    :: windUx, Vcmxx, eJmxx
      real    :: QLair, QLleaf ! to Tsoil_simu
      real    :: QLsoil ! same in soil module
      ! calculate in other subroutines
      real    :: emair, Qabs(3,2), Rnstar(2), grdn
      real    :: Tleaf(2)
      
      data Gaussx/0.0469101, 0.2307534, 0.5, 0.7692465, 0.9530899/        ! 5-point
      data Gaussw/0.1184635, 0.2393144, 0.2844444, 0.2393144, 0.1184635/
      data Gaussw_cum/0.11846, 0.35777, 0.64222, 0.88153, 1.0/

      ! print *, "LAI", spec%LAI
      flait  = spec%LAI                 ! Jian: LAI relavant variable from vegetable to energy
      wind   = iforcing%WS
      Dair   = iforcing%VPD    ! air water vapour defficit? Unit Pa
      RH     = AMAX1(0.01,AMIN1(99.99,iforcing%RH))                ! relative humidity
      esat1  = 610.78*exp(17.27*iforcing%Tair/(iforcing%Tair + 237.3))      ! intermediate parameter
      eairP  = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure
      Dair   = esat1-eairP                                ! Jian: confused that SPRUCE has the VPD data, why calculate it again?
      if (wind .lt. 0.01) wind = 0.01   ! set windspeed to the minimum speed to avoid zero Gb
      ! soil water conditions
      WILTPT = st%wsmin/100.
      FILDCP = st%wsmax/100.
      ! reset the vairables
      Rnst1  = 0.0        !net rad, sunlit
      Rnst2  = 0.0        !net rad, shaded
      Qcan1  = 0.0        !vis rad
      Qcan2  = 0.0
      Rcan1  = 0.0        !NIR rad
      Rcan2  = 0.0
      Acan1  = 0.0        !CO2
      Acan2  = 0.0
      Ecan1  = 0.0        !Evap
      Ecan2  = 0.0
      Hcan1  = 0.0        !Sens heat
      Hcan2  = 0.0
      Gbwc1  = 0.0        !Boundary layer conductance
      Gbwc2  = 0.0
      Gswc1  = 0.0        !Canopy conductance
      Gswc2  = 0.0
      Tleaf1 = 0.0       !Leaf Temp
      Tleaf2 = 0.0  
      ! aerodynamic resistance
      raero = 50./wind
      ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*spec%xfang - 0.33*spec%xfang*spec%xfang
      xphi2 = 0.877*(1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*coszen      ! G-function: Projection of unit leaf area in direction of beam
      if (coszen .gt. 0) then          ! check if day or night
         extKb = funG/coszen           ! beam extinction coeff - black leaves
      else
         extKb = 100.
      end if
      ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      ! Effective extinction coefficient for diffuse radiation Goudriaan & van Laar Eq 6.6)
      pi180   = 3.1416/180.
      cozen15 = cos(pi180*15)
      cozen45 = cos(pi180*45)
      cozen75 = cos(pi180*75)
      xK15    = xphi1/cozen15 + xphi2
      xK45    = xphi1/cozen45 + xphi2
      xK75    = xphi1/cozen75 + xphi2
      !--------------------------------------
      transd  = 0.308*exp(-xK15*flait) + 0.514*exp(-xK45*flait) +  &
                &  0.178*exp(-xK75*flait)
      extkd   = (-1./flait)*alog(transd)
      extkn   = extkd                        !N distribution coeff

      ! canopy reflection coefficients (Array indices: first -> 1=VIS,  2=NIR; second -> 1=beam, 2=diffuse)
      do nw = 1, 2  ! nw:1=VIS, 2=NIR
         scatt(nw)   = tauL(nw) + rhoL(nw)                                               ! scattering coeff
         if ((1.-scatt(nw)) < 0.0) scatt(nw) = 0.9999                                    ! Weng 10/31/2008
         kpr(nw, 1)  = extKb*sqrt(1.-scatt(nw))                                          ! modified k beam scattered (6.20)
         kpr(nw, 2)  = extkd*sqrt(1.-scatt(nw))                                          ! modified k diffuse (6.20)
         rhoch       = (1.-sqrt(1.-scatt(nw)))/(1.+sqrt(1.-scatt(nw)))                   ! canopy reflection black horizontal leaves (6.19)
         rhoc15      = 2.*xK15*rhoch/(xK15 + extkd)                                      ! canopy reflection (6.21) diffuse
         rhoc45      = 2.*xK45*rhoch/(xK45 + extkd)
         rhoc75      = 2.*xK75*rhoch/(xK75 + extkd)
         rhoc(nw, 2) = 0.308*rhoc15 + 0.514*rhoc45 + 0.178*rhoc75
         rhoc(nw, 1) = 2.*extKb/(extKb + extkd)*rhoch                                    ! canopy reflection (6.21) beam
         reff(nw, 1) = rhoc(nw, 1) + (rhoS(nw) - rhoc(nw, 1))*exp(-2.*kpr(nw, 1)*flait)  ! effective canopy-soil reflection coeff - beam (6.27)
         reff(nw, 2) = rhoc(nw, 2) + (rhoS(nw) - rhoc(nw, 2))*exp(-2.*kpr(nw, 2)*flait)  ! effective canopy-soil reflection coeff - diffuse (6.27)
      end do
      ! isothermal net radiation & radiation conductance at canopy top - needed to calc emair
      call Radiso(iforcing, fbeam, extkd, flait, flai, Qabs, emair, Rnstar, grdn)           ! Jian: some parameters not initialization.
      TairK = iforcing%Tair + 273.2
      Tleaf1 = 0.
      do ng = 1, 5
         flai = gaussx(ng)*flait
         ! radiation absorption for visible and near infra-red
         call goudriaan(spec, fbeam, coszen, Radabv, kpr, reff, flai, scatt, Qabs)
         ! isothermal net radiation & radiation conductance at canopy top
         call Radiso(iforcing, fbeam, extkd, flait, flai, Qabs, emair, Rnstar, grdn)
         windUx = wind*exp(- st%extkU*flai)     ! windspeed at depth xi
         scalex = exp(-extkn*flai)              ! scale Vcmx0 & Jmax0
         Vcmxx  = spec%Vcmx0*scalex             ! Vcmx0 ---> Vcmax0
         ! print*, "vcmx0: ",Vcmxx, spec%Vcmx0, scalex, extkn, flai, gaussx(ng), flait,transd,spec%LAI, &
         ! spec%bmleaf, spec%SLA, spec%QC(1)
         eJmxx  = spec%eJmx0*scalex
         if (radabv(1) .ge. 10.0) then                          !check solar Radiation > 10 W/m2
            ! leaf stomata-photosynthesis-transpiration model - daytime
            call agsean_day(spec, iforcing, windUx, grdn, Qabs, Rnstar, Vcmxx, eJmxx, &
               Tleaf, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf)
         else
            call agsean_ngt(spec, iforcing, windUx, grdn, Qabs, Rnstar, Vcmxx, &
               Tleaf, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf)
         end if
         fslt      = exp(-extKb*flai)                        !fraction of sunlit leaves
         fshd      = 1.0 - fslt                                !fraction of shaded leaves
         Rnst1     = Rnst1 + fslt*Rnstar(1)*Gaussw(ng)*flait  !Isothermal net rad`
         Rnst2     = Rnst2 + fshd*Rnstar(2)*Gaussw(ng)*flait
         RnstL(ng) = Rnst1 + Rnst2

         Qcan1     = Qcan1 + fslt*Qabs(1, 1)*Gaussw(ng)*flait  !visible
         Qcan2     = Qcan2 + fshd*Qabs(1, 2)*Gaussw(ng)*flait
         QcanL(ng) = Qcan1 + Qcan2

         Rcan1     = Rcan1 + fslt*Qabs(2, 1)*Gaussw(ng)*flait  !NIR
         Rcan2     = Rcan2 + fshd*Qabs(2, 2)*Gaussw(ng)*flait
         RcanL(ng) = Rcan1 + Rcan2

         if (Aleaf(1) .lt. 0.0) Aleaf(1) = 0.0      !Weng 2/16/2006
         if (Aleaf(2) .lt. 0.0) Aleaf(2) = 0.0      !Weng 2/16/2006

         Acan1      = Acan1 + fslt*Aleaf(1)*Gaussw(ng)*flait*spec%stom_n    !amphi/hypostomatous
         Acan2      = Acan2 + fshd*Aleaf(2)*Gaussw(ng)*flait*spec%stom_n

         AcanL(ng)  = Acan1 + Acan2

         layer1(ng) = Aleaf(1)
         layer2(ng) = Aleaf(2)

         Ecan1      = Ecan1 + fslt*Eleaf(1)*Gaussw(ng)*flait
         Ecan2      = Ecan2 + fshd*Eleaf(2)*Gaussw(ng)*flait
         EcanL(ng)  = Ecan1 + Ecan2
         ! print *, "test_transp1: ", Ecan1, fslt, Eleaf(1), Gaussw(ng), flait

         Hcan1     = Hcan1 + fslt*Hleaf(1)*Gaussw(ng)*flait
         Hcan2     = Hcan2 + fshd*Hleaf(2)*Gaussw(ng)*flait
         HcanL(ng) = Hcan1 + Hcan2

         Gbwc1     = Gbwc1 + fslt*gbleaf(1)*Gaussw(ng)*flait*spec%stom_n
         Gbwc2     = Gbwc2 + fshd*gbleaf(2)*Gaussw(ng)*flait*spec%stom_n

         Gswc1     = Gswc1 + fslt*gsleaf(1)*Gaussw(ng)*flait*spec%stom_n
         Gswc2     = Gswc2 + fshd*gsleaf(2)*Gaussw(ng)*flait*spec%stom_n

         Tleaf1    = Tleaf1 + fslt*Tleaf(1)*Gaussw(ng)*flait
         Tleaf2    = Tleaf2 + fshd*Tleaf(2)*Gaussw(ng)*flait
      end do  ! 5 layers

      FLAIT1 = (1.0 - exp(-extKb*flait))/extkb
      Tleaf1 = Tleaf1/FLAIT1
      Tleaf2 = Tleaf2/(flait - FLAIT1)
      ! Soil surface energy and water fluxes
      ! Radiation absorbed by soil
      Rsoilab1 = fbeam*(1.-reff(1, 1))*exp(-kpr(1, 1)*flait)        &
          &         + (1.-fbeam)*(1.-reff(1, 2))*exp(-kpr(1, 2)*flait)          !visible
      Rsoilab2 = fbeam*(1.-reff(2, 1))*exp(-kpr(2, 1)*flait)        &
          &         + (1.-fbeam)*(1.-reff(2, 2))*exp(-kpr(2, 2)*flait)          !NIR
      Rsoilab1 = Rsoilab1*Radabv(1)
      Rsoilab2 = Rsoilab2*Radabv(2)
      Tlk1     = Tleaf1 + 273.2
      Tlk2     = Tleaf2 + 273.2
      QLair    = emair*sigma*(TairK**4)*exp(-extkd*flait)
      QLleaf   = emleaf*sigma*(Tlk1**4)*exp(-extkb*flait)           &
                 &      + emleaf*sigma*(Tlk2**4)*(1.0 - exp(-extkb*flait))
      if(ISNAN(QLleaf)) then
         ! write(*,*) extKb, flai, exp(-extKb*flai)
         ! write(*,*) fslt, Tleaf(1), Gaussw, flait
         write(*,*)"QLleaf is NAN1111: ",QLleaf, emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd, Tleaf1, FLAIT1
         stop
      endif
      QLleaf   = QLleaf*(1.0 - exp(-extkd*flait))
      QLsoil   = emsoil*sigma*(TairK**4)
      Rsoilab3 = (QLair + QLleaf)*(1.0 - rhoS(3)) - QLsoil
      if(ISNAN(QLleaf)) then
         write(*,*)"QLleaf is NAN: ", emleaf, sigma, Tlk1, extKb, flait, Tlk2, extkd
         stop
      endif
      ! Net radiation absorbed by soil
      ! the old version of net long-wave radiation absorbed by soils
      ! (with isothermal assumption)
      ! Rsoil3=(sigma*TairK**4)*(emair-emleaf)*exp(-extkd*flait)         !Longwave
      ! Rsoilab3=(1-rhoS(3))*Rsoil3

      ! Total radiation absorbed by soil
      Rsoilabs = Rsoilab1 + Rsoilab2 + Rsoilab3       ! Jian: Rsoilabs is re-calcualted in soil module

      ! thermodynamic parameters for air
      TairK  = iforcing%Tair + 273.2
      rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv  = H2oLv0 - 2.365e3*iforcing%Tair
      slope  = (esat(iforcing%Tair + 0.1) - esat(iforcing%Tair))/0.1
      psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      fw1    = AMIN1(AMAX1((FILDCP - st%wcl(1))/(FILDCP - WILTPT), 0.05), 1.0)
      Rsoil  = 30.*exp(0.2/fw1)
      rLAI   = exp(flait)
      ! latent heat flux into air from soil
      ! Eleaf(ileaf)=1.0*
      ! &     (slope*Y*Rnstar(ileaf)+rhocp*Dair/(rbH_L+raero))/    !2* Weng 0215
      ! &     (slope*Y+psyc*(rswv+rbw+raero)/(rbH_L+raero))
      ! Jian: seem have no effect in this module, Esoil and Hsoil will be calculated in soil module
      spec%Esoil  = (slope*(Rsoilabs - st%G) + rhocp*Dair/(raero + rLAI))/       &
               &      (slope + psyc*(rsoil/(raero + rLAI) + 1.))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
      ! sensible heat flux into air from soil
      spec%Hsoil = Rsoilabs - spec%Esoil - st%G  ! Jian: it is also calculated in soil module?
      return
   end subroutine xlayers

   subroutine Radiso(iforcing, fbeam, extkd, flait, flai, Qabs, emair, Rnstar, grdn)   
      ! Rnstar(type): type=1 for sunlit; =2 for shaded leaves (W/m2)
      ! 23 Dec 1994
      ! calculates isothermal net radiation for sunlit and shaded leaves under clear skies
      ! -----------------------------------------------------------------------------------
      type(forcing_data_type), intent(in) :: iforcing
      real, intent(in) :: fbeam, extkd, flait, flai, Qabs(3,2)
      ! local vars
      real emsky, ep8z, tau8, emcloud, Bn0, Bnxi
      real TairK, rhocp, RH, esat1, eairP
      ! update here and use it here and other subroutine
      real emair, Rnstar(2), grdn

      TairK = iforcing%Tair + 273.2
      RH    = AMAX1(0.01,AMIN1(99.99,iforcing%RH))                ! relative humidity
      esat1 = 610.78*exp(17.27*iforcing%Tair/(iforcing%Tair + 237.3))      ! intermediate parameter
      eairP = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure

      ! thermodynamic properties of air
      rhocp   = cpair*Patm*airMa/(Rconst*TairK)   ! volumetric heat capacity (J/m3/K)
      ! apparent atmospheric emissivity for clear skies (Brutsaert, 1975)
      emsky   = 0.642*(eairP/Tairk)**(1./7)       ! note eair in Pa
      ! apparent emissivity from clouds (Kimball et al 1982)
      ep8z    = 0.24 + 2.98e-12*eairP*eairP*exp(3000/TairK)
      tau8    = amin1(1.0, 1.0 - ep8z*(1.4 - 0.4*ep8z))            !ensure tau8<1
      emcloud = 0.36*tau8*(1.-fbeam)*(1 - 10./TairK)**4      !10 from Tcloud = Tair-10
      ! apparent emissivity from sky plus clouds
      !      emair=emsky+emcloud
      ! 20/06/96
      emair   = emsky
      if (emair .gt. 1.0) emair = 1.0
      ! net isothermal outgoing longwave radiation per unit leaf area at canopy
      ! top & thin layer at flai (Note Rn* = Sn + Bn is used rather than Rn* = Sn - Bn in Leuning et al 1985)
      Bn0  = sigma*(TairK**4.)
      Bnxi = Bn0*extkd*(exp(-extkd*flai)*(emair - emleaf)       &
          &    + exp(-extkd*(flait - flai))*(emsoil - emleaf))
      ! isothermal net radiation per unit leaf area for thin layer of sunlit and
      ! shaded leaves
      Rnstar(1) = Qabs(1, 1) + Qabs(2, 1) + Bnxi
      Rnstar(2) = Qabs(1, 2) + Qabs(2, 2) + Bnxi
      ! radiation conductance (m/s) @ flai
      grdn = 4.*sigma*(TairK**3.)*extkd*emleaf*               &       ! corrected by Jiang Jiang 2015/9/29
          &    (exp(-extkd*flai) + exp(-extkd*(flait - flai)))       &
          &    /rhocp
      return
   end subroutine Radiso

   subroutine goudriaan(spec, fbeam, coszen, Radabv, kpr, reff, flai, scatt, Qabs)
      ! for spheric leaf angle distribution only
      ! compute within canopy radiation (PAR and near infra-red bands)
      ! using two-stream approximation (Goudriaan & vanLaar 1994)
      ! tauL: leaf transmittance
      ! rhoL: leaf reflectance
      ! rhoS: soil reflectance
      ! sfang XiL function of Ross (1975) - allows for departure from spherical LAD
      !     (-1 vertical, +1 horizontal leaves, 0 spherical)
      ! FLAI: canopy leaf area index
      ! funG: Ross' G function
      ! scatB: upscatter parameter for direct beam
      ! scatD: upscatter parameter for diffuse
      ! albedo: single scattering albedo
      ! output:
      ! Qabs(nwave,type), nwave=1 for visible; =2 for NIR,
      !                   type =1 for sunlit;  =2 for shaded (W/m2)
      !------------------------------------------------------------------------------
      type(spec_data_type), intent(in) :: spec
      real, intent(in) :: fbeam, coszen, Radabv(2)    !cos zenith angle
      real, intent(in) :: kpr(3,2), reff(3,2), scatt(2), flai
      real :: Qabs(3,2)

      real :: extkb ! the same calculation process in xlayer
      real xphi1, xphi2, funG, Qd0, Qb0
      integer nw

      ! Ross-Goudriaan function for G(u) (see Sellers 1985, Eq 13)
      xphi1 = 0.5 - 0.633*spec%xfang - 0.33*spec%xfang*spec%xfang
      xphi2 = 0.877*(1.0 - 2.0*xphi1)
      funG = xphi1 + xphi2*coszen                        !G-function: Projection of unit leaf area in direction of beam
      if (coszen .gt. 0) then                            !check if day or night
         extKb = funG/coszen                             !beam extinction coeff - black leaves
      else
         extKb = 100.
      end if
      ! Goudriaan theory as used in Leuning et al 1995 (Eq Nos from Goudriaan & van Laar, 1994)
      do nw = 1, 2
         Qd0 = (1.-fbeam)*radabv(nw)                    !diffuse incident radiation
         Qb0 = fbeam*radabv(nw)                         !beam incident radiation
         Qabs(nw, 2) = Qd0*(kpr(nw, 2)*(1.-reff(nw, 2))*exp(-kpr(nw, 2)*FLAI)) +  & ! absorbed radiation - shaded leaves, diffuse
                &      Qb0*(kpr(nw, 1)*(1.-reff(nw, 1))*exp(-kpr(nw, 1)*FLAI)  -  & ! beam scattered
                &      extKb*(1.-scatt(nw))*exp(-extKb*FLAI))
         Qabs(nw, 1) = Qabs(nw, 2) + extKb*Qb0*(1.-scatt(nw))                     !absorbed radiation - sunlit leaves
      end do
      return
   end subroutine goudriaan

   subroutine agsean_day(spec, iforcing, windUx, grdn, Qabs, Rnstar, Vcmxx, eJmxx, &
      Tleaf, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf)
      implicit none
      type(spec_data_type), intent(inout) :: spec
      type(forcing_data_type), intent(in) :: iforcing
      real, intent(in) :: windUx, grdn, Qabs(3,2), Rnstar(2), Vcmxx, eJmxx
      real, intent(inout) :: Eleaf(2), Tleaf(2),  Hleaf(2), gbleaf(2), gsleaf(2), Aleaf(2)
      ! real :: Tleaf(2),  Hleaf(2), gbleaf(2), gsleaf(2)  ! output?
      ! local variables
      integer :: kr1, ileaf
      real    :: Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn, Y
      real    :: gsw, gswv, rswv
      real    :: TairK, Tlk1, esat1, eairP, RH, Dair, wind 
      real    :: rhocp, H2OLv, Cmolar, raero
      real    :: gbHu, Tlk, Dleaf, slope, psyc
      real    :: gbc, gsc0, co2cs, Qapar
      real    :: co2ci(2)
      ! variables from other subroutine
      real    :: Aleafx, Gscx

      ! thermodynamic parameters for air
      TairK  = iforcing%Tair + 273.2
      esat1  = 610.78*exp(17.27*iforcing%Tair/(iforcing%Tair + 237.3))      ! intermediate parameter
      RH     = AMAX1(0.01,AMIN1(99.99,iforcing%RH))                ! relative humidity
      eairP  = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure
      Dair   = esat1-eairP
      wind   = iforcing%WS
      if (wind .lt. 0.01) wind = 0.01   ! set windspeed to the minimum speed to avoid zero Gb

      rhocp  = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv  = H2oLv0 - 2.365e3*iforcing%Tair
      slope  = (esat(iforcing%Tair + 0.1) - esat(iforcing%Tair))/0.1
      psyc   = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      ! boundary layer conductance for heat - single sided, forced convection
      ! (Monteith 1973, P106 & notes dated 23/12/94)
      if (windUx/wleaf >= 0.0) then
         gbHu = 0.003*sqrt(windUx/wleaf)    !m/s
      else
         gbHu = 0.003 !*sqrt(-windUx/wleaf)
      end if         ! Weng 10/31/2008
      raero = 50./wind                 ! aerodynamic resistance s/m
      do ileaf = 1, 2                  ! loop over sunlit and shaded leaves
         ! first estimate of leaf temperature - assume air temp
         Tleaf(ileaf) = iforcing%Tair
         Tlk = Tleaf(ileaf) + 273.2    ! Tleaf to deg K
         ! first estimate of deficit at leaf surface - assume Da
         Dleaf = Dair                  ! Pa
         ! first estimate for co2cs
         co2cs = co2ca                 ! mol/mol
         Qapar = (4.6e-6)*Qabs(1, ileaf)
         ! --------------------------------------------------------------
         kr1 = 0                     !iteration counter for LE
         ! return point for evaporation iteration
         do ! iteration for leaf temperature
            ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras   = 1.595e8*ABS(Tleaf(ileaf) - iforcing%Tair)*(wleaf**3.)     !Grashof
            gbHf   = 0.5*Dheat*(Gras**0.25)/wleaf
            gbH    = gbHu + gbHf                         !m/s
            rbH    = 1./gbH                            !b/l resistance to heat transfer
            rbw    = 0.93*rbH                          !b/l resistance to water vapour
            ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L  = rbH*spec%stom_n/2.                   !final b/l resistance for heat
            rrdn   = 1./grdn
            Y      = 1./(1.+(rbH_L + raero)/rrdn)
            ! boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc    = Cmolar*gbH/1.32            !mol/m2/s
            gsc0   = gsw0/1.57                 !convert conductance for H2O to that for CO2
            ! varQc  = 0.0
            ! -------------------------------------------
            call photosyn(spec, co2cs, Qapar, Vcmxx, Tlk, gsc0, eJmxx, Dleaf, Gbc, &
                           Aleafx, Gscx)  !outputs: Aleafx
            ! choose smaller of Ac, Aq
               
            Aleaf(ileaf) = Aleafx      !0.7 Weng 3/22/2006          !mol CO2/m2/s
            ! calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca - Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs - Aleaf(ileaf)/gscx
            ! scale variables
            gsw  = gscx*1.56       !gsw in mol/m2/s, oreginal:gsw=gscx*1.56,Weng20090226
            gswv = gsw/Cmolar                           !gsw in m/s
            rswv = 1./gswv
            ! calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf) = 1.0*(slope*Y*Rnstar(ileaf) + rhocp*Dair/(rbH_L + raero))/   &   !2* Weng 0215
                           & (slope*Y + psyc*(rswv + rbw + raero)/(rbH_L + raero))   !*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))!*(0.08+(1.-0.08)* Amin1(1.0, 0.3*omega))
            ! calculate sensible heat flux
            Hleaf(ileaf) = Y*(Rnstar(ileaf) - Eleaf(ileaf))
            ! calculate new leaf temperature (K)
            Tlk1 = 273.2 + iforcing%Tair + Hleaf(ileaf)*(rbH/2.+raero)/rhocp
            ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf = psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf) = gbc*1.32*1.075
            gsleaf(ileaf) = gsw
            ! compare current and previous leaf temperatures
            if (abs(Tlk1 - Tlk) .le. 0.1) exit ! original is 0.05 C Weng 10/31/2008
            ! update leaf temperature  ! leaf temperature calculation has many problems! Weng 10/31/2008
            Tlk = Tlk1
            Tleaf(ileaf) = Tlk1 - 273.2
            kr1 = kr1 + 1
            if (kr1 > 500) then
               Tlk = TairK
               exit
            end if
            if (Tlk < 200.) then
               Tlk = TairK
               exit
            end if                     ! Weng 10/31/2008
            ! goto 100                          !solution not found yet
         end do
         ! 10  continue
      end do
      return
   end subroutine agsean_day

   ! -------------------------------------------------------------------------
   subroutine agsean_ngt(spec, iforcing, windUx, grdn, Qabs, Rnstar, Vcmxx,  &
      Tleaf, Aleaf, Eleaf, Hleaf, gbleaf, gsleaf)
      implicit none
      type(spec_data_type), intent(inout) :: spec
      type(forcing_data_type), intent(in) :: iforcing
      real, intent(in) :: windUx, grdn, Qabs(3,2), Rnstar(2), Vcmxx
      real :: Tleaf(2), Eleaf(2), Hleaf(2), gbleaf(2), gsleaf(2)  ! output?
      ! local variables
      integer :: kr1, ileaf
      real    :: Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn, Y
      real    :: gsw, gswv, rswv
      real    :: TairK, Tlk1, esat1, eairP, RH, Dair, wind
      real    :: rhocp, H2OLv, Cmolar, raero
      real    :: gbHu, Tlk, Dleaf, slope, psyc
      real    :: gbc, gsc0, co2cs, Qapar, gsc
      real    :: Aleaf(2), co2ci(2)
      ! variables from other subroutine
      real    :: Aleafx
      
      ! integer kr1, ileaf
      ! real Gras, gbHf, gbH, rbH, rbw, rbH_L, rrdn, Y
      ! real gsw, gswv, rswv
      ! real gsc
      ! thermodynamic parameters for air

      TairK = iforcing%Tair + 273.2
      esat1  = 610.78*exp(17.27*iforcing%Tair/(iforcing%Tair + 237.3))      ! intermediate parameter
      RH     = AMAX1(0.01,AMIN1(99.99,iforcing%RH))                ! relative humidity
      eairP  = esat1*RH/100.                              ! Added for SPRUCE, due to lack of VPD data. Jian: ? SPRUCE has the data? !air water vapour pressure
      Dair   = esat1-eairP
      wind   = iforcing%WS
      if (wind .lt. 0.01) wind = 0.01   ! set windspeed to the minimum speed to avoid zero Gb

      rhocp = cpair*Patm*AirMa/(Rconst*TairK)
      H2OLv = H2oLv0 - 2.365e3*iforcing%Tair
      slope = (esat(iforcing%Tair + 0.1) - esat(iforcing%Tair))/0.1
      psyc = Patm*cpair*AirMa/(H2OLv*H2OMw)
      Cmolar = Patm/(Rconst*TairK)
      ! weighJ = 1.0
      ! boundary layer conductance for heat - single sided, forced convection
      ! (Monteith 1973, P106 & notes dated 23/12/94)
      gbHu = 0.003*sqrt(windUx/wleaf)    !m/s
      ! raero=0.0                        !aerodynamic resistance s/m
      raero = 50./wind                 ! aerodynamic resistance s/m
      do ileaf = 1, 2                  ! loop over sunlit and shaded leaves
         ! first estimate of leaf temperature - assume air temp
         Tleaf(ileaf) = iforcing%Tair
         Tlk = Tleaf(ileaf) + 273.2    !Tleaf to deg K
         ! first estimate of deficit at leaf surface - assume Da
         Dleaf = Dair                !Pa
         ! first estimate for co2cs
         co2cs = co2ca               !mol/mol
         Qapar = (4.6e-6)*Qabs(1, ileaf)
         ! ********************************************************************
         kr1 = 0                     !iteration counter for LE
         do
            !100        continue !    return point for evaporation iteration
            ! single-sided boundary layer conductance - free convection (see notes 23/12/94)
            Gras = 1.595e8*abs(Tleaf(ileaf) - iforcing%Tair)*(wleaf**3)     !Grashof
            gbHf = 0.5*Dheat*(Gras**0.25)/wleaf
            gbH = gbHu + gbHf                         !m/s
            rbH = 1./gbH                            !b/l resistance to heat transfer
            rbw = 0.93*rbH                          !b/l resistance to water vapour
            ! Y factor for leaf: stom_n = 1.0 for hypostomatous leaf;  stom_n = 2.0 for amphistomatous leaf
            rbH_L = rbH*spec%stom_n/2.                   !final b/l resistance for heat
            rrdn = 1./grdn
            Y = 1./(1.+(rbH_L + raero)/rrdn)
            ! boundary layer conductance for CO2 - single side only (mol/m2/s)
            gbc = Cmolar*gbH/1.32            !mol/m2/s
            gsc0 = gsw0/1.57                        !convert conductance for H2O to that for CO2
            ! varQc = 0.0
            ! weighR = 1.0
            ! respiration
            Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2))
            gsc = gsc0
            ! choose smaller of Ac, Aq
            if (ISNAN(Aleafx)) then
               print*, "Aleafx is nan", Aleafx, Vcmxx, Tlk
               stop
            endif
            Aleaf(ileaf) = Aleafx                     !mol CO2/m2/s
            ! calculate new values for gsc, cs (Lohammer model)
            co2cs = co2ca - Aleaf(ileaf)/gbc
            co2Ci(ileaf) = co2cs - Aleaf(ileaf)/gsc
            ! scale variables
            gsw = gsc*1.56                              !gsw in mol/m2/s
            gswv = gsw/Cmolar                           !gsw in m/s
            rswv = 1./gswv
            ! calculate evap'n using combination equation with current estimate of gsw
            Eleaf(ileaf) = (slope*Y*Rnstar(ileaf) + rhocp*Dair/(rbH_L + raero))/   &
                &      (slope*Y + psyc*(rswv + rbw + raero)/(rbH_L + raero))!*omega!*AMIN1(0.5, AMAX1(0.333, 0.333 + omega))
            ! calculate sensible heat flux
            Hleaf(ileaf) = Y*(Rnstar(ileaf) - Eleaf(ileaf))
            ! calculate new leaf temperature (K)
            Tlk1 = 273.2 + iforcing%Tair + Hleaf(ileaf)*(rbH/2.+raero)/rhocp
            ! write (*,*) "Tlk1",Tlk1,"Hleaf(ileaf)",Hleaf(ileaf),"rbH",rbH,"raero",raero,"rhocp",rhocp
            ! calculate Dleaf use LE=(rhocp/psyc)*gsw*Ds
            Dleaf = psyc*Eleaf(ileaf)/(rhocp*gswv)
            gbleaf(ileaf) = gbc*1.32*1.075
            gsleaf(ileaf) = gsw

            ! compare current and previous leaf temperatures
            if (abs(Tlk1 - Tlk) .le. 0.1) exit
            if (kr1 .gt. 500) exit
            ! update leaf temperature
            Tlk = Tlk1
            Tleaf(ileaf) = Tlk1 - 273.2
            kr1 = kr1 + 1
         end do                          !solution not found yet
! 10       continue
      end do
      return
   end subroutine agsean_ngt

   subroutine photosyn(spec, co2cs, Qapar, Vcmxx, Tlk, gsc0, eJmxx, Dleaf, Gbc, &
      Aleafx, Gscx)
      implicit none
      type(spec_data_type), intent(inout) :: spec
      real, intent(in) :: co2cs, Qapar, Vcmxx, Tlk, gsc0, eJmxx
      real, intent(in) :: Dleaf, Gbc

      ! outputs
      real :: Aleafx, Gscx
      ! calculate Vcmax, Jmax at leaf temp (Eq 9, Harley et al 1992)
      ! turned on by Weng, 2012-03-13

      ! -------------------------------------------------------

      ! local variables
      real :: weighJ, weighR, a1
      real :: CO2Csx
      real :: TminV, TmaxV, ToptV, TminJ, TmaxJ, ToptJ, Tlf
      real :: VcmxT, eJmxT, eJ
      real :: conKcT, conKoT
      real :: Rd, Tdiff, gammas, gamma
      real :: X, Gma, Bta
      real :: Acx, Aqx ! from ciandA

      VcmxT = Vjmax(Tlk,Trefk,Vcmxx,Eavm,Edvm,Rconst,spec%Entrpy)
      eJmxT = Vjmax(Tlk,Trefk,eJmxx,Eajm,Edjm,Rconst,spec%Entrpy)
      CO2Csx = AMAX1(CO2Cs, 0.6*CO2Ca)
      ! check if it is dark - if so calculate respiration and g0 to assign conductance
      if (Qapar .le. 0.) then                            !night, umol quanta/m2/s
         Aleafx = -0.0089*Vcmxx*exp(0.069*(Tlk - 293.2))   ! original: 0.0089 Weng 3/22/2006
         Gscx = gsc0
      end if
      ! calculate  Vcmax, Jmax at leaf temp using Reed et al (1976) function J appl Ecol 13:925
      TminV = spec%gddonset/10.  ! original -5.        !-Jiang Jiang 2015/10/13
      TmaxV = 50.
      ToptV = 35.

      TminJ = TminV
      TmaxJ = TmaxV
      ToptJ = ToptV

      Tlf = Tlk - 273.2
      ! VcmxT = VJtemp(Tlf, TminV, TmaxV, ToptV, Vcmxx)
      ! eJmxT = VJtemp(Tlf, TminJ, TmaxJ, ToptJ, eJmxx)
      ! calculate J, the asymptote for RuBP regeneration rate at given Q
      weighJ = 1.0
      weighR = 1.0
      eJ = weighJ*fJQres(eJmxT, spec%alpha, Qapar, theta)
      ! calculate Kc, Ko, Rd gamma*  & gamma at leaf temp
      conKcT = EnzK(Tlk, Trefk, conKc0, Rconst, Ekc)
      conKoT = EnzK(Tlk, Trefk, conKo0, Rconst, Eko)
      ! following de Pury 1994, eq 7, make light respiration a fixed proportion of Vcmax
      Rd     = 0.0089*VcmxT*weighR                              !de Pury 1994, Eq7
      Tdiff  = Tlk - Trefk
      gammas = gam0*(1.+gam1*Tdiff + gam2*Tdiff*Tdiff)       !gamma*
      ! gamma = (gammas+conKcT*(1.+O2ci/conKoT)*Rd/VcmxT)/(1.-Rd/VcmxT)
      gamma = 0.0
      ! ***********************************************************************
      ! Analytical solution for ci. This is the ci which satisfies supply and demand
      ! functions simultaneously
      ! calculate X using Lohammer model, and scale for soil moisture
      a1 = 1./(1.-0.7)
      X = a1*st%fwsoil/((co2csx - gamma)*(1.0 + Dleaf/spec%Ds0))
      ! calculate solution for ci when Rubisco activity limits A
      Gma = VcmxT
      Bta = conKcT*(1.0 + o2ci/conKoT)

      call ciandA(Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Acx)
      ! calculate +ve root for ci when RuBP regeneration limits A
      Gma = eJ/4.
      Bta = 2.*gammas
      ! calculate coefficients for quadratic equation for ci
      call ciandA(Gma, Bta, gsc0, X, Rd, co2Csx, gammas, Aqx)
      ! choose smaller of Ac, Aq
      ! sps = AMAX1(0.001, sps)                  !Weng, 3/30/2006
      Aleafx = (amin1(Acx, Aqx) - Rd) !*sps     ! Weng 4/4/2006

      ! if(Aleafx.lt.0.0) Aleafx=0.0    ! by Weng 3/21/2006
      ! calculate new values for gsc, cs (Lohammer model)
      CO2csx = co2ca - Aleafx/Gbc
      Gscx = gsc0 + X*Aleafx  ! revised by Weng
      return
   end subroutine photosyn

   subroutine ciandA(Gma,Bta,g0,X,Rd,co2Csx,gammas,Aquad)      ! Gma,Bta,g0,X,Rd,co2Cs,gammas,ciquad,Aquad
      real Gma,Bta,g0,X,Rd,co2Csx,gammas,ciquad,Aquad
      real b2, b1, b0, bx
      ! calculate coefficients for quadratic equation for ci
      b2 =  g0 + X*(Gma - Rd)
      b1 =  (1.-co2Csx*X)*(Gma - Rd) + g0*(Bta - co2Csx) - X*(Gma*gammas + Bta*Rd)
      b0 = -(1.-co2Csx*X)*(Gma*gammas + Bta*Rd) - g0*Bta*co2Csx

      bx = b1*b1 - 4.*b2*b0
      if (bx .gt. 0.0) then
         ! calculate larger root of quadratic
         ciquad = (-b1 + sqrt(bx))/(2.*b2)
      end if

      IF (ciquad .lt. 0 .or. bx .lt. 0.) THEN
         Aquad = 0.0
         ciquad = 0.7*co2Csx
      ELSE
         Aquad = Gma*(ciquad - gammas)/(ciquad + Bta)
      END IF
      return
   end

   ! functions 
   real function sinbet(doy, hour, lat)
      integer, intent(in) :: doy, hour
      real, intent(in) :: lat
      real rad, sinlat, coslat, sindec, cosdec, A, B
      ! sin(bet), bet = elevation angle of sun
      ! calculations according to Goudriaan & van Laar 1994 P30
      rad = pi/180.
      ! sine and cosine of latitude
      sinlat = sin(rad*lat)
      coslat = cos(rad*lat)
      ! sine of maximum declination
      sindec = -sin(23.45*rad)*cos(2.0*pi*(doy + 10.0)/365.0)
      cosdec = sqrt(1.-sindec*sindec)
      ! terms A & B in Eq 3.3
      A = sinlat*sindec
      B = coslat*cosdec
      sinbet = A + B*cos(pi*(hour - 12.)/12.)
      return
   end function sinbet

   ! real function esat(T)
   !    real T
   !    ! returns saturation vapour pressure in Pa
   !    esat = 610.78*exp(17.27*T/(T + 237.3))
   !    return
   ! end

   ! ****************************************************************************
   ! Reed et al (1976, J appl Ecol 13:925) equation for temperature response
   ! used for Vcmax and Jmax
   real function VJtemp(Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0)
      real Tlf,TminVJ,TmaxVJ,ToptVJ,VJmax0
      real pwr
      if (Tlf .lt. TminVJ) Tlf = TminVJ   !constrain leaf temperatures between min and max
      if (Tlf .gt. TmaxVJ) Tlf = TmaxVJ
      pwr    = (TmaxVJ - ToptVJ)/(ToptVj - TminVj)
      VJtemp = VJmax0*((Tlf - TminVJ)/(ToptVJ - TminVJ))*     &
               &       ((TmaxVJ - Tlf)/(TmaxVJ - ToptVJ))**pwr
      return
   end
   !****************************************************************************
   real function Vjmax(Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop)
      real :: Tk,Trefk,Vjmax0,Eactiv,Edeact,Rconst,Entrop, aden, anum
      anum = Vjmax0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      aden = 1. + EXP((Entrop*Tk-Edeact)/(Rconst*Tk))
      Vjmax = anum/aden
      return
   end

   ! ****************************************************************************
   real function fJQres(eJmx,alpha,Q,theta)
      real eJmx,alpha,Q,theta
      real AX, BX, CX
      AX = theta                                 !a term in J fn
      BX = alpha*Q + eJmx                          !b term in J fn
      CX = alpha*Q*eJmx                          !c term in J fn
      if ((BX*BX - 4.*AX*CX) >= 0.0) then
         fJQres = (BX - SQRT(BX*BX - 4.*AX*CX))/(2*AX)
      else
         fJQres = (BX)/(2*AX)                   !Weng 10/31/2008
      end if
      return
   end
   ! *************************************************************************
   real function EnzK(Tk,Trefk,EnzK0,Rconst,Eactiv)
      real Tk,Trefk,EnzK0,Rconst,Eactiv
      real temp1 
      temp1 = (Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk)
      ! if (temp1<50.)then
      EnzK = EnzK0*EXP((Eactiv/(Rconst*Trefk))*(1.-Trefk/Tk))
      ! else
      ! EnzK = EnzK0*EXP(50.)                                          ! Weng 10/31/2008
      ! endif
      return
   end

end module vegetation