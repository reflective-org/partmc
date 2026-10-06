! Copyright (C) 2026 Ali Akherati
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.
!
! The Maattanen et al. (2018) neutral channel below is ported from the
! authors' reference code newbinapara_v1.0.f90 (Copyright (C) 2018
! Maattanen et al.; Zenodo doi:10.5281/zenodo.1217782; please cite the
! article). Its coefficients are kept to the 17 digits the authors wrote
! (double-precision literals; the reference's default-REAL literals would
! round them to single precision) and its print warnings are removed.

!> \file
!> The pmc_nucleate_binary module.

!> Neutral binary H2SO4-H2O nucleation: formation rate and new particle.
!!
!! Each scheme returns the formation rate together with the particle it
!! forms (sulfuric acid molecules, total molecules and wet radius), as
!! defined by the shared nucleation schemes of the aerosol box-model
!! intercomparison (\c harness/nucleation_schemes.py). Units follow the
!! Maattanen et al. (2018) reference code: concentrations in cm^{-3},
!! rates in cm^{-3} s^{-1}, radii in m. Neutral channels only.
!!
!! References:
!!   - A. Maattanen, J. Merikanto, H. Henschel, J. Duplissy,
!!     R. Makkonen, I. K. Ortega and H. Vehkamaki (2018), New
!!     parameterizations for neutral and ion-induced sulfuric acid-water
!!     particle formation in nucleation and kinetic regimes,
!!     <i>J. Geophys. Res. Atmos.</i>, 123, 1269-1296,
!!     doi:10.1002/2017JD027429.
!!   - E. M. Dunne et al. (2016), Global atmospheric particle formation
!!     from CERN CLOUD measurements, <i>Science</i>, 354, 1119-1124,
!!     doi:10.1126/science.aaf2649.
!!   - A. Tabazadeh, O. B. Toon, S. L. Clegg and P. Hamill (1997), A new
!!     parameterization of H2SO4/H2O aerosol composition: Atmospheric
!!     implications, <i>Geophys. Res. Lett.</i>, 24, 1931-1934.
!!   - H. Vehkamaki et al. (2002), An improved parameterization for
!!     sulfuric acid-water nucleation rates for tropospheric and
!!     stratospheric conditions, <i>J. Geophys. Res.</i>, 107, 4622.
module pmc_nucleate_binary

  use pmc_constants

  !> H2SO4 molar mass of the shared schemes (g mol^{-1}).
  real(kind=dp), parameter :: NUCLEATE_BINARY_MW_H2SO4 = 98.0784d0
  !> H2O molar mass of the shared schemes (g mol^{-1}).
  real(kind=dp), parameter :: NUCLEATE_BINARY_MW_H2O = 18.01528d0
  !> Avogadro constant of the shared schemes (mol^{-1}).
  real(kind=dp), parameter :: NUCLEATE_BINARY_AVOGADRO = 6.02214076d23
  !> Pi.
  real(kind=dp), parameter :: NUCLEATE_BINARY_PI = 3.14159265358979323846d0

  !> Dunne et al. (2016) binary-neutral exponent \f$p\f$.
  real(kind=dp), parameter :: NUCLEATE_DUNNE_P_BN = 3.95d0
  !> Dunne et al. (2016) binary-neutral coefficient \f$u\f$.
  real(kind=dp), parameter :: NUCLEATE_DUNNE_U_BN = 9.70d0
  !> Dunne et al. (2016) binary-neutral coefficient \f$v\f$.
  real(kind=dp), parameter :: NUCLEATE_DUNNE_V_BN = 12.6d0
  !> Dunne et al. (2016) binary-neutral coefficient \f$w\f$.
  real(kind=dp), parameter :: NUCLEATE_DUNNE_W_BN = -0.00707d0
  !> Geometric diameter of the Dunne et al. (2016) particle (m): 1.7 nm
  !> mobility diameter minus 0.3 nm.
  real(kind=dp), parameter :: NUCLEATE_DUNNE_DIAM = 1.4d-9

  !> Monomer radius of the Maattanen et al. (2018) kinetic regime (m).
  real(kind=dp), parameter :: NUCLEATE_MAATTANEN_MONOMER_RADIUS = 0.3d-9

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> H2SO4 dimer collision rate (cm^{-3} s^{-1}): the kinetic-regime rate
  !> of the Maattanen et al. (2018) reference code (hard spheres of radius
  !> 0.3 nm and 98.07 g/mol), at the unbounded inputs as in the reference.
  real(kind=dp) function nucleate_binary_dimer_collision_rate(t, rhoa)

    !> Temperature (K).
    real(kind=dp), intent(in) :: t
    !> Sulfuric acid concentration (cm^{-3}).
    real(kind=dp), intent(in) :: rhoa

    nucleate_binary_dimer_collision_rate = 1.d6*(2.d0*0.3d-9)**2 &
         * sqrt(8.d0*3.141593d0*1.38d-23*(1.d0/(1.661d-27*98.07d0) &
         + 1.d0/(1.661d-27*98.07d0)))/2.d0*sqrt(t)*rhoa**2

  end function nucleate_binary_dimer_collision_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Neutral channel of \c newbinapara (Maattanen et al. 2018): rate of
  !> formation of critical clusters and the critical cluster itself.
  !!
  !! Inputs are bounded exactly as in the reference code (T 165-400 K,
  !! saturation ratio 1e-5-1, [H2SO4] 1e4-1e13 cm^{-3}). In the kinetic
  !! regime the rate is the dimer collision rate and the critical cluster
  !! is one acid molecule of radius 0.3 nm. Rates below 1e-7 cm^{-3}
  !! s^{-1} are set to zero.
  subroutine nucleate_binary_maattanen2018_cluster(t, satrat, rhoa, &
       jnuc_n, ntot_n, x_n, na_n, rc_n, kinetic_n)

    !> Temperature (K).
    real(kind=dp), intent(in) :: t
    !> Saturation ratio of water (between zero and 1).
    real(kind=dp), intent(in) :: satrat
    !> Sulfuric acid concentration (cm^{-3}).
    real(kind=dp), intent(in) :: rhoa
    !> Neutral nucleation rate (cm^{-3} s^{-1}).
    real(kind=dp), intent(out) :: jnuc_n
    !> Total number of molecules in the critical cluster.
    real(kind=dp), intent(out) :: ntot_n
    !> Mole fraction of H2SO4 in the critical cluster.
    real(kind=dp), intent(out) :: x_n
    !> Sulfuric acid molecules in the critical cluster.
    real(kind=dp), intent(out) :: na_n
    !> Radius of the critical cluster (m).
    real(kind=dp), intent(out) :: rc_n
    !> Whether the kinetic regime applies.
    logical, intent(out) :: kinetic_n

    real(kind=dp) :: satratln, rhoaln, tln, kinrhotresn

    kinetic_n = .false.
    jnuc_n = 0d0
    ntot_n = 0d0
    na_n = 0d0
    rc_n = 0d0
    x_n = 0d0
    satratln = satrat
    rhoaln = rhoa
    tln = t

    ! bounds of the parameterization (reference code, neutral case)
    if (t .le. 165d0) tln = 165d0
    if (t .ge. 400d0) tln = 400d0
    if (satrat .lt. 1d-5) satratln = 1d-5
    if (satrat .gt. 1d0) satratln = 1d0
    if (rhoa .le. 1d4) rhoaln = 1d4
    if (rhoa .gt. 1d13) rhoaln = 1d13

    ! critical cluster composition
    x_n=  7.9036365428891719d-1 - 2.8414059650092153d-3*tln + 1.4976802556584141d-2*LOG(satratln) &
         & - 2.4511581740839115d-4*tln*LOG(satratln) + 3.4319869471066424d-3 *LOG(satratln)**2     &
         & - 2.8799393617748428d-5*tln*LOG(satratln)**2 + 3.0174314126331765d-4*LOG(satratln)**3 &
         & - 2.2673492408841294d-6*tln*LOG(satratln)**3 - 4.3948464567032377d-3*LOG(rhoaln)&
         & + 5.3305314722492146d-5*tln*LOG(rhoaln)
    x_n=MIN(MAX(x_n,1.d-30),1.d0)

    ! kinetic limit check
    if (satratln .ge. 1.d-2 .and. satratln .le. 1.d0) then
       kinrhotresn=exp(7.8920778706888086d+1 + 7.3665492897447082d0*satratln - 1.2420166571163805d+4/tln &
            & + (-6.1831234251470971d+2*satratln)/tln - 2.4501159970109945d-2*tln &
            & -1.3463066443605762d-2*satratln*tln + 8.3736373989909194d-06*tln**2   &
            & -1.4673887785408892d0*Log(satratln) + (-3.2141890006517094d+1*Log(satratln))/tln &
            & + 2.7137429081917556d-3*tln*Log(satratln)) !1/cm3
       if(kinrhotresn.lt.rhoaln) kinetic_n=.true.
    endif
    if (satratln .ge. 1.d-4  .and. satratln .lt. 1.d-2) then
       kinrhotresn=exp(7.9074383049843647d+1 - 2.8746005462158347d+1*satratln - 1.2070272068458380d+4/tln &
            & + (-5.9205040320056632d+3*satratln)/tln - 2.4800372593452726d-2*tln &
            & -4.3983007681295948d-2*satratln*tln + 2.5943854791342071d-5*tln**2   &
            & -2.3141363245211317d0*Log(satratln) + (9.9186787997857735d+1*Log(satratln))/tln &
            & + 5.6819382556144681d-3*tln*Log(satratln)) !1/cm3
       if(kinrhotresn.lt.rhoaln) kinetic_n=.true.
    endif
    if (satratln .ge. 5.d-6  .and. satratln .lt. 1.d-4) then
       kinrhotresn=exp(8.5599712000361677d+1 + 2.7335119660796581d+3*satratln - 1.1842350246291651d+4/tln &
            & + (-1.2439843468881438d+6*satratln)/tln - 5.4536964974944230d-2*tln &
            & + 5.0886987425326087d0*satratln*tln + 7.1964722655507067d-5*tln**2   &
            & -2.4472627526306372d0*Log(satratln) + (1.7561478001423779d+2*Log(satratln))/tln &
            & + 6.2640132818141811d-3*tln*Log(satratln)) !1/cm3
       if(kinrhotresn.lt.rhoaln) kinetic_n=.true.
    endif

    if (kinetic_n) then
       ! dimer formation rate at the unbounded T and [H2SO4]
       jnuc_n = nucleate_binary_dimer_collision_rate(t, rhoa)
       ntot_n = 1d0
       ! the critical cluster contains one molecule, but the produced
       ! cluster contains two
       na_n = 1d0
       x_n = na_n / ntot_n
       rc_n = NUCLEATE_MAATTANEN_MONOMER_RADIUS
    else
       jnuc_n= 2.1361182605986115d-1 + 3.3827029855551838d0 *tln -3.2423555796175563d-2*tln**2 +  &
            &  7.0120069477221989d-5*tln**3 +8.0286874752695141d0/x_n +  &
            &  -2.6939840579762231d-1*LOG(satratln) +1.6079879299099518d0*tln*LOG(satratln) +  &
            &  -1.9667486968141933d-2*tln**2*LOG(satratln) +  &
            &  5.5244755979770844d-5*tln**3*LOG(satratln) + (7.8884704837892468d0*LOG(satratln))/x_n +  &
            &  4.6374659198909596d0*LOG(satratln)**2 - 8.2002809894792153d-2*tln*LOG(satratln)**2 +  &
            &  8.5077424451172196d-4*tln**2*LOG(satratln)**2 +  &
            &  -2.6518510168987462d-6*tln**3*LOG(satratln)**2 +  &
            &  (-1.4625482500575278d0*LOG(satratln)**2)/x_n - 5.2413002989192037d-1*LOG(satratln)**3 +  &
            &  5.2755117653715865d-3*tln*LOG(satratln)**3 +  &
            &  -2.9491061332113830d-6*tln**2*LOG(satratln)**3 +  &
            &  -2.4815454194486752d-8*tln**3*LOG(satratln)**3 +  &
            &  (-5.2663760117394626d-2*LOG(satratln)**3)/x_n +  &
            &  1.6496664658266762d0*LOG(rhoaln) +  &
            &  -8.0809397859218401d-1*tln*LOG(rhoaln) +  &
            &  8.9302927091946642d-3*tln**2*LOG(rhoaln) +  &
            &  -1.9583649496497497d-5*tln**3*LOG(rhoaln) +  &
            &  (-8.9505572676891685d0*LOG(rhoaln))/x_n +  &
            &  -3.0025283601622881d+1*LOG(satratln)*LOG(rhoaln) +  &
            &  3.0783365644763633d-1*tln*LOG(satratln)*LOG(rhoaln) +  &
            &  -7.4521756337984706d-4*tln**2*LOG(satratln)*LOG(rhoaln) +  &
            &  -5.7651433870681853d-7*tln**3*LOG(satratln)*LOG(rhoaln) +  &
            &  (1.2872868529673207d0*LOG(satratln)*LOG(rhoaln))/x_n +  &
            &  -6.1739867501526535d-1*LOG(satratln)**2*LOG(rhoaln) +  &
            &  7.2347385705333975d-3*tln*LOG(satratln)**2*LOG(rhoaln) +  &
            &  -3.0640494530822439d-5*tln**2*LOG(satratln)**2*LOG(rhoaln) +  &
            &  6.5944609194346214d-8*tln**3*LOG(satratln)**2*LOG(rhoaln) +  &
            &  (-2.8681650332461055d-2*LOG(satratln)**2*LOG(rhoaln))/x_n +  &
            &  6.5213802375160306d0*LOG(rhoaln)**2 +  &
            &  -4.7907162004793016d-2*tln*LOG(rhoaln)**2 +  &
            &  -1.0727890114215117d-4*tln**2*LOG(rhoaln)**2 +  &
            &  5.6401818280534507d-7*tln**3*LOG(rhoaln)**2 +  &
            &  (5.4113070888923009d-1*LOG(rhoaln)**2)/x_n +  &
            &  5.2062808476476330d-1*LOG(satratln)*LOG(rhoaln)**2 +  &
            &  -6.0696882500824584d-3*tln*LOG(satratln)*LOG(rhoaln)**2 +  &
            &  2.3851383302608477d-5*tln**2*LOG(satratln)*LOG(rhoaln)**2 +  &
            &  -1.5243837103067096d-8*tln**3*LOG(satratln)*LOG(rhoaln)**2 +  &
            &  (-5.6543192378015687d-2*LOG(satratln)*LOG(rhoaln)**2)/x_n +  &
            &  -1.1630806410696815d-1*LOG(rhoaln)**3 +  &
            &  1.3806404273119610d-3*tln*LOG(rhoaln)**3 +  &
            &  -2.0199865087650833d-6*tln**2*LOG(rhoaln)**3 +  &
            &  -3.0200284885763192d-9*tln**3*LOG(rhoaln)**3 +  &
            &  (-6.9425267104126316d-3*LOG(rhoaln)**3)/x_n
       jnuc_n=EXP(jnuc_n)
       ntot_n =-3.5863435141979573d-3 - 1.0098670235841110d-1 *tln + 8.9741268319259721d-4 *tln**2 - 1.4855098605195757d-6*tln**3 &
            &   - 1.2080330016937095d-1/x_n + 1.1902674923928015d-3*LOG(satratln) - 1.9211358507172177d-2*tln*LOG(satratln) +  &
            &   2.4648094311204255d-4*tln**2*LOG(satratln) - 7.5641448594711666d-7*tln**3*LOG(satratln) +  &
            &   (-2.0668639384228818d-02*LOG(satratln))/x_n - 3.7593072011595188d-2*LOG(satratln)**2 + &
            &   9.0993182774415718d-4 *tln*LOG(satratln)**2 +&
            &   -9.5698412164297149d-6*tln**2*LOG(satratln)**2 + 3.7163166416110421d-8*tln**3*LOG(satratln)**2 +  &
            &   (1.1026579525210847d-2*LOG(satratln)**2)/x_n + 1.1530844115561925d-2 *LOG(satratln)**3 +  &
            &   - 1.8083253906466668d-4 *tln*LOG(satratln)**3 + 8.0213604053330654d-7*tln**2*LOG(satratln)**3 +  &
            &   -8.5797885383051337d-10*tln**3*LOG(satratln)**3 + (1.0243693899717402d-3*LOG(satratln)**3)/x_n +  &
            &   -1.7248695296299649d-2*LOG(rhoaln) + 1.1294004162437157d-2*tln*LOG(rhoaln) +  &
            &   -1.2283640163189278d-4*tln**2*LOG(rhoaln) + 2.7391732258259009d-7*tln**3*LOG(rhoaln) +  &
            &   (6.8505583974029602d-2*LOG(rhoaln))/x_n +2.9750968179523635d-1*LOG(satratln)*LOG(rhoaln) +  &
            &   -3.6681154503992296d-3 *tln*LOG(satratln)*LOG(rhoaln) + 1.0636473034653114d-5*tln**2*LOG(satratln)*LOG(rhoaln) +  &
            &   5.8687098466515866d-9*tln**3*LOG(satratln)*LOG(rhoaln) + (-5.2028866094191509d-3*LOG(satratln)*LOG(rhoaln))/x_n +  &
            &   7.6971988880587231d-4*LOG(satratln)**2*LOG(rhoaln) - 2.4605575820433763d-5*tln*LOG(satratln)**2*LOG(rhoaln) +  &
            &   2.3818484400893008d-7*tln**2*LOG(satratln)**2*LOG(rhoaln) +  &
            &   -8.8474102392445200d-10*tln**3*LOG(satratln)**2*LOG(rhoaln) +  &
            &   (-1.6640566678168968d-4*LOG(satratln)**2*LOG(rhoaln))/x_n - 7.7390093776705471d-2*LOG(rhoaln)**2 +  &
            &   5.8220163188828482d-4*tln*LOG(rhoaln)**2 + 1.2291679321523287d-6*tln**2*LOG(rhoaln)**2 +  &
            &   -7.4690997508075749d-9*tln**3*LOG(rhoaln)**2 + (-5.6357941220497648d-3*LOG(rhoaln)**2)/x_n +  &
            &   -4.7170109625089768d-3*LOG(satratln)*LOG(rhoaln)**2 + 6.9828868534370193d-5*tln*LOG(satratln)*LOG(rhoaln)**2 +  &
            &   -3.1738912157036403d-7*tln**2*LOG(satratln)*LOG(rhoaln)**2 +  &
            &   2.3975538706787416d-10*tln**3*LOG(satratln)*LOG(rhoaln)**2 +  &
            &   (4.2304213386288567d-4*LOG(satratln)*LOG(rhoaln)**2)/x_n + 1.3696520973423231d-3*LOG(rhoaln)**3 +  &
            &   -1.6863387574788199d-5*tln*LOG(rhoaln)**3 + 2.7959499278844516d-8*tln**2*LOG(rhoaln)**3 +  &
            &   3.9423927013227455d-11*tln**3*LOG(rhoaln)**3 + (8.6136359966337272d-5*LOG(rhoaln)**3)/x_n
       ntot_n=EXP(ntot_n)
       rc_n=EXP(-22.378268374023630d0+0.44462953606125100d0*x_n+0.33499495707849131d0*LOG(ntot_n)) !in meters
       na_n=x_n*ntot_n
       if (na_n .lt. 1d0) then
          na_n = 1d0
       end if
    end if

    ! set the neutral nucleation rate to zero if less than 1e-7 cm^{-3} s^{-1}
    if (jnuc_n .lt. 1d-7) then
       jnuc_n = 0d0
    end if

  end subroutine nucleate_binary_maattanen2018_cluster

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Maattanen et al. (2018) neutral formation rate and new particle.
  !!
  !! Nucleation regime: the critical cluster (at least one acid
  !! molecule). Kinetic regime: the dimer formed by collisions, i.e. two
  !! acid molecules and radius \f$0.3\ {\rm nm} \cdot 2^{1/3}\f$.
  subroutine nucleate_binary_maattanen2018(temp, rh, h2so4, rate, n_acid, &
       n_total, radius, kinetic)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Relative humidity over liquid water (1).
    real(kind=dp), intent(in) :: rh
    !> Sulfuric acid concentration (cm^{-3}).
    real(kind=dp), intent(in) :: h2so4
    !> Formation rate (cm^{-3} s^{-1}).
    real(kind=dp), intent(out) :: rate
    !> H2SO4 molecules per new particle.
    real(kind=dp), intent(out) :: n_acid
    !> H2SO4 + H2O molecules per new particle.
    real(kind=dp), intent(out) :: n_total
    !> Wet radius of the new particle (m).
    real(kind=dp), intent(out) :: radius
    !> Whether the kinetic regime applies.
    logical, intent(out) :: kinetic

    real(kind=dp) :: ntot, x, na, rc

    call nucleate_binary_maattanen2018_cluster(temp, rh, h2so4, rate, ntot, &
         x, na, rc, kinetic)
    if (kinetic) then
       n_acid = 2d0
       n_total = 2d0
       radius = NUCLEATE_MAATTANEN_MONOMER_RADIUS * 2d0**(1d0 / 3d0)
    else
       n_acid = na
       n_total = max(ntot, na)
       radius = rc
    end if

  end subroutine nucleate_binary_maattanen2018

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dunne et al. (2016) binary-neutral formation rate
  !> \f$J_{1.7} = k(T) ([{\rm H_2SO_4}]/10^6)^p\f$ with
  !> \f$k = \exp(u - \exp(v (T/1000 - w)))\f$ (cm^{-3} s^{-1}), uncapped.
  real(kind=dp) function nucleate_binary_dunne2016_rate(temp, h2so4)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Sulfuric acid concentration (cm^{-3}).
    real(kind=dp), intent(in) :: h2so4

    real(kind=dp) :: c, k

    c = max(h2so4, 0d0) * 1d-6
    k = exp(NUCLEATE_DUNNE_U_BN &
         - exp(NUCLEATE_DUNNE_V_BN * (temp / 1000d0 - NUCLEATE_DUNNE_W_BN)))
    nucleate_binary_dunne2016_rate = k * c**NUCLEATE_DUNNE_P_BN

  end function nucleate_binary_dunne2016_rate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> H2SO4 mass fraction of a binary solution droplet in equilibrium with
  !> water activity \c rh (Tabazadeh et al. 1997), bounded to 15-99 wt%.
  real(kind=dp) function nucleate_binary_tabazadeh1997_mass_frac(temp, rh)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Water activity, i.e. relative humidity (1).
    real(kind=dp), intent(in) :: rh

    real(kind=dp) :: aw, y1, y2, molal, wp

    aw = min(max(rh, 1d-12), 1d0)
    if (aw <= 0.05d0) then
       y1 = 12.372089320d0 * aw**(-0.16125516114d0) - 30.490657554d0 * aw &
            - 2.1133114241d0
       y2 = 13.455394705d0 * aw**(-0.19213122550d0) - 34.285174607d0 * aw &
            - 1.7620073078d0
    elseif (aw <= 0.85d0) then
       y1 = 11.820654354d0 * aw**(-0.20786404244d0) - 4.8073063730d0 * aw &
            - 5.1727540348d0
       y2 = 12.891938068d0 * aw**(-0.23233847708d0) - 6.4261237757d0 * aw &
            - 4.9005471319d0
    else
       y1 = -180.06541028d0 * aw**(-0.38601102592d0) - 93.317846778d0 * aw &
            + 273.88132245d0
       y2 = -176.95814097d0 * aw**(-0.36257048154d0) - 90.469744201d0 * aw &
            + 267.45509988d0
    end if
    molal = y1 + (temp - 190d0) * (y2 - y1) / 70d0
    wp = 9800d0 * molal / (98d0 * molal + 1000d0)
    nucleate_binary_tabazadeh1997_mass_frac = min(max(wp, 15d0), 99d0) * 0.01d0

  end function nucleate_binary_tabazadeh1997_mass_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> H2SO4-H2O solution density (g cm^{-3}) of Vehkamaki et al. (2002)
  !> at acid mass fraction \c w.
  real(kind=dp) function nucleate_binary_vehkamaki2002_density(temp, w)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> H2SO4 mass fraction (1).
    real(kind=dp), intent(in) :: w

    real(kind=dp) :: a, b, c

    a = 0.7681724d0 + 2.184714d0 * w + 7.163002d0 * w**2 - 44.31447d0 * w**3 &
         + 88.74606d0 * w**4 - 75.73729d0 * w**5 + 23.43228d0 * w**6
    b = 1.808225d-3 - 9.294656d-3 * w - 3.742148d-2 * w**2 &
         + 2.565321d-1 * w**3 - 5.362872d-1 * w**4 + 4.857736d-1 * w**5 &
         - 1.629592d-1 * w**6
    c = -3.478524d-6 + 1.335867d-5 * w + 5.195706d-5 * w**2 &
         - 3.717636d-4 * w**3 + 7.990811d-4 * w**4 - 7.458060d-4 * w**5 &
         + 2.581390d-4 * w**6
    nucleate_binary_vehkamaki2002_density = a + b * temp + c * temp * temp

  end function nucleate_binary_vehkamaki2002_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Dunne et al. (2016) binary-neutral formation rate and its new
  !> particle: 1.4 nm geometric diameter, Tabazadeh et al. (1997) acid mass
  !> fraction at the ambient T and RH, Vehkamaki et al. (2002) density.
  subroutine nucleate_binary_dunne2016(temp, rh, h2so4, rate, n_acid, &
       n_total, radius, kinetic)

    !> Temperature (K).
    real(kind=dp), intent(in) :: temp
    !> Relative humidity over liquid water (1).
    real(kind=dp), intent(in) :: rh
    !> Sulfuric acid concentration (cm^{-3}).
    real(kind=dp), intent(in) :: h2so4
    !> Formation rate (cm^{-3} s^{-1}).
    real(kind=dp), intent(out) :: rate
    !> H2SO4 molecules per new particle.
    real(kind=dp), intent(out) :: n_acid
    !> H2SO4 + H2O molecules per new particle.
    real(kind=dp), intent(out) :: n_total
    !> Wet radius of the new particle (m).
    real(kind=dp), intent(out) :: radius
    !> Always false (no kinetic branch).
    logical, intent(out) :: kinetic

    real(kind=dp) :: w, rho, mass_g, n_water

    rate = nucleate_binary_dunne2016_rate(temp, h2so4)
    radius = 0.5d0 * NUCLEATE_DUNNE_DIAM
    w = nucleate_binary_tabazadeh1997_mass_frac(temp, rh)
    rho = nucleate_binary_vehkamaki2002_density(temp, w)
    mass_g = (4d0 / 3d0) * NUCLEATE_BINARY_PI * (radius * 100d0)**3 * rho
    n_acid = w * mass_g / NUCLEATE_BINARY_MW_H2SO4 * NUCLEATE_BINARY_AVOGADRO
    n_water = (1d0 - w) * mass_g / NUCLEATE_BINARY_MW_H2O &
         * NUCLEATE_BINARY_AVOGADRO
    n_total = n_acid + n_water
    kinetic = .false.

  end subroutine nucleate_binary_dunne2016

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_nucleate_binary
