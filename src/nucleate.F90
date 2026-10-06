! Copyright (C) 2010-2015 Matthew West
! Licensed under the GNU General Public License version 2 or (at your
! option) any later version. See the file COPYING for details.

!> \file
!> The pmc_nucleate module.

!> Aerosol nucleation functions.
module pmc_nucleate

  use pmc_env_state
  use pmc_aero_state
  use pmc_aero_data
  use pmc_gas_data
  use pmc_gas_state
  use pmc_nucleate_binary

  !> Type code for unknown or invalid nucleation type.
  integer, parameter :: NUCLEATE_TYPE_INVALID   = 0
  !> Type code for H2SO4 to SO4 nucleation with quadratic rate.
  integer, parameter :: NUCLEATE_TYPE_SULF_ACID = 1
  !> Type code for Maattanen et al. (2018) neutral H2SO4-H2O nucleation.
  integer, parameter :: NUCLEATE_TYPE_MAATTANEN2018 = 2
  !> Type code for Dunne et al. (2016) binary-neutral H2SO4-H2O nucleation.
  integer, parameter :: NUCLEATE_TYPE_DUNNE2016 = 3

  !> Source name for nucleated particles.
  character(len=AERO_SOURCE_NAME_LEN), parameter :: NUCLEATE_SOURCE_NAME &
       = "nucleate"

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Do nucleation of the type given by the first argument.
  subroutine nucleate(nucleate_type, nucleate_source, nucleate_weight_class, &
       nucleate_h2so4_fixed_conc, nucleate_trace_file, env_state, gas_data, &
       aero_data, aero_state, gas_state, del_t, allow_doubling, allow_halving)

    !> Type of nucleation.
    integer, intent(in) :: nucleate_type
    !> Nucleate source number.
    integer, intent(in) :: nucleate_source
    !> Nucleate weight class.
    integer, intent(in) :: nucleate_weight_class
    !> Fixed H2SO4 concentration for the binary H2SO4-H2O types
    !> (molecules m^{-3}); negative to use the gas state.
    real(kind=dp), intent(in) :: nucleate_h2so4_fixed_conc
    !> Nucleation trace file for the binary H2SO4-H2O types ("none" for
    !> no trace).
    character(len=*), intent(in) :: nucleate_trace_file
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    if (nucleate_type == NUCLEATE_TYPE_SULF_ACID) then
       call nucleate_sulf_acid(nucleate_source, nucleate_weight_class, &
            env_state, gas_data, aero_data, aero_state, gas_state, del_t, &
            allow_doubling, allow_halving)
    elseif ((nucleate_type == NUCLEATE_TYPE_MAATTANEN2018) &
         .or. (nucleate_type == NUCLEATE_TYPE_DUNNE2016)) then
       call nucleate_binary_h2so4(nucleate_type, nucleate_source, &
            nucleate_h2so4_fixed_conc, nucleate_trace_file, env_state, &
            gas_data, aero_data, aero_state, gas_state, del_t, &
            allow_doubling, allow_halving)
    else
       call die_msg(983831728, &
            "unknown nucleation type: " &
            // trim(integer_to_string(nucleate_type)))
    end if

  end subroutine nucleate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Nucleate sulfuric acid into aerosol particles, using a power-law
  !> dependence, for time \c del_t.
  !!
  !! The modeled emission rate is \f$ J = K H^2 \f$, where \f$H\f$ is
  !! the concentration of \f$ \rm H_2SO_4 \f$ and \f$K\f$ is a
  !! constant coefficient.
  !!
  !! The reference is:
  !!
  !! C. Kuang, P. H. McMurry, A. V. McCormick, and F. L. Eisele
  !! (2008), Dependence of nucleation rates on sulfuric acid vapor
  !! concentration in diverse atmospheric locations,
  !! <i>J. Geophys. Res.</i>, 113, D10209, doi:<a
  !! href="http://dx.doi.org/10.1029/2007JD009253">10.1029/2007JD009253</a>.
  subroutine nucleate_sulf_acid(nucleate_source, nucleate_weight_class, &
       env_state, gas_data, aero_data, aero_state, gas_state, del_t, &
       allow_doubling, allow_halving)

    !> Nucleate source number.
    integer, intent(in) :: nucleate_source
    !> Nucleate weight class.
    integer, intent(in) :: nucleate_weight_class
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    real(kind=dp), parameter :: nucleate_coeff = 1d-18 ! K (m^3 s^{-1})
    real(kind=dp), parameter :: nucleate_diam = 1d-9   ! diameter of new
                                                       ! particles (m)

    integer :: i_gas_h2so4, i_aero_so4, n_samp, i_samp, i_bin, i_group, n_group
    integer :: i_class, i_source
    real(kind=dp) :: sulf_acid_conc, nucleate_rate, n_samp_avg
    real(kind=dp) :: total_so4_vol, so4_vol, h2so4_removed_conc
    type(aero_particle_t) :: aero_particle

    ! look up the species numbers
    i_gas_h2so4 = gas_data_spec_by_name(gas_data, "H2SO4")
    call assert_msg(886839228, i_gas_h2so4 > 0, &
         "nucleate_sulf_acid requires H2SO4 as a gas species")
    i_aero_so4 = aero_data_spec_by_name(aero_data, "SO4")
    call assert_msg(551966998, i_aero_so4 > 0, &
         "nucleate_sulf_acid requires SO4 as an aerosol species")

    ! H2SO4 concentration in molecules / m^3
    sulf_acid_conc = env_state_ppb_to_conc(env_state, &
         gas_state%mix_rat(i_gas_h2so4))

    ! particle nucleation rate in (particles m^{-3} s^{-1})
    nucleate_rate = nucleate_coeff * sulf_acid_conc**2

    ! weight class to nucleate into
    i_class = aero_state_weight_class_for_source(aero_state, nucleate_source)
    i_source = nucleate_source

    ! add particles to each weight group
    total_so4_vol = 0d0
    do i_group = 1,aero_weight_array_n_group(aero_state%awa)
       ! adjust weight if necessary
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), diam2rad(nucleate_diam))
       call aero_state_prepare_weight_for_add(aero_state, aero_data, &
            i_group, i_class, n_samp_avg, allow_doubling, allow_halving)

       ! determine number of nucleated particles
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), diam2rad(nucleate_diam))
       n_samp = rand_poisson(n_samp_avg)

       ! create the particles
       do i_samp = 1,n_samp
          so4_vol = aero_data_diam2vol(aero_data, nucleate_diam)
          total_so4_vol = total_so4_vol + so4_vol

          call aero_particle_zero(aero_particle, aero_data)
          call aero_particle_set_component(aero_particle, i_source, &
               env_state%elapsed_time)
          aero_particle%vol(i_aero_so4) = so4_vol
          call aero_particle_new_id(aero_particle)
          call aero_particle_set_weight(aero_particle, i_group, i_class)
          call aero_state_add_particle(aero_state, aero_particle, aero_data)
       end do
    end do

    ! remove gases that formed new particles
    h2so4_removed_conc = &
         total_so4_vol &                        ! volume of SO4
         * aero_weight_array_num_conc_at_radius(aero_state%awa, i_class, &
         diam2rad(nucleate_diam)) &             ! volume conc of SO4
         * aero_data%density(i_aero_so4) &      ! mass conc of SO4
         / aero_data%molec_weight(i_aero_so4) & ! mole conc of SO4
         * const%avagadro                       ! molecule conc of SO4
    gas_state%mix_rat(i_gas_h2so4) = gas_state%mix_rat(i_gas_h2so4) &
         - env_state_conc_to_ppb(env_state, h2so4_removed_conc)
    if (gas_state%mix_rat(i_gas_h2so4) < 0d0) then
       gas_state%mix_rat(i_gas_h2so4) = 0d0
    end if

  end subroutine nucleate_sulf_acid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Nucleate neutral binary H2SO4-H2O particles for time \c del_t, with
  !> the Maattanen et al. (2018) or Dunne et al. (2016) scheme of
  !> pmc_nucleate_binary.
  !!
  !! The scheme gives the formation rate \f$J\f$ together with the
  !! particle it forms (acid molecules \f$n_a\f$ and total molecules
  !! \f$n_t\f$ at the ambient temperature, relative humidity and
  !! H2SO4 concentration). Each new particle holds \f$n_a\f$ molecules
  !! of SO4 (one SO4 per H2SO4 molecule, at the SO4 molecular weight and
  !! density of \c aero_data) and \f$n_t - n_a\f$ molecules of H2O (at
  !! the H2O molecular weight and density of \c aero_data); its PartMC
  !! radius follows from these volumes and is used for the particle
  !! weighting. New particles are Poisson-sampled into the "nucleate"
  !! weight class and gas-phase H2SO4 is depleted by the acid they
  !! contain, as in nucleate_sulf_acid(). No survival factor is applied:
  !! the particle is inserted at the rate the scheme gives for it.
  !!
  !! If \c h2so4_fixed_conc is non-negative the gas-phase H2SO4 is reset
  !! to it before the rate is evaluated (an external reservoir holding
  !! the vapor fixed across timesteps); the depletion by the new
  !! particles is still applied afterwards. Unless \c trace_file is
  !! "none", one line per call is appended to it: elapsed time (s),
  !! timestep (s), temperature (K), relative humidity (1), H2SO4
  !! concentration used (m^{-3}), formation rate (m^{-3} s^{-1}), acid
  !! molecules, total molecules, scheme radius (m), kinetic flag (1 or 0),
  !! PartMC radius of the new particle (m), number of computational
  !! particles added, number concentration added (m^{-3}), H2SO4 removed
  !! from the gas (molecules m^{-3}) and gas-phase H2SO4 afterwards
  !! (molecules m^{-3}).
  subroutine nucleate_binary_h2so4(nucleate_type, nucleate_source, &
       h2so4_fixed_conc, trace_file, env_state, gas_data, aero_data, &
       aero_state, gas_state, del_t, allow_doubling, allow_halving)

    !> Type of nucleation (NUCLEATE_TYPE_MAATTANEN2018 or
    !> NUCLEATE_TYPE_DUNNE2016).
    integer, intent(in) :: nucleate_type
    !> Nucleate source number.
    integer, intent(in) :: nucleate_source
    !> Fixed H2SO4 concentration (molecules m^{-3}); negative to use the
    !> gas state.
    real(kind=dp), intent(in) :: h2so4_fixed_conc
    !> Trace file ("none" for no trace).
    character(len=*), intent(in) :: trace_file
    !> Environment state.
    type(env_state_t), intent(in) :: env_state
    !> Gas data.
    type(gas_data_t), intent(in) :: gas_data
    !> Aerosol data.
    type(aero_data_t), intent(in) :: aero_data
    !> Aerosol state.
    type(aero_state_t), intent(inout) :: aero_state
    !> Gas state.
    type(gas_state_t), intent(inout) :: gas_state
    !> Time to perform nucleation for.
    real(kind=dp), intent(in) :: del_t
    !> Whether to allow doubling of the population.
    logical, intent(in) :: allow_doubling
    !> Whether to allow halving of the population.
    logical, intent(in) :: allow_halving

    integer :: i_gas_h2so4, i_aero_so4, i_aero_h2o, i_group, i_class
    integer :: n_samp, i_samp, n_samp_total, unit, ios
    real(kind=dp) :: h2so4_conc, rate_cm3, nucleate_rate, n_acid, n_total
    real(kind=dp) :: scheme_radius, so4_vol, h2o_vol, radius, n_samp_avg
    real(kind=dp) :: num_conc, added_num_conc, h2so4_removed_conc
    real(kind=dp) :: kinetic_flag
    logical :: kinetic
    type(aero_particle_t) :: aero_particle

    ! look up the species numbers
    i_gas_h2so4 = gas_data_spec_by_name(gas_data, "H2SO4")
    call assert_msg(301928374, i_gas_h2so4 > 0, &
         "nucleate_binary_h2so4 requires H2SO4 as a gas species")
    i_aero_so4 = aero_data_spec_by_name(aero_data, "SO4")
    call assert_msg(301928375, i_aero_so4 > 0, &
         "nucleate_binary_h2so4 requires SO4 as an aerosol species")
    i_aero_h2o = aero_data%i_water
    call assert_msg(301928376, i_aero_h2o > 0, &
         "nucleate_binary_h2so4 requires H2O as an aerosol species")

    ! H2SO4 concentration in molecules / m^3
    if (h2so4_fixed_conc >= 0d0) then
       h2so4_conc = h2so4_fixed_conc
       gas_state%mix_rat(i_gas_h2so4) = env_state_conc_to_ppb(env_state, &
            h2so4_fixed_conc)
    else
       h2so4_conc = env_state_ppb_to_conc(env_state, &
            gas_state%mix_rat(i_gas_h2so4))
    end if

    ! formation rate and new particle (cm^{-3} units in the schemes)
    if (nucleate_type == NUCLEATE_TYPE_MAATTANEN2018) then
       call nucleate_binary_maattanen2018(env_state%temp, &
            env_state%rel_humid, h2so4_conc * 1d-6, rate_cm3, n_acid, &
            n_total, scheme_radius, kinetic)
    elseif (nucleate_type == NUCLEATE_TYPE_DUNNE2016) then
       call nucleate_binary_dunne2016(env_state%temp, env_state%rel_humid, &
            h2so4_conc * 1d-6, rate_cm3, n_acid, n_total, scheme_radius, &
            kinetic)
    else
       call die_msg(301928377, "unknown binary nucleation type: " &
            // trim(integer_to_string(nucleate_type)))
    end if
    ! particle nucleation rate in (particles m^{-3} s^{-1})
    nucleate_rate = rate_cm3 * 1d6

    ! species volumes of one new particle (m^3)
    so4_vol = n_acid / const%avagadro * aero_data%molec_weight(i_aero_so4) &
         / aero_data%density(i_aero_so4)
    h2o_vol = (n_total - n_acid) / const%avagadro &
         * aero_data%molec_weight(i_aero_h2o) / aero_data%density(i_aero_h2o)
    radius = aero_data_vol2rad(aero_data, so4_vol + h2o_vol)

    ! weight class to nucleate into
    i_class = aero_state_weight_class_for_source(aero_state, nucleate_source)

    ! add particles to each weight group
    n_samp_total = 0
    do i_group = 1,aero_weight_array_n_group(aero_state%awa)
       ! adjust weight if necessary
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), radius)
       call aero_state_prepare_weight_for_add(aero_state, aero_data, &
            i_group, i_class, n_samp_avg, allow_doubling, allow_halving)

       ! determine number of nucleated particles
       n_samp_avg = nucleate_rate * del_t / aero_weight_num_conc_at_radius( &
            aero_state%awa%weight(i_group, i_class), radius)
       n_samp = rand_poisson(n_samp_avg)
       n_samp_total = n_samp_total + n_samp

       ! create the particles
       do i_samp = 1,n_samp
          call aero_particle_zero(aero_particle, aero_data)
          call aero_particle_set_component(aero_particle, nucleate_source, &
               env_state%elapsed_time)
          aero_particle%vol(i_aero_so4) = so4_vol
          aero_particle%vol(i_aero_h2o) = h2o_vol
          call aero_particle_new_id(aero_particle)
          call aero_particle_set_weight(aero_particle, i_group, i_class)
          call aero_state_add_particle(aero_state, aero_particle, aero_data)
       end do
    end do

    ! remove gases that formed new particles
    num_conc = aero_weight_array_num_conc_at_radius(aero_state%awa, i_class, &
         radius)
    added_num_conc = real(n_samp_total, kind=dp) * num_conc
    h2so4_removed_conc = &
         added_num_conc * so4_vol &               ! volume conc of SO4
         * aero_data%density(i_aero_so4) &        ! mass conc of SO4
         / aero_data%molec_weight(i_aero_so4) &   ! mole conc of SO4
         * const%avagadro                         ! molecule conc of SO4
    gas_state%mix_rat(i_gas_h2so4) = gas_state%mix_rat(i_gas_h2so4) &
         - env_state_conc_to_ppb(env_state, h2so4_removed_conc)
    if (gas_state%mix_rat(i_gas_h2so4) < 0d0) then
       gas_state%mix_rat(i_gas_h2so4) = 0d0
    end if

    if (trim(trace_file) /= "none") then
       kinetic_flag = 0d0
       if (kinetic) kinetic_flag = 1d0
       unit = get_unit()
       open(unit=unit, file=trace_file, status='unknown', &
            position='append', action='write', iostat=ios)
       call assert_msg(301928378, ios == 0, &
            'unable to open nucleation trace file ' // trim(trace_file))
       write(unit, '(15ES25.16E3)') env_state%elapsed_time, del_t, &
            env_state%temp, env_state%rel_humid, h2so4_conc, nucleate_rate, &
            n_acid, n_total, scheme_radius, kinetic_flag, radius, &
            real(n_samp_total, kind=dp), added_num_conc, &
            h2so4_removed_conc, &
            env_state_ppb_to_conc(env_state, gas_state%mix_rat(i_gas_h2so4))
       call close_file(unit)
    end if

  end subroutine nucleate_binary_h2so4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine spec_file_read_nucleate_type(file, aero_data, nucleate_type, &
       nucleate_source, nucleate_weight_class, nucleate_h2so4_fixed_conc, &
       nucleate_trace_file)

    !> Spec file.
    type(spec_file_t), intent(inout) :: file
    !> Aerosol data.
    type(aero_data_t), intent(inout) :: aero_data
    !> Nucleate type.
    integer, intent(out) :: nucleate_type
    !> Nucleate source number.
    integer, intent(out) :: nucleate_source
    !> Nucleate weight class.
    integer, intent(out) :: nucleate_weight_class
    !> Fixed H2SO4 concentration for the binary H2SO4-H2O types
    !> (molecules m^{-3}); negative to use the gas state.
    real(kind=dp), intent(out) :: nucleate_h2so4_fixed_conc
    !> Nucleation trace file for the binary H2SO4-H2O types ("none" for
    !> no trace).
    character(len=*), intent(out) :: nucleate_trace_file

    character(len=SPEC_LINE_MAX_VAR_LEN) :: nucleate_type_name

    !> \page input_format_nucleate Input File Format: Nucleation Parameterization
    !!
    !! The nucleation parameterization is specified by the parameter:
    !!   - \b nucleate (string): the type of nucleation
    !!     parameterization --- must be one of: "none" for no
    !!     nucleation; "sulf_acid" for the nucleate_sulf_acid()
    !!     parameterization; or "maattanen2018" or "dunne2016" for the
    !!     neutral binary H2SO4-H2O schemes of nucleate_binary_h2so4()
    !!     (Maattanen et al. 2018 neutral channel; Dunne et al. 2016
    !!     binary-neutral channel), which require H2SO4 as a gas species
    !!     and SO4 and H2O as aerosol species
    !!
    !! If \c nucleate is "maattanen2018" or "dunne2016" then the
    !! following parameters must also be provided:
    !!   - \b nucleate_h2so4_fixed_conc (real, unit molecules m^{-3}):
    !!     if non-negative, the gas-phase H2SO4 concentration is reset to
    !!     this value before each nucleation step (an external reservoir);
    !!     if negative, the evolving gas state is used
    !!   - \b nucleate_trace (string): name of a file to which one line
    !!     per nucleation step is appended (see nucleate_binary_h2so4()),
    !!     or "none"
    !!
    !! See also:
    !!   - \ref spec_file_format --- the input file text format

    call spec_file_read_string(file, 'nucleate', nucleate_type_name)
    nucleate_h2so4_fixed_conc = -1d0
    nucleate_trace_file = "none"
    if (nucleate_type_name == 'sulf_acid') then
       nucleate_type = NUCLEATE_TYPE_SULF_ACID
    elseif (nucleate_type_name == 'maattanen2018') then
       nucleate_type = NUCLEATE_TYPE_MAATTANEN2018
    elseif (nucleate_type_name == 'dunne2016') then
       nucleate_type = NUCLEATE_TYPE_DUNNE2016
    else
       call spec_file_die_msg(707263678, file, "unknown nucleate type: " &
            // trim(nucleate_type_name))
    end if
    nucleate_source = aero_data_source_by_name(aero_data, &
         NUCLEATE_SOURCE_NAME)
    nucleate_weight_class = aero_data_weight_class_by_name(aero_data, &
         NUCLEATE_SOURCE_NAME)
    if (nucleate_type /= NUCLEATE_TYPE_SULF_ACID) then
       call spec_file_read_real(file, 'nucleate_h2so4_fixed_conc', &
            nucleate_h2so4_fixed_conc)
       call spec_file_read_string(file, 'nucleate_trace', nucleate_trace_file)
    end if

  end subroutine spec_file_read_nucleate_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pmc_nucleate
