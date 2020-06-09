!> A column-wise toolbox for implementing neutral diffusion
module MOM_neutral_diffusion

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,             only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,             only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_domains,               only : pass_var
use MOM_diag_mediator,         only : diag_ctrl, time_type
use MOM_diag_mediator,         only : post_data, register_diag_field
use MOM_EOS,                   only : EOS_type, EOS_manual_init, EOS_domain
use MOM_EOS,                   only : calculate_density, calculate_density_derivs
use MOM_EOS,                   only : extract_member_EOS, EOS_LINEAR, EOS_TEOS10, EOS_WRIGHT
use MOM_error_handler,         only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,           only : get_param, log_version, param_file_type
use MOM_file_parser,           only : openParameterBlock, closeParameterBlock
use MOM_grid,                  only : ocean_grid_type
use MOM_remapping,             only : remapping_CS, initialize_remapping
use MOM_remapping,             only : extract_member_remapping_CS, build_reconstructions_1d
use MOM_remapping,             only : average_value_ppoly, remappingSchemesDoc, remappingDefaultScheme
use MOM_tracer_registry,       only : tracer_registry_type, tracer_type
use MOM_unit_scaling,          only : unit_scale_type
use MOM_verticalGrid,          only : verticalGrid_type
use polynomial_functions,      only : evaluation_polynomial, first_derivative_polynomial, analytic_derivative
use PPM_functions,             only : PPM_reconstruction, PPM_boundary_extrapolation
use regrid_edge_values,        only : edge_values_implicit_h4
use MOM_CVMix_KPP,             only : KPP_get_BLD, KPP_CS
use MOM_energetic_PBL,         only : energetic_PBL_get_MLD, energetic_PBL_CS
use MOM_diabatic_driver,       only : diabatic_CS, extract_diabatic_member
use MOM_lateral_boundary_diffusion, only : boundary_k_range, SURFACE, BOTTOM

use iso_fortran_env, only : stdout=>output_unit, stderr=>error_unit

implicit none ; private

#include <MOM_memory.h>

public neutral_diffusion, neutral_diffusion_init, neutral_diffusion_end
public neutral_diffusion_calc_coeffs
public neutral_diffusion_unit_tests

type :: column_properties 
  real, allocatable, dimension(:,:) :: T_at_interface !< Temperature at top (1) and bottom (2) interfaces
  real, allocatable, dimension(:,:) :: S_at_interface !< Salinity at top (1) and bottom (2) interfaces
  real, allocatable, dimension(:,:) :: P_at_interface !< Pressure at top (1) and bottom (2) interfaces
  real, allocatable, dimension(:,:) :: T_poly         !< Polynomial coefficients for temperature 
  real, allocatable, dimension(:,:) :: S_poly         !< Polynomial coefficients for salinity
  logical, allocatable, dimension(:) :: stable_cell    !< Polynomial coefficients for salinity
  real, dimension(2)                :: dRdT           !< dRho_dT at the top (1) and bottom (2) interface
  real, dimension(2)                :: dRdS           !< dRho_dS at the top (1) and bottom (2) interface
  integer                           :: first_stable   !< Index of the first stable layer
  integer                           :: last_stable    !< Index of the last stable layer
end type column_properties

!> The control structure for the MOM_neutral_diffusion module
type, public :: neutral_diffusion_CS ; private
  integer :: nkp1     !< Number of interfaces for a column = nk + 1
  integer :: nsurf    !< Number of neutral surfaces
  integer :: deg = 2  !< Degree of polynomial used for reconstructions
  logical :: continuous_reconstruction = .true. !< True if using continuous PPM reconstruction at interfaces
  logical :: debug = .false. !< If true, write verbose debugging messages
  integer :: max_iter !< Maximum number of iterations if refine_position is defined
  real :: drho_tol    !< Convergence criterion representing density difference from true neutrality [R ~> kg m-3]
  real :: x_tol       !< Convergence criterion for how small an update of the position can be
  real :: ref_pres    !< Reference pressure, negative if using locally referenced neutral
                      !! density [R L2 T-2 ~> Pa]
  logical :: interior_only !< If true, only applies neutral diffusion in the ocean interior.
                      !! That is, the algorithm will exclude the surface and bottom boundary layers.
  ! Positions of neutral surfaces in both the u, v directions
  real,    allocatable, dimension(:,:,:) :: uPoL  !< Non-dimensional position with left layer uKoL-1, u-point
  real,    allocatable, dimension(:,:,:) :: uPoR  !< Non-dimensional position with right layer uKoR-1, u-point
  integer, allocatable, dimension(:,:,:) :: uKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a u-point
  integer, allocatable, dimension(:,:,:) :: uKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a u-point
  real,    allocatable, dimension(:,:,:) :: uHeff !< Effective thickness at u-point [H ~> m or kg m-2]
  real,    allocatable, dimension(:,:,:) :: vPoL  !< Non-dimensional position with left layer uKoL-1, v-point
  real,    allocatable, dimension(:,:,:) :: vPoR  !< Non-dimensional position with right layer uKoR-1, v-point
  integer, allocatable, dimension(:,:,:) :: vKoL  !< Index of left interface corresponding to neutral surface,
                                                  !! at a v-point
  integer, allocatable, dimension(:,:,:) :: vKoR  !< Index of right interface corresponding to neutral surface,
                                                  !! at a v-point
  real,    allocatable, dimension(:,:,:) :: vHeff !< Effective thickness at v-point [H ~> m or kg m-2]
  ! Coefficients of polynomial reconstructions for temperature and salinity
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_T !< Polynomial coefficients for temperature
  real,    allocatable, dimension(:,:,:,:) :: ppoly_coeffs_S !< Polynomial coefficients for salinity
  ! Variables needed for continuous reconstructions
  real,    allocatable, dimension(:,:,:) :: dRdT !< dRho/dT [R degC-1 ~> kg m-3 degC-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: dRdS !< dRho/dS [R ppt-1 ~> kg m-3 ppt-1] at interfaces
  real,    allocatable, dimension(:,:,:) :: Tint !< Interface T [degC]
  real,    allocatable, dimension(:,:,:) :: Sint !< Interface S [ppt]
  real,    allocatable, dimension(:,:,:) :: Pint !< Interface pressure [R L2 T-2 ~> Pa]
  ! Variables needed for discontinuous reconstructions
  real,    allocatable, dimension(:,:,:,:) :: T_i    !< Top edge reconstruction of temperature [degC]
  real,    allocatable, dimension(:,:,:,:) :: S_i    !< Top edge reconstruction of salinity [ppt]
  real,    allocatable, dimension(:,:,:,:) :: P_i    !< Interface pressures [R L2 T-2 ~> Pa]
  real,    allocatable, dimension(:,:,:,:) :: dRdT_i !< dRho/dT [R degC-1 ~> kg m-3 degC-1] at top edge
  real,    allocatable, dimension(:,:,:,:) :: dRdS_i !< dRho/dS [R ppt-1 ~> kg m-3 ppt-1] at top edge
  integer, allocatable, dimension(:,:)     :: ns     !< Number of interfacs in a column
  logical, allocatable, dimension(:,:,:) :: stable_cell !< True if the cell is stably stratified wrt to the next cell
  real :: R_to_kg_m3 = 1.0                   !< A rescaling factor translating density to kg m-3 for
                                             !! use in diagnostic messages [kg m-3 R-1 ~> 1].
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                             !! regulate the timing of diagnostic output.
  integer :: neutral_pos_method              !< Method to find the position of a neutral surface within the layer
  character(len=40)  :: delta_rho_form       !< Determine which (if any) approximation is made to the
                                             !! equation describing the
                                             !difference in density
  
  integer :: drho_degree = -1 !< Degree of the polynomial when T and S are multiplied by the equation of state
                                  !! This applies only when the alpha and beta vary linearly 

  integer :: id_uhEff_2d = -1 !< Diagnostic IDs
  integer :: id_vhEff_2d = -1 !< Diagnostic IDs

  type(EOS_type), pointer :: EOS   !< Equation of state parameters
  type(remapping_CS) :: remap_CS   !< Remapping control structure used to create sublayers
  logical :: remap_answers_2018    !< If true, use the order of arithmetic and expressions that
                                   !! recover the answers for remapping from the end of 2018.
                                   !! Otherwise, use more robust forms of the same expressions.
  type(KPP_CS),           pointer :: KPP_CSp => NULL()          !< KPP control structure needed to get BLD
  type(energetic_PBL_CS), pointer :: energetic_PBL_CSp => NULL()!< ePBL control structure needed to get MLD

  type(column_properties), dimension (:,:), allocatable :: column !< Store properties of a column needed
                                                                  !! to calculate the position of neutral surfaces
end type neutral_diffusion_CS

! This include declares and sets the variable "version".
#include "version_variable.h"
character(len=40)  :: mdl = "MOM_neutral_diffusion" !< module name

contains

!> Read parameters and allocate control str/mark_unstableucture for neutral_diffusion module.
logical function neutral_diffusion_init(Time, G, US, param_file, diag, EOS, diabatic_CSp, CS)
  type(time_type), target,    intent(in)    :: Time       !< Time structure
  type(ocean_grid_type),      intent(in)    :: G          !< Grid structure
  type(unit_scale_type),      intent(in)    :: US         !< A dimensional unit scaling type
  type(diag_ctrl), target,    intent(inout) :: diag       !< Diagnostics control structure
  type(param_file_type),      intent(in)    :: param_file !< Parameter file structure
  type(EOS_type),  target,    intent(in)    :: EOS        !< Equation of state
  type(diabatic_CS),          pointer       :: diabatic_CSp!< KPP control structure needed to get BLD
  type(neutral_diffusion_CS), pointer       :: CS         !< Neutral diffusion control structure

  ! Local variables
  character(len=256) :: mesg    ! Message for error messages.
  character(len=80)  :: string  ! Temporary strings
  integer :: remap_degree
  logical :: default_2018_answers
  logical :: boundary_extrap
  integer :: i, j

  if (associated(CS)) then
    call MOM_error(FATAL, "neutral_diffusion_init called with associated control structure.")
    return
  endif


  ! Log this module and master switch for turning it on/off
  call log_version(param_file, mdl, version, &
       "This module implements neutral diffusion of tracers")
  call get_param(param_file, mdl, "USE_NEUTRAL_DIFFUSION", neutral_diffusion_init, &
                 "If true, enables the neutral diffusion module.", &
                 default=.false.)

  if (.not.neutral_diffusion_init) then
    return
  endif

  allocate(CS)
  CS%diag => diag
  CS%EOS => EOS
 ! call openParameterBlock(param_file,'NEUTRAL_DIFF')

  ! Read all relevant parameters and write them to the model log.
  call get_param(param_file, mdl, "NDIFF_CONTINUOUS", CS%continuous_reconstruction, &
                 "If true, uses a continuous reconstruction of T and S when "//&
                 "finding neutral surfaces along which diffusion will happen. "//&
                 "If false, a PPM discontinuous reconstruction of T and S "//&
                 "is done which results in a higher order routine but exacts "//&
                 "a higher computational cost.", default=.true.)
  call get_param(param_file, mdl, "NDIFF_REF_PRES", CS%ref_pres,                    &
                 "The reference pressure (Pa) used for the derivatives of "//&
                 "the equation of state. If negative (default), local pressure is used.", &
                 units="Pa", default = -1., scale=US%kg_m3_to_R*US%m_s_to_L_T**2)
  call get_param(param_file, mdl, "NDIFF_INTERIOR_ONLY", CS%interior_only, &
                 "If true, only applies neutral diffusion in the ocean interior."//&
                 "That is, the algorithm will exclude the surface and bottom"//&
                 "boundary layers.", default = .false.)

  ! Initialize and configure remapping
  if ( .not.CS%continuous_reconstruction ) then
    call get_param(param_file, mdl, "NDIFF_BOUNDARY_EXTRAP", boundary_extrap, &
                   "Extrapolate at the top and bottommost cells, otherwise   \n"//  &
                   "assume boundaries are piecewise constant",                      &
                   default=.false.)
    call get_param(param_file, mdl, "NDIFF_REMAPPING_SCHEME", string, &
                   "This sets the reconstruction scheme used "//&
                   "for vertical remapping for all variables. "//&
                   "It can be one of the following schemes: "//&
                   trim(remappingSchemesDoc), default=remappingDefaultScheme)
    call get_param(param_file, mdl, "DEFAULT_2018_ANSWERS", default_2018_answers, &
                 "This sets the default value for the various _2018_ANSWERS parameters.", &
                 default=.true.)
    call get_param(param_file, mdl, "REMAPPING_2018_ANSWERS", CS%remap_answers_2018, &
                 "If true, use the order of arithmetic and expressions that recover the "//&
                 "answers from the end of 2018.  Otherwise, use updated and more robust "//&
                 "forms of the same expressions.", default=default_2018_answers)
    call initialize_remapping( CS%remap_CS, string, boundary_extrapolation=boundary_extrap, &
                               answers_2018=CS%remap_answers_2018 )
    call extract_member_remapping_CS(CS%remap_CS, degree=CS%deg)
    call get_param(param_file, mdl, "NEUTRAL_POS_METHOD", CS%neutral_pos_method,   &
                   "Method used to find the neutral position                 \n"// &
                   "1. Delta_rho varies linearly, find 0 crossing            \n"// &
                   "2. Alpha and beta vary linearly from top to bottom,      \n"// &
                   "   Newton's method for neutral position                  \n"// &
                   "3. Full nonlinear equation of state, use regula falsi    \n"// &
                   "   for neutral position", default=3)
    if (CS%neutral_pos_method > 4 .or. CS%neutral_pos_method < 0) then
      call MOM_error(FATAL,"Invalid option for NEUTRAL_POS_METHOD")
    endif
    if (CS%neutral_pos_method == 2) then
      CS%drho_degree = CS%deg + 1
    endif
    call get_param(param_file, mdl, "DELTA_RHO_FORM", CS%delta_rho_form,           &
                   "Determine how the difference in density is calculated    \n"// &
                   "  full       : Difference of in-situ densities           \n"// &
                   "  no_pressure: Calculated from dRdT, dRdS, but no        \n"// &
                   "               pressure dependence",                           &
                   default="mid_pressure")
    if (CS%neutral_pos_method > 1) then
      call get_param(param_file, mdl, "NDIFF_DRHO_TOL", CS%drho_tol,            &
                     "Sets the convergence criterion for finding the neutral\n"// &
                     "position within a layer in kg m-3.",                        &
                     default=1.e-10, scale=US%kg_m3_to_R)
      call get_param(param_file, mdl, "NDIFF_X_TOL", CS%x_tol,            &
                     "Sets the convergence criterion for a change in nondim\n"// &
                     "position within a layer.",                        &
                     default=0.)
      call get_param(param_file, mdl, "NDIFF_MAX_ITER", CS%max_iter,              &
                    "The maximum number of iterations to be done before \n"//     &
                     "exiting the iterative loop to find the neutral surface",    &
                     default=10)
    endif
    call get_param(param_file, mdl, "NDIFF_DEBUG", CS%debug,             &
                   "Turns on verbose output for discontinuous neutral "//&
                   "diffusion routines.", &
                   default = .false.)
  endif

  ! Store a rescaling factor for use in diagnostic messages.
  CS%R_to_kg_m3 = US%R_to_kg_m3

  if (CS%interior_only) then
    call extract_diabatic_member(diabatic_CSp, KPP_CSp=CS%KPP_CSp)
    call extract_diabatic_member(diabatic_CSp, energetic_PBL_CSp=CS%energetic_PBL_CSp)
    if ( .not. ASSOCIATED(CS%energetic_PBL_CSp) .and. .not. ASSOCIATED(CS%KPP_CSp) ) then
      call MOM_error(FATAL,"NDIFF_INTERIOR_ONLY is true, but no valid boundary layer scheme was found")
    endif
  endif

! call get_param(param_file, mdl, "KHTR", CS%KhTr, &
!                "The background along-isopycnal tracer diffusivity.", &
!                units="m2 s-1", default=0.0)
!  call closeParameterBlock(param_file)
  if (CS%continuous_reconstruction) then
    CS%nsurf = 2*G%ke+2 ! Continuous reconstruction means that every interface has two connections
    allocate(CS%dRdT(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%dRdT(:,:,:) = 0.
    allocate(CS%dRdS(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%dRdS(:,:,:) = 0.
  else
    CS%nsurf = 4*G%ke   ! Discontinuous means that every interface has four connections
    allocate(CS%column(SZI_(G),SZJ_(G)))
    do i = G%isd,G%ied; do j=G%jsd, G%jed
      allocate(CS%column(i,j)%T_at_interface(G%ke,2)); CS%column(i,j)%T_at_interface(:,:) = 0.
      allocate(CS%column(i,j)%S_at_interface(G%ke,2)); CS%column(i,j)%S_at_interface(:,:) = 0.
      allocate(CS%column(i,j)%P_at_interface(G%ke,2)); CS%column(i,j)%P_at_interface(:,:) = 0.
      allocate(CS%column(i,j)%T_poly(G%ke,CS%deg+1)); CS%column(i,j)%T_poly(:,:) = 0.
      allocate(CS%column(i,j)%S_poly(G%ke,CS%deg+1)); CS%column(i,j)%S_poly(:,:) = 0.
      allocate(CS%column(i,j)%stable_cell(G%ke)); CS%column(i,j)%stable_cell(:) = .true.
      CS%column(i,j)%dRdT(:) = 0.
      CS%column(i,j)%dRdS(:) = 0.
    enddo; enddo

    allocate(CS%T_i(SZI_(G),SZJ_(G),SZK_(G),2))    ; CS%T_i(:,:,:,:) = 0.
    allocate(CS%S_i(SZI_(G),SZJ_(G),SZK_(G),2))    ; CS%S_i(:,:,:,:) = 0.
    allocate(CS%P_i(SZI_(G),SZJ_(G),SZK_(G),2))    ; CS%P_i(:,:,:,:) = 0.
    allocate(CS%dRdT_i(SZI_(G),SZJ_(G),SZK_(G),2)) ; CS%dRdT_i(:,:,:,:) = 0.
    allocate(CS%dRdS_i(SZI_(G),SZJ_(G),SZK_(G),2)) ; CS%dRdS_i(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_T(SZI_(G),SZJ_(G),SZK_(G),CS%deg+1)) ; CS%ppoly_coeffs_T(:,:,:,:) = 0.
    allocate(CS%ppoly_coeffs_S(SZI_(G),SZJ_(G),SZK_(G),CS%deg+1)) ; CS%ppoly_coeffs_S(:,:,:,:) = 0.
    allocate(CS%ns(SZI_(G),SZJ_(G)))    ; CS%ns(:,:) = 0.
  endif
  ! T-points
  allocate(CS%Tint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Tint(:,:,:) = 0.
  allocate(CS%Sint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Sint(:,:,:) = 0.
  allocate(CS%Pint(SZI_(G),SZJ_(G),SZK_(G)+1)) ; CS%Pint(:,:,:) = 0.
  allocate(CS%stable_cell(SZI_(G),SZJ_(G),SZK_(G))) ; CS%stable_cell(:,:,:) = .true.
  ! U-points
  allocate(CS%uPoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uPoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uPoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uPoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0.
  allocate(CS%uKoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uKoL(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uKoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%uKoR(G%isc-1:G%iec,G%jsc:G%jec,:)   = 0
  allocate(CS%uHeff(G%isd:G%ied,G%jsd:G%jed,CS%nsurf-1)); CS%uHeff(G%isc-1:G%iec,G%jsc:G%jec,:) = 0
  ! V-points
  allocate(CS%vPoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vPoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vPoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vPoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0.
  allocate(CS%vKoL(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vKoL(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vKoR(G%isd:G%ied,G%jsd:G%jed, CS%nsurf)); CS%vKoR(G%isc:G%iec,G%jsc-1:G%jec,:)   = 0
  allocate(CS%vHeff(G%isd:G%ied,G%jsd:G%jed,CS%nsurf-1)); CS%vHeff(G%isc:G%iec,G%jsc-1:G%jec,:) = 0

end function neutral_diffusion_init

!> Calculate remapping factors for u/v columns used to map adjoining columns to
!! a shared coordinate space.
subroutine neutral_diffusion_calc_coeffs(G, GV, US, h, T, S, CS, p_surf)
  type(ocean_grid_type),                    intent(in) :: G   !< Ocean grid structure
  type(verticalGrid_type),                  intent(in) :: GV  !< ocean vertical grid structure
  type(unit_scale_type),                    intent(in) :: US  !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: h   !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: T   !< Potential temperature [degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in) :: S   !< Salinity [ppt]
  type(neutral_diffusion_CS),               pointer    :: CS  !< Neutral diffusion control structure
  real, dimension(SZI_(G),SZJ_(G)), optional, intent(in) :: p_surf !< Surface pressure to include in pressures used
                                                              !! for equation of state calculations [R L2 T-2 ~> Pa]

  ! Local variables
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k
  ! Variables used for reconstructions
  real, dimension(SZK_(G),2) :: ppoly_r_S       ! Reconstruction slopes
  real, dimension(SZI_(G), SZJ_(G)) :: hEff_sum ! Summed effective face thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G))  :: hbl      ! Boundary layer depth [H ~> m or kg m-2]
  integer :: iMethod
  real, dimension(SZI_(G)) :: ref_pres ! Reference pressure used to calculate alpha/beta [R L2 T-2 ~> Pa]
  real, dimension(SZI_(G)) :: rho_tmp  ! Routine to calculate drho_dp, returns density which is not used
  real :: h_neglect, h_neglect_edge    ! Negligible thicknesses [H ~> m or kg m-2]
  integer, dimension(SZI_(G), SZJ_(G)) :: k_top  ! Index of the first layer within the boundary
  real,    dimension(SZI_(G), SZJ_(G)) :: zeta_top ! Distance from the top of a layer to the intersection of the
                                                   ! top extent of the boundary layer (0 at top, 1 at bottom) [nondim]
  integer, dimension(SZI_(G), SZJ_(G)) :: k_bot    ! Index of the last layer within the boundary
  real,    dimension(SZI_(G), SZJ_(G)) :: zeta_bot ! Distance of the lower layer to the boundary layer depth
  real :: pa_to_H                      ! A conversion factor from pressure to H units [H T2 R-1 Z-2 ~> m Pa-1 or s2 m-2]

  pa_to_H = 1. / (GV%H_to_RZ * GV%g_Earth)

  k_top(:,:) = 1     ; k_bot(:,:) = 1
  zeta_top(:,:) = 0. ; zeta_bot(:,:) = 1.

  ! Check if hbl needs to be extracted
  if (CS%interior_only) then
    hbl(:,:) = 0.
    if (ASSOCIATED(CS%KPP_CSp)) call KPP_get_BLD(CS%KPP_CSp, hbl, G, US, m_to_BLD_units=GV%m_to_H)
    if (ASSOCIATED(CS%energetic_PBL_CSp)) &
      call energetic_PBL_get_MLD(CS%energetic_PBL_CSp, hbl, G, US, m_to_MLD_units=GV%m_to_H)
    call pass_var(hbl, G%Domain)
    ! get k-indices and zeta
    do j=G%jsc-1, G%jec+1 ; do i=G%isc-1,G%iec+1
      call boundary_k_range(SURFACE, G%ke, h(i,j,:), hbl(i,j), k_top(i,j), zeta_top(i,j), k_bot(i,j), zeta_bot(i,j))
    enddo; enddo
    ! TODO: add similar code for BOTTOM boundary layer
  endif

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  elseif (GV%Boussinesq) then
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  else
    h_neglect = GV%kg_m2_to_H*1.0e-30 ; h_neglect_edge = GV%kg_m2_to_H*1.0e-10
  endif

  ! If doing along isopycnal diffusion (as opposed to neutral diffusion, set the reference pressure)
  if (CS%ref_pres>=0.) then
    ref_pres(:) = CS%ref_pres
  endif

  if (CS%continuous_reconstruction) then
    CS%dRdT(:,:,:) = 0.
    CS%dRdS(:,:,:) = 0.
  else
    CS%T_i(:,:,:,:) = 0.
    CS%S_i(:,:,:,:) = 0.
    CS%dRdT_i(:,:,:,:) = 0.
    CS%dRdS_i(:,:,:,:) = 0.
    CS%ns(:,:) = 0.
    CS%stable_cell(:,:,:) = .true.
  endif

  ! Calculate pressure at interfaces and layer averaged alpha/beta
  if (present(p_surf)) then
     do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
       CS%Pint(i,j,1) = p_surf(i,j)
     enddo ; enddo
  else
    CS%Pint(:,:,1) = 0.
  endif
  do k=1,G%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
    CS%Pint(i,j,k+1) = CS%Pint(i,j,k) + h(i,j,k)*(GV%g_Earth*GV%H_to_RZ)
  enddo ; enddo ; enddo

  ! Pressures at the interfaces, this is redundant as P_i(k,1) = P_i(k-1,2) however retain this
  ! for now to ensure consitency of indexing for diiscontinuous reconstructions
  if (.not. CS%continuous_reconstruction) then
    if (present(p_surf)) then
      do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
        CS%column(i,j)%P_at_interface(1,1) = p_surf(i,j)
        CS%column(i,j)%P_at_interface(1,2) = p_surf(i,j) + h(i,j,1)*(GV%H_to_RZ*GV%g_Earth)
      enddo ; enddo
    else
      do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
        CS%column(i,j)%P_at_interface(1,1) = 0.
        CS%column(i,j)%P_at_interface(1,2) = h(i,j,1)*(GV%H_to_RZ*GV%g_Earth)
      enddo ; enddo
    endif
    do k=2,G%ke ; do j=G%jsc-1,G%jec+1 ; do i=G%isc-1,G%iec+1
      CS%column(i,j)%P_at_interface(k,1) = CS%column(i,j)%P_at_interface(k-1,2)
      CS%column(i,j)%P_at_interface(k,2) = CS%column(i,j)%P_at_interface(k-1,2) + h(i,j,k)*(GV%H_to_RZ*GV%g_Earth)
    enddo ; enddo ; enddo
  endif

  EOSdom(:) = EOS_domain(G%HI, halo=1)

  do j = G%jsc-1, G%jec+1
    ! Interpolate state to interface
    do i = G%isc-1, G%iec+1
      if (CS%continuous_reconstruction) then
        call interface_scalar(G%ke, h(i,j,:), T(i,j,:), CS%Tint(i,j,:), 2, h_neglect)
        call interface_scalar(G%ke, h(i,j,:), S(i,j,:), CS%Sint(i,j,:), 2, h_neglect)
      else
        call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), T(i,j,:), CS%column(i,j)%T_poly(:,:), &
                                       CS%column(i,j)%T_at_interface(:,:), ppoly_r_S,                     &
                                       iMethod, h_neglect, h_neglect_edge )
        call build_reconstructions_1d( CS%remap_CS, G%ke, h(i,j,:), S(i,j,:), CS%column(i,j)%S_poly(:,:), &
                                       CS%column(i,j)%S_at_interface(:,:), ppoly_r_S,                      &
                                       iMethod, h_neglect, h_neglect_edge )
        ! In the current ALE formulation, interface values are not exactly at the 0. or 1. of the
        ! polynomial reconstructions
        do k=1,G%ke
          CS%column(i,j)%T_at_interface(k,1) = evaluation_polynomial( CS%column(i,j)%T_poly(k,:), CS%deg+1, 0. )
          CS%column(i,j)%T_at_interface(k,2) = evaluation_polynomial( CS%column(i,j)%T_poly(k,:), CS%deg+1, 1. )
          CS%column(i,j)%S_at_interface(k,1) = evaluation_polynomial( CS%column(i,j)%S_poly(k,:), CS%deg+1, 0. )
          CS%column(i,j)%S_at_interface(k,2) = evaluation_polynomial( CS%column(i,j)%S_poly(k,:), CS%deg+1, 1. )
        enddo
      endif
    enddo

    ! Continuous reconstruction
    if (CS%continuous_reconstruction) then
      do k = 1, G%ke+1
        if (CS%ref_pres<0) ref_pres(:) = CS%Pint(:,j,k)
        call calculate_density_derivs(CS%Tint(:,j,k), CS%Sint(:,j,k), ref_pres, CS%dRdT(:,j,k), &
                                      CS%dRdS(:,j,k), CS%EOS, EOSdom)
      enddo
    endif
  enddo

  if (.not. CS%continuous_reconstruction) then
    do j = G%jsc-1, G%jec+1 ; do i = G%isc-1, G%iec+1
      call column_mark_unstable( CS, G%ke, CS%column(i,j) ) 
    enddo ; enddo
  endif

  CS%uhEff(:,:,:) = 0.
  CS%vhEff(:,:,:) = 0.
  CS%uPoL(:,:,:) = 0.
  CS%vPoL(:,:,:) = 0.
  CS%uPoR(:,:,:) = 0.
  CS%vPoR(:,:,:) = 0.
  CS%uKoL(:,:,:) = 1
  CS%vKoL(:,:,:) = 1
  CS%uKoR(:,:,:) = 1
  CS%vKoR(:,:,:) = 1

  ! Neutral surface factors at U points
  do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
    if (G%mask2dCu(I,j) > 0.) then
      if (CS%continuous_reconstruction) then
        call find_neutral_surface_positions_continuous(G%ke,                                               &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),            &
                CS%Pint(i+1,j,:), CS%Tint(i+1,j,:), CS%Sint(i+1,j,:), CS%dRdT(i+1,j,:), CS%dRdS(i+1,j,:),  &
                CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:),           &
                k_bot(I,j), k_bot(I+1,j), 1.-zeta_bot(I,j), 1.-zeta_bot(I+1,j))
      else
        call find_neutral_surface_positions_discontinuous(CS, G%ke, CS%column(i,j), CS%column(i+1,j),                &
                h(i,j,:), h(i+1,j,:), CS%uPoL(I,j,:), CS%uPoR(I,j,:), CS%uKoL(I,j,:), CS%uKoR(I,j,:), CS%uhEff(I,j,:))
      endif
    endif
  enddo ; enddo

  ! Neutral surface factors at V points
  do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
    if (G%mask2dCv(i,J) > 0.) then
      if (CS%continuous_reconstruction) then
        call find_neutral_surface_positions_continuous(G%ke,                                              &
                CS%Pint(i,j,:), CS%Tint(i,j,:), CS%Sint(i,j,:), CS%dRdT(i,j,:), CS%dRdS(i,j,:),           &
                CS%Pint(i,j+1,:), CS%Tint(i,j+1,:), CS%Sint(i,j+1,:), CS%dRdT(i,j+1,:), CS%dRdS(i,j+1,:), &
                CS%vPoL(i,J,:), CS%vPoR(i,J,:), CS%vKoL(i,J,:), CS%vKoR(i,J,:), CS%vhEff(i,J,:), &
                k_bot(i,J), k_bot(i,J+1), 1.-zeta_bot(i,J), 1.-zeta_bot(i,J+1))
      else
        call find_neutral_surface_positions_discontinuous(CS, G%ke, CS%column(i,j), CS%column(i,j+1), &
            h(i,j,:), h(i,j+1,:), CS%vPoL(I,j,:), CS%vPoR(I,j,:), CS%vKoL(I,j,:), CS%vKoR(I,j,:), CS%vhEff(I,j,:))
      endif
    endif
  enddo ; enddo

  ! Continuous reconstructions calculate hEff as the difference between the pressures of the
  ! neutral surfaces which need to be reconverted to thickness units. The discontinuous version
  ! calculates hEff from the nondimensional fraction of the layer spanned by adjacent neutral
  ! surfaces, so hEff is already in thickness units.
  if (CS%continuous_reconstruction) then
    do k = 1, CS%nsurf-1 ; do j = G%jsc, G%jec ; do I = G%isc-1, G%iec
      if (G%mask2dCu(I,j) > 0.) CS%uhEff(I,j,k) = CS%uhEff(I,j,k) * pa_to_H
    enddo ; enddo ; enddo
    do k = 1, CS%nsurf-1 ; do J = G%jsc-1, G%jec ; do i = G%isc, G%iec
      if (G%mask2dCv(i,J) > 0.) CS%vhEff(i,J,k) = CS%vhEff(i,J,k) * pa_to_H
    enddo ; enddo ; enddo
  endif

  if (CS%id_uhEff_2d>0) then
    hEff_sum(:,:) = 0.
    do k = 1,CS%nsurf-1 ; do j=G%jsc,G%jec ; do i=G%isc-1,G%iec
      hEff_sum(i,j) = hEff_sum(i,j) + CS%uhEff(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_uhEff_2d, hEff_sum, CS%diag)
  endif
  if (CS%id_vhEff_2d>0) then
    hEff_sum(:,:) = 0.
    do k = 1,CS%nsurf-1 ; do j=G%jsc-1,G%jec ; do i=G%isc,G%iec
      hEff_sum(i,j) = hEff_sum(i,j) + CS%vhEff(i,j,k)
    enddo ; enddo ; enddo
    call post_data(CS%id_vhEff_2d, hEff_sum, CS%diag)
  endif

end subroutine neutral_diffusion_calc_coeffs

!> Update tracer concentration due to neutral diffusion; layer thickness unchanged by this update.
subroutine neutral_diffusion(G, GV, h, Coef_x, Coef_y, dt, Reg, US, CS)
  type(ocean_grid_type),                     intent(in)    :: G      !< Ocean grid structure
  type(verticalGrid_type),                   intent(in)    :: GV     !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)),  intent(in)    :: h      !< Layer thickness [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: Coef_x !< dt * Kh * dy / dx at u-points [L2 ~> m2]
  real, dimension(SZI_(G),SZJB_(G)),         intent(in)    :: Coef_y !< dt * Kh * dx / dy at v-points [L2 ~> m2]
  real,                                      intent(in)    :: dt     !< Tracer time step * I_numitts [T ~> s]
                                                                     !! (I_numitts in tracer_hordiff)
  type(tracer_registry_type),                pointer       :: Reg    !< Tracer registry
  type(unit_scale_type),                     intent(in)    :: US     !< A dimensional unit scaling type
  type(neutral_diffusion_CS),                pointer       :: CS     !< Neutral diffusion control structure

  ! Local variables
  real, dimension(SZIB_(G),SZJ_(G),CS%nsurf-1) :: uFlx        ! Zonal flux of tracer [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJB_(G),CS%nsurf-1) :: vFlx        ! Meridional flux of tracer
                                                              ! [H conc ~> m conc or conc kg m-2]
  real, dimension(SZI_(G),SZJ_(G),G%ke)        :: tendency    ! tendency array for diagn
  real, dimension(SZI_(G),SZJ_(G))             :: tendency_2d ! depth integrated content tendency for diagn
  real, dimension(SZIB_(G),SZJ_(G))            :: trans_x_2d  ! depth integrated diffusive tracer x-transport diagn
  real, dimension(SZI_(G),SZJB_(G))            :: trans_y_2d  ! depth integrated diffusive tracer y-transport diagn
  real, dimension(G%ke)                        :: dTracer     ! change in tracer concentration due to ndiffusion

  type(tracer_type), pointer                   :: Tracer => NULL() ! Pointer to the current tracer

  integer :: i, j, k, m, ks, nk
  real :: Idt  ! The inverse of the time step [T-1 ~> s-1]
  real :: h_neglect, h_neglect_edge

  if (.not.CS%remap_answers_2018) then
    h_neglect = GV%H_subroundoff ; h_neglect_edge = GV%H_subroundoff
  else
    h_neglect = GV%m_to_H*1.0e-30 ; h_neglect_edge = GV%m_to_H*1.0e-10
  endif

  nk = GV%ke

  do m = 1,Reg%ntr ! Loop over tracer registry

    tracer => Reg%Tr(m)

    ! for diagnostics
    if (tracer%id_dfxy_conc > 0 .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 .or. &
        tracer%id_dfx_2d > 0 .or. tracer%id_dfy_2d > 0) then
      Idt = 1.0 / dt
      tendency(:,:,:)  = 0.0
    endif

    uFlx(:,:,:) = 0.
    vFlx(:,:,:) = 0.

    ! x-flux
    do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
      if (G%mask2dCu(I,j)>0.) then
        call neutral_surface_flux(nk, CS%nsurf, CS%deg, h(i,j,:), h(i+1,j,:),       &
                                  tracer%t(i,j,:), tracer%t(i+1,j,:), &
                                  CS%uPoL(I,j,:), CS%uPoR(I,j,:), &
                                  CS%uKoL(I,j,:), CS%uKoR(I,j,:), &
                                  CS%uhEff(I,j,:), uFlx(I,j,:), &
                                  CS%continuous_reconstruction, h_neglect, CS%remap_CS, h_neglect_edge)
      endif
    enddo ; enddo

    ! y-flux
    do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
      if (G%mask2dCv(i,J)>0.) then
        call neutral_surface_flux(nk, CS%nsurf, CS%deg, h(i,j,:), h(i,j+1,:),       &
                                  tracer%t(i,j,:), tracer%t(i,j+1,:), &
                                  CS%vPoL(i,J,:), CS%vPoR(i,J,:), &
                                  CS%vKoL(i,J,:), CS%vKoR(i,J,:), &
                                  CS%vhEff(i,J,:), vFlx(i,J,:),   &
                                  CS%continuous_reconstruction, h_neglect, CS%remap_CS, h_neglect_edge)
      endif
    enddo ; enddo

    ! Update the tracer concentration from divergence of neutral diffusive flux components
    do j = G%jsc,G%jec ; do i = G%isc,G%iec
      if (G%mask2dT(i,j)>0.) then

        dTracer(:) = 0.
        do ks = 1,CS%nsurf-1
          k = CS%uKoL(I,j,ks)
          if (k>0) dTracer(k) = dTracer(k) + Coef_x(I,j)   * uFlx(I,j,ks)
          k = CS%uKoR(I-1,j,ks)
          if (k>0) dTracer(k) = dTracer(k) - Coef_x(I-1,j) * uFlx(I-1,j,ks)
          k = CS%vKoL(i,J,ks)
          if (k>0) dTracer(k) = dTracer(k) + Coef_y(i,J)   * vFlx(i,J,ks)
          k = CS%vKoR(i,J-1,ks)
          if (k>0) dTracer(k) = dTracer(k) - Coef_y(i,J-1) * vFlx(i,J-1,ks)
        enddo
        do k = 1, GV%ke
          tracer%t(i,j,k) = tracer%t(i,j,k) + dTracer(k) * &
                          ( G%IareaT(i,j) / ( h(i,j,k) + GV%H_subroundoff ) )
        enddo

        if (tracer%id_dfxy_conc > 0  .or. tracer%id_dfxy_cont > 0 .or. tracer%id_dfxy_cont_2d > 0 ) then
          do k = 1, GV%ke
            tendency(i,j,k) = dTracer(k) * G%IareaT(i,j) * Idt
          enddo
        endif

      endif
    enddo ; enddo

    ! Diagnose vertically summed zonal flux, giving zonal tracer transport from ndiff.
    ! Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfx_2d > 0) then
      do j = G%jsc,G%jec ; do I = G%isc-1,G%iec
        trans_x_2d(I,j) = 0.
        if (G%mask2dCu(I,j)>0.) then
          do ks = 1,CS%nsurf-1
            trans_x_2d(I,j) = trans_x_2d(I,j) - Coef_x(I,j) * uFlx(I,j,ks)
          enddo
          trans_x_2d(I,j) = trans_x_2d(I,j) * Idt
        endif
      enddo ; enddo
      call post_data(tracer%id_dfx_2d, trans_x_2d(:,:), CS%diag)
    endif

    ! Diagnose vertically summed merid flux, giving meridional tracer transport from ndiff.
    ! Note sign corresponds to downgradient flux convention.
    if (tracer%id_dfy_2d > 0) then
      do J = G%jsc-1,G%jec ; do i = G%isc,G%iec
        trans_y_2d(i,J) = 0.
        if (G%mask2dCv(i,J)>0.) then
          do ks = 1,CS%nsurf-1
            trans_y_2d(i,J) = trans_y_2d(i,J) - Coef_y(i,J) * vFlx(i,J,ks)
          enddo
          trans_y_2d(i,J) = trans_y_2d(i,J) * Idt
        endif
      enddo ; enddo
      call post_data(tracer%id_dfy_2d, trans_y_2d(:,:), CS%diag)
    endif

    ! post tendency of tracer content
    if (tracer%id_dfxy_cont > 0) then
      call post_data(tracer%id_dfxy_cont, tendency(:,:,:), CS%diag)
    endif

    ! post depth summed tendency for tracer content
    if (tracer%id_dfxy_cont_2d > 0) then
      tendency_2d(:,:) = 0.
      do j = G%jsc,G%jec ; do i = G%isc,G%iec
        do k = 1, GV%ke
          tendency_2d(i,j) = tendency_2d(i,j) + tendency(i,j,k)
        enddo
      enddo ; enddo
      call post_data(tracer%id_dfxy_cont_2d, tendency_2d(:,:), CS%diag)
    endif

    ! post tendency of tracer concentration; this step must be
    ! done after posting tracer content tendency, since we alter
    ! the tendency array.
    if (tracer%id_dfxy_conc > 0) then
      do k = 1, GV%ke ; do j = G%jsc,G%jec ; do i = G%isc,G%iec
        tendency(i,j,k) =  tendency(i,j,k) / ( h(i,j,k) + GV%H_subroundoff )
      enddo ; enddo ; enddo
      call post_data(tracer%id_dfxy_conc, tendency, CS%diag)
    endif
  enddo ! Loop over tracer registry

end subroutine neutral_diffusion

!> Returns interface scalar, Si, for a column of layer values, S.
subroutine interface_scalar(nk, h, S, Si, i_method, h_neglect)
  integer,               intent(in)    :: nk       !< Number of levels
  real, dimension(nk),   intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2]
  real, dimension(nk),   intent(in)    :: S        !< Layer scalar (conc, e.g. ppt)
  real, dimension(nk+1), intent(inout) :: Si       !< Interface scalar (conc, e.g. ppt)
  integer,               intent(in)    :: i_method !< =1 use average of PLM edges
                                                   !! =2 use continuous PPM edge interpolation
  real,                  intent(in)    :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]
  ! Local variables
  integer :: k, km2, kp1
  real, dimension(nk) :: diff
  real :: Sb, Sa

  call PLM_diff(nk, h, S, 2, 1, diff)
  Si(1) = S(1) - 0.5 * diff(1)
  if (i_method==1) then
    do k = 2, nk
      ! Average of the two edge values (will be bounded and,
      ! when slopes are unlimited, notionally second-order accurate)
      Sa = S(k-1) + 0.5 * diff(k-1) ! Lower edge value of a PLM reconstruction for layer above
      Sb = S(k) - 0.5 * diff(k) ! Upper edge value of a PLM reconstruction for layer below
      Si(k) = 0.5 * ( Sa + Sb )
    enddo
  elseif (i_method==2) then
    do k = 2, nk
      ! PPM quasi-fourth order interpolation for edge values following
      ! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
      km2 = max(1, k-2)
      kp1 = min(nk, k+1)
      Si(k) = ppm_edge(h(km2), h(k-1), h(k), h(kp1),  S(k-1), S(k), diff(k-1), diff(k), h_neglect)
    enddo
  endif
  Si(nk+1) = S(nk) + 0.5 * diff(nk)

end subroutine interface_scalar

!> Returns the PPM quasi-fourth order edge value at k+1/2 following
!! equation 1.6 in Colella & Woodward, 1984: JCP 54, 174-201.
real function ppm_edge(hkm1, hk, hkp1, hkp2,  Ak, Akp1, Pk, Pkp1, h_neglect)
  real, intent(in) :: hkm1 !< Width of cell k-1
  real, intent(in) :: hk   !< Width of cell k
  real, intent(in) :: hkp1 !< Width of cell k+1
  real, intent(in) :: hkp2 !< Width of cell k+2
  real, intent(in) :: Ak   !< Average scalar value of cell k
  real, intent(in) :: Akp1 !< Average scalar value of cell k+1
  real, intent(in) :: Pk   !< PLM slope for cell k
  real, intent(in) :: Pkp1 !< PLM slope for cell k+1
  real, intent(in) :: h_neglect !< A negligibly small thickness [H ~> m or kg m-2]

  ! Local variables
  real :: R_hk_hkp1, R_2hk_hkp1, R_hk_2hkp1, f1, f2, f3, f4

  R_hk_hkp1 = hk + hkp1
  if (R_hk_hkp1 <= 0.) then
    ppm_edge = 0.5 * ( Ak + Akp1 )
    return
  endif
  R_hk_hkp1 = 1. / R_hk_hkp1
  if (hk<hkp1) then
    ppm_edge = Ak + ( hk * R_hk_hkp1 ) * ( Akp1 - Ak )
  else
    ppm_edge = Akp1 + ( hkp1 * R_hk_hkp1 ) * ( Ak - Akp1 )
  endif

  R_2hk_hkp1 = 1. / ( ( 2. * hk + hkp1 ) + h_neglect )
  R_hk_2hkp1 = 1. / ( ( hk + 2. * hkp1 ) + h_neglect )
  f1 = 1./ ( ( hk + hkp1) + ( hkm1 + hkp2 ) )
  f2 = 2. * ( hkp1 * hk ) * R_hk_hkp1 * &
            ( ( hkm1 + hk ) * R_2hk_hkp1  - ( hkp2 + hkp1 ) * R_hk_2hkp1 )
  f3 = hk * ( hkm1 + hk ) * R_2hk_hkp1
  f4 = hkp1 * ( hkp1 + hkp2 ) * R_hk_2hkp1

  ppm_edge = ppm_edge + f1 * ( f2 * ( Akp1 - Ak ) - ( f3 * Pkp1 - f4 * Pk ) )

end function ppm_edge

!> Returns the average of a PPM reconstruction between two
!! fractional positions.
real function ppm_ave(xL, xR, aL, aR, aMean)
  real, intent(in) :: xL    !< Fraction position of left bound (0,1)
  real, intent(in) :: xR    !< Fraction position of right bound (0,1)
  real, intent(in) :: aL    !< Left edge scalar value, at x=0
  real, intent(in) :: aR    !< Right edge scalar value, at x=1
  real, intent(in) :: aMean !< Average scalar value of cell

  ! Local variables
  real :: dx, xave, a6, a6o3

  dx = xR - xL
  xave = 0.5 * ( xR + xL )
  a6o3 = 2. * aMean - ( aL + aR ) ! a6 / 3.
  a6 = 3. * a6o3

  if (dx<0.) then
    stop 'ppm_ave: dx<0 should not happend!'
  elseif (dx>1.) then
    stop 'ppm_ave: dx>1 should not happend!'
  elseif (dx==0.) then
    ppm_ave = aL + ( aR - aL ) * xR + a6 * xR * ( 1. - xR )
  else
    ppm_ave = ( aL + xave * ( ( aR - aL ) + a6 ) )  - a6o3 * ( xR**2 + xR * xL + xL**2 )
  endif
end function ppm_ave

!> A true signum function that returns either -abs(a), when x<0; or abs(a) when x>0; or 0 when x=0.
real function signum(a,x)
  real, intent(in) :: a !< The magnitude argument
  real, intent(in) :: x !< The sign (or zero) argument

  signum = sign(a,x)
  if (x==0.) signum = 0.

end function signum

!> Returns PLM slopes for a column where the slopes are the difference in value across each cell.
!! The limiting follows equation 1.8 in Colella & Woodward, 1984: JCP 54, 174-201.
subroutine PLM_diff(nk, h, S, c_method, b_method, diff)
  integer,             intent(in)    :: nk       !< Number of levels
  real, dimension(nk), intent(in)    :: h        !< Layer thickness [H ~> m or kg m-2]
  real, dimension(nk), intent(in)    :: S        !< Layer salinity (conc, e.g. ppt)
  integer,             intent(in)    :: c_method !< Method to use for the centered difference
  integer,             intent(in)    :: b_method !< =1, use PCM in first/last cell, =2 uses linear extrapolation
  real, dimension(nk), intent(inout) :: diff     !< Scalar difference across layer (conc, e.g. ppt)
                                                 !! determined by the following values for c_method:
                                                 !!   1. Second order finite difference (not recommended)
                                                 !!   2. Second order finite volume (used in original PPM)
                                                 !!   3. Finite-volume weighted least squares linear fit
                                                 !! \todo  The use of c_method to choose a scheme is inefficient
                                                 !! and should eventually be moved up the call tree.

  ! Local variables
  integer :: k
  real :: hkm1, hk, hkp1, Skm1, Sk, Skp1, diff_l, diff_r, diff_c

  do k = 2, nk-1
    hkm1 = h(k-1)
    hk = h(k)
    hkp1 = h(k+1)

    if ( ( hkp1 + hk ) * ( hkm1 + hk ) > 0.) then
      Skm1 = S(k-1)
      Sk = S(k)
      Skp1 = S(k+1)
      if (c_method==1) then
        ! Simple centered diff (from White)
        if ( hk + 0.5 * (hkm1 + hkp1) /= 0. ) then
          diff_c = ( Skp1 - Skm1 ) * ( hk / ( hk + 0.5 * (hkm1 + hkp1) ) )
        else
          diff_c = 0.
        endif
      elseif (c_method==2) then
        ! Second order accurate centered FV slope (from Colella and Woodward, JCP 1984)
        diff_c = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      elseif (c_method==3) then
        ! Second order accurate finite-volume least squares slope
        diff_c = hk * fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
      endif
      ! Limit centered slope by twice the side differenced slopes
      diff_l = 2. * ( Sk - Skm1 )
      diff_r = 2. * ( Skp1 - Sk )
      if ( signum(1., diff_l) * signum(1., diff_r) <= 0. ) then
        diff(k) = 0. ! PCM for local extrema
      else
        diff(k) = sign( min( abs(diff_l), abs(diff_c), abs(diff_r) ), diff_c )
      endif
    else
      diff(k) = 0. ! PCM next to vanished layers
    endif
  enddo
  if (b_method==1) then ! PCM for top and bottom layer
    diff(1) = 0.
    diff(nk) = 0.
  elseif (b_method==2) then ! Linear extrapolation for top and bottom interfaces
    diff(1) = ( S(2) - S(1) ) * 2. * ( h(1) / ( h(1) + h(2) ) )
    diff(nk) = S(nk) - S(nk-1) * 2. * ( h(nk) / ( h(nk-1) + h(nk) ) )
  endif

end subroutine PLM_diff

!> Returns the cell-centered second-order finite volume (unlimited PLM) slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a difference across the central cell (i.e. units of scalar S).
!! Discretization follows equation 1.7 in Colella & Woodward, 1984: JCP 54, 174-201.
real function fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: h_sum, hp, hm

  h_sum = ( hkm1 + hkp1 ) + hk
  if (h_sum /= 0.) h_sum = 1./ h_sum
  hm =  hkm1 + hk
  if (hm /= 0.) hm = 1./ hm
  hp =  hkp1 + hk
  if (hp /= 0.) hp = 1./ hp
  fv_diff = ( hk * h_sum ) * &
            (   ( 2. * hkm1 + hk ) * hp * ( Skp1 - Sk ) &
              + ( 2. * hkp1 + hk ) * hm * ( Sk - Skm1 ) )
end function fv_diff


!> Returns the cell-centered second-order weighted least squares slope
!! using three consecutive cell widths and average values. Slope is returned
!! as a gradient (i.e. units of scalar S over width units).
real function fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  real, intent(in) :: hkm1 !< Left cell width
  real, intent(in) :: hk   !< Center cell width
  real, intent(in) :: hkp1 !< Right cell width
  real, intent(in) :: Skm1 !< Left cell average value
  real, intent(in) :: Sk   !< Center cell average value
  real, intent(in) :: Skp1 !< Right cell average value

  ! Local variables
  real :: xkm1, xkp1
  real :: h_sum, hx_sum, hxsq_sum, hxy_sum, hy_sum, det

  xkm1 = -0.5 * ( hk + hkm1 )
  xkp1 = 0.5 * ( hk + hkp1 )
  h_sum = ( hkm1 + hkp1 ) + hk
  hx_sum = hkm1*xkm1 + hkp1*xkp1
  hxsq_sum = hkm1*(xkm1**2) + hkp1*(xkp1**2)
  hxy_sum = hkm1*xkm1*Skm1 + hkp1*xkp1*Skp1
  hy_sum = ( hkm1*Skm1 + hkp1*Skp1 ) + hk*Sk
  det = h_sum * hxsq_sum - hx_sum**2
  if (det /= 0.) then
    !a = ( hxsq_sum * hy_sum - hx_sum*hxy_sum ) / det ! a would be mean of straight line fit
    fvlsq_slope = ( h_sum * hxy_sum - hx_sum*hy_sum ) / det ! Gradient of straight line fit
  else
    fvlsq_slope = 0. ! Adcroft's reciprocal rule
  endif
end function fvlsq_slope


!> Returns positions within left/right columns of combined interfaces using continuous reconstructions of T/S
subroutine find_neutral_surface_positions_continuous(nk, Pl, Tl, Sl, dRdTl, dRdSl, Pr, Tr, Sr, &
                                                     dRdTr, dRdSr, PoL, PoR, KoL, KoR, hEff, bl_kl, bl_kr, bl_zl, bl_zr)
  integer,                    intent(in)    :: nk    !< Number of levels
  real, dimension(nk+1),      intent(in)    :: Pl    !< Left-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tl    !< Left-column interface potential temperature [degC]
  real, dimension(nk+1),      intent(in)    :: Sl    !< Left-column interface salinity [ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTl !< Left-column dRho/dT [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSl !< Left-column dRho/dS [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(nk+1),      intent(in)    :: Pr    !< Right-column interface pressure [R L2 T-2 ~> Pa] or other units
  real, dimension(nk+1),      intent(in)    :: Tr    !< Right-column interface potential temperature [degC]
  real, dimension(nk+1),      intent(in)    :: Sr    !< Right-column interface salinity [ppt]
  real, dimension(nk+1),      intent(in)    :: dRdTr !< Left-column dRho/dT [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk+1),      intent(in)    :: dRdSr !< Left-column dRho/dS [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(2*nk+2),    intent(inout) :: PoL   !< Fractional position of neutral surface within
                                                     !! layer KoL of left column
  real, dimension(2*nk+2),    intent(inout) :: PoR   !< Fractional position of neutral surface within
                                                     !! layer KoR of right column
  integer, dimension(2*nk+2), intent(inout) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(2*nk+2), intent(inout) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(2*nk+1),    intent(inout) :: hEff  !< Effective thickness between two neutral surfaces
                                                     !! [R L2 T-2 ~> Pa] or other units following Pl and Pr.
  integer, optional,          intent(in)    :: bl_kl !< Layer index of the boundary layer (left)
  integer, optional,          intent(in)    :: bl_kr !< Layer index of the boundary layer (right)
  real, optional,             intent(in)    :: bl_zl !< Nondimensional position of the boundary layer (left)
  real, optional,             intent(in)    :: bl_zr !< Nondimensional position of the boundary layer (right)

  ! Local variables
  integer :: ns                     ! Number of neutral surfaces
  integer :: k_surface              ! Index of neutral surface
  integer :: kl                     ! Index of left interface
  integer :: kr                     ! Index of right interface
  real    :: dRdT, dRdS             ! dRho/dT [kg m-3 degC-1] and dRho/dS [kg m-3 ppt-1] for the neutral surface
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  logical :: reached_bottom         ! True if one of the bottom-most interfaces has been used as the target
  integer :: krm1, klm1
  real    :: dRho, dRhoTop, dRhoBot ! Potential density differences at various points [R ~> kg m-3]
  real    :: hL, hR                 ! Pressure thicknesses [R L2 T-2 ~> Pa]
  integer :: lastK_left, lastK_right ! Layers used during the last iteration
  real    :: lastP_left, lastP_right ! Fractional positions during the last iteration [nondim]
  logical :: interior_limit

  ns = 2*nk+2

  ! Initialize variables for the search
  kr = 1 ;
  kl = 1 ;
  lastP_right = 0.
  lastP_left = 0.
  lastK_right = 1
  lastK_left  = 1
  reached_bottom = .false.

  ! Check to see if we should limit the diffusion to the interior
  interior_limit = PRESENT(bl_kl) .and. PRESENT(bl_kr) .and. PRESENT(bl_zr) .and. PRESENT(bl_zl)

  ! Loop over each neutral surface, working from top to bottom
  neutral_surfaces: do k_surface = 1, ns
    klm1 = max(kl-1, 1)
    if (klm1>nk) stop 'find_neutral_surface_positions(): klm1 went out of bounds!'
    krm1 = max(kr-1, 1)
    if (krm1>nk) stop 'find_neutral_surface_positions(): krm1 went out of bounds!'

    ! Potential density difference, rho(kr) - rho(kl)
    dRho = 0.5 * ( ( dRdTr(kr) + dRdTl(kl) ) * ( Tr(kr) - Tl(kl) ) &
                 + ( dRdSr(kr) + dRdSl(kl) ) * ( Sr(kr) - Sl(kl) ) )
    ! Which column has the lighter surface for the current indexes, kr and kl
    if (.not. reached_bottom) then
      if (dRho < 0.) then
        searching_left_column = .true.
        searching_right_column = .false.
      elseif (dRho > 0.) then
        searching_right_column = .true.
        searching_left_column = .false.
      else ! dRho == 0.
        if (kl + kr == 2) then ! Still at surface
          searching_left_column = .true.
          searching_right_column = .false.
        else ! Not the surface so we simply change direction
          searching_left_column = .not.  searching_left_column
          searching_right_column = .not.  searching_right_column
        endif
      endif
    endif

    if (searching_left_column) then
      ! Interpolate for the neutral surface position within the left column, layer klm1
      ! Potential density difference, rho(kl-1) - rho(kr) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTl(klm1) + dRdTr(kr) ) * ( Tl(klm1) - Tr(kr) ) &
                     + ( dRdSl(klm1) + dRdSr(kr) ) * ( Sl(klm1) - Sr(kr) ) )
      ! Potential density difference, rho(kl) - rho(kr) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTl(klm1+1) + dRdTr(kr) ) * ( Tl(klm1+1) - Tr(kr) ) &
                      + ( dRdSl(klm1+1) + dRdSr(kr) ) * ( Sl(klm1+1) - Sr(kr) ) )

      ! Because we are looking left, the right surface, kr, is lighter than klm1+1 and should be denser than klm1
      ! unless we are still at the top of the left column (kl=1)
      if (dRhoTop > 0. .or. kr+kl==2) then
        PoL(k_surface) = 0. ! The right surface is lighter than anything in layer klm1
      elseif (dRhoTop >= dRhoBot) then ! Left layer is unstratified
        PoL(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pl(kl-1) and Pl(kl) where the density difference
        ! between right and left is zero. The Pl here are only used to handle massless layers.
        PoL(k_surface) = interpolate_for_nondim_position( dRhoTop, Pl(klm1), dRhoBot, Pl(klm1+1) )
      endif
      if (PoL(k_surface)>=1. .and. klm1<nk) then ! >= is really ==, when PoL==1 we point to the bottom of the cell
        klm1 = klm1 + 1
        PoL(k_surface) = PoL(k_surface) - 1.
      endif
      if (real(klm1-lastK_left)+(PoL(k_surface)-lastP_left)<0.) then
        PoL(k_surface) = lastP_left
        klm1 = lastK_left
      endif
      KoL(k_surface) = klm1
      if (kr <= nk) then
        PoR(k_surface) = 0.
        KoR(k_surface) = kr
      else
        PoR(k_surface) = 1.
        KoR(k_surface) = nk
      endif
      if (kr <= nk) then
        kr = kr + 1
      else
        reached_bottom = .true.
        searching_right_column = .true.
        searching_left_column = .false.
      endif
    elseif (searching_right_column) then
      ! Interpolate for the neutral surface position within the right column, layer krm1
      ! Potential density difference, rho(kr-1) - rho(kl) (should be negative)
      dRhoTop = 0.5 * ( ( dRdTr(krm1) + dRdTl(kl) ) * ( Tr(krm1) - Tl(kl) ) + &
                        ( dRdSr(krm1) + dRdSl(kl) ) * ( Sr(krm1) - Sl(kl) ) )
      ! Potential density difference, rho(kr) - rho(kl) (will be positive)
      dRhoBot = 0.5 * ( ( dRdTr(krm1+1) + dRdTl(kl) ) * ( Tr(krm1+1) - Tl(kl) ) + &
                        ( dRdSr(krm1+1) + dRdSl(kl) ) * ( Sr(krm1+1) - Sl(kl) ) )

      ! Because we are looking right, the left surface, kl, is lighter than krm1+1 and should be denser than krm1
      ! unless we are still at the top of the right column (kr=1)
      if (dRhoTop >= 0. .or. kr+kl==2) then
        PoR(k_surface) = 0. ! The left surface is lighter than anything in layer krm1
      elseif (dRhoTop >= dRhoBot) then ! Right layer is unstratified
        PoR(k_surface) = 1.
      else
        ! Linearly interpolate for the position between Pr(kr-1) and Pr(kr) where the density difference
        ! between right and left is zero. The Pr here are only used to handle massless layers.
        PoR(k_surface) = interpolate_for_nondim_position( dRhoTop, Pr(krm1), dRhoBot, Pr(krm1+1) )
      endif
      if (PoR(k_surface)>=1. .and. krm1<nk) then ! >= is really ==, when PoR==1 we point to the bottom of the cell
        krm1 = krm1 + 1
        PoR(k_surface) = PoR(k_surface) - 1.
      endif
      if (real(krm1-lastK_right)+(PoR(k_surface)-lastP_right)<0.) then
        PoR(k_surface) = lastP_right
        krm1 = lastK_right
      endif
      KoR(k_surface) = krm1
      if (kl <= nk) then
        PoL(k_surface) = 0.
        KoL(k_surface) = kl
      else
        PoL(k_surface) = 1.
        KoL(k_surface) = nk
      endif
      if (kl <= nk) then
        kl = kl + 1
      else
        reached_bottom = .true.
        searching_right_column = .false.
        searching_left_column = .true.
      endif
    else
      stop 'Else what?'
    endif
    if (interior_limit) then
      if (KoL(k_surface)<=bl_kl) then
        KoL(k_surface) = bl_kl
        if (PoL(k_surface)<bl_zl) then
          PoL(k_surface) = bl_zl
        endif
      endif
      if (KoR(k_surface)<=bl_kr) then
        KoR(k_surface) = bl_kr
        if (PoR(k_surface)<bl_zr) then
          PoR(k_surface) = bl_zr
        endif
      endif
    endif

    lastK_left = KoL(k_surface) ; lastP_left = PoL(k_surface)
    lastK_right = KoR(k_surface) ; lastP_right = PoR(k_surface)
    ! Effective thickness
    ! NOTE: This would be better expressed in terms of the layers thicknesses rather
    ! than as differences of position - AJA
    if (k_surface>1) then
      hL = absolute_position(nk,ns,Pl,KoL,PoL,k_surface) - absolute_position(nk,ns,Pl,KoL,PoL,k_surface-1)
      hR = absolute_position(nk,ns,Pr,KoR,PoR,k_surface) - absolute_position(nk,ns,Pr,KoR,PoR,k_surface-1)
      if ( hL + hR > 0.) then
        hEff(k_surface-1) = 2. * hL * hR / ( hL + hR ) ! Harmonic mean of layer thicknesses
      else
        hEff(k_surface-1) = 0.
      endif
    endif

  enddo neutral_surfaces

end subroutine find_neutral_surface_positions_continuous

!> Returns the non-dimensional position between Pneg and Ppos where the
!! interpolated density difference equals zero.
!! The result is always bounded to be between 0 and 1.
real function interpolate_for_nondim_position(dRhoNeg, Pneg, dRhoPos, Ppos)
  real, intent(in) :: dRhoNeg !< Negative density difference [R ~> kg m-3]
  real, intent(in) :: Pneg    !< Position of negative density difference [R L2 T-2 ~> Pa] or [nondim]
  real, intent(in) :: dRhoPos !< Positive density difference [R ~> kg m-3]
  real, intent(in) :: Ppos    !< Position of positive density difference [R L2 T-2 ~> Pa] or [nondim]

  character(len=120) :: mesg

  if (Ppos < Pneg) then
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Ppos<Pneg')
  elseif (dRhoNeg>dRhoPos) then
    write(stderr,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=',dRhoNeg, Pneg, dRhoPos, Ppos
    write(mesg,*) 'dRhoNeg, Pneg, dRhoPos, Ppos=', dRhoNeg, Pneg, dRhoPos, Ppos
    call MOM_error(WARNING, 'interpolate_for_nondim_position: '//trim(mesg))
  elseif (dRhoNeg>dRhoPos) then !### Does this duplicated test belong here?
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! dRhoNeg>dRhoPos')
  endif
  if (Ppos<=Pneg) then ! Handle vanished or inverted layers
    interpolate_for_nondim_position = 0.5
  elseif ( dRhoPos - dRhoNeg > 0. ) then
    interpolate_for_nondim_position = min( 1., max( 0., -dRhoNeg / ( dRhoPos - dRhoNeg ) ) )
  elseif ( dRhoPos - dRhoNeg == 0) then
    if (dRhoNeg>0.) then
      interpolate_for_nondim_position = 0.
    elseif (dRhoNeg<0.) then
      interpolate_for_nondim_position = 1.
    else ! dRhoPos = dRhoNeg = 0
      interpolate_for_nondim_position = 0.5
    endif
  else ! dRhoPos - dRhoNeg < 0
    interpolate_for_nondim_position = 0.5
  endif
  if ( interpolate_for_nondim_position < 0. ) &
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Pint < Pneg')
  if ( interpolate_for_nondim_position > 1. ) &
    call MOM_error(FATAL, 'interpolate_for_nondim_position: Houston, we have a problem! Pint > Ppos')
end function interpolate_for_nondim_position

!> Higher order version of find_neutral_surface_positions. Returns positions within left/right columns
!! of combined interfaces using intracell reconstructions of T/S. Note that the polynomial reconstrcutions
!! of T and S are optional to aid with unit testing, but will always be passed otherwise
subroutine find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hl, hr, &
                   PoL, PoR, KoL, KoR, hEff, zeta_sbl_L, zeta_sbl_R, k_sbl_L, k_sbl_R)

  type(neutral_diffusion_CS),     intent(inout) :: CS        !< Neutral diffusion control structure
  integer,                        intent(in   ) :: nk        !< Number of levels
  type(column_properties),        intent(inout) :: left_column  !< Properties of the lefthand column
  type(column_properties),        intent(inout) :: right_column !< Properties of the righthand column
  real, dimension(nk),            intent(in   ) :: hl        !< Thicknesses of left column
  real, dimension(nk),            intent(in   ) :: hr        !< Thicknesses of right column
  real, dimension(4*nk),          intent(inout) :: PoL       !< Fractional position of neutral surface within
                                                             !! layer KoL of left column [nondim]
  real, dimension(4*nk),          intent(inout) :: PoR       !< Fractional position of neutral surface within
                                                             !! layer KoR of right column [nondim]
  integer, dimension(4*nk),       intent(inout) :: KoL       !< Index of first left interface above neutral surface
  integer, dimension(4*nk),       intent(inout) :: KoR       !< Index of first right interface above neutral surface
  real, dimension(4*nk-1),        intent(inout) :: hEff      !< Effective thickness between two neutral surfaces
                                                             !! [H ~> m or kg m-2] or other units taken from hcol_l
  real, optional,                 intent(in)    :: zeta_sbl_L!< Non-dimensional distance to where the boundary layer
                                                             !! intersects the cell (left) [nondim]
  real, optional,                 intent(in)    :: zeta_sbl_R!< Non-dimensional distance to where the boundary layer
                                                             !! intersects the cell (right) [nondim]

  integer, optional,              intent(in)    :: k_sbl_L   !< k-index for the boundary layer (left) [nondim]
  integer, optional,              intent(in)    :: k_sbl_R   !< k-index for the boundary layer (right) [nondim]
  ! Local variables
  integer :: kl_left, kl_right      ! Index of layers on the left/right
  integer :: ki_left, ki_right      ! Index of interfaces on the left/right
  logical :: searching_left_column  ! True if searching for the position of a right interface in the left column
  logical :: searching_right_column ! True if searching for the position of a left interface in the right column
  integer :: k_sub                  ! Index of the current sublayer with finite thickness
  integer :: k                      ! Loop variable
  real    :: dRho, dRhoTop, dRhoBot ! A density difference between columns [R ~> kg m-3]
  real    :: pos
  ! Indices referring to the top (1) or bottom interface (2). Note these are the same for the left (l) or 
  ! right (r) but are used as an aid to understanding subroutine calls
  integer, parameter :: top_left = 1, top_right = 1
  integer, parameter :: bot_left = 2, bot_right = 2
  real :: dRdT_top_l, dRdT_top_r
  real :: dRdS_top_l, dRdS_top_r
  real :: p_top_l, p_top_r, p_bot_l, p_bot_r
  integer :: k_top_l, k_top_r, k_bot_l, k_bot_r
  real :: hsub_l, hsub_r
  logical :: interior_only, left_finished, right_finished
  integer :: n_stable_left, n_stable_right
  
  interior_only = PRESENT(zeta_sbl_L) .and. PRESENT(zeta_sbl_R) .and. PRESENT(k_sbl_L) .and. PRESENT(k_sbl_R)
 
  ! Initialize all outputs
  KoL(:) = -1; KoR(:) = -1
  PoL(:) = 0.; PoR(:) = 0.
  hEff(:) = 0.
  ! Return immediately nothing if one column is unstable
  if ( ALL(.not. left_column%stable_cell(:)) ) return
  if ( ALL(.not. right_column%stable_cell(:)) ) return
  kl_left = left_column%first_stable
  kl_right = right_column%first_stable  
  ! Initialize variables for first iteration
  k_sub = 1
  left_finished  = .false. 
  right_finished = .false.
  searching_left_column  = .true.
  searching_right_column = .true.

  ! Sweep down layer by layer to find neutral surfaces
  do k = 1,2*nk
    if (left_finished .or. right_finished) return

    p_top_l = -1; p_top_r = -1
    p_bot_l = -1; p_bot_r = -1
    k_top_l = kl_left; k_top_r = kl_right 
    k_bot_l = kl_left; k_bot_r = kl_right 
    ! 
    call interface_drho_and_derivs(CS, left_column, right_column, kl_left, kl_right, top_left, top_right, dRhoTop)
    if (dRhoTop<0.) then
      searching_left_column  = .true.
      searching_right_column = .false.
      ! No need to switch sign because we will compare in the same direction
      ! dRhoTop = dRhoTop
    elseif (dRhoTop>0.) then
      searching_left_column  = .false. 
      searching_right_column = .true.
      ! Need to switch the sign of dRhoTop because we will compare left to right
    else
      searching_left_column  = .not. searching_left_column
      searching_right_column = .not. searching_right_column
    endif
    if (searching_left_column) then
      call interface_drho_and_derivs(CS, left_column, right_column,                 &
                                     kl_left, kl_right, bot_left, top_right, dRhoBot)
      if (dRhoTop<0. .and. dRhoBot<0.) then ! Right cell is lighter than the left layer
        call next_layer(nk, left_column, kl_left, left_finished)
        cycle
      else
        p_top_l = search_other_column(CS, dRhoTop, dRhoBot, right_column, left_column, &
                                      kl_right, top_right, kl_left)
        p_top_r = 0.
        ! Now check the bot interfaces
        call interface_drho_and_derivs(CS, left_column, right_column,                       & 
                                       kl_left, kl_right, bot_left, bot_right, dRhoBot)
        if (dRhoBot == 0.) then ! Same density
          p_bot_l = 1.
          p_bot_r = 1.
          call next_layer(nk, left_column,  kl_left,  left_finished)
          call next_layer(nk, right_column, kl_right, right_finished)
        elseif (dRhoBot > 0.) then
          ! Right interface is lighter than left interface, so search left
          call interface_drho_and_derivs(CS, left_column, right_column,                    &
                                         kl_left, kl_right, top_left, bot_right, dRhoTop)
          ! Search was the same direction so dRhoBot = dRhoBot
          p_bot_l = search_other_column(CS, dRhoTop, dRhoBot, right_column, left_column, &
                                        kl_right, bot_right, kl_left)
          p_bot_r = 1.
          call next_layer(nk, right_column, kl_right, right_finished)
        elseif (dRhoBot < 0.) then
          ! Left interface is lighter, so search right
          call interface_drho_and_derivs(CS, right_column, left_column, &
                                         kl_right, kl_left, top_right, bot_left, dRhoTop)
          dRhoBot = - dRhoBot
          ! Searched in a different direction
          p_bot_r = search_other_column(CS, dRhoTop, dRhoBot, left_column, right_column, &
                                        kl_left, bot_left, kl_right)
          p_bot_l = 1.
          call next_layer(nk, left_column, kl_left, left_finished)
        endif
      endif
    elseif (searching_right_column) then
      if (dRhoTop /= 0.) dRhoTop = -dRhoTop
      call interface_drho_and_derivs(CS, right_column, left_column,                 &
                                     kl_right, kl_left, bot_right, top_left, dRhobot)
      if (dRhoTop<0. .and. dRhoBot<0.) then ! Left cell is lighter than right layer 
        call next_layer(nk, right_column, kl_right, right_finished)
        cycle
      else
        p_top_r = search_other_column(CS, dRhoTop, dRhoBot, left_column, right_column, &
                                      kl_left, top_left, kl_right)
        p_top_l = 0.
        ! Now check the bot interfaces
        call interface_drho_and_derivs(CS, left_column, right_column,                       & 
                                       kl_left, kl_right, bot_left, bot_right, dRhoBot)
        if (dRhoBot == 0.) then ! Same density
          p_bot_l = 1.
          p_bot_r = 1.
          call next_layer(nk, left_column,  kl_left,  left_finished)
          call next_layer(nk, right_column, kl_right, right_finished)
        elseif (dRhoBot > 0.) then
          ! Right interface is lighter than left interface, so search left
          call interface_drho_and_derivs(CS, left_column, right_column,                    &
                                         kl_left, kl_right, top_left, bot_right, dRhoTop)
          ! Search was the same direction so dRhoBot = dRhoBot
          p_bot_l = search_other_column(CS, dRhoTop, dRhoBot, right_column, left_column, &
                                        kl_right, bot_right, kl_left)
          p_bot_r = 1.
          call next_layer(nk, right_column, kl_right, right_finished)
        elseif (dRhoBot < 0.) then
          ! Left interface is lighter, so search right
          call interface_drho_and_derivs(CS, right_column, left_column, &
                                         kl_right, kl_left, top_right, bot_left, dRhoTop)
          dRhoBot = - dRhoBot
          ! Searched in a different direction
          p_bot_r = search_other_column(CS, dRhoTop, dRhoBot, left_column, right_column, &
                                        kl_left, bot_left, kl_right)
          p_bot_l = 1.
          call next_layer(nk, left_column, kl_left, left_finished)
        endif
      endif
    endif
       
    ! Store the positions of the sublayers
    if (p_bot_l>=0. .and. p_bot_r>=0. .and. p_top_l>=0. .and. p_top_r>=0.) then
      ! Do not do neutral diffusion within the boundary layer by requiring that the
      ! interfaces be no shallower than the boundary layer. This will mean that the
      ! sublayer may not have neutral fluxes
      if (interior_only) then
        if (k_top_l <= k_sbl_l .and. p_top_l <= zeta_sbl_l) then
          k_top_l = k_sbl_l
          p_top_l = zeta_sbl_l
        endif
        if (k_bot_l <= k_sbl_l .and. p_bot_l <= zeta_sbl_l) then
          k_bot_l = k_sbl_l
          p_bot_l = zeta_sbl_l
        endif
      endif
      ! Calculate the harmonic mean 
      hsub_l = hl(k_top_l)*(p_bot_l-p_top_l)
      hsub_r = hr(k_top_r)*(p_bot_r-p_top_r)
      if (hsub_l /= 0. .and. hsub_r /= 0.) then
        hEff(k_sub) = 2.*((hsub_l*hsub_r)/(hsub_l+hsub_r))
      else
        hEff(k_sub) = 0.
      endif
      if (hEff(k_sub) > 0.) then
        KoL(k_sub) = k_top_l; KoR(k_sub) = k_top_r
        PoL(k_sub) = p_top_l; PoR(k_sub) = p_top_r
        k_sub = k_sub + 1
        KoL(k_sub) = k_top_l; KoR(k_sub) = k_top_r
        PoL(k_sub) = p_bot_l; PoR(k_sub) = p_bot_r
        k_sub = k_sub + 1
      endif
    endif
  enddo

end subroutine find_neutral_surface_positions_discontinuous

!> Wrapper for mark_unstable_cells that works for a column_properties type
subroutine column_mark_unstable(CS, nk, column)
  type(neutral_diffusion_CS), intent(in   ) :: CS !< Neutral diffusion control structure
  integer,                    intent(in   ) :: nk !< Nubmer of levels in the column
  type(column_properties),    intent(inout) :: column !< Column to be marked

  call mark_unstable_cells(CS, nk, column%T_at_interface(:,:), column%S_at_interface(:,:), column%P_at_interface,&
                           column%stable_cell(:), column%first_stable, column%last_stable)

end subroutine column_mark_unstable

!> Sweep down through the column and mark as stable if the bottom interface of a cell is denser than the top
subroutine mark_unstable_cells(CS, nk, T, S, P, stable_cell, first_stable, last_stable)
  type(neutral_diffusion_CS), intent(in) :: CS      !< Neutral diffusion control structure
  integer,                intent(in)    :: nk          !< Number of levels in a column
  real, dimension(nk,2),  intent(in)    :: T           !< Temperature at interfaces [degC]
  real, dimension(nk,2),  intent(in)    :: S           !< Salinity at interfaces [ppt]
  real, dimension(nk,2),  intent(in)    :: P           !< Pressure at interfaces [R L2 T-2 ~> Pa]
  logical, dimension(nk), intent(  out) :: stable_cell !< True if this cell is unstably stratified
  integer,                intent(  out) :: first_stable !< Index of the first stable cell
  integer,                intent(  out) :: last_stable  !< Index of the last stable cell

  integer :: k
  real :: delta_rho ! A density difference [R ~> kg m-3]
  first_stable = nk+1
  last_stable  = 0

  do k = 1,nk
    call calc_delta_rho_and_derivs( CS, T(k,2), S(k,2), P(k,2), T(k,1), S(k,1), P(k,1), delta_rho )
    stable_cell(k) = (delta_rho > 0.)
    if (stable_cell(k)) then
      first_stable = MIN(k,first_stable)
      last_stable  = MAX(k,last_stable)
    endif
  enddo
end subroutine mark_unstable_cells

!> Searches the "other" (searched) column for the position of the neutral surface
real function search_other_column(CS, dRhoTop, dRhoBot, from_column, other_column, kl_from, ki_from, kl_other) &
                                  result(pos)
  type(neutral_diffusion_CS), intent(in   ) :: CS       !< Neutral diffusion control structure
  real,                       intent(in   ) :: dRhoTop     !< Difference in density between 'from' interface and
                                                           !! the 'other' column's top interface
  real,                       intent(in   ) :: dRhoBot     !< Difference in density between 'from' interface and
                                                           !! the 'other' column's top interface
  type(column_properties),    intent(in   ) :: from_column !< Properties of column being searched from 
  type(column_properties),    intent(in   ) :: other_column!< Properties of column being searched
  integer,                    intent(in   ) :: kl_from     !< The layer index of the column being searched from
  integer,                    intent(in   ) :: ki_from     !< The interface of the column being searched from
  integer,                    intent(in   ) :: kl_other    !< The layer index of the column to be searched

  real :: dRdT_top, dRdT_bot, dRdS_top, dRdS_bot, dRdT1, dRdT2, dRdS1, dRdS2, drho
  real :: T_other, S_other, P_other, T_from, S_from, P_from

  if ( dRhoBot == 0. ) then      ! Matches perfectly at the bottom
    pos = 1.
  elseif ( dRhoTop == 0. ) then  ! Matches perfectly at the top
    pos = 0.
  elseif (dRhoTop <0. .and. dRhoBot > 0) then ! Connects within a layer
    if (CS%neutral_pos_method==1) then
      pos = interpolate_for_nondim_position( dRhoTop, 0., dRhoBot, 1. )
    ! For the 'Linear' case of finding the neutral position, the fromerence pressure to use is the average
    ! of the midpoint of the layer being searched and the interface being searched from
    elseif (CS%neutral_pos_method == 2) then
      if ( (CS%drho_degree == 2) .or. (CS%drho_degree == 3) ) then
        T_from = from_column%T_at_interface(kl_from, ki_from)
        S_from = from_column%S_at_interface(kl_from, ki_from)
        P_from = from_column%P_at_interface(kl_from, ki_from)
        T_other = from_column%T_at_interface(kl_other, 1)
        S_other = from_column%S_at_interface(kl_other, 1)
        P_other = from_column%P_at_interface(kl_other, 1)
        call calc_delta_rho_and_derivs( CS, T_from, S_from, P_from, T_other, S_other, P_other, drho, &
                                        dRdT1, dRdS1, dRdT2, dRdS2)
        dRdT_top = dRdT1 + dRdT2
        dRdS_top = dRdS1 + dRdS1
        T_other = from_column%T_at_interface(kl_other, 2)
        S_other = from_column%S_at_interface(kl_other, 2)
        P_other = from_column%P_at_interface(kl_other, 2)
        call calc_delta_rho_and_derivs( CS, T_from, S_from, P_from, T_other, S_other, P_other, drho, &
                                        dRdT1, dRdS1, dRdT2, dRdS2)
        dRdT_bot = dRdT1 + dRdT2
        dRdS_bot = dRdS1 + dRdS1
        pos = find_neutral_pos_linear_by_poly( CS, from_column%T_at_interface(kl_from, ki_from),                   &
                                                   from_column%S_at_interface(kl_from, ki_from),                   &
                                                   dRdT_top, dRdS_top, dRdT_bot, dRdS_bot,                         &
                                                   other_column%T_poly(kl_other,:), other_column%S_poly(kl_other,:))
      else
        pos = find_neutral_pos_linear( CS, from_column%T_at_interface(kl_from, ki_from),                   &
                                           from_column%S_at_interface(kl_from, ki_from),                   &
                                           from_column%dRdT(ki_from), from_column%dRdS(ki_from),           &
                                           other_column%dRdT(1), other_column%dRdS(1),                     &
                                           other_column%dRdT(2), other_column%dRdS(2),                     &
                                           other_column%T_poly(kl_other,:), other_column%S_poly(kl_other,:))
      endif
    elseif (CS%neutral_pos_method == 3) then
      pos = find_neutral_pos_full( CS, from_column%T_at_interface(kl_from, ki_from), &
                                       from_column%S_at_interface(kl_from, ki_from), &
                                       from_column%P_at_interface(kl_from, ki_from), &
                                       other_column%P_at_interface(kl_other, 1),     &
                                       other_column%P_at_interface(kl_other, 2),     &
                                       other_column%T_poly(kl_other,:),              &
                                       other_column%S_poly(kl_other,:) )
    endif
  else ! Not within the layer 
    pos = -1
  endif

end function search_other_column

!> Find the next stable layer and whether it would be below the bottom
subroutine next_layer(nk, column, kl, reached_bottom)
  integer,                 intent(in   ) :: nk             !< Number of layers
  type(column_properties), intent(in   ) :: column         !< Column to find the next layer
  integer,                 intent(inout) :: kl             !< Current layer index
  logical,                 intent(  out) :: reached_bottom !< True if the next index is below the bottom

  if (kl == column%last_stable .or. kl == nk) then
    reached_bottom = .true.
  else
    kl = kl + 1
    reached_bottom = .false.
  endif

end subroutine next_layer

!> Similar to find_neutral_pos_linear, except that a root-finding algorithm is applied to the N+1 degree polynomial
!! resulting from using N-degree polynomial representations of T and S in delta rho
function find_neutral_pos_linear_by_poly( CS, T_ref, S_ref, dRdT_top, dRdS_top, dRdT_bot, dRdS_bot, ppoly_T, ppoly_S )&
          result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [ppt]
  real,                      intent(in) :: dRdT_top  !< Sum of alpha at the top interface and the searched from interface
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_top  !< Sum of beta at the searched from interface
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_bot  !< Sum of dRho/dT at bottom interface and searched from interface
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_bot  !< Sum of dRho/dS at bottom interface and searched from interface
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [degC].
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of S within
                                                     !! the layer to be searched [ppt].
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
  real(kind=8), dimension(CS%drho_degree+1) :: drho_poly, drho_poly_reversed
  real(kind=8), dimension(CS%drho_degree)   :: roots
  logical :: failed, degenerate
  integer :: m, drho_degree, nonzero_idx
  real :: p, q, asquared, a, b, c, d, bsquared, discriminant, sub_root_calc, coef, pi

  ! The delta rho polynomial coeffs follow the remapping convention: a1 + a2*x + a3*x^2 ...
  ! Note that a factor of 1/2 can (and has been) be dropped because 
  ! (a1 + a2*x +a3*x^2) = 0 and (0.5)*(a1+a2*x+a3*x^2) = 0 have the same roots
  drho_poly(1) = (ppoly_T(1) - T_ref)*dRdT_top + (ppoly_S(1) - S_ref)*dRdS_top

  drho_poly(2) = ((ppoly_T(1) - T_ref)*(dRdT_bot-dRdT_top) + ppoly_T(2)*dRdT_top) + & 
                 ((ppoly_S(1) - S_ref)*(dRdS_bot-dRdS_top) + ppoly_S(2)*dRdS_top)
  if (CS%drho_degree == 2) then
    drho_poly(3) = ppoly_T(2)*(dRdT_bot - dRdT_top) + ppoly_S(2)*(dRdS_bot - dRdS_top) 
  else if (CS%drho_degree == 3) then
    drho_poly(3) = (ppoly_T(2)*(dRdT_bot-dRdT_top) + ppoly_T(3)*dRdT_top) + &
                   (ppoly_S(2)*(dRdS_bot-dRdS_top) + ppoly_S(3)*dRdS_top)
    drho_poly(4) = ppoly_T(3)*(dRdT_bot-dRdT_top) + ppoly_S(3)*(dRdS_bot-dRdS_top)
  else
    call MOM_error(FATAL,"Combination of reconstruction and EOS approximation not supported")
  endif

  ! Check the coeffs to find the actual degree of the polynomial. For example, if using PPM and a linear EOS
  ! drho_poly(4) would be 0.
  nonzero_idx = CS%drho_degree + 1
  do while( drho_poly(nonzero_idx)==0. )
    nonzero_idx = nonzero_idx - 1
    if (nonzero_idx == 1) exit
  enddo
  drho_degree = nonzero_idx - 1

  if (drho_degree == 0) then
    ! Undefined position
    z = -1
  elseif (drho_degree == 1) then
    if (drho_poly(1) == 0.) then
      roots(1) = 0
    else
      roots(1) = -drho_poly(2)/drho_poly(1)
    endif
  elseif (drho_degree == 2) then    
    roots(:) = -1.
    ! Convert to a monic quadratic
    p = drho_poly(2)/drho_poly(3)
    q = drho_poly(1)/drho_poly(3)
    discriminant = p*p-4*q
    ! Two real roots
    if ( discriminant > 0.) then
      roots(1) = 0.5*(-p + SQRT(discriminant))
      roots(2) = 0.5*(-p - SQRT(discriminant))
      ! Following should be more accurate if one root is much larger than the other
      if ( ABS(roots(1))>ABS(roots(2)) ) then
        roots(2) = q / roots(1)
      else
        roots(1) = q / roots(2)
      endif
      ! No real roots
    elseif ( discriminant < 0. ) then
      roots(:) = -1
    else
      roots(:) = -0.5*p
    endif
  elseif (drho_degree == 3) then
    roots(:) = -1.
    a = drho_poly(4)
    b = drho_poly(3)
    c = drho_poly(2)
    d = drho_poly(1)
    ! Convert to depressed cubic t**3 + p*t + q
    asquared = a*a
    bsquared = b*b
    p = (3.*(a*c) - bsquared)/(3.*asquared)
    q = (2.*(bsquared*b) - 9.*(a*(b*c)) + 27.*(asquared*d))/(27.*(asquared*a))
    discriminant = -(4.*(p**3) + 27.*q*q)
    if (discriminant == 0.) then ! One double root and one simple root
      roots(1) = 3.*(q/p)
      roots(1) = roots(1) - b/(3.*a)
      roots(2) = -1.5*(q/p)
      roots(2) = roots(2) - b/(3.*a)
    elseif (discriminant < 0.) then ! One real root
      sub_root_calc = SQRT(0.25*q*q + p**3/27.)
      roots(1) = cbrt(-0.5*q + sub_root_calc) + cbrt(-0.5*q - sub_root_calc)
      roots(1) = roots(1) - b/(3.*a)
    else ! Three real roots
      pi = 4.*ATAN(1.)
      coef = 2.*SQRT(-p/3.)
      sub_root_calc = ACOS(1.5*((q/p)*SQRT(-3./p)))/3.
      do m = 0,2
        roots(m+1) = coef*COS(sub_root_calc - (2./3.)*(pi*m))
        roots(m+1) = roots(m+1) - b/(3.*a)
      enddo
    endif
  endif

  if (drho_degree > 0) then
    do m = 1,drho_degree
      if (roots(m)<0. .or. roots(m)>1.) then
        roots(m) = -1
      else
        roots(m) = polish_root(roots(m), drho_poly(1:drho_degree+1), drho_degree)
      endif
    enddo
    z = MAXVAL(roots)
  endif

end function find_neutral_pos_linear_by_poly

! Compute the cube root of a value that works for both positive and negative numbers
real function cbrt(x)
  real, intent(in) :: x !< Number for which the cube root will be calculated

  cbrt = SIGN( EXP(LOG(ABS(X))/3.), x )

end function cbrt

!> Use Halley's method to polish off the root
real function polish_root( z, poly_coeffs, degree )
  real, intent(in)          :: z           !< Original root
  real, dimension(degree+1) :: poly_coeffs !< Coefficients of the polynomial 
  integer                   :: degree      !< Degree of the polynomial

  real, dimension(degree)   :: deriv1 !< First derivative
  real, dimension(degree-1) :: deriv2 !< Second derivative
  real :: f, f1, f2, diff_z
  integer :: m

  polish_root = z
  if (degree > 1) then
    call analytic_derivative(poly_coeffs, degree, deriv1)
    call analytic_derivative(deriv1, degree-1, deriv2)
    f  = evaluation_polynomial( poly_coeffs, degree+1, polish_root)
    f1 = evaluation_polynomial( deriv1, degree, polish_root)
    f2 = evaluation_polynomial( deriv2, degree-1, polish_root)
    diff_z = polish_root
    polish_root = polish_root - 2.*((f*f1)/(2.*(f1*f1) - f*f2))
  endif

end function polish_root

!> Search a layer to find where delta_rho = 0 based on a linear interpolation of alpha and beta of the top and bottom
!! being searched and polynomial reconstructions of T and S. Compressibility is not needed because either, we are
!! assuming incompressibility in the equation of state for this module or alpha and beta are calculated having been
!! displaced to the average pressures of the two pressures We need Newton's method because the T and S reconstructions
!! make delta_rho a polynomial function of z if using PPM or higher. If Newton's method would search fall out of the
!! interval [0,1], a bisection step would be taken instead. Also this linearization of alpha, beta means that second
!! derivatives of the EOS are not needed. Note that delta in variable names below refers to horizontal differences and
!! 'd' refers to vertical differences
function find_neutral_pos_linear( CS, T_ref, S_ref, dRdT_ref, dRdS_ref, &
                                  dRdT_top, dRdS_top, dRdT_bot, dRdS_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [ppt]
  real,                      intent(in) :: dRdT_ref  !< dRho/dT at the searched from interface
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_ref  !< dRho/dS at the searched from interface
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_top  !< dRho/dT at top of layer being searched
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_top  !< dRho/dS at top of layer being searched
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real,                      intent(in) :: dRdT_bot  !< dRho/dT at bottom of layer being searched
                                                     !! [R degC-1 ~> kg m-3 degC-1]
  real,                      intent(in) :: dRdS_bot  !< dRho/dS at bottom of layer being searched
                                                     !! [R ppt-1 ~> kg m-3 ppt-1]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [degC].
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of S within
                                                     !! the layer to be searched [ppt].
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
  real :: dRdT_diff  ! Difference in the partial derivative of density with temperature across the
                     ! layer [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS_diff  ! Difference in the partial derivative of density with salinity across the
                     ! layer [R ppt-1 ~> kg m-3 ppt-1]
  real :: drho, drho_dz ! Density anomaly and its derivative with fracitonal position [R ~> kg m-3]
  real :: dRdT_z     ! Partial derivative of density with temperature at a point [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS_z     ! Partial derivative of density with salinity at a point [R ppt-1 ~> kg m-3 ppt-1]
  real :: T_z, dT_dz ! Temperature at a point and its derivative with fractional position [degC]
  real :: S_z, dS_dz ! Salinity at a point and its derivative with fractional position [ppt]
  real :: drho_min, drho_max ! Bounds on density differences [R ~> kg m-3]
  real :: ztest, zmin, zmax ! Fractional positions in the cell [nondim]
  real :: dz         ! Change in position in the cell [nondim]
  real :: a1, a2     ! Fractional weights of the top and bottom values [nondim]
  integer :: iter
  integer :: nterm

  nterm = SIZE(ppoly_T)

  ! Position independent quantities
  dRdT_diff = dRdT_bot - dRdT_top
  dRdS_diff = dRdS_bot - dRdS_top
  ! Initial starting drho (used for bisection)
  zmin = 0        ! Lower bounding interval
  zmax = 1.        ! Maximum bounding interval (bottom of layer)
  a1 = 1. - zmin
  a2 = zmin
  T_z = evaluation_polynomial( ppoly_T, nterm, zmin )
  S_z = evaluation_polynomial( ppoly_S, nterm, zmin )
  dRdT_z = a1*dRdT_top + a2*dRdT_bot
  dRdS_z = a1*dRdS_top + a2*dRdS_bot
  drho_min = 0.5*((dRdT_z+dRdT_ref)*(T_z-T_ref) + (dRdS_z+dRdS_ref)*(S_z-S_ref))

  T_z = evaluation_polynomial( ppoly_T, nterm, 1. )
  S_z = evaluation_polynomial( ppoly_S, nterm, 1. )
  drho_max = 0.5*((dRdT_bot+dRdT_ref)*(T_z-T_ref) + (dRdS_bot+dRdS_ref)*(S_z-S_ref))

  if (drho_min >= 0.) then
    z = zmin
    return
  elseif (drho_max == 0.) then
    z = 1.
    return
  endif
  if ( SIGN(1.,drho_min) == SIGN(1.,drho_max) ) then
    call MOM_error(FATAL, "drho_min is the same sign as dhro_max")
  endif

  z = zmin
  ztest = zmin
  do iter = 1, CS%max_iter
    ! Calculate quantities at the current nondimensional position
    a1 = 1.-z
    a2 = z
    dRdT_z    = a1*dRdT_top + a2*dRdT_bot
    dRdS_z    = a1*dRdS_top + a2*dRdS_bot
    T_z       = evaluation_polynomial( ppoly_T, nterm, z )
    S_z       = evaluation_polynomial( ppoly_S, nterm, z )
    drho = 0.5*((dRdT_z+dRdT_ref)*(T_z-T_ref) + (dRdS_z+dRdS_ref)*(S_z-S_ref))

    ! Check for convergence
    if (ABS(drho) <= CS%drho_tol) exit
    ! Update bisection bracketing intervals
    if (drho < 0. .and. drho > drho_min) then
      drho_min = drho
      zmin = z
    elseif (drho > 0. .and. drho < drho_max) then
      drho_max = drho
      zmax = z
    endif

    ! Calculate a Newton step
    dT_dz = first_derivative_polynomial( ppoly_T, nterm, z )
    dS_dz = first_derivative_polynomial( ppoly_S, nterm, z )
    drho_dz = 0.5*( (dRdT_diff*(T_z - T_ref) + (dRdT_ref+dRdT_z)*dT_dz) + &
                    (dRdS_diff*(S_z - S_ref) + (dRdS_ref+dRdS_z)*dS_dz) )

    ztest = z - drho/drho_dz
    ! Take a bisection if z falls out of [zmin,zmax]
    if (ztest < zmin .or. ztest > zmax) then
      if ( drho < 0. ) then
        ztest = 0.5*(z + zmax)
      else
        ztest = 0.5*(zmin + z)
      endif
    endif

    ! Test to ensure we haven't stalled out
    if ( abs(z-ztest) <= CS%x_tol ) exit
    ! Reset for next iteration
    z = ztest
  enddo

end function find_neutral_pos_linear

!> Use the full equation of state to calculate the difference in locally referenced potential density. The derivatives
!! in this case are not trivial to calculate, so instead we use a regula falsi method
function find_neutral_pos_full( CS, T_ref, S_ref, P_ref, P_top, P_bot, ppoly_T, ppoly_S ) result( z )
  type(neutral_diffusion_CS),intent(in) :: CS        !< Control structure with parameters for this module
  real,                      intent(in) :: T_ref     !< Temperature at the searched from interface [degC]
  real,                      intent(in) :: S_ref     !< Salinity at the searched from interface [ppt]
  real,                      intent(in) :: P_ref     !< Pressure at the searched from interface [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_top     !< Pressure at top of layer being searched [R L2 T-2 ~> Pa]
  real,                      intent(in) :: P_bot     !< Pressure at bottom of layer being searched [R L2 T-2 ~> Pa]
  real, dimension(:),        intent(in) :: ppoly_T   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [degC]
  real, dimension(:),        intent(in) :: ppoly_S   !< Coefficients of the polynomial reconstruction of T within
                                                     !! the layer to be searched [ppt]
  real                                  :: z         !< Position where drho = 0 [nondim]
  ! Local variables
  integer :: iter
  integer :: nterm

  real :: drho_a, drho_b, drho_c ! Density differences [R ~> kg m-3]
  real :: a, b, c     ! Fractional positions [nondim]
  real :: Ta, Tb, Tc  ! Temperatures [degC]
  real :: Sa, Sb, Sc  ! Salinities [ppt]
  real :: Pa, Pb, Pc  ! Pressures [R L2 T-2 ~> Pa]
  integer :: side

  side = 0
  ! Set the first two evaluation to the endpoints of the interval
  b = 0 ; c = 1
  nterm = SIZE(ppoly_T)

  ! Calculate drho at the minimum bound
  Tb = evaluation_polynomial( ppoly_T, nterm, b )
  Sb = evaluation_polynomial( ppoly_S, nterm, b )
  Pb = P_top*(1.-b) + P_bot*b
  call calc_delta_rho_and_derivs(CS, Tb, Sb, Pb, T_ref, S_ref, P_ref, drho_b)

  ! Calculate drho at the maximum bound
  Tc = evaluation_polynomial( ppoly_T, nterm, 1. )
  Sc = evaluation_polynomial( ppoly_S, nterm, 1. )
  Pc = P_Bot
  call calc_delta_rho_and_derivs(CS, Tc, Sc, Pc, T_ref, S_ref, P_ref, drho_c)

  if (drho_b >= 0.) then
    z = 0
    return
  elseif (drho_c == 0.) then
    z = 1.
    return
  endif
  if ( SIGN(1.,drho_b) == SIGN(1.,drho_c) ) then
    z = 0 
    return
  endif

  do iter = 1, CS%max_iter
    ! Calculate new position and evaluate if we have converged
    a = (drho_b*c - drho_c*b)/(drho_b-drho_c)
    Ta = evaluation_polynomial( ppoly_T, nterm, a )
    Sa = evaluation_polynomial( ppoly_S, nterm, a )
    Pa = P_top*(1.-a) + P_bot*a
    call calc_delta_rho_and_derivs(CS, Ta, Sa, Pa, T_ref, S_ref, P_ref, drho_a)
    if (ABS(drho_a) < CS%drho_tol) then
      z = a
      return
    endif

    if (drho_a*drho_c > 0.) then
      if ( ABS(a-c)<CS%x_tol) then
        z = a
        return
      endif
      c = a ; drho_c = drho_a;
      if (side == -1) drho_b = 0.5*drho_b
      side = -1
    elseif ( drho_b*drho_a > 0 ) then
      if ( ABS(a-b)<CS%x_tol) then
        z = a
        return
      endif
      b = a ; drho_b = drho_a
      if (side == 1) drho_c = 0.5*drho_c
      side = 1
    else
      z = a
      return
    endif
  enddo

  z = a

end function find_neutral_pos_full

!> Calculate difference in density between interfaces of two columns and update derivates
subroutine interface_drho_and_derivs(CS, col1, col2, kl_1, kl_2, ki_1, ki_2, delta_rho)
  type(neutral_diffusion_CS), intent(in   ) :: CS        !< Neutral diffusion control structure
  type(column_properties)   , intent(inout) :: col1      !< First column
  type(column_properties)   , intent(inout) :: col2      !< Second column
  integer                   , intent(in   ) :: kl_1      !< Layer index of first column
  integer                   , intent(in   ) :: kl_2      !< Layer index of second column
  integer                   , intent(in   ) :: ki_1      !< Top (1) or bottom (2) interface of column 1
  integer                   , intent(in   ) :: ki_2      !< Top (1) or bottom (2) interface of column 2
  real                      , intent(  out) :: delta_rho
  
  real :: T1, S1 , P1, T2, S2, P2

  T1 = col1%T_at_interface(kl_1, ki_1)
  S1 = col1%S_at_interface(kl_1, ki_1)
  P1 = col1%P_at_interface(kl_1, ki_1)
  T2 = col2%T_at_interface(kl_2, ki_2)
  S2 = col2%S_at_interface(kl_2, ki_2)
  P2 = col2%P_at_interface(kl_2, ki_2)

  call calc_delta_rho_and_derivs(CS, T1, S1, P1, T2, S2, P2, delta_rho,                            &
                                 col1%dRdT(ki_1), col1%dRdS(ki_1), col2%dRdT(ki_2), col2%dRdS(ki_2))
end

!> Calculate the difference in density between two points in a variety of ways
subroutine calc_delta_rho_and_derivs(CS, T1, S1, p1_in, T2, S2, p2_in, drho, &
                                     drdt1_out, drds1_out, drdt2_out, drds2_out )
  type(neutral_diffusion_CS)    :: CS        !< Neutral diffusion control structure
  real,           intent(in   ) :: T1        !< Temperature at point 1 [degC]
  real,           intent(in   ) :: S1        !< Salinity at point 1 [ppt]
  real,           intent(in   ) :: p1_in     !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real,           intent(in   ) :: T2        !< Temperature at point 2 [degC]
  real,           intent(in   ) :: S2        !< Salinity at point 2 [ppt]
  real,           intent(in   ) :: p2_in     !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real,           intent(  out) :: drho      !< Difference in density between the two points [R ~> kg m-3]
  real, optional, intent(  out) :: dRdT1_out !< drho_dt at point 1 [R degC-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS1_out !< drho_ds at point 1 [R ppt-1 ~> kg m-3 ppt-1]
  real, optional, intent(  out) :: dRdT2_out !< drho_dt at point 2 [R degC-1 ~> kg m-3 degC-1]
  real, optional, intent(  out) :: dRdS2_out !< drho_ds at point 2 [R ppt-1 ~> kg m-3 ppt-1]
  ! Local variables
  real :: rho1, rho2   ! Densities [R ~> kg m-3]
  real :: p1, p2, pmid ! Pressures [R L2 T-2 ~> Pa]
  real :: drdt1, drdt2 ! Partial derivatives of density with temperature [R degC-1 ~> kg m-3 degC-1]
  real :: drds1, drds2 ! Partial derivatives of density with salinity [R ppt-1 ~> kg m-3 ppt-1]
  real :: drdp1, drdp2 ! Partial derivatives of density with pressure [T2 L-2 ~> s2 m-2]

  ! Use the same reference pressure or the in-situ pressure
  if (CS%ref_pres > 0.) then
    p1 = CS%ref_pres
    p2 = CS%ref_pres
  else
    p1 = p1_in
    p2 = p2_in
  endif
  ! Use the full linear equation of state to calculate the difference in density (expensive!)
  if     (TRIM(CS%delta_rho_form) == 'full') then
    pmid = 0.5 * (p1 + p2)
    call calculate_density( T1, S1, pmid, rho1, CS%EOS)
    call calculate_density( T2, S2, pmid, rho2, CS%EOS)
    drho = rho1 - rho2
  ! Use the density derivatives at the average of pressures and the differentces int temperature
  elseif (TRIM(CS%delta_rho_form) == 'mid_pressure') then
    pmid = 0.5 * (p1 + p2)
    if (CS%ref_pres>=0) pmid = CS%ref_pres
    call calculate_density_derivs(T1, S1, pmid, drdt1, drds1, CS%EOS)
    call calculate_density_derivs(T2, S2, pmid, drdt2, drds2, CS%EOS)
    drho = delta_rho_from_derivs( T1, S1, p1, drdt1, drds1, T2, S2, p2, drdt2, drds2)
  elseif (TRIM(CS%delta_rho_form) == 'local_pressure') then
    call calculate_density_derivs(T1, S1, p1, drdt1, drds1, CS%EOS)
    call calculate_density_derivs(T2, S2, p2, drdt2, drds2, CS%EOS)
    drho = delta_rho_from_derivs( T1, S1, p1, drdt1, drds1, T2, S2, p2, drdt2, drds2)
  else
    call MOM_error(FATAL, "delta_rho_form is not recognized")
  endif

  if (PRESENT(drdt1_out)) drdt1_out = drdt1
  if (PRESENT(drds1_out)) drds1_out = drds1
  if (PRESENT(drdt2_out)) drdt2_out = drdt2
  if (PRESENT(drds2_out)) drds2_out = drds2

end subroutine calc_delta_rho_and_derivs

!> Calculate delta rho from derivatives and gradients of properties
!! drho < 0 implies that parcel 2 is lighter than parcel 1
function delta_rho_from_derivs( T1, S1, P1, dRdT1, dRdS1, &
                                T2, S2, P2, dRdT2, dRdS2  ) result (drho)
  real :: T1    !< Temperature at point 1 [degC]
  real :: S1    !< Salinity at point 1 [ppt]
  real :: P1    !< Pressure at point 1 [R L2 T-2 ~> Pa]
  real :: dRdT1 !< The partial derivative of density with temperature at point 1 [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS1 !< The partial derivative of density with salinity at point 1 [R ppt-1 ~> kg m-3 ppt-1]
  real :: T2    !< Temperature at point 2 [degC]
  real :: S2    !< Salinity at point 2 [ppt]
  real :: P2    !< Pressure at point 2 [R L2 T-2 ~> Pa]
  real :: dRdT2 !< The partial derivative of density with temperature at point 2 [R degC-1 ~> kg m-3 degC-1]
  real :: dRdS2 !< The partial derivative of density with salinity at point 2 [R ppt-1 ~> kg m-3 ppt-1]
  ! Local variables
  real :: drho  ! The density difference [R ~> kg m-3]

  drho = 0.5 * ( (dRdT1+dRdT2)*(T1-T2) + (dRdS1+dRdS2)*(S1-S2))

end function delta_rho_from_derivs

!> Converts non-dimensional position within a layer to absolute position (for debugging)
function absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  integer, intent(in) :: n            !< Number of levels
  integer, intent(in) :: ns           !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1)    !< Position of interfaces [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)     !< Index of interface above position
  real,    intent(in) :: NParr(ns)    !< Non-dimensional position within layer Karr(:) [nondim]
  integer, intent(in) :: k_surface    !< k-interface to query
  real                :: absolute_position !< The absolute position of a location [R L2 T-2 ~> Pa]
                                      !! or other units following Pint
  ! Local variables
  integer :: k

  k = Karr(k_surface)
  if (k>n) stop 'absolute_position: k>nk is out of bounds!'
  absolute_position = Pint(k) + NParr(k_surface) * ( Pint(k+1) - Pint(k) )

end function absolute_position

!> Converts non-dimensional positions within layers to absolute positions (for debugging)
function absolute_positions(n,ns,Pint,Karr,NParr)
  integer, intent(in) :: n         !< Number of levels
  integer, intent(in) :: ns        !< Number of neutral surfaces
  real,    intent(in) :: Pint(n+1) !< Position of interface [R L2 T-2 ~> Pa] or other units
  integer, intent(in) :: Karr(ns)  !< Indexes of interfaces about positions
  real,    intent(in) :: NParr(ns) !< Non-dimensional positions within layers Karr(:)

  real,  dimension(ns) :: absolute_positions !< Absolute positions [R L2 T-2 ~> Pa]
                                   !! or other units following Pint

  ! Local variables
  integer :: k_surface, k

  do k_surface = 1, ns
    absolute_positions(k_surface) = absolute_position(n,ns,Pint,Karr,NParr,k_surface)
  enddo

end function absolute_positions

!> Returns a single column of neutral diffusion fluxes of a tracer.
subroutine neutral_surface_flux(nk, nsurf, deg, hl, hr, Tl, Tr, PiL, PiR, KoL, KoR, &
                                hEff, Flx, continuous, h_neglect, remap_CS, h_neglect_edge)
  integer,                      intent(in)    :: nk    !< Number of levels
  integer,                      intent(in)    :: nsurf !< Number of neutral surfaces
  integer,                      intent(in)    :: deg   !< Degree of polynomial reconstructions
  real, dimension(nk),          intent(in)    :: hl    !< Left-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: hr    !< Right-column layer thickness [H ~> m or kg m-2]
  real, dimension(nk),          intent(in)    :: Tl    !< Left-column layer tracer (conc, e.g. degC)
  real, dimension(nk),          intent(in)    :: Tr    !< Right-column layer tracer (conc, e.g. degC)
  real, dimension(nsurf),       intent(in)    :: PiL   !< Fractional position of neutral surface
                                                       !! within layer KoL of left column
  real, dimension(nsurf),       intent(in)    :: PiR   !< Fractional position of neutral surface
                                                       !! within layer KoR of right column
  integer, dimension(nsurf),    intent(in)    :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(nsurf),    intent(in)    :: KoR   !< Index of first right interface above neutral surface
  real, dimension(nsurf-1),     intent(in)    :: hEff  !< Effective thickness between two neutral
                                                       !! surfaces [H ~> m or kg m-2]
  real, dimension(nsurf-1),     intent(inout) :: Flx   !< Flux of tracer between pairs of neutral layers (conc H)
  logical,                      intent(in)    :: continuous !< True if using continuous reconstruction
  real,                         intent(in)    :: h_neglect !< A negligibly small width for the
                                             !! purpose of cell reconstructions [H ~> m or kg m-2]
  type(remapping_CS), optional, intent(in)    :: remap_CS !< Remapping control structure used
                                             !! to create sublayers
  real,               optional, intent(in)    :: h_neglect_edge !< A negligibly small width used for
                                             !! edge value calculations if continuous is false [H ~> m or kg m-2]
  ! Local variables
  integer :: k_sublayer, klb, klt, krb, krt, k
  real :: T_right_top, T_right_bottom, T_right_layer, T_right_sub, T_right_top_int, T_right_bot_int
  real :: T_left_top, T_left_bottom, T_left_layer, T_left_sub, T_left_top_int, T_left_bot_int
  real :: dT_top, dT_bottom, dT_layer, dT_ave, dT_sublayer, dT_top_int, dT_bot_int
  real, dimension(nk+1) :: Til !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk+1) :: Tir !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_l !< Left-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_l !< Left-column right edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aL_r !< Right-column left edge value of tracer (conc, e.g. degC)
  real, dimension(nk) :: aR_r !< Right-column right edge value of tracer (conc, e.g. degC)
  ! Discontinuous reconstruction
  integer               :: iMethod
  real, dimension(nk,2) :: Tid_l !< Left-column interface tracer (conc, e.g. degC)
  real, dimension(nk,2) :: Tid_r !< Right-column interface tracer (conc, e.g. degC)
  real, dimension(nk,deg+1) :: ppoly_r_coeffs_l
  real, dimension(nk,deg+1) :: ppoly_r_coeffs_r
  real, dimension(nk,deg+1) :: ppoly_r_S_l
  real, dimension(nk,deg+1) :: ppoly_r_S_r
  logical :: down_flux
  ! Setup reconstruction edge values
  if (continuous) then
    call interface_scalar(nk, hl, Tl, Til, 2, h_neglect)
    call interface_scalar(nk, hr, Tr, Tir, 2, h_neglect)
    call ppm_left_right_edge_values(nk, Tl, Til, aL_l, aR_l)
    call ppm_left_right_edge_values(nk, Tr, Tir, aL_r, aR_r)
  else
    ppoly_r_coeffs_l(:,:) = 0.
    ppoly_r_coeffs_r(:,:) = 0.
    Tid_l(:,:) = 0.
    Tid_r(:,:) = 0.

    call build_reconstructions_1d( remap_CS, nk, hl, Tl, ppoly_r_coeffs_l, Tid_l, &
                                   ppoly_r_S_l, iMethod, h_neglect, h_neglect_edge )
    call build_reconstructions_1d( remap_CS, nk, hr, Tr, ppoly_r_coeffs_r, Tid_r, &
                                   ppoly_r_S_r, iMethod, h_neglect, h_neglect_edge )
  endif

  do k_sublayer = 1, nsurf-1
    if ( (KoL(k_sublayer) < 0) .or. (KoR(k_sublayer)<0) ) exit
    if (hEff(k_sublayer) == 0.) then
      Flx(k_sublayer) = 0.
    else
      if (continuous) then
        klb = KoL(k_sublayer+1)
        T_left_bottom = ( 1. - PiL(k_sublayer+1) ) * Til(klb) + PiL(k_sublayer+1) * Til(klb+1)
        klt = KoL(k_sublayer)
        T_left_top = ( 1. - PiL(k_sublayer) ) * Til(klt) + PiL(k_sublayer) * Til(klt+1)
        T_left_layer = ppm_ave(PiL(k_sublayer), PiL(k_sublayer+1) + real(klb-klt), &
                               aL_l(klt), aR_l(klt), Tl(klt))

        krb = KoR(k_sublayer+1)
        T_right_bottom = ( 1. - PiR(k_sublayer+1) ) * Tir(krb) + PiR(k_sublayer+1) * Tir(krb+1)
        krt = KoR(k_sublayer)
        T_right_top = ( 1. - PiR(k_sublayer) ) * Tir(krt) + PiR(k_sublayer) * Tir(krt+1)
        T_right_layer = ppm_ave(PiR(k_sublayer), PiR(k_sublayer+1) + real(krb-krt), &
                                aL_r(krt), aR_r(krt), Tr(krt))
        dT_top = T_right_top - T_left_top
        dT_bottom = T_right_bottom - T_left_bottom
        dT_ave = 0.5 * ( dT_top + dT_bottom )
        dT_layer = T_right_layer - T_left_layer
        if (signum(1.,dT_top) * signum(1.,dT_bottom) <= 0. .or. signum(1.,dT_ave) * signum(1.,dT_layer) <= 0.) then
          dT_ave = 0.
        else
          dT_ave = dT_layer
        endif
        Flx(k_sublayer) = dT_ave * hEff(k_sublayer)
      else ! Discontinuous reconstruction
        ! Calculate tracer values on left and right side of the neutral surface
        call neutral_surface_T_eval(nk, nsurf, k_sublayer, KoL, PiL, Tl, Tid_l, deg, iMethod, &
                                    ppoly_r_coeffs_l, T_left_top, T_left_bottom, T_left_sub, &
                                    T_left_top_int, T_left_bot_int, T_left_layer)
        call neutral_surface_T_eval(nk, nsurf, k_sublayer, KoR, PiR, Tr, Tid_r, deg, iMethod, &
                                    ppoly_r_coeffs_r, T_right_top, T_right_bottom, T_right_sub, &
                                    T_right_top_int, T_right_bot_int, T_right_layer)

        dT_top      = T_right_top     - T_left_top
        dT_bottom   = T_right_bottom  - T_left_bottom
        dT_sublayer = T_right_sub     - T_left_sub
        dT_top_int  = T_right_top_int - T_left_top_int
        dT_bot_int  = T_right_bot_int - T_left_bot_int
        ! Enforcing the below criterion incorrectly zero out fluxes
        !dT_layer = T_right_layer - T_left_layer

        down_flux = dT_top <= 0. .and. dT_bottom <= 0. .and.       &
                    dT_sublayer <= 0. .and. dT_top_int <= 0. .and. &
                    dT_bot_int <= 0.
        down_flux = down_flux .or.                                 &
                    (dT_top >= 0. .and. dT_bottom >= 0. .and.      &
                    dT_sublayer >= 0. .and. dT_top_int >= 0. .and. &
                    dT_bot_int >= 0.)
        if (down_flux) then
          Flx(k_sublayer) = dT_sublayer * hEff(k_sublayer)
        else
          Flx(k_sublayer) = 0.
        endif
      endif
    endif
  enddo

end subroutine neutral_surface_flux

!> Evaluate various parts of the reconstructions to calculate gradient-based flux limter
subroutine neutral_surface_T_eval(nk, ns, k_sub, Ks, Ps, T_mean, T_int, deg, iMethod, T_poly, &
                                  T_top, T_bot, T_sub, T_top_int, T_bot_int, T_layer)
  integer,                   intent(in   ) :: nk        !< Number of cell everages
  integer,                   intent(in   ) :: ns        !< Number of neutral surfaces
  integer,                   intent(in   ) :: k_sub     !< Index of current neutral layer
  integer, dimension(ns),    intent(in   ) :: Ks        !< List of the layers associated with each neutral surface
  real, dimension(ns),       intent(in   ) :: Ps        !< List of the positions within a layer of each surface
  real, dimension(nk),       intent(in   ) :: T_mean    !< Cell average of tracer
  real, dimension(nk,2),     intent(in   ) :: T_int     !< Cell interface values of tracer from reconstruction
  integer,                   intent(in   ) :: deg       !< Degree of reconstruction polynomial (e.g. 1 is linear)
  integer,                   intent(in   ) :: iMethod   !< Method of integration to use
  real, dimension(nk,deg+1), intent(in   ) :: T_poly    !< Coefficients of polynomial reconstructions
  real,                      intent(  out) :: T_top     !< Tracer value at top (across discontinuity if necessary)
  real,                      intent(  out) :: T_bot     !< Tracer value at bottom (across discontinuity if necessary)
  real,                      intent(  out) :: T_sub     !< Average of the tracer value over the sublayer
  real,                      intent(  out) :: T_top_int !< Tracer value at top interface of neutral layer
  real,                      intent(  out) :: T_bot_int !< Tracer value at bottom interface of neutral layer
  real,                      intent(  out) :: T_layer   !< Cell-average that the the reconstruction belongs to

  integer :: kl, ks_top, ks_bot

  ks_top = k_sub
  ks_bot = k_sub + 1
  if ( Ks(ks_top) /= Ks(ks_bot) ) then
    call MOM_error(FATAL, "Neutral surfaces span more than one layer")
  endif
  kl = Ks(k_sub)
  ! First if the neutral surfaces spans the entirety of a cell, then do not search across the discontinuity
  if ( (Ps(ks_top) == 0.) .and. (Ps(ks_bot) == 1.)) then
    T_top = T_int(kl,1)
    T_bot = T_int(kl,2)
  else
    ! Search across potential discontinuity at top
    if ( (kl > 1) .and. (Ps(ks_top) == 0.)  ) then
      T_top = T_int(kl-1,2)
    else
      T_top = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_top) )
    endif
    ! Search across potential discontinuity at bottom
    if ( (kl < nk) .and. (Ps(ks_bot) == 1.) ) then
      T_bot = T_int(kl+1,1)
    else
      T_bot = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_bot) )
    endif
  endif
  T_sub = average_value_ppoly(nk, T_mean, T_int, T_poly, iMethod, kl, Ps(ks_top), Ps(ks_bot))
  T_top_int = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_top))
  T_bot_int = evaluation_polynomial( T_poly(kl,:), deg+1, Ps(ks_bot))
  T_layer = T_mean(kl)

end subroutine neutral_surface_T_eval

!> Discontinuous PPM reconstructions of the left/right edge values within a cell
subroutine ppm_left_right_edge_values(nk, Tl, Ti, aL, aR)
  integer,                    intent(in)    :: nk !< Number of levels
  real, dimension(nk),        intent(in)    :: Tl !< Layer tracer (conc, e.g. degC)
  real, dimension(nk+1),      intent(in)    :: Ti !< Interface tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aL !< Left edge value of tracer (conc, e.g. degC)
  real, dimension(nk),        intent(inout) :: aR !< Right edge value of tracer (conc, e.g. degC)

  integer :: k
  ! Setup reconstruction edge values
  do k = 1, nk
    aL(k) = Ti(k)
    aR(k) = Ti(k+1)
    if ( signum(1., aR(k) - Tl(k))*signum(1., Tl(k) - aL(k)) <= 0.0 ) then
      aL(k) = Tl(k)
      aR(k) = Tl(k)
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) > abs(aR(k) - aL(k)) ) then
      aL(k) = Tl(k) + 2.0 * ( Tl(k) - aR(k) )
    elseif ( sign(3., aR(k) - aL(k)) * ( (Tl(k) - aL(k)) + (Tl(k) - aR(k))) < -abs(aR(k) - aL(k)) ) then
      aR(k) = Tl(k) + 2.0 * ( Tl(k) - aL(k) )
    endif
  enddo
end subroutine ppm_left_right_edge_values

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function neutral_diffusion_unit_tests(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout

  neutral_diffusion_unit_tests = .false. .or. &
    ndiff_unit_tests_continuous(verbose) .or. ndiff_unit_tests_discontinuous(verbose)

end function neutral_diffusion_unit_tests

!> Returns true if unit tests of neutral_diffusion functions fail. Otherwise returns false.
logical function ndiff_unit_tests_continuous(verbose)
  logical, intent(in) :: verbose !< If true, write results to stdout
  ! Local variables
  integer, parameter         :: nk = 4
  real, dimension(nk+1)      :: TiL, TiR1, TiR2, TiR4, Tio ! Test interface temperatures
  real, dimension(nk)        :: TL                         ! Test layer temperatures
  real, dimension(nk+1)      :: SiL                        ! Test interface salinities
  real, dimension(nk+1)      :: PiL, PiR4                  ! Test interface positions
  real, dimension(2*nk+2)    :: PiLRo, PiRLo               ! Test positions
  integer, dimension(2*nk+2) :: KoL, KoR                   ! Test indexes
  real, dimension(2*nk+1)    :: hEff                       ! Test positions
  real, dimension(2*nk+1)    :: Flx                        ! Test flux
  integer :: k
  logical :: v
  real :: h_neglect

  h_neglect = 1.0e-30

  v = verbose

  ndiff_unit_tests_continuous = .false. ! Normally return false
  write(stdout,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_continuous ='

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,1.,1., 0.,1.,2., 1., 'FV: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,1.,0., 0.,4.,8., 7., 'FV: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,1.,1., 0.,4.,8., 7., 'FV: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,2.,4., 0.,3.,9., 4., 'FV: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,2.,0.,2., 0.,1.,2., 0., 'FV: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,1.,0., 0.,1.,2., 2., 'FV: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,1.,0.,0., 0.,1.,2., 0., 'FV: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fv_diff(v,0.,0.,0., 0.,1.,2., 0., 'FV: All vanished cells')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,1.,1., 0.,1.,2., 1., 'LSQ: Straight line on uniform grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,1.,0., 0.,1.,2., 1., 'LSQ: Vanished right cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,1.,1., 0.,1.,2., 1., 'LSQ: Vanished left cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,2.,4., 0.,3.,9., 2., 'LSQ: Stretched grid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,0.,1., 0.,1.,2., 2., 'LSQ: Vanished middle cell')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,1.,0., 0.,1.,2., 0., 'LSQ: Vanished on both sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,1.,0.,0., 0.,1.,2., 0., 'LSQ: Two vanished cell sides')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_fvlsq_slope(v,0.,0.,0., 0.,1.,2., 0., 'LSQ: All vanished cells')

  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 1, h_neglect)
  !ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
  !  test_data1d(5, Tio, (/27.,21.,15.,9.,3./), 'Linear profile, interface temperatures')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_data1d(v,5, Tio, (/24.,22.5,15.,7.5,6./), 'Linear profile, linear interface temperatures')
  call interface_scalar(4, (/10.,10.,10.,10./), (/24.,18.,12.,6./), Tio, 2, h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_data1d(v,5, Tio, (/24.,22.,15.,8.,6./), 'Linear profile, PPM interface temperatures')

  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  1.0, 1.0, 0.5, 'Check mid-point')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.0, 0.,  1.0, 1.0, 0.0, 'Check bottom')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.1, 0.,  1.1, 1.0, 0.0, 'Check below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  0.0, 1.0, 1.0, 'Check top')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0., -0.1, 1.0, 1.0, 'Check above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0.,  3.0, 1.0, 0.25, 'Check 1/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-3.0, 0.,  1.0, 1.0, 0.75, 'Check 3/4')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 1.0, 0.,  1.0, 1.0, 0.0, 'Check dRho=0 below')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-1.0, 0., -1.0, 1.0, 1.0, 'Check dRho=0 above')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v, 0.0, 0.,  0.0, 1.0, 0.5, 'Check dRho=0 mid')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. &
    test_ifndp(v,-2.0, .5,  5.0, 0.5, 0.5, 'Check dP=0')

  ! Identical columns
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! KoL
                                   (/1,1,2,2,3,3,3,3/), & ! KoR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Identical columns')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... left positions')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,10.,10.,20.,20.,30.,30./), '... right positions')
  call neutral_surface_flux(3, 2*3+2, 2, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/20.,16.,12./), (/20.,16.,12./), & ! Tl, Tr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true., h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 7, Flx, &
              (/0.,0.,0.,0.,0.,0.,0./), 'Identical columns, rho flux (=0)')
  call neutral_surface_flux(3, 2*3+2, 2, (/10.,10.,10./), (/10.,10.,10./), & ! nk, hL, hR
                               (/-1.,-1.,-1./), (/1.,1.,1./), & ! Sl, Sr
                               PiLRo, PiRLo, KoL, KoR, hEff, Flx, .true., h_neglect)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 7, Flx, &
              (/0.,20.,0.,20.,0.,20.,0./), 'Identical columns, S flux')

  ! Right column slightly cooler than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/20.,16.,12.,8./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,1,2,2,3,3,3/), & ! kR
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pL
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly cooler')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoL, PiLRo), &
                                   (/0.,5.,10.,15.,20.,25.,30.,30./), '... left positions')
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or. test_data1d(v, 8, &
                                   absolute_positions(3, 8, (/0.,10.,20.,30./), KoR, PiRLo), &
                                   (/0.,0.,5.,10.,15.,20.,25.,30./), '... right positions')

  ! Right column slightly warmer than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/24.,20.,16.,12./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,2,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,0.,0.5,1./), & ! pL
                                   (/0.,0.5,0.,0.5,0.,0.5,1.,1./), & ! pR
                                   (/0.,5.,5.,5.,5.,5.,0./), & ! hEff
                                   'Right column slightly warmer')

  ! Right column somewhat cooler than left
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/16.,12.,8.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,2,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,2,3,3/), & ! kR
                                   (/0.,0.,0.5,0.,0.5,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.,0.5,0.,1./), & ! pR
                                   (/0.,0.,5.,5.,5.,0.,0./), & ! hEff
                                   'Right column somewhat cooler')

  ! Right column much colder than left with no overlap
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/9.,7.,5.,3./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,1,2,3,3/), & ! kR
                                   (/0.,0.,0.,1.,1.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,0.,0./), & ! hEff
                                   'Right column much cooler')

  ! Right column with mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/22.,18.,14.,10./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,1,2,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,1.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,0.,1./), & ! pR
                                   (/0.,0.,0.,0.,10.,0.,0./), & ! hEff
                                   'Right column with mixed layer')

  ! Identical columns with mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,2,2,3,3,3,3/), & ! kL
                                   (/1,1,2,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,1.,1./), & ! pR
                                   (/0.,10.,0.,10.,0.,10.,0./), & ! hEff
                                   'Identical columns with mixed layer')

  ! Right column with unstable mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,2,3,3,3,3,3,3/), & ! kL
                                   (/1,1,1,2,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Right column with unstable mixed layer')

  ! Left column with unstable mixed layer
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/14.,14.,10.,2./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,2,3,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.25,1.,1./), & ! pL
                                   (/0.,0.,0.,0.,0.,0.,.75,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,7.5,0./), & ! hEff
                                   'Left column with unstable mixed layer')

  ! Two unstable mixed layers
  call find_neutral_surface_positions_continuous(3, &
             (/0.,10.,20.,30./), (/8.,12.,10.,2./), (/0.,0.,0.,0./), & ! Left positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Left dRdT and dRdS
             (/0.,10.,20.,30./), (/10.,14.,12.,4./), (/0.,0.,0.,0./), & ! Right positions, T and S
             (/-1.,-1.,-1.,-1./), (/1.,1.,1.,1./), &! Right dRdT and dRdS
             PiLRo, PiRLo, KoL, KoR, hEff)
  ndiff_unit_tests_continuous = ndiff_unit_tests_continuous .or.  test_nsp(v, 8, KoL, KoR, PiLRo, PiRLo, hEff, &
                                   (/1,1,1,1,2,3,3,3/), & ! kL
                                   (/1,2,3,3,3,3,3,3/), & ! kR
                                   (/0.,0.,0.,0.,0.,0.,0.75,1./), & ! pL
                                   (/0.,0.,0.,0.5,0.5,0.5,1.,1./), & ! pR
                                   (/0.,0.,0.,0.,0.,6.,0./), & ! hEff
                                   'Two unstable mixed layers')

  if (.not. ndiff_unit_tests_continuous) write(stdout,*) 'Pass'

end function ndiff_unit_tests_continuous

logical function ndiff_unit_tests_discontinuous(verbose)
  logical, intent(in) :: verbose !< It true, write results to stdout
  ! Local variables
  integer, parameter          :: nk = 3
  integer, parameter          :: ns = nk*4
  real, dimension(nk)         :: Sl, Sr, Tl, Tr ! Salinities [ppt] and temperatures [degC]
  real, dimension(nk)         :: hl, hr    ! Thicknesses in pressure units [R L2 T-2 ~> Pa]
  real, dimension(nk,2)       :: TiL, SiL, TiR, SiR ! Cell edge salinities [ppt] and temperatures [degC]
  real, dimension(nk,2)       :: Pres_l, Pres_r ! Interface pressures [R L2 T-2 ~> Pa]
  integer, dimension(ns)      :: KoL, KoR
  real, dimension(ns)         :: PoL, PoR
  real, dimension(ns-1)       :: hEff, Flx
  type(neutral_diffusion_CS)  :: CS        !< Neutral diffusion control structure
  type(EOS_type),     pointer :: EOS       !< Structure for linear equation of state
  type(remapping_CS), pointer :: remap_CS  !< Remapping control structure (PLM)
  real, dimension(nk,2)       :: ppoly_T_l, ppoly_T_r ! Linear reconstruction for T
  real, dimension(nk,2)       :: ppoly_S_l, ppoly_S_r ! Linear reconstruction for S
  real, dimension(nk,2)       :: dRdT      !< Partial derivative of density with temperature at
                                           !! cell edges [R degC-1 ~> kg m-3 degC-1]
  real, dimension(nk,2)       :: dRdS      !< Partial derivative of density with salinity at
                                           !! cell edges [R ppt-1 ~> kg m-3 ppt-1]
  logical, dimension(nk)      :: stable_l, stable_r
  integer                     :: iMethod
  integer                     :: ns_l, ns_r
  integer :: k
  logical :: v
  type(column_properties) :: left_column, right_column

  ndiff_unit_tests_discontinuous = .false. ! Normally return false
  v = verbose
  write(stdout,*) '==== MOM_neutral_diffusion: ndiff_unit_tests_discontinuous ='
 
  ! Initialize the 'column' types
  allocate(left_column%T_at_interface(nk,2))
  allocate(left_column%S_at_interface(nk,2))
  allocate(left_column%P_at_interface(nk,2))
  allocate(left_column%T_poly(nk,2))
  allocate(left_column%stable_cell(nk))
  allocate(right_column%T_at_interface(nk,2))
  allocate(right_column%S_at_interface(nk,2))
  allocate(right_column%P_at_interface(nk,2))
  allocate(right_column%T_poly(nk,2))
  allocate(right_column%stable_cell(nk))
 
  ! Unit tests for find_neutral_surface_positions_discontinuous
  ! Salinity is 0 for all these tests
  allocate(CS%EOS)
  call EOS_manual_init(CS%EOS, form_of_EOS=EOS_LINEAR, dRho_dT=-1., dRho_dS=0.)
  Sl(:) = 0. ; Sr(:) = 0. ; ; SiL(:,:) = 0. ; SiR(:,:) = 0.
  ppoly_T_l(:,:) = 0.; ppoly_T_r(:,:) = 0.
  ppoly_S_l(:,:) = 0.; ppoly_S_r(:,:) = 0.
  ! Intialize any control structures needed for unit tests
  CS%ref_pres = -1.
 
  hL = (/10.,10.,10./) ; hR = (/10.,10.,10./)
 
  left_column%P_at_interface(1,1)  = 0.; left_column%P_at_interface(1,2)  = hL(1)
  right_column%P_at_interface(1,1) = 0.; right_column%P_at_interface(1,2) = hR(1)
  do k = 2,nk
    left_column%P_at_interface(k,1) = left_column%P_at_interface(k-1,2)
    left_column%P_at_interface(k,2) = left_column%P_at_interface(k,1) + hL(k)
    right_column%P_at_interface(k,1) = right_column%P_at_interface(k-1,2)
    right_column%P_at_interface(k,2) = right_column%P_at_interface(k,1) + hR(k)
  enddo
  CS%delta_rho_form = 'mid_pressure'
  CS%neutral_pos_method = 1
  
  left_column%T_at_interface(1,:)  = (/ 22.00, 18.00 /)
  left_column%T_at_interface(2,:)  = (/ 18.00, 14.00 /)
  left_column%T_at_interface(3,:)  = (/ 14.00, 10.00 /)
  right_column%T_at_interface(1,:) = (/ 22.00, 18.00 /)
  right_column%T_at_interface(2,:) = (/ 18.00, 14.00 /)
  right_column%T_at_interface(3,:) = (/ 14.00, 10.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 2, 2, 3, 3, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 1, 1, 2, 2, 3, 3, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 10.00, 0.00, 10.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Identical Columns')
 
  left_column%T_at_interface(1,:)  = (/ 22.00, 18.00 /)
  left_column%T_at_interface(2,:)  = (/ 18.00, 14.00 /)
  left_column%T_at_interface(3,:)  = (/ 14.00, 10.00 /)
  right_column%T_at_interface(1,:) = (/ 20.00, 16.00 /)
  right_column%T_at_interface(2,:) = (/ 16.00, 12.00 /)
  right_column%T_at_interface(3,:) = (/ 12.00, 8.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, -1, -1 /),  & ! KoL
    (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1 /),  & ! KoR
    (/ 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.00, 0.00 /),  & ! PoR
    (/ 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 0.00 /),  & ! hEff
    'Right slightly cooler')
 
  left_column%T_at_interface(1,:)  = (/ 20.00, 16.00 /)
  left_column%T_at_interface(2,:)  = (/ 16.00, 12.00 /)
  left_column%T_at_interface(3,:)  = (/ 12.00, 8.00 /)
  right_column%T_at_interface(1,:) = (/ 22.00, 18.00 /)
  right_column%T_at_interface(2,:) = (/ 18.00, 14.00 /)
  right_column%T_at_interface(3,:) = (/ 14.00, 10.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, -1, -1 /),  & ! KoL
    (/ 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, -1, -1 /),  & ! KoR
    (/ 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.00, 0.00 /),  & ! PoL
    (/ 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.50, 0.50, 1.00, 0.00, 0.00 /),  & ! PoR
    (/ 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 5.00, 0.00, 0.00 /),  & ! hEff
    'Right slightly cooler')
 
  left_column%T_at_interface(1,:)  = (/ 20.00, 16.00 /)
  left_column%T_at_interface(2,:)  = (/ 15.00, 12.00 /)
  left_column%T_at_interface(3,:)  = (/ 12.00, 9.00 /)
  right_column%T_at_interface(1,:) = (/ 20.00, 19.00 /)
  right_column%T_at_interface(2,:) = (/ 19.00, 18.00 /)
  right_column%T_at_interface(3,:) = (/ 18.00, 17.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 1, 1, 2, 2, 3, 3, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 0.25, 0.25, 0.50, 0.50, 0.75, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 4.00, 0.00, 4.00, 0.00, 4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Entire right column to one layer')

  left_column%T_at_interface(1,:)  = (/ 22.00, 20.00 /)
  left_column%T_at_interface(2,:)  = (/ 18.00, 16.00 /)
  left_column%T_at_interface(3,:)  = (/ 14.00, 12.00 /)
  right_column%T_at_interface(1,:) = (/ 32.00, 24.00 /)
  right_column%T_at_interface(2,:) = (/ 22.00, 14.00 /)
  right_column%T_at_interface(3,:) = (/ 12.00, 4.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 2, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 0.25, 0.50, 0.75, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 4.00, 0.00, 4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Right more strongly stratified')

  left_column%T_at_interface(1,:)  = (/ 22.00, 18.00 /)
  left_column%T_at_interface(2,:)  = (/ 18.00, 14.00 /)
  left_column%T_at_interface(3,:)  = (/ 14.00, 10.00 /)
  right_column%T_at_interface(1,:) = (/ 14.00, 14.00 /)
  right_column%T_at_interface(2,:) = (/ 14.00, 14.00 /)
  right_column%T_at_interface(3,:) = (/ 12.00, 8.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.50, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 5.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Deep Mixed layer on the right')

  left_column%T_at_interface(1,:)  = (/ 14.00, 14.00 /)
  left_column%T_at_interface(2,:)  = (/ 14.00, 12.00 /)
  left_column%T_at_interface(3,:)  = (/ 10.00, 8.00 /)
  right_column%T_at_interface(1,:) = (/ 14.00, 14.00 /)
  right_column%T_at_interface(2,:) = (/ 14.00, 14.00 /)
  right_column%T_at_interface(3,:) = (/ 14.00, 14.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Right unstratified column')

  left_column%T_at_interface(1,:)  = (/ 14.00, 14.00 /)
  left_column%T_at_interface(2,:)  = (/ 14.00, 10.00 /)
  left_column%T_at_interface(3,:)  = (/ 10.00, 2.00 /)
  right_column%T_at_interface(1,:) = (/ 14.00, 14.00 /)
  right_column%T_at_interface(2,:) = (/ 14.00, 10.00 /)
  right_column%T_at_interface(3,:) = (/ 10.00, 2.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 2, 2, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 2, 2, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 10.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Identical columns with mixed layer')

  left_column%T_at_interface(1,:)  = (/ 14.00, 12.00 /)
  left_column%T_at_interface(2,:)  = (/ 10.00, 10.00 /)
  left_column%T_at_interface(3,:)  = (/ 8.00, 2.00 /)
  right_column%T_at_interface(1,:) = (/ 14.00, 12.00 /)
  right_column%T_at_interface(2,:) = (/ 12.00, 8.00 /)
  right_column%T_at_interface(3,:) = (/ 8.00, 2.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 1, 1, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 1, 1, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 10.00, 0.00, 10.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Left interior unstratified')

  left_column%T_at_interface(1,:)  = (/ 12.00, 12.00 /)
  left_column%T_at_interface(2,:)  = (/ 12.00, 10.00 /)
  left_column%T_at_interface(3,:)  = (/ 10.00, 6.00 /)
  right_column%T_at_interface(1,:) = (/ 12.00, 10.00 /)
  right_column%T_at_interface(2,:) = (/ 10.00, 12.00 /)
  right_column%T_at_interface(3,:) = (/ 8.00, 4.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 2, 2, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 1, 1, 3, 3, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.50, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.00, 1.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 10.00, 0.00, 5.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Left mixed layer, Right unstable interior')

  left_column%T_at_interface(1,:)  = (/ 14.00, 14.00 /)
  left_column%T_at_interface(2,:)  = (/ 10.00, 10.00 /)
  left_column%T_at_interface(3,:)  = (/ 8.00, 6.00 /)
  right_column%T_at_interface(1,:) = (/ 10.00, 14.00 /)
  right_column%T_at_interface(2,:) = (/ 16.00, 16.00 /)
  right_column%T_at_interface(3,:) = (/ 12.00, 4.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.50, 0.75, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 4.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Left thick mixed layer, Right unstable mixed')

  left_column%T_at_interface(1,:)  = (/ 8.00, 12.00 /)
  left_column%T_at_interface(2,:)  = (/ 12.00, 10.00 /)
  left_column%T_at_interface(3,:)  = (/ 8.00, 4.00 /)
  right_column%T_at_interface(1,:) = (/ 10.00, 14.00 /)
  right_column%T_at_interface(2,:) = (/ 14.00, 12.00 /)
  right_column%T_at_interface(3,:) = (/ 10.00, 6.00 /)
  call column_mark_unstable( CS, nk, left_column )
  call column_mark_unstable( CS, nk, right_column )
  call find_neutral_surface_positions_discontinuous(CS, nk, left_column, right_column, hL, hR, PoL, PoR, KoL, KoR, hEff)
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or.  test_nsp(v, 12, KoL, KoR, PoL, PoR, hEff, &
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoL
    (/ 3, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 /),  & ! KoR
    (/ 0.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoL
    (/ 0.50, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! PoR
    (/ 5.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /),  & ! hEff
    'Unstable mixed layers, left cooler')


  call EOS_manual_init(CS%EOS, form_of_EOS = EOS_LINEAR, dRho_dT = -1., dRho_dS = 2.)
  ! Tests for linearized version of searching the layer for neutral surface position
  ! EOS linear in T, uniform alpha
  CS%max_iter = 10
  ! Unit tests require explicit initialization of tolerance
  CS%Drho_tol = 0.
  CS%x_tol = 0.
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 10., 35., -0.2, 0., &
                                     -0.2, 0., -0.2, 0.,                     &
                                     (/12.,-4./), (/34.,0./)), "Temp Uniform Linearized Alpha/Beta"))
  ! EOS linear in S, uniform beta
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 10., 35., 0., 0.8, &
                                     0., 0.8, 0., 0.8,                &
                                    (/12.,0./), (/34.,2./)), "Salt Uniform Linearized Alpha/Beta"))
  ! EOS linear in T/S, uniform alpha/beta
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5,   &
             find_neutral_pos_linear(CS, 10., 35., -0.5, 0.5,                &
                                     -0.5, 0.5, -0.5, 0.5,  &
                                     (/12.,-4./), (/34.,2./)), "Temp/salt Uniform Linearized Alpha/Beta"))
  ! EOS linear in T, insensitive to So
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 10., 35., -0.2, 0., &
                                     -0.4, 0., -0.6, 0.,  &
                                     (/12.,-4./), (/34.,0./)), "Temp stratified Linearized Alpha/Beta"))
  ! EOS linear in S, insensitive to T
  ndiff_unit_tests_discontinuous = ndiff_unit_tests_discontinuous .or. (test_rnp(0.5, &
             find_neutral_pos_linear(CS, 10., 35., 0., 0.8,  &
                                      0., 1.0,  0., 0.5,  &
                                     (/12.,0./), (/34.,2./)), "Salt stratified Linearized Alpha/Beta"))
  if (.not. ndiff_unit_tests_discontinuous) write(stdout,*) 'Pass'

end function ndiff_unit_tests_discontinuous

!> Returns true if a test of fv_diff() fails, and conditionally writes results to stream
logical function test_fv_diff(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width [nondim]
  real,             intent(in) :: hk    !< Center cell width [nondim]
  real,             intent(in) :: hkp1  !< Right cell width [nondim]
  real,             intent(in) :: Skm1  !< Left cell average value
  real,             intent(in) :: Sk    !< Center cell average value
  real,             intent(in) :: Skp1  !< Right cell average value
  real,             intent(in) :: Ptrue !< True answer [nondim]
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fv_diff(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fv_diff = (Pret /= Ptrue)

  if (test_fv_diff .or. verbose) then
    stdunit = stdout
    if (test_fv_diff) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_fv_diff) then
      write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_fv_diff

!> Returns true if a test of fvlsq_slope() fails, and conditionally writes results to stream
logical function test_fvlsq_slope(verbose, hkm1, hk, hkp1, Skm1, Sk, Skp1, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: hkm1  !< Left cell width
  real,             intent(in) :: hk    !< Center cell width
  real,             intent(in) :: hkp1  !< Right cell width
  real,             intent(in) :: Skm1  !< Left cell average value
  real,             intent(in) :: Sk    !< Center cell average value
  real,             intent(in) :: Skp1  !< Right cell average value
  real,             intent(in) :: Ptrue !< True answer
  character(len=*), intent(in) :: title !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = fvlsq_slope(hkm1, hk, hkp1, Skm1, Sk, Skp1)
  test_fvlsq_slope = (Pret /= Ptrue)

  if (test_fvlsq_slope .or. verbose) then
    stdunit = stdout
    if (test_fvlsq_slope) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_fvlsq_slope) then
      write(stdunit,'(2(x,a,f20.16),x,a)') 'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(2(x,a,f20.16))') 'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_fvlsq_slope

!> Returns true if a test of interpolate_for_nondim_position() fails, and conditionally writes results to stream
logical function test_ifndp(verbose, rhoNeg, Pneg, rhoPos, Ppos, Ptrue, title)
  logical,          intent(in) :: verbose !< If true, write results to stdout
  real,             intent(in) :: rhoNeg !< Lighter density [R ~> kg m-3]
  real,             intent(in) :: Pneg   !< Interface position of lighter density [nondim]
  real,             intent(in) :: rhoPos !< Heavier density [R ~> kg m-3]
  real,             intent(in) :: Ppos   !< Interface position of heavier density [nondim]
  real,             intent(in) :: Ptrue  !< True answer [nondim]
  character(len=*), intent(in) :: title  !< Title for messages

  ! Local variables
  integer :: stdunit
  real :: Pret

  Pret = interpolate_for_nondim_position(rhoNeg, Pneg, rhoPos, Ppos)
  test_ifndp = (Pret /= Ptrue)

  if (test_ifndp .or. verbose) then
    stdunit = stdout
    if (test_ifndp) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    if (test_ifndp) then
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15),x,a)') &
            'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue,'WRONG!'
    else
      write(stdunit,'(4(x,a,f20.16),2(x,a,1pe22.15))') &
            'r1=',rhoNeg,'p1=',Pneg,'r2=',rhoPos,'p2=',Ppos,'pRet=',Pret,'pTrue=',Ptrue
    endif
  endif

end function test_ifndp

!> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
logical function test_data1d(verbose, nk, Po, Ptrue, title)
  logical,             intent(in) :: verbose !< If true, write results to stdout
  integer,             intent(in) :: nk    !< Number of layers
  real, dimension(nk), intent(in) :: Po    !< Calculated answer
  real, dimension(nk), intent(in) :: Ptrue !< True answer
  character(len=*),    intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit

  test_data1d = .false.
  do k = 1,nk
    if (Po(k) /= Ptrue(k)) test_data1d = .true.
  enddo

  if (test_data1d .or. verbose) then
    stdunit = stdout
    if (test_data1d) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) then
        test_data1d = .true.
        write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15,x,a)') &
              'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k),'WRONG!'
      else
        if (verbose) &
          write(stdunit,'(a,i2,2(x,a,f20.16),x,a,1pe22.15)') &
                'k=',k,'Po=',Po(k),'Ptrue=',Ptrue(k),'err=',Po(k)-Ptrue(k)
      endif
    enddo
  endif

end function test_data1d

!> Returns true if comparison of Po and Ptrue fails, and conditionally writes results to stream
logical function test_data1di(verbose, nk, Po, Ptrue, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: nk    !< Number of layers
  integer, dimension(nk), intent(in) :: Po    !< Calculated answer
  integer, dimension(nk), intent(in) :: Ptrue !< True answer
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit

  test_data1di = .false.
  do k = 1,nk
    if (Po(k) /= Ptrue(k)) test_data1di = .true.
  enddo

  if (test_data1di .or. verbose) then
    stdunit = stdout
    if (test_data1di) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,nk
      if (Po(k) /= Ptrue(k)) then
        test_data1di = .true.
        write(stdunit,'(a,i2,2(x,a,i5),x,a)') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k),'WRONG!'
      else
        if (verbose) &
          write(stdunit,'(a,i2,2(x,a,i5))') 'k=',k,'Io=',Po(k),'Itrue=',Ptrue(k)
      endif
    enddo
  endif

end function test_data1di

!> Returns true if output of find_neutral_surface_positions() does not match correct values,
!! and conditionally writes results to stream
logical function test_nsp(verbose, ns, KoL, KoR, pL, pR, hEff, KoL0, KoR0, pL0, pR0, hEff0, title)
  logical,                intent(in) :: verbose !< If true, write results to stdout
  integer,                intent(in) :: ns    !< Number of surfaces
  integer, dimension(ns), intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer, dimension(ns), intent(in) :: KoR   !< Index of first right interface above neutral surface
  real, dimension(ns),    intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real, dimension(ns),    intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  real, dimension(ns-1),  intent(in) :: hEff  !< Effective thickness between two neutral surfaces [R L2 T-2 ~> Pa]
  integer, dimension(ns), intent(in) :: KoL0  !< Correct value for KoL
  integer, dimension(ns), intent(in) :: KoR0  !< Correct value for KoR
  real, dimension(ns),    intent(in) :: pL0   !< Correct value for pL
  real, dimension(ns),    intent(in) :: pR0   !< Correct value for pR
  real, dimension(ns-1),  intent(in) :: hEff0 !< Correct value for hEff
  character(len=*),       intent(in) :: title !< Title for messages

  ! Local variables
  integer :: k, stdunit
  logical :: this_row_failed

  test_nsp = .false.
  do k = 1,ns
    test_nsp = test_nsp .or. compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
    if (k < ns) then
      if (hEff(k) /= hEff0(k)) test_nsp = .true.
    endif
  enddo

  if (test_nsp .or. verbose) then
    stdunit = stdout
    if (test_nsp) stdunit = stderr ! In case of wrong results, write to error stream
    write(stdunit,'(a)') title
    do k = 1,ns
      this_row_failed = compare_nsp_row(KoL(k), KoR(k), pL(k), pR(k), KoL0(k), KoR0(k), pL0(k), pR0(k))
      if (this_row_failed) then
        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k),' <-- WRONG!'
        write(stdunit,10) k,KoL0(k),pL0(k),KoR0(k),pR0(k),' <-- should be this'
      else
        write(stdunit,10) k,KoL(k),pL(k),KoR(k),pR(k)
      endif
      if (k < ns) then
        if (hEff(k) /= hEff0(k)) then
          write(stdunit,'(i3,8x,"layer hEff =",2(f20.16,a))') k,hEff(k)," .neq. ",hEff0(k),' <-- WRONG!'
        else
          write(stdunit,'(i3,8x,"layer hEff =",f20.16)') k,hEff(k)
        endif
      endif
    enddo
  endif
  if (test_nsp) call MOM_error(FATAL,"test_nsp failed")

10 format("ks=",i3," kL=",i3," pL=",f20.16," kR=",i3," pR=",f20.16,a)
end function test_nsp

!> Compares a single row, k, of output from find_neutral_surface_positions()
logical function compare_nsp_row(KoL, KoR, pL, pR, KoL0, KoR0, pL0, pR0)
  integer,  intent(in) :: KoL   !< Index of first left interface above neutral surface
  integer,  intent(in) :: KoR   !< Index of first right interface above neutral surface
  real,     intent(in) :: pL    !< Fractional position of neutral surface within layer KoL of left column
  real,     intent(in) :: pR    !< Fractional position of neutral surface within layer KoR of right column
  integer,  intent(in) :: KoL0  !< Correct value for KoL
  integer,  intent(in) :: KoR0  !< Correct value for KoR
  real,     intent(in) :: pL0   !< Correct value for pL
  real,     intent(in) :: pR0   !< Correct value for pR

  compare_nsp_row = .false.
  if (KoL /= KoL0) compare_nsp_row = .true.
  if (KoR /= KoR0) compare_nsp_row = .true.
  if (pL /= pL0) compare_nsp_row = .true.
  if (pR /= pR0) compare_nsp_row = .true.
end function compare_nsp_row

!> Compares output position from refine_nondim_position with an expected value
logical function test_rnp(expected_pos, test_pos, title)
  real,             intent(in) :: expected_pos !< The expected position
  real,             intent(in) :: test_pos !< The position returned by the code
  character(len=*), intent(in) :: title    !< A label for this test
  ! Local variables
  integer :: stdunit

  stdunit = stdout ! Output to standard error
  test_rnp = ABS(expected_pos - test_pos) > 2*EPSILON(test_pos)
  if (test_rnp) then
    write(stdunit,'(A, f20.16, " .neq. ", f20.16, " <-- WRONG")') title, expected_pos, test_pos
  else
    write(stdunit,'(A, f20.16, " ==  ", f20.16)') title, expected_pos, test_pos
  endif
end function test_rnp
!> Deallocates neutral_diffusion control structure
subroutine neutral_diffusion_end(CS)
  type(neutral_diffusion_CS), pointer :: CS  !< Neutral diffusion control structure

  if (associated(CS)) deallocate(CS)

end subroutine neutral_diffusion_end

end module MOM_neutral_diffusion
