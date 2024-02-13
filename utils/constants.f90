module constants
  implicit none

  ! Set precision
  ! Single precision real numbers, 6 digits, range 10⁻³⁷ to 10³⁷-1; 32 bits
  ! integer, parameter :: wp = selected_real_kind(6, 37)
  !> Double precision real numbers, 15 digits, range 10⁻³⁰⁷ to 10³⁰⁷-1; 64 bits
  integer, parameter :: wp = selected_real_kind(15, 307)
  ! Quadruple precision real numbers, 33 digits, range 10⁻⁴⁹³¹ to 10⁴⁹³¹-1; 128 bits
  ! integer, parameter :: wp = selected_real_kind(33, 4931)
  integer, parameter :: li = selected_int_kind(16)

   ! Physical constants
  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp), parameter :: kb = 3.166811563e-06_wp !> Boltzmann constant in au i.e. Hartree/K
  real(wp), parameter :: Da2au = 1822.88849_wp !> 1 amu to mass of electron (au)
  real(wp), parameter :: ang2au = 1.8897261246_wp !> 1 Ang to Bohr (au)
  real(wp), parameter :: au2cm_1 = 219474.63136320_wp !> 1 Ha in cm^-1
  real(wp), parameter :: fine_struct_alpha = 0.0072973525628_wp !> Fine structure constant (~1/137), dimensionless
  real(wp), parameter :: Ha2kcalpermol =  627.509469_wp !> 1 Ha in kcal/mol

end module constants
