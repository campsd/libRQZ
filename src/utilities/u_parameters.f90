! libRQZ
! daan. camps@cs.kuleuven.be
! Description:
!   Commonly used parameters
! ___________________________________________________________________
module u_parameters
  implicit none

  integer, parameter  ::  dp = SELECTED_REAL_KIND(15)
  real(kind=dp), parameter  ::  mp = epsilon(1.0_dp)
  real(kind=dp)  ::  sfmin

  ! For complex immplementation (z)
  complex(kind=dp), parameter :: czero = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter :: cone = (1.0_dp,0.0_dp)
  integer, parameter          :: ssth = 80

  ! For double implementation (d)
  real(kind=dp), parameter :: dzero = 0.0_dp
  real(kind=dp), parameter :: done = 1.0_dp

end module
