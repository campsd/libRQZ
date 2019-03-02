! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Computes elementary rotations
! ___________________________________________________________________
module d_ctransformations

use u_parameters
implicit none
private
public d_compute_ct, d_apply_ct_r, d_apply_ct_l

contains

  subroutine d_compute_ct(x,y,c,s,r)
  ! Computes an elementary ctransformation that introduces a zero
  ! in a real 2x1 vector :
  ! [ c  s ]    [ x  ]        [ r   ]
  ! [      ] *  [    ]    =   [     ]
  ! [-s  c ]    [ y  ]        [ 0   ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! x         double [IN]
  !             First element of 2x1 vector
  ! y         double [IN]
  !             Second element of 2x1 vector
  ! c         double [OUT]
  !             Cosine of the rotation matrix
  ! s         double [OUT]
  !             Sine of the rotation matrix
  ! r         double [OUT]
  !             Remainder of the vector after rotating (2-norm)
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp), intent(in)  :: x,y
    real(kind=dp), intent(out) :: c,s,r

    ! real(kind=dp)              :: t ! Only required in own implementation

    call DLARTG( x, y, c, s, r )

    ! if (abs(y) .lt. sfmin) then
    !   c = done
    !   s = dzero
    !   r = x
    ! else
    !   if ( abs(x) .gt. abs(y) ) then
    !     t = y/x
    !     r = sqrt(done + abs(t)**2)
    !     c = done / r
    !     if ( x .lt. dzero ) then
    !       c = -c
    !     end if
    !     s = t * c
    !     r = abs(x) * r
    !   else
    !     t = x/y
    !     r = sqrt(done + abs(t)**2)
    !     s = done / r
    !     if (y .lt. dzero ) then
    !       s = -s
    !     end if
    !     c = t * s
    !     r = abs(y) * r
    !   end if
    ! end if
  end subroutine

  subroutine d_apply_ct_l(v1,v2,c,s)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Applies an elementary rotation from the left to two row
  ! vectors v1 and v2 of equal lenght :
  ! [ c  s ]    [ v1  ]
  ! [      ] *  [     ]
  ! [-s  c ]    [ v2  ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! v1        double array [INOUT]
  !             First row vector, on output equal to c * v1 + s * v2
  ! v2        double array [INOUT]
  !             Second row vector, on output equal to -s * v1 + c * v2
  ! c         double [IN]
  !             Cosine of rotation
  ! s         double [IN]
  !             Sine of rotation
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp), intent(inout)  :: v1(:),v2(:)
    real(kind=dp), intent(in)     :: c, s

    integer                       :: i, n
    real(kind=dp)                 :: t

    n = size(v1)
    do i=1, n
      t = c * v1(i) + s * v2(i)
      v2(i) = -s * v1(i) + c * v2(i)
      v1(i) = t
    end do
  end subroutine

  subroutine d_apply_ct_r(v1,v2,c,s)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Applies an elementary rotation from the left to two column
  ! vectors v1 and v2 of equal lenght :
  ! [     ]   [ c  s ]
  ! [v1 v2] * [      ]
  ! [     ]   [-s  c ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! v1        double array [INOUT]
  !             First column vector, on output equal to c * v1 - s * v2
  ! v2        double array [INOUT]
  !             Second row vector, on output equal to s * v1 + c * v2
  ! c         double [IN]
  !             Cosine of rotation
  ! s         double [IN]
  !             Sine of rotation
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp), intent(inout)  :: v1(:),v2(:)
    real(kind=dp), intent(in)     :: c, s

    integer                       :: i, n
    real(kind=dp)                 :: t

    n = size(v1)
    do i=1, n
      t = c * v1(i) - s * v2(i)
      v2(i) = s * v1(i) + c * v2(i)
      v1(i) = t
    end do
  end subroutine
end module
