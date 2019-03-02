! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Computes elementary rotations
! ___________________________________________________________________
module z_ctransformations

use u_parameters
implicit none
private
public z_compute_ct, z_apply_ct_r, z_apply_ct_l, &
       z_compute_ct_complexc, z_apply_ct_r_complexc, z_apply_ct_l_complexc

contains

  subroutine z_compute_ct(x,y,c,s,r)
    ! Computes an elementary ctransformation that introduces a zero
    ! in a 2x1 vector :
    ! [ c  s ]    [ x  ]        [ r   ]
    ! [ _  _ ] *  [    ]    =   [     ]
    ! [-s  c ]    [ y  ]        [ 0   ]
    ! Wrapper around LAPACKs ZLARTG
    ! Our own implementation uses a complex variable for c. This is
    ! more efficient.
    !___________________________________________________________________________
    complex(kind=dp), intent(in)  :: x, y
    complex(kind=dp), intent(out) :: s, r
    real(kind=dp), intent(out)    :: c

    call ZLARTG( x, y, c, s, r )
  end subroutine

  subroutine z_apply_ct_l(v1,v2,c,s)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Applies an elementary rotation from the left to two row
  ! vectors v1 and v2 of equal lenght :
  ! [ c  s ]    [ v1  ]
  ! [ _    ] *  [     ]
  ! [-s  c ]    [ v2  ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! v1        complex double array [INOUT]
  !             First row vector, on output equal to c * v1 + s * v2
  ! v2        complex double array [INOUT]             _        _
  !             Second row vector, on output equal to -s * v1 + c * v2
  ! c         double [IN]
  !             Cosine of rotation
  ! s         complex double [IN]
  !             Sine of rotation
  !___________________________________________________________________________
  ! last edit: June 7, 2018
    complex(kind=dp), intent(inout)   ::  v1(:), v2(:)
    complex(kind=dp), intent(in)      ::  s
    real(kind=dp), intent(in)         ::  c

    integer i, n
    complex(kind=dp) t

    n = size(v1)

    do i = 1, n
        t = c * v1(i) + s * v2(i)
        !v2(i) = -dconjg(s) * v1(i) + dconjg(c) * v2(i) ! For own implementation
        v2(i) = -dconjg(s) * v1(i) + c * v2(i)
        v1(i) = t
    end do

  end subroutine

  subroutine z_apply_ct_r(v1,v2,c,s)
    ! DESCRIPTION
    !___________________________________________________________________________
    ! Applies an elementary rotation from the right to
    ! two column vectors v1 and v2 of equal length:
    ! [     ]   [ c  s ]
    ! [v1 v2] * [ _  _ ]
    ! [     ]   [-s  c ]
    !
    ! ARGUMENTS
    !___________________________________________________________________________
    ! v1        complex double array [INOUT]                       _
    !             First column vector, on output equal to c * v1 - s * v2
    ! v2        complex double array [INOUT]                     _
    !             Second column vector, on output equal to s * v1 + c * v2
    ! c         double [IN]
    !             Cosine of rotation
    ! s         complex double [IN]
    !             Sine of rotation
    !___________________________________________________________________________
    ! last edit: June 7, 2018
    complex(kind=dp), intent(inout) ::  v1(:), v2(:)
    complex(kind=dp), intent(in)    ::  s
    real(kind=dp), intent(in)    ::  c

    integer i, n
    complex(kind=dp)  t

    n = size(v1)

    do i = 1, n
        t = c * v1(i) - dconjg(s) * v2(i)
        !v2(i) = s * v1(i) + dconjg(c) * v2(i) ! For own implementation
        v2(i) = s * v1(i) + c * v2(i)
        v1(i) = t
    end do
  end subroutine

  subroutine z_compute_ct_complexc(x,y,c,s,r)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes an elementary ctransformation that introduces a zero
  ! in a 2x1 vector :
  ! [ c  s ]    [ x  ]        [ r   ]
  ! [ _  _ ] *  [    ]    =   [     ]
  ! [-s  c ]    [ y  ]        [ 0   ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! x         complex double [IN]
  !             First element of 2x1 vector
  ! y         complex double [IN]
  !             Second element of 2x1 vector
  ! c         complex double [OUT]
  !             Cosine of the rotation matrix
  ! s         complex double [OUT]
  !             Sine of the rotation matrix
  ! r         complex double [OUT]
  !             Remainder of the vector after rotating (2-norm)
  !
  !___________________________________________________________________________
  ! last edit: June 7, 2018
    complex(kind=dp), intent(in) :: x,y
    complex(kind=dp), intent(out) :: c,s,r

    complex(kind=dp)  :: theta, t

    if (abs(y) .lt. sfmin) then 
      c = cone
      s = czero
      r = x
    else
        if( abs(x) .gt. abs(y) ) then
          theta = dconjg(x/abs(x))
          t = y/x
          r = sqrt(1.0_dp + abs(t)**2)
          c = theta / r
          s = dconjg(t) * c
          r = theta * x * r
        else
          theta = dconjg(y/abs(y))
          t = x/y
          r = sqrt(1.0_dp + abs(t)**2)
          s = theta / r
          c = dconjg(t) * s
          r = theta * y * r
        end if
    end if
  end subroutine

  subroutine z_apply_ct_l_complexc(v1,v2,c,s)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Applies an elementary rotation from the left to two row
  ! vectors v1 and v2 of equal lenght :
  ! [ c  s ]    [ v1  ]
  ! [ _  _ ] *  [     ]
  ! [-s  c ]    [ v2  ]
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! v1        complex double array [INOUT]
  !             First row vector, on output equal to c * v1 + s * v2
  ! v2        complex double array [INOUT]             _        _
  !             Second row vector, on output equal to -s * v1 + c * v2
  ! c         double [IN]
  !             Cosine of rotation
  ! s         complex double [IN]
  !             Sine of rotation
  !___________________________________________________________________________
  ! last edit: June 7, 2018
    complex(kind=dp), intent(inout)   ::  v1(:), v2(:)
    complex(kind=dp), intent(in)      ::  c, s

    integer i, n
    complex(kind=dp) t

    n = size(v1)

    do i = 1, n
        t = c * v1(i) + s * v2(i)
        v2(i) = -dconjg(s) * v1(i) + dconjg(c) * v2(i)
        v1(i) = t
    end do

  end subroutine

  subroutine z_apply_ct_r_complexc(v1,v2,c,s)
    ! DESCRIPTION
    !___________________________________________________________________________
    ! Applies an elementary rotation from the right to
    ! two column vectors v1 and v2 of equal length:
    ! [     ]   [ c  s ]
    ! [v1 v2] * [ _  _ ]
    ! [     ]   [-s  c ]
    !
    ! ARGUMENTS
    !___________________________________________________________________________
    ! v1        complex double array [INOUT]                       _
    !             First column vector, on output equal to c * v1 - s * v2
    ! v2        complex double array [INOUT]                     _
    !             Second column vector, on output equal to s * v1 + c * v2
    ! c         double [IN]
    !             Cosine of rotation
    ! s         complex double [IN]
    !             Sine of rotation
    !___________________________________________________________________________
    ! last edit: June 7, 2018
    complex(kind=dp), intent(inout) ::  v1(:), v2(:)
    complex(kind=dp), intent(in)    ::  c, s

    integer i, n
    complex(kind=dp)  t

    n = size(v1)

    do i = 1, n
        t = c * v1(i) - dconjg(s) * v2(i)
        v2(i) = s * v1(i) + dconjg(c) * v2(i)
        v1(i) = t
    end do
  end subroutine
end module
