! libRQZ
! author: daan.camps@cs.kuleuven.be
! Description:
!   Implements the routines for selecting the shifts and poles during
!   the iteration.
! ___________________________________________________________________
module z_computepoles
 use u_parameters

 implicit none
 private
 public z_get_shift, z_get_pole

contains

 subroutine z_get_shift(A, B, idx, alpha, beta)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes the eigenvalues of the trailing two by two block
  ! and returns the one closest to the last diagonal element
  ! based on the input.
  ! The user can overwrite this function to use a different shift
  ! selection procedure.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       complex double array [IN]
  !           Upper Hessenberg matrix A
  ! B       complex double array [IN]
  !           Upper Hessenberg matrix B
  ! idx     integer [IN]
  !           Index indicating the two-by-two block,
  !             idx-1:idx,idx-1:idx
  ! alpha   complex double [OUT]
  !           Defines the eigenvalue of the
  !           block closest to element idx,idx
  ! beta    complex double [OUT]
  !           Defines the eigenvalue of the
  !           block closest to element idx,idx
  !___________________________________________________________________________
  ! last edit: June 7, 2018
  complex(kind=dp),intent(in) :: A(:,:), B(:,:)
  integer, intent(in) :: idx
  complex(kind=dp),intent(out) :: alpha, beta

  call z_get_eigenvalues2x2(A(idx-1:idx,idx-1:idx),B(idx-1:idx,idx-1:idx),&
                            alpha,beta)

  if (abs(B(idx,idx)*alpha-A(idx,idx)) .lt. &
      abs(B(idx,idx)*beta-A(idx,idx))) then
   beta = cone
  else
   alpha = beta
   beta = cone
  endif

 end subroutine

 subroutine z_get_pole(A, B, idx, alpha, beta)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes the eigenvalues of the leading two by two block
  ! and returns the one closest to the first diagonal element
  ! based on the input.
  ! The user can overwrite this function to use a different shift
  ! selection procedure.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       complex double array [IN]
  !           Upper Hessenberg matrix A
  ! B       complex double array [IN]
  !           Upper Hessenberg matrix B
  ! idx     integer [IN]
  !           Index indicating the two-by-two block,
  !             idx:idx+1,idx:idx+1
  ! alpha   complex double [OUT]
  !           Defines the eigenvalue of the
  !           block closest to element idx,idx
  ! beta    complex double [OUT]
  !           Defines the eigenvalue of the
  !           block closest to element idx,idx
  !___________________________________________________________________________
  ! last edit: June 7, 2018
  complex(kind=dp),intent(in) :: A(:,:), B(:,:)
  integer, intent(in) :: idx
  complex(kind=dp),intent(out) :: alpha, beta

  call z_get_eigenvalues2x2(A(idx:idx+1,idx:idx+1),B(idx:idx+1,idx:idx+1),&
                            alpha,beta)

  if (abs(B(idx,idx)*alpha-A(idx,idx)) .lt. &
      abs(B(idx,idx)*beta-A(idx,idx))) then
   beta = cone
  else
   alpha = beta
   beta = cone
  endif
 end subroutine

 subroutine z_get_eigenvalues2x2(A,B,l1,l2)
 ! DESCRIPTION
 !___________________________________________________________________________
 ! Computes the eigenvalues of a 2-by-2 block
 !
 ! ARGUMENTS
 !___________________________________________________________________________
 ! A       complex double array [IN]
 !           2 x 2 complex array
 ! B       complex double array [IN]
 !           2 x 2 complex array
 ! l1      complex double [OUT]
 !          eigenvalue of (A,B)
 ! l2      complex double [OUT]
 !          eigenvalue of (A,B)
 !
 !___________________________________________________________________________
 ! last edit: June 7, 2018
  complex(kind=dp),intent(in) :: A(:,:), B(:,:)
  complex(kind=dp),intent(out) :: l1, l2

  complex(kind=dp) :: DET, AB11, AB12, AB21, AB22, T1, RTDISC

  DET = B(1,1)*B(2,2) - B(1,2)* B(2,1)

  if (abs(DET) .gt. mp) then
   AB11 = A(1,1)*B(2,2)-A(2,1)*B(1,2)
   AB12 = A(1,2)*B(2,2)-A(2,2)*B(1,2)
   AB21 = A(2,1)*B(1,1)-A(1,1)*B(2,1)
   AB22 = A(2,2)*B(1,1)-A(1,2)*B(2,1)

   T1 =  (1.0_dp/2.0_dp)*(AB11+AB22)
   RTDISC = sqrt(T1**2 + AB12*AB21-AB11*AB22)

   l1 = (T1 + RTDISC)/DET
   l2 = (T1 - RTDISC)/DET
  else
    write(*,*) 'SINGULAR B - TODO'
    write(*,*) A,B
    stop
  end if
 end subroutine
end module
