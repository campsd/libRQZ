! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Single shift RQZ for complex Hessenberg pairs
! ___________________________________________________________________
module z_RQZ

use u_parameters
use u_activeparts
use z_setpoles
use z_swappoles
use z_computepoles
use z_deflations

implicit none
private
public z_RQZ1

contains
  subroutine z_RQZ1(A, B, compAB, strtidx, stpidx, Q, compQ, qcol, &
                    Z, compZ, zcol, alpha, beta)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Implements the single shift complex-valued RQZ method
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       complex double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with A
  !           if compAB is true.
  !
  ! B       complex double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with B
  !           if compAB is true.
  !
  ! compAB  boolean       [IN]
  !           If true, the equivalences are applied to (A,B) such
  !           that the output triangular pair is equivalent to
  !           the input pair. If not, the output is not
  !           equivalent but the eigenvalues of the original pair
  !           are available on the diagonal.
  !
  ! strtidx integer     [IN]
  !           Start index of part of the matrix pencil that needs
  !           to be reduced to Schur form. Provide '0' if start is
  !           the beginning of the pencil.
  !
  ! stpidx  integer     [IN]
  !           Stop index of the part of the matrix pencil that needs
  !           to be reduced to Schur form. Provide N (matrix size) if
  !           stop is the end of the pencil. The strt and stp parameters
  !           need to satisfy:
  !               0 <= strtidx < stpidx <= N
  !
  ! Q       complex double array [INOUT]
  !           Unitary matrix containing the left Schur
  !           vectors. The number of columns of Q should be at least
  !           stpidx-strtidx. The Q array is not initialized in this routine.
  !
  ! compQ   boolean       [IN]
  !           If true, the left Schur vectors are computed.
  !
  ! qcol    integer       [IN]
  !           Column index from whereon the equivalence transformations
  !           are applied. If they should be applied from the beginning
  !           of the Q array, qcol needs to be set to 1. If qcol > 1, then
  !           the Q array should have at least qcol+stpidx-strtidx-1 columns.
  !
  ! Z       complex double array [INOUT]
  !           Unitary matrix containing the right Schur
  !           vectors. The number of columns of Z should be at least
  !           stp-strt. The Z array is not initialized in this routine.
  !
  ! compZ   boolean       [IN]
  !           If true, the right Schur vectors are computed.
  !
  ! zcol    integer       [IN]
  !           Column index from whereon the equivalence transformations
  !           are applied. If they should be applied from the beginning
  !           of the Z array, zcol needs to be set to 1. If zcol > 1, then
  !           the Z array should have at least zcol+stpidx-strtidx-1 columns.
  !
  ! alpha   complex double array [OUT]
  !           Array of length stpidx-strtidx containing the numerators for
  !           the generalized eigenvalues
  !
  ! beta    complex double array [OUT]
  !           Array of length stpidx-strtidx containing the denominators
  !           for the generalized eigenvalues
  !
  ! On successful completion of the routine the following equivalences hold:
  !   A(out) = Q^* A(in) Z
  !   B(out) = Q^* B(in) Z
  !___________________________________________________________________________
  ! last edit: July 31, 2019

  complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
  logical, intent(in)               ::  compAB, compQ, compZ
  integer, intent(in)               ::  strtidx, stpidx, qcol, zcol
  complex(kind=dp), intent(out)     ::  alpha(:), beta(:)

  ! Internal variables
  type(tAp), pointer                :: apHead, apCurrent, apTotal
  integer                           :: i, N, nit, nswaps, ld, qci, zci
  integer,pointer                   :: cstrt, cstp
  logical                           :: problemsolved, defl
  complex(kind=dp)                  :: mu, nu


  ! Initialize variables
  N = stpidx-strtidx
  ! apHead stores the active parts throughout the iteration
  call tApInit(apHead,strt=strtidx,stp=stpidx)
  ! apTotal stores the active part at the start and is never modified
  call tApInit(apTotal,strt=strtidx,stp=stpidx)
  ! Q and Z indices
  qci = qcol-strtidx
  zci = zcol-strtidx

  ! Performance tracking
  nit = 0
  nswaps = 0
  ld = 0
  problemsolved = .false.

  ! Start of main program
  do while (.not. problemsolved)
  ! Main loop is repeated until no active parts are left
  apCurrent => apHead
  do
  ! Loop over all active parts
    call tApGet(apCurrent,cstrt,cstp)
    if (cstp - cstrt .le. 1) then
      ! Not an active part anymore, remove from list
      if (.not. (associated(tApPrevious(apCurrent)))) then
        ! apCurrent is the head, update head
        apHead => tApNext(apCurrent) ! can be null()
      end if
      ! Delete and restart
      call tApDelete(apCurrent)
      apCurrent => apHead
    else
      ! Handle the current active part
      ! check for deflations
      ! at front:
      call z_check_start_deflation(defl,A,B,compAB,apCurrent,apTotal,&
                                   Q,compQ,qci+cstrt)
      ! apCurrent is updated; consequently cstrt, and cstp are !
      if (defl) then
        ld = 0
      end if
      ! along subdiagonal:
      call z_check_interior_deflations(defl,A, B, apCurrent)
      ! apCurrent is updated; consequently cstrt, and cstp
      if (defl) then
        ! Make sure that we're handling the correct range
        call tApGet(apCurrent,cstrt,cstp)
        ld = 0
      end if

      if (cstp - cstrt .gt. 1) then
        nit = nit + 1
        ld = ld + 1
        ! set the shift
        if ((ld .gt. 0) .and. (mod(ld,10) .eq. 0)) then
          ! exceptional shift
          mu = dconjg( A(cstp-1,cstp) )
          nu = dconjg( B(cstp-1,cstp) )
        else
          call z_get_shift(A, B, cstp, mu, nu)
        end if
        call z_set_first_pole(mu, nu, A, B, compAB, apCurrent, apTotal,&
                              Q, compQ, qci+cstrt)
        ! chase
        do i = cstrt+1, cstp-2
          ! swap
          call z_swap(i+1,i, A, B, compAB, apCurrent, apTotal, &
                      Q, compQ, qci+i, Z, compZ, zci+i-1)
          nswaps = nswaps + 1
        end do
        ! check for deflations at end
        call z_check_stop_deflation(defl, A, B, compAB,apCurrent, apTotal,&
                                    Z, compZ, zci+cstp-2)
        ! apCurrent is updated; consequently cstrt, and cstp
        if (defl) then
          ld = 0
        end if
        if (cstp - cstrt .gt. 1) then ! Check again
          ! set the pole
          call z_get_pole(A, B, cstrt+1, mu, nu)
          call z_set_last_pole(mu, nu, A, B, compAB, apCurrent,apTotal,&
                               Z, compZ, zci+cstp-2)
        end if
      end if
    end if

    ! Check for exit condition
    if (associated(tApNext(apCurrent))) then
      apCurrent => tApNext(apCurrent) ! Set the next one
    else
      exit ! Exit this loop over the active parts
    end if
  ! end of loop over active parts
  end do

  ! Check if entire problem is solved
  call tApGet(apHead,cstrt,cstp)
  if (.not. (associated(cstrt) .or. associated(cstp))) then
    problemsolved = .true.
  end if

  ! end of main loop
  end do

  ! return results
  do i = 1, N
    alpha(i) = A(strtidx+i, strtidx+i)
    beta(i) = B(strtidx+i,strtidx+i)
  end do

  ! free everything
  apCurrent => null()
  call tApFree(apHead)
  call tApFree(apTotal)

  end subroutine

end module z_RQZ
