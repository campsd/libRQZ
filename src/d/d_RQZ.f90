! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Single shift RQZ for real Hessenberg pairs
! ___________________________________________________________________
module d_RQZ

use u_parameters
use u_activeparts
use d_computepoles
use d_setpoles
use d_deflations
use d_swappoles12
use d_swappoles22
use d_memorymgmt

implicit none
private
public d_RQZ2, d_getEigenvalues

contains
  subroutine d_RQZ2(A, B, compAB, strtidx, stpidx, Q, compQ, qcol, &
                    Z, compZ, zcol, sa)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Implements the combination shift real-valued RQZ method
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [INOUT]
  !           Matrix in upper block-Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper blocktriangular form and is equivalent with A
  !           if compAB is true.
  !
  ! B       double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with B
  !           if compAB is true.
  !
  ! compAB  boolean       [IN]
  !           If true, the equivalences are applied to (A,B) such
  !           that the output (block)triangular pair is equivalent to
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
  ! Q       double array [INOUT]
  !           Orthonormal matrix containing the left Schur
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
  ! Z       double array [INOUT]
  !           Orthonormal matrix containing the right Schur
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
  ! sa   boolean [OPTIONAL - IN]
  !           Optional boolean, defaults to true. If true, the routine is
  !           assumed to be called 'standalone'. This means that:
  !             * auxilary arrays allocated in this subroutine are deallocated
  !           A standalone call of d_RQZ2 should omit this optional
  !           argument, or provide true. Otherwise memory can leak.
  !
  ! On successful completion of the routine the following equivalences hold:
  !   A(out) = Q^T A(in) Z
  !   B(out) = Q^T B(in) Z
  !___________________________________________________________________________
  ! last edit: January 15, 2019

  real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
  logical, intent(in)               ::  compAB, compQ, compZ
  integer, intent(in)               ::  strtidx, stpidx, qcol, zcol
  logical, intent(in), optional     ::  sa

  ! Internal variables
  type(tAp), pointer                :: apHead, apCurrent, apTotal
  integer                           :: i, N, nit, nswaps, ld, qci, zci
  integer,pointer                   :: cstrt, cstp
  logical                           :: problemsolved, defl, isreal, ds, sa2
  real(kind=dp)                     :: mu, nu

  external                          :: DGEMM
  real(kind=dp)                     :: DLAMCH

  ! Initialize variables
  sfmin = DLAMCH('S')
  N = stpidx-strtidx
  if (present(sa)) then
    sa2 = sa
  else
    sa2 = .true.
  end if

  ! apHead stores the active parts throughout the iteration
  call tApInit(apHead,strt=strtidx,stp=stpidx)
  ! apTotal stores the active part at the start and is never modified
  call tApInit(apTotal,strt=strtidx,stp=stpidx)
  ! Q and Z indices
  qci = qcol-strtidx
  zci = zcol-strtidx

  ! For larger swaps
  call d_allocateSwap4Mem
  call d_allocateRCMem(4,N)

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
    if (cstp - cstrt .le. 3) then
      if (cstp - cstrt .eq. 2) then
        ! 2 x 2 Block
        if (abs(A(cstrt+2,cstrt+1)) .gt. dzero) then
          ! Check if cc
          call d_eigenvalues2x2(A,B,cstrt+1,cstrt+1,mu,nu,isreal)
          if (isreal) then
            call d_make_hess(mu,cstrt+1,cstrt+1,A,B,compAB,apCurrent,apTotal,&
                   Q,compQ,qci+cstrt,Z,compZ,zci+cstrt)
          end if
        end if
      elseif (cstp - cstrt .eq. 3) then
        ! 3 x 3 Block
        ! Handle direct to avoid getting stuck
        if (abs(A(cstrt+3,cstrt+1)) .gt. dzero) then
          call d_first_double_to_inf(A,B,compAB,apCurrent,apTotal,&
                  Q,compQ,qci+cstrt)
        end if
        defl = .false.
        ld = 0
        do while (.not. defl)
          call d_start_deflation_single(defl,A,B,compAB,apCurrent,apTotal,&
                  Q,compQ,qci+cstrt)
          if (defl) exit

          if ((ld .gt. 0) .and. (mod(ld,10) .eq. 0)) then
            ! exceptional shift
            isreal = .true.
            mu = A(cstp-1,cstp)
            nu = B(cstp-1,cstp)
          else
            call d_get_shift(A, B, cstp, mu, nu, isreal)
          end if

          if (isreal) then
            call d_first_single_pole(mu,nu,A,B,compAB,apCurrent,apTotal,&
                     Q,compQ,qci+cstrt)
            call d_swap11(cstrt+2,cstrt+1,A,B,compAB,apCurrent,apTotal,&
                     Q,compQ,qci+cstrt+1,Z,compZ,zci+cstrt)
            call d_last_single_pole(done,dzero,A,B,compAB,apCurrent,apTotal,&
                     Z,compZ,zci+cstrt+1)
          else
            call d_first_double_pole(cmplx(mu,nu,kind=dp),A,B,compAB,&
                     apCurrent,apTotal,Q,compQ,qci+cstrt,Z,compZ,zci+cstrt)
            call d_last_double_to_inf(A,B,compAB,apCurrent,apTotal,&
                     Z,compZ,zci+cstrt)
          end if
          call d_stop_deflation_single(defl,A,B,compAB,apCurrent,apTotal,&
                  Z,compZ,zci+cstrt+1)
          if (defl) exit
          ld = ld + 1
        end do
      end if
      ! Not an active part anymore, remove from list
      if (.not. (associated(tApPrevious(apCurrent)))) then
        ! apCurrent is the head, update head
        apHead => tApNext(apCurrent) ! can be null()
      end if
      ! Delete and restart
      call tApDelete(apCurrent)
      apCurrent => apHead
    else
      ! Handle the current active part with an RQZ sweep
      ! Check for deflations
      ! at front:
      if (abs(A(cstrt+3,cstrt+1)) .gt. dzero) then
        call d_start_deflation_double(defl,A,B,compAB,apCurrent,apTotal,&
                                   Q,compQ,qci+cstrt,Z,compZ,zci+cstrt)
      else
        call d_start_deflation_single(defl,A,B,compAB,apCurrent,apTotal,&
                                   Q,compQ,qci+cstrt)
      end if
      ! apCurrent is updated; consequently cstrt, and cstp are !
      if (defl) then
        ! Make sure that we're handling the correct range
        call tApGet(apCurrent,cstrt,cstp)
        ld = 0
      end if
      ! along subdiagonal:
      call d_check_interior_deflations(defl,A, B, compAB, apCurrent, apTotal, &
             Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)
      ! apCurrent is updated; consequently cstrt, and cstp
      if (defl) then
        ! Make sure that we're handling the correct range
        call tApGet(apCurrent,cstrt,cstp)
        ld = 0
      end if

      if (cstp - cstrt .gt. 3) then
        ! Get the shift
        if ((ld .gt. 0) .and. (mod(ld,10) .eq. 0)) then
          ! exceptional shift
          mu = A(cstp-1,cstp)
          nu = B(cstp-1,cstp)
          isreal = .true.
        else
          call d_get_shift(A, B, cstp, mu, nu, isreal)
        endif
        ld = ld + 1

        if (isreal) then
          ! Will chase a real shift
          if (abs(A(cstrt+3,cstrt+1)) .gt. dzero) then
            ! We need to shift it to infinity
            call d_first_double_to_inf(A, B, compAB, apCurrent, apTotal, &
                    Q, compQ, qci+cstrt)
          end if
          call d_first_single_pole(mu, nu, A, B, compAB, apCurrent, apTotal, &
                                         Q, compQ, qci+cstrt)
        else
          ! Will chase a double shift
          if (abs(A(cstrt+4,cstrt+2)) .gt. dzero) then
            ! We place the double block first
            call d_swap12(cstrt+2,cstrt+1, A, B, compAB, apCurrent, apTotal, &
                        Q, compQ, qci+cstrt+1, Z, compZ, zci+cstrt)
          end if
          call d_first_double_pole(cmplx(mu,nu,kind=dp), A, B, compAB, &
                  apCurrent, apTotal, Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)
        end if

        !Chase
        nit = nit + 1
        i = cstrt+1
        if (isreal) then
          ! Chase a single shift
          do while(i .lt. cstp-2)
            if (abs(A(i+3,i+1)) .gt. dzero) then
              ! Next is a double pole
              call d_swap12(i+1,i, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+i, Z, compZ, zci+i-1)
              i = i + 2
              nswaps = nswaps + 2
            else
              ! Next is a single pole
              call d_swap11(i+1,i, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+i, Z, compZ, zci+i-1)
              i = i + 1
              nswaps = nswaps + 1
            end if
          end do
          if (i .eq. cstp-2) then
            ! can take one more single swap
            call d_swap11(i+1,i, A, B, compAB, apCurrent, apTotal, &
                        Q, compQ, qci+i, Z, compZ, zci+i-1)
            i = i + 1
            nswaps = nswaps + 1
          end if
        else
          ! Chase a double shift
          do while (i .lt. cstp-3)
            if (abs(A(i+4,i+2)) .gt. dzero) then
              ! Next is a double pole
              call d_swap22(i+1,i, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+i, Z, compZ, zci+i-1)
              i = i + 2
              nswaps = nswaps + 4

            else
              ! Next is a single pole
              call d_swap21(i+1,i, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+i, Z, compZ, zci+i-1)
              i = i + 1
              nswaps = nswaps + 2
            end if
          end do
          if (i .eq. cstp-3) then
            ! Can take one more single swap
            call d_swap21(i+1,i, A, B, compAB, apCurrent, apTotal, &
                        Q, compQ, qci+i, Z, compZ, zci+i-1)
            i = i + 1
            nswaps = nswaps + 1
          end if
        end if

        ! Check for deflations at end
        if (.not.(abs(A(cstp,cstp-2)) .gt. dzero)) then
          call d_stop_deflation_single(defl, A, B, compAB,apCurrent, apTotal,&
                                    Z, compZ, zci+cstp-2)
        else
          call d_stop_deflation_double(defl, A, B, compAB,apCurrent, apTotal,&
                                    Q, compQ, qci+cstp-2,Z, compZ, zci+cstp-3)
        end if
        ! apCurrent is updated; consequently cstrt, and cstp
        if (defl) then
          ld = 0
        end if

        if (cstp - cstrt .gt. 3) then ! Check again
          ! set the pole
          ! Francis poles
          call d_get_pole(A, B, cstrt+1, mu, nu,isreal)
          ! infinity poles
          ! isreal = .true.
          ! mu = done
          ! nu = dzero
          if (isreal) then
            if (abs(A(cstp,cstp-2)) .gt. dzero) then
              ! Shift to Inf
              call d_last_double_to_inf( A, B, compAB, apCurrent, apTotal, &
                      Z, compZ, zci+cstp-3)
              ! Introduce real pole
              call d_last_single_pole(mu, nu, A, B, compAB, apCurrent,apTotal,&
                                   Z, compZ, zci+cstp-2)
              call d_swap11(cstp-1,cstp-2, A, B, compAB, apCurrent, apTotal, &
                             Q, compQ, qci+cstp-2, Z, compZ, zci+cstp-3)
            else
              ! Introduce real pole
              call d_last_single_pole(mu, nu, A, B, compAB, apCurrent,apTotal,&
                                   Z, compZ, zci+cstp-2)
            end if
          else
            ! Need to introduce a C.C. pole
            if (abs(A(cstp-1,cstp-3)) .gt. dzero) then
              call d_swap21(cstp-2,cstp-3, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+cstp-3, Z, compZ, zci+cstp-4)
              ds = .true.
            else
              ds = .false.
            end if
            call d_last_double_pole(cmplx(mu,nu,kind=dp), A, B, compAB, &
                   apCurrent,apTotal, Q, compQ,qci+cstp-2, Z, compZ, zci+cstp-3)
            if (ds) then
              ! Avoids stagnation
              call d_swap12(cstp-2,cstp-3, A, B, compAB, apCurrent, apTotal, &
                          Q, compQ, qci+cstp-3, Z, compZ, zci+cstp-4)
            end if
          end if
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

  ! Normalize all 2x2 blocks
  apCurrent => apTotal ! This allows us to run make_hess if needed
  call d_normalizeSchur(A,B,compAB,apCurrent,apTotal,&
                            Q,compQ,qci,Z,compZ,zci,sa2)

  if (sa2) then
    call d_deallocateSwap4Mem
    call d_deallocateRCMem
  end if

  ! free everything
  apCurrent => null()
  call tApFree(apHead)
  call tApFree(apTotal)

  end subroutine

  subroutine d_getEigenvalues(A,B,strtidx,stpidx,realev,k,ccevr,ccevc,l)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Retrieve the eigenvalues for a pencil that has been (partially)
  ! reduced to generalized Schur form.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [IN]
  !         block upper triangular matrix
  ! B       double array [IN]
  !         upper triangular matrix
  ! strtidx integer [IN]
  !         pencil (A,B) is assumed to be in generalized Schur form from
  !         strtidx+1:stpidx
  ! stpidx  integer [IN]
  !         pencil (A,B) is assumed to be in generalized Schur form from
  !         strtidx+1:stpidx
  ! realev  double array [OUT]
  !         array that contains the k real eigenvalues in the range
  !         strtidx+1:stpidx
  ! k       integer [OUT]
  !         number of real eigenvalues in range strtidx+1:stpidx
  ! ccevr   double array [OUT]
  !         real parts of the l pairs of complex conjugate eigenvalues
  !         in range strtidx+1:stpidx
  ! ccevc   double array [OUT]
  !         complex parts of the l pairs of complex conjugate eigenvalues
  !         in range strtidx+1:stpidx
  ! l       integer [OUT]
  !         number of pairs of complex conjugate eigenvalues
  !         l and k satisfy : k + 2*l = stpidx - strtidx
  !___________________________________________________________________________
  ! last edit: January 15, 2019
    real(kind=dp), intent(in)   :: A(:,:), B(:,:)
    integer, intent(in)         :: strtidx, stpidx
    real(kind=dp), intent(out)  :: realev(:), ccevr(:), ccevc(:)
    integer, intent(out)        :: k, l

    integer :: i
    logical :: isreal
    real(kind=dp) :: mu, nu

    i = strtidx+1
    k = 0
    l = 0
    do while (i .lt. stpidx)
      if (abs(A(i+1,i)) .gt. dzero) then
        call d_eigenvalues2x2(A,B,i,i,mu,nu,isreal)
        if (isreal) then
          realev(k+1) = mu
          realev(k+2) = nu
          k = k + 2
        else
          ccevr(l+1) = mu
          ccevc(l+1) = nu
          l = l + 1
        end if
        i = i + 2
      else
        realev(k+1) = A(i,i)/B(i,i)
        k = k + 1
        i = i + 1
      end if
    end do
    if (i .eq. stpidx) then
      realev(k+1) = A(i,i)/B(i,i)
      k = k + 1
    end if
  end subroutine
end module d_RQZ
