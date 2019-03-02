! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Swap poles in real block-Hessenberg pairs:
!     2x2 with 2x2
! ___________________________________________________________________
module d_swappoles22

use u_parameters
use u_activeparts
use d_ctransformations
use d_memorymgmt
use d_computepoles
use d_swappoles12

implicit none
private
public d_swap22

contains
  subroutine d_swap22(i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
  !   DESCRIPTION
  !___________________________________________________________________________
  ! Swaps complex conjugate pole/ev at  (i:i+1,j:j+1)
  ! with complex conjugate pole/ev at (i+2:i+3,j+2:j+3)
  ! The pencil (A,B) and equivalences Q and Z are updated as requested
  ! It is assumed that both are complex conjugate pairs but this
  ! is still explicitly checked in this subroutine and if one or both
  ! are actually real, alternative swapping methods are used.
  !
  ! This procedure is based on:
  ! "Swapping 2x2 blocks in the Schur and generalized Schur form"
  ! -- D. Camps, N. Mastronardi, R. Vandebril, and P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (row in A,B) of the first pole/ev
  ! j       integer [IN]
  !           Index (column in A,B) of the first pole/ev
  !           It is assumed that (i:i+3,j:j+3) is block upper triangular
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [IN]
  !           Current active part.
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  !           Normally, the following order,
  !             tstrt <= cstrt < i < i+3 < cstp <= tstp,
  !           holds when calling d_swap. c = current, t = total
  ! Q       double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1, qcol+2, qcol+3
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1, zcol+2, zcol+3
  !___________________________________________________________________________
  ! last edit: January 16, 2019
    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)             :: compAB, compQ, compZ
    integer, intent(in)             :: i, j, qcol, zcol
    type(tAp), pointer, intent(in)  :: apC, apT


    real(kind=dp)                   :: l1,l2,l3,l4,M(8,8),TAU(2),WORK(138),&
                                       c,s,r, RHS(8), Ac(4,4), Bc(4,4)
    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDAB, LDC, LDR, LDQ, LDZ, ERR, k, &
                                       kmax, IPIV(8)
    logical                         :: isreal1, isreal2, nconv
    real(kind=dp) DLANGE
    external                        :: DGEMM, DGEQRF, DGERQF, DORMQR, DORGQR,&
                                         DORGRQ, DTGEX2, DGESVD, ZGESVD, DGETRF

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)

    call d_eigenvalues2x2(A,B,i,j,l1,l2,isreal1)
    call d_eigenvalues2x2(A,B,i+2,j+2,l3,l4,isreal2)

    ! Check if all eigenvalues are complex conjugates
    if (isreal1) then
      call d_make_hess(l1,i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
      if (isreal2) then
        call d_make_hess(l3,i+2,j+2,A,B,compAB,apC,apT,&
                 Q,compQ,qcol+2,Z,compZ,zcol+2)
        call d_swap11(i+1,j+1,A,B,compAB,apC,apT,&
                 Q,compQ,qcol+1,Z,compZ,zcol+1)
        call d_swap11(i+2,j+2,A,B,compAB,apC,apT,&
                 Q,compQ,qcol+2,Z,compZ,zcol+2)
        call d_swap11(i,j,A,B,compAB,apC,apT,&
                 Q,compQ,qcol,Z,compZ,zcol)
        call d_swap11(i+1,j+1,A,B,compAB,apC,apT,&
                 Q,compQ,qcol+1,Z,compZ,zcol+1)
      else
        call d_swap12(i+1,j+1,A,B,compAB,apC,apT,&
                 Q,compQ,qcol+1,Z,compZ,zcol+1)
        call d_swap12(i,j,A,B,compAB,apC,apT,&
                 Q,compQ,qcol,Z,compZ,zcol)
      end if
      return
    elseif (isreal2) then
      call d_make_hess(l3,i+2,j+2,A,B,compAB,apC,apT,&
               Q,compQ,qcol+2,Z,compZ,zcol+2)
      call d_swap21(i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
      call d_swap21(i+1,j+1,A,B,compAB,apC,apT,Q,compQ,qcol+1,Z,compZ,zcol+1)
      return
    end if

    ! Check the gap between the eigenvalues
    if (sqrt((l1-l3)**2 + (l2-l4)**2) .lt. 2 * mp * &
          (DLANGE('F',4,4,A(i,j),LDAB,WORK) + &
           DLANGE('F',4,4,B(i,j),LDAB,WORK))) then
      return
    end if
    Ac = A(i:i+3,j:j+3)
    Bc = B(i:i+3,j:j+3)
    ! True 2x2 swap with large enough gap
    ! Create Kronecker system by copying data from (A,B)
    M = dzero
    ! I (x) A11
    M(1:2,1:2) = A(i:i+1,j:j+1)
    M(3:4,3:4) = A(i:i+1,j:j+1)
    ! I (x) B11
    M(5:6,1:2) = B(i:i+1,j:j+1)
    M(7:8,3:4) = B(i:i+1,j:j+1)
    ! A22T (x) I
    M(1,5) = A(i+2,j+2)
    M(2,6) = A(i+2,j+2)
    M(3,5) = A(i+2,j+3)
    M(4,6) = A(i+2,j+3)
    M(1,7) = A(i+3,j+2)
    M(2,8) = A(i+3,j+2)
    M(3,7) = A(i+3,j+3)
    M(4,8) = A(i+3,j+3)
    ! B22T (x) I
    M(5,5) = B(i+2,j+2)
    M(6,6) = B(i+2,j+2)
    M(7,5) = B(i+2,j+3)
    M(8,6) = B(i+2,j+3)
    M(7,7) = B(i+3,j+3)
    M(8,8) = B(i+3,j+3)

    ! Create RHS
    RHS(1:2) = -A(i:i+1,j+2)
    RHS(3:4) = -A(i:i+1,j+3)
    RHS(5:6) = -B(i:i+1,j+2)
    RHS(7:8) = -B(i:i+1,j+3)

    ! LU factorization of M with partial pivoting
    call DGETRF(8,8,M,8,IPIV,ERR)
    if (ERR .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
                '[0m factorization in 2 with 2 swap failed. element: ', ERR, &
                 'value: ', M(ERR,ERR)
      ! We try with a small perturbation
      M(ERR,ERR) = 2*mp
    end if
    ! Solve system
    call DGETRS('N',8,1,M,8,IPIV,RHS,8,ERR)
    if (ERR .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
      '[0m solving 2 with 2 swap failed. '
      ! We undo it ...
      A(i:i+3,j:j+3) = Ac
      B(i:i+3,j:j+3) = Bc
      return
    end if
    ! TODO iterative refinement when required

    ! Create L in Q4(3:4,:)
    Q4 = dzero
    Q4(3,1) = done
    Q4(4,2) = done
    Q4(3:4,3) = RHS(5:6)
    Q4(3:4,4) = RHS(7:8)
    ! Create Q
    call DGERQF(2,4,Q4(3,1),4,TAU,WORK,128,ERR)
    call DORGRQ(4,4,2,Q4(1,1),4,TAU,WORK,128,ERR)

    ! Row update within window
    call DGEMM('N','N',4,4,4,done,Q4(1,1),4,A(i,j),LDAB,dzero,Ru(1,1),LDR)
    A(i:i+3,j:j+3) = Ru(1:4,1:4)
    call DGEMM('N','N',4,4,4,done,Q4(1,1),4,B(i,j),LDAB,dzero,Ru(1,1),LDR)
    B(i:i+3,j:j+3) = Ru(1:4,1:4)

    ! Create R in Z4(:,1:2)
    Z4 = dzero
    Z4(1:2,1) = RHS(1:2)
    Z4(1:2,2) = RHS(3:4)
    Z4(3,1) = done
    Z4(4,2) = done
    ! Create Z
    call DGEQRF(4,2,Z4(1,1),4,TAU,WORK,128,ERR)
    call DORGQR(4,4,2,Z4(1,1),4,TAU,WORK,128,ERR)

    ! Column update within window
    call DGEMM('N','N',4,4,4,done,A(i,j),LDAB,Z4(1,1),4,dzero,Cu(1,1),LDC)
    A(i:i+3,j:j+3) = Cu(1:4,1:4)
    call DGEMM('N','N',4,4,4,done,B(i,j),LDAB,Z4(1,1),4,dzero,Cu(1,1),LDC)
    B(i:i+3,j:j+3) = Cu(1:4,1:4)

    ! Check error and update if required
    if ((DLANGE('F',2,2,A(i+2,j),LDAB,WORK) .gt. &
         mp*(DLANGE('F',4,4,A(i,j),LDAB,WORK))) .or. &
         (DLANGE('F',2,2,B(i+2,j),LDAB,WORK) .gt. &
         mp*(DLANGE('F',4,4,B(i,j),LDAB,WORK)))) then
      nconv = .true.
      k = 1
      kmax = 4
      do while (nconv .and. (k .lt. kmax))
        call updateTransformation(A,B,i,j)
        k = k + 1
        if ((DLANGE('F',2,2,A(i+2,j),LDAB,WORK) .lt. &
            mp*(DLANGE('F',4,4,A(i,j),LDAB,WORK))) .and. &
            (DLANGE('F',2,2,B(i+2,j),LDAB,WORK) .lt. &
            mp*(DLANGE('F',4,4,B(i,j),LDAB,WORK)))) then
          nconv = .false.
        end if
      end do
      if ((k .eq. kmax) .and. nconv) then
        ! Swap did not converge within reasonable number of iterations
        write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//'[0m 2x2 with &
                  &2x2 swap not converged'
        ! We undo it ...
        A(i:i+3,j:j+3) = Ac
        B(i:i+3,j:j+3) = Bc
        return
      end if
    end if

    ! Make B triangular if required
    if (abs(B(i+1,j)) .gt. dzero) then
      call d_compute_ct(B(i+1,j+1),B(i+1,j),c,s,r)
      call d_apply_ct_r(A(i:i+1,j),A(i:i+1,j+1),c,s)
      call d_apply_ct_r(B(i:i+1,j),B(i:i+1,j+1),c,s)
      call d_apply_ct_r(Z4(:,1),Z4(:,2),c,s)
    end if
    if (abs(B(i+3,j+2)) .gt. dzero) then
      call d_compute_ct(B(i+3,j+3),B(i+3,j+2),c,s,r)
      call d_apply_ct_r(A(i:i+3,j+2),A(i:i+3,j+3),c,s)
      call d_apply_ct_r(B(i:i+3,j+2),B(i:i+3,j+3),c,s)
      call d_apply_ct_r(Z4(:,3),Z4(:,4),c,s)
    end if

    ! Update outside window and equivalence
    if (compAB) then
      ! Row
      if (j+3 .lt. tstp) then
        call DGEMM('N','N',4,tstp-j-3,4,done,Q4(1,1),4,A(i,j+4),LDAB,&
               dzero,Ru(1,1),LDR)
        A(i:i+3,j+4:tstp) = Ru(1:4,1:tstp-j-3)
        call DGEMM('N','N',4,tstp-j-3,4,done,Q4(1,1),4,B(i,j+4),LDAB,&
               dzero,Ru(1,1),LDR)
        B(i:i+3,j+4:tstp) = Ru(1:4,1:tstp-j-3)
      end if
      ! Column
      if (i .gt. tstrt+1) then
        call DGEMM('N','N',i-tstrt-1,4,4,done,A(tstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(tstrt+1:i-1,j:j+3) = Cu(1:i-tstrt-1,1:4)
        call DGEMM('N','N',i-tstrt-1,4,4,done,B(tstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(tstrt+1:i-1,j:j+3) = Cu(1:i-tstrt-1,1:4)
      end if
    else
      ! Row
      if (j+3 .lt. cstp) then
        call DGEMM('N','N',4,cstp-j-3,4,done,Q4(1,1),4,A(i,j+4),LDAB,&
              dzero,Ru(1,1),LDR)
        A(i:i+3,j+4:cstp) = Ru(1:4,1:cstp-j-3)
        call DGEMM('N','N',4,cstp-j-3,4,done,Q4(1,1),4,B(i,j+4),LDAB,&
              dzero,Ru(1,1),LDR)
        B(i:i+3,j+4:cstp) = Ru(1:4,1:cstp-j-3)
      end if
      ! Column
      if (i .gt. cstrt+1) then
        call DGEMM('N','N',i-cstrt-1,4,4,done,A(cstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(cstrt+1:i-1,j:j+3) = Cu(1:i-cstrt-1,1:4)
        call DGEMM('N','N',i-cstrt-1,4,4,done,B(cstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(cstrt+1:i-1,j:j+3) = Cu(1:i-cstrt-1,1:4)
      end if
    end if

    if (compQ) then
      call DGEMM('N','T',LDQ,4,4,done,Q(1,qcol),LDQ,Q4(1,1),4,dzero,Cu(1,1),LDC)
      Q(:,qcol:qcol+3) = Cu(1:LDQ,1:4)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,4,4,done,Z(1,zcol),LDZ,Z4(1,1),4,dzero,Cu(1,1),LDC)
      Z(:,zcol:zcol+3) = Cu(1:LDZ,1:4)
    end if

    ! Set to zero
    A(i+2:i+3,j:j+1) = dzero
    B(i+2:i+3,j:j+1) = dzero
    B(i+1,j) = dzero
    B(i+3,j+2) = dzero
  end subroutine

  subroutine updateTransformation(A,B,i,j)
  !   DESCRIPTION
  !___________________________________________________________________________
  ! Update the swapping transformation to reduce magnitude of subdiagonal block
  ! in 2-by-2 with 2-by-2 swap, after invoking this function, the
  ! norm of the A_{21} and B_{21} blocks should decrease.
  ! B is not made upper triangular in this routine.
  !
  ! This procedure is based on:
  ! "Swapping 2x2 blocks in the Schur and generalized Schur form"
  ! -- D. Camps, N. Mastronardi, R. Vandebril, and P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper Hessenberg matrix B
  ! i       integer [IN]
  !           Index (row in A,B) of the first pole/ev
  ! j       integer [IN]
  !           Index (column in A,B) of the first pole/ev
  !           It is assumed that (i:i+3,j:j+3) is block upper triangular
  !
  ! NOTA
  !___________________________________________________________________________
  ! Internal subroutine, not public
  !___________________________________________________________________________
  ! last edit: January 3, 2019
   real(kind=dp), intent(inout)    :: A(:,:), B(:,:)
   integer,intent(in)              :: i,j

   integer       :: LDC, LDR, LDAB, ERR, IPIV(8)
   real(kind=dp) :: At(2,2), Bt(2,2),D(2,2),E(2,2),X(2,2),&
                    Y(2,2),Qup(4,4),Zup(4,4), SVAL(2), CS(2), SN(2), WORK(134),&
                    Z(8,8), RHS(8)
   external DGEMM, DLASV2, DGESVD, DGESV

   call d_getLeadingDim(LDAB,A)
   call d_getLeadingDim(LDC,Cu)
   call d_getLeadingDim(LDR,Ru)

   ! Create Kronecker product system
   Z = dzero
   Z(1:2,1:2) = A(i+2:i+3,j+2:j+3)
   Z(3:4,3:4) = A(i+2:i+3,j+2:j+3)
   Z(5:6,1:2) = B(i+2:i+3,j+2:j+3)
   Z(7:8,3:4) = B(i+2:i+3,j+2:j+3)
   Z(1,5) = -A(i,j)
   Z(2,6) = -A(i,j)
   Z(1,7) = -A(i+1,j)
   Z(2,8) = -A(i+1,j)
   Z(3,5) = -A(i,j+1)
   Z(4,6) = -A(i,j+1)
   Z(3,7) = -A(i+1,j+1)
   Z(4,8) = -A(i+1,j+1)

   Z(5,5) = -B(i,j)
   Z(6,6) = -B(i,j)
   Z(5,7) = -B(i+1,j)
   Z(6,8) = -B(i+1,j)
   Z(7,5) = -B(i,j+1)
   Z(8,6) = -B(i,j+1)
   Z(7,7) = -B(i+1,j+1)
   Z(8,8) = -B(i+1,j+1)

   RHS(1:2) = A(i+2:i+3,j)
   RHS(3:4) = A(i+2:i+3,j+1)
   RHS(5:6) = B(i+2:i+3,j)
   RHS(7:8) = B(i+2:i+3,j+1)

   ! Solve it
   call DGESV(8,1,Z,8,IPIV,RHS,8,ERR)
   X = reshape(RHS(1:4),(/2,2/))
   Y = reshape(RHS(5:8),(/2,2/))

   ! TWO SIDED FACTORIZED APPROACH
   ! --------------------------------------------------------------------------
   ! X
   call DGESVD('A','A',2,2,X(1,1),2,SVAL(1),D(1,1),2,E(1,1),2,WORK,134,ERR)
   call cosine_sv(SVAL,CS,SN)

   ! Zup
   ! 11 block
   At = transpose(E)
   At(:,1) = At(:,1) * CS(1)
   At(:,2) = At(:,2) * CS(2)
   Bt = matmul(At,E)
   Zup(1:2,1:2) = Bt
   ! 12 block
   At = transpose(E)
   At(:,1) = At(:,1) * SN(1)
   At(:,2) = At(:,2) * SN(2)
   Bt = matmul(At,transpose(D))
   Zup(1:2,3:4) = Bt
   ! 21 block
   At = -D
   At(:,1) = At(:,1) * SN(1)
   At(:,2) = At(:,2) * SN(2)
   Bt = matmul(At,E)
   Zup(3:4,1:2) = Bt
   ! 22 block
   At = D
   At(:,1) = At(:,1) * CS(1)
   At(:,2) = At(:,2) * CS(2)
   Bt = matmul(At,transpose(D))
   Zup(3:4,3:4) = Bt

   ! Y
   call DGESVD('A','A',2,2,Y(1,1),2,SVAL(1),D(1,1),2,E(1,1),2,WORK,134,ERR)
   call cosine_sv(SVAL,CS,SN)

   ! Qup^T
   ! 11 block
   At = transpose(E)
   At(:,1) = At(:,1) * CS(1)
   At(:,2) = At(:,2) * CS(2)
   Bt = matmul(At,E)
   Qup(1:2,1:2) = Bt
   ! 12 block
   At = transpose(E)
   At(:,1) = At(:,1) * SN(1)
   At(:,2) = At(:,2) * SN(2)
   Bt = matmul(At,transpose(D))
   Qup(1:2,3:4) = Bt
   ! 21 block
   At = -D
   At(:,1) = At(:,1) * SN(1)
   At(:,2) = At(:,2) * SN(2)
   Bt = matmul(At,E)
   Qup(3:4,1:2) = Bt
   ! 22 block
   At = D
   At(:,1) = At(:,1) * CS(1)
   At(:,2) = At(:,2) * CS(2)
   Bt = matmul(At,transpose(D))
   Qup(3:4,3:4) = Bt

   ! Update within window
   ! Row update within window
   call DGEMM('T','N',4,4,4,done,Qup(1,1),4,A(i,j),LDAB,dzero,Ru(1,1),LDR)
   A(i:i+3,j:j+3) = Ru(1:4,1:4)
   call DGEMM('T','N',4,4,4,done,Qup(1,1),4,B(i,j),LDAB,dzero,Ru(1,1),LDR)
   B(i:i+3,j:j+3) = Ru(1:4,1:4)
   ! Column update within window
   call DGEMM('N','N',4,4,4,done,A(i,j),LDAB,Zup(1,1),4,dzero,Cu(1,1),LDC)
   A(i:i+3,j:j+3) = Cu(1:4,1:4)
   call DGEMM('N','N',4,4,4,done,B(i,j),LDAB,Zup(1,1),4,dzero,Cu(1,1),LDC)
   B(i:i+3,j:j+3) = Cu(1:4,1:4)

   ! Update Q4,Z4
   call DGEMM('T','N',4,4,4,done,Qup(1,1),4,Q4(1,1),4,dzero,Cu(1,1),LDC)
   Q4(:,:) = Cu(1:4,1:4)
   call DGEMM('N','N',4,4,4,done,Z4(1,1),4,Zup(1,1),4,dzero,Cu(1,1),LDC)
   Z4(:,:) = Cu(1:4,1:4)
  end subroutine

  subroutine cosine_sv(sv,cs,sn)
  !   DESCRIPTION
  !___________________________________________________________________________
  ! Computes the cosine and sine of the singular values up to highest precision
  !
  ! This procedure is based on:
  ! "Swapping 2x2 blocks in the Schur and generalized Schur form"
  ! -- D. Camps, N. Mastronardi, R. Vandebril, and P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! sv       double array [IN]
  !           length 2 array with singular values
  ! cs       double array [OUT]
  !           length 2 array with cosine of singular values
  ! sn       double array [OUT]
  !           length 2 array with sine of singular values
  !
  ! NOTA
  !___________________________________________________________________________
  ! Internal subroutine, not public
  !___________________________________________________________________________
  ! last edit: January 3, 2019
    real(kind=dp), intent(in)   :: sv(2)
    real(kind=dp), intent(out)  :: cs(2), sn(2)

    integer k
    do k = 1,2
     if (sv(k) .gt. 1) then
       cs(k) = done/sqrt(done+sv(k)**2)
       sn(k) = cs(k) * sv(k)
     else
       sn(k) = done/sqrt(done+done/sv(k)**2)
       cs(k) = sn(k) / sv(k)
     end if
    end do
  end subroutine
end module
