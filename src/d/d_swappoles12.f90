! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Swap poles in real block-Hessenberg pairs:
!     1x1 with 1x1, 1x1 with 2x2, and 2x2 with 1x1
! ___________________________________________________________________
module d_swappoles12

use u_parameters
use u_activeparts
use d_ctransformations
use d_memorymgmt
use d_computepoles
use z_ctransformations

implicit none
private
public d_swap11, d_swap12, d_swap21

contains

  subroutine d_swap11(i, j, A, B, compAB, apC, apT, Q, compQ, qcol, &
                    Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! swaps real poles (bHbH) or eigenvalues (bTbT) i and i+1 in (A,B)
  ! The pencil (A,B) and Q and Z are updated as requested
  ! It is assumed that both i and i+1 are single, real poles/evs
  !
  ! This procedure is based on:
  ! "A generalized eigenvalue approach for solving Ricatti equations"
  ! -- P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (row in A,B) of the first pole/ev
  ! j       integer [IN]
  !           Index (column in A,B) of the first pole/ev
  !           It is assumed that (i:i+1,j:j+1) is upper triangular
  ! A       double array [INOUT]
  !           Upper Hessenberg matrix A
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
  !             tstrt <= cstrt < i < i+1 < cstp <= tstp,
  !           holds when calling d_swap. c = current, t = total
  ! Q       double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)             :: compAB, compQ, compZ
    integer, intent(in)             :: i, j, qcol, zcol
    type(tAp), pointer, intent(in)  :: apC, apT


    real(kind=dp)                   :: c, s, r
    real(kind=dp)                   :: an, bn
    integer,pointer                 :: cstrt, cstp, tstrt, tstp

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    an = abs(A(i+1,j+1))
    bn = abs(B(i+1,j+1))

    ! Compute Z
    call d_compute_ct(A(i+1,j+1)*B(i,j+1) - B(i+1,j+1)*A(i,j+1), &
                      A(i+1,j+1)*B(i,j) - B(i+1,j+1)*A(i,j), &
                     c,s,r)

    if (compAB) then
     call d_apply_ct_r(A(tstrt+1:i+1,j),A(tstrt+1:i+1,j+1),c,s)
     call d_apply_ct_r(B(tstrt+1:i+1,j),B(tstrt+1:i+1,j+1),c,s)
    else
     call d_apply_ct_r(A(cstrt+1:i+1,j),A(cstrt+1:i+1,j+1),c,s)
     call d_apply_ct_r(B(cstrt+1:i+1,j),B(cstrt+1:i+1,j+1),c,s)
    endif

    if (compZ) then
      call d_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c,s)
    end if

    ! Compute Q
    !if (abs(B(i+2,i+1)) .ge. abs(A(i+2,i+1))) then
    if (bn .ge. an) then
    !if (abs(B(i+2,i)) .ge. abs(A(i+2,i))) then
      call d_compute_ct(B(i,j),B(i+1,j),c,s,r)
    else
      call d_compute_ct(A(i,j),A(i+1,j),c,s,r)
    end if

    if (compAB) then
     call d_apply_ct_l(A(i,j:tstp),A(i+1,j:tstp),c,s)
     call d_apply_ct_l(B(i,j:tstp),B(i+1,j:tstp),c,s)
    else
     call d_apply_ct_l(A(i,j:cstp),A(i+1,j:cstp),c,s)
     call d_apply_ct_l(B(i,j:cstp),B(i+1,j:cstp),c,s)
    end if
    A(i+1,j) = dzero
    B(i+1,j) = dzero
    if (compQ) then
     call d_apply_ct_r(Q(:,qcol),Q(:,qcol+1),c,-s)
    end if
  end subroutine

  subroutine d_swap21(i, j, A, B, compAB, apC, apT, Q, compQ, qcol, &
                    Z, compZ, zcol)
  !   DESCRIPTION
  !___________________________________________________________________________
  ! swaps complex conjugate poles (bHbH) or eigenvalues (bTbT) i/i+1
  ! with real pole/ev i+2 in (A,B)
  ! The pencil (A,B) and Q and Z are updated as requested
  ! It is assumed that i/i+1 is a complex conjugate pair of poles in
  ! standard form, and i+2 is a real pole.
  !
  ! This procedure is based on:
  ! "A generalized eigenvalue approach for solving Ricatti equations"
  ! -- P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (row in A,B) of the first pole/ev
  ! j       integer [IN]
  !           Index (column in A,B) of the first pole/ev
  !           It is assumed that (i:i+2,j:j+2) is block upper triangular
  ! A       double array [INOUT]
  !           Upper Hessenberg matrix A
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
  !             tstrt <= cstrt < i < i+1 < cstp <= tstp,
  !           holds when calling d_swap. c = current, t = total
  ! Q       double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: September 11, 2018
    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)             :: compAB, compQ, compZ
    integer, intent(in)             :: i, j, qcol, zcol
    type(tAp), pointer, intent(in)  :: apC, apT


    real(kind=dp)                   :: c, c2, s, s2, x1, x2, r, H(2,3)
    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDAB, LDC, LDR, LDQ, LDZ

    external                        :: DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    if (abs(A(i+1,j)) .gt. dzero) then

    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)

    ! Initialize Q4 and Z4
    Q4(1:3,1:3) = dzero
    Q4(1,1) = done
    Q4(2,2) = done
    Q4(3,3) = done
    Z4(1:3,1:3) = dzero
    Z4(1,1) = done
    Z4(2,2) = done
    Z4(3,3) = done

    if (abs(B(i+2,j+2)) .ge. abs(A(i+2,j+2))) then
      H = A(i+2,j+2) * B(i:i+1,j:j+2) - B(i+2,j+2) * A(i:i+1,j:j+2)
      call d_compute_ct(H(1,1),H(2,1),c,s,x2)
      call d_apply_ct_l(H(1,2:3),H(2,2:3),c,s)
      ! Compute Z
      call d_compute_ct(H(2,3),H(2,2),c,s,x1)
      x1 = c * H(1,2) -s * H(1,3)
      call d_compute_ct(x1,x2,c2,s2,r)
      ! Form Z
      call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c,s)
      call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c2,s2)
      ! Apply to (A,B) and compute Q to make B upper triangular
      ! Z1
      call d_apply_ct_r(B(i:i+2,j+1),B(i:i+2,j+2),c,s)
      call d_apply_ct_r(A(i:i+2,j+1),A(i:i+2,j+2),c,s)
      ! Q1
      call d_compute_ct(B(i+1,j+1),B(i+2,j+1),c,s,r)
      call d_apply_ct_l(B(i+1,j:j+2),B(i+2,j:j+2),c,s)
      call d_apply_ct_l(A(i+1,j:j+2),A(i+2,j:j+2),c,s)
      B(i+2,j+1) = dzero
      ! Z2
      call d_apply_ct_r(B(i:i+2,j),B(i:i+2,j+1),c2,s2)
      call d_apply_ct_r(A(i:i+2,j),A(i:i+2,j+1),c2,s2)
      ! Q2
      call d_compute_ct(B(i,j),B(i+1,j),c2,s2,r)
      call d_apply_ct_l(B(i,j:j+2),B(i+1,j:j+2),c2,s2)
      call d_apply_ct_l(A(i,j:j+2),A(i+1,j:j+2),c2,s2)
      B(i+1,j) = dzero
      ! Form Q
      call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c,s)
      call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c2,s2)
      A(i+1:i+2,j) = dzero
    else
      ! Make A upper triangular and reverse roles of A and B
      call d_compute_ct(A(i+1,j+1),A(i+1,j),c,s,r)
      call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c,s)
      call d_apply_ct_r(A(i:i+1,j),A(i:i+1,j+1),c,s)
      call d_apply_ct_r(B(i:i+1,j),B(i:i+1,j+1),c,s)
      A(i+1,j) = dzero
      ! Continue
      H = A(i+2,j+2) * B(i:i+1,j:j+2) - B(i+2,j+2) * A(i:i+1,j:j+2)
      call d_compute_ct(H(1,1),H(2,1),c,s,x2)
      call d_apply_ct_l(H(1,2:3),H(2,2:3),c,s)
      ! Compute Z
      call d_compute_ct(H(2,3),H(2,2),c,s,x1)
      x1 = c * H(1,2) -s * H(1,3)
      call d_compute_ct(x1,x2,c2,s2,r)
      ! Form Z
      call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c,s)
      call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c2,s2)
      ! Apply to (A,B) and compute Q to make A upper triangular
      ! Z1
      call d_apply_ct_r(B(i:i+2,j+1),B(i:i+2,j+2),c,s)
      call d_apply_ct_r(A(i:i+2,j+1),A(i:i+2,j+2),c,s)
      ! Q1
      call d_compute_ct(A(i+1,j+1),A(i+2,j+1),c,s,r)
      call d_apply_ct_l(B(i+1,j:j+2),B(i+2,j:j+2),c,s)
      call d_apply_ct_l(A(i+1,j:j+2),A(i+2,j:j+2),c,s)
      ! Z2
      call d_apply_ct_r(B(i:i+2,j),B(i:i+2,j+1),c2,s2)
      call d_apply_ct_r(A(i:i+2,j),A(i:i+2,j+1),c2,s2)
      ! Q2
      call d_compute_ct(A(i,j),A(i+1,j),c2,s2,r)
      call d_apply_ct_l(B(i,j:j+2),B(i+1,j:j+2),c2,s2)
      call d_apply_ct_l(A(i,j:j+2),A(i+1,j:j+2),c2,s2)
      ! Form Q
      call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c,s)
      call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c2,s2)
      ! Finish off by making B triangular
      call d_compute_ct(B(i+2,j+2),B(i+2,j+1),c,s,r)
      call d_apply_ct_r(B(i:i+2,j+1),B(i:i+2,j+2),c,s)
      call d_apply_ct_r(A(i:i+2,j+1),A(i:i+2,j+2),c,s)
      call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c,s)
      B(i+1:i+2,j) = dzero
      B(i+2,j+1) = dzero
      A(i+1:i+2,j) = dzero
    end if

    if (compAB) then
      ! Row update
      if (j+3 .le. tstp) then
        call DGEMM('N','N',3,tstp-j-2,3,done,Q4(1,1),4,A(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        A(i:i+2,j+3:tstp) = Ru(1:3,1:tstp-j-2)
        call DGEMM('N','N',3,tstp-j-2,3,done,Q4(1,1),4,B(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        B(i:i+2,j+3:tstp) = Ru(1:3,1:tstp-j-2)
      end if
      ! Column update
      if (i-1 .ge. tstrt+1) then
        call DGEMM('N','N',i-tstrt-1,3,3,done,A(tstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(tstrt+1:i-1,j:j+2) = Cu(1:i-tstrt-1,1:3)
        call DGEMM('N','N',i-tstrt-1,3,3,done,B(tstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(tstrt+1:i-1,j:j+2) = Cu(1:i-tstrt-1,1:3)
      end if
    else
      ! Row update
      if (j+3 .le. cstp) then
        call DGEMM('N','N',3,cstp-j-2,3,done,Q4(1,1),4,A(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        A(i:i+2,j+3:cstp) = Ru(1:3,1:cstp-j-2)
        call DGEMM('N','N',3,cstp-j-2,3,done,Q4(1,1),4,B(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        B(i:i+2,j+3:cstp) = Ru(1:3,1:cstp-j-2)
      end if
      ! Column update
      if (i-1 .ge. cstrt+1) then
        call DGEMM('N','N',i-cstrt-1,3,3,done,A(cstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(cstrt+1:i-1,j:j+2) = Cu(1:i-cstrt-1,1:3)
        call DGEMM('N','N',i-cstrt-1,3,3,done,B(cstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(cstrt+1:i-1,j:j+2) = Cu(1:i-cstrt-1,1:3)
      end if
    end if

    if (compQ) then
      call DGEMM('N','T',LDQ,3,3,done,Q(1,qcol),LDQ,Q4(1,1),4,dzero,Cu(1,1),LDC)
      Q(1:LDQ,qcol:qcol+2) = Cu(1:LDQ,1:3)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,3,3,done,Z(1,zcol),LDZ,Z4(1,1),4,dzero,Cu(1,1),LDC)
      Z(1:LDZ,zcol:zcol+2) = Cu(1:LDZ,1:3)
    end if

    else
      ! We use the single swap twice
      call d_swap11(i+1,j+1,A,B,compAB,apC,apT,Q,compQ,qcol+1,Z,compZ,zcol+1)
      call d_swap11(i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
    end if
  end subroutine

  subroutine d_swap12(i, j, A, B, compAB, apC, apT, Q, compQ, qcol, &
                    Z, compZ, zcol)
  !   DESCRIPTION
  !___________________________________________________________________________
  ! swaps real pole (bHbH) or eigenvalue (bTbT) i
  ! with complex conjugate pole/ev i+1/i+2 in (A,B)
  ! The pencil (A,B) and Q and Z are updated as requested
  ! It is assumed that i is a real pole and i+1/i+2 is a complex conjugate
  ! pair of poles in standard form
  !
  ! This procedure is based on:
  ! "A generalized eigenvalue approach for solving Ricatti equations"
  ! -- P. Van Dooren
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (row in A,B) of the first pole/ev
  ! j       integer [IN]
  !           Index (column in A,B) of the first pole/ev
  !           It is assumed that (i:i+2,j:j+2) is block upper triangular
  ! A       double array [INOUT]
  !           Upper Hessenberg matrix A
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
  !             tstrt <= cstrt < i < i+1 < cstp <= tstp,
  !           holds when calling d_swap. c = current, t = total
  ! Q       double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: September 12, 2018
    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)             :: compAB, compQ, compZ
    integer, intent(in)             :: i, j, qcol, zcol
    type(tAp), pointer, intent(in)  :: apC, apT

    real(kind=dp)                   :: c, c2, s, s2, x1, x2, r, H(3,2)
    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDAB, LDC, LDR, LDQ, LDZ

    external                        :: DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)

    ! Initialize Q4 and Z4
    Q4(1:3,1:3) = dzero
    Q4(1,1) = done
    Q4(2,2) = done
    Q4(3,3) = done
    Z4(1:3,1:3) = dzero
    Z4(1,1) = done
    Z4(2,2) = done
    Z4(3,3) = done

    if (abs(B(i,j)) .ge. abs(A(i,j))) then
      H = A(i,j) * B(i:i+2,j+1:j+2) - B(i,j) * A(i:i+2,j+1:j+2)
      call d_compute_ct(H(3,2),H(3,1),c,s,x2)
      call d_apply_ct_r(H(1:2,1),H(1:2,2),c,s)
      ! Compute Q
      call d_compute_ct(H(1,1),H(2,1),c,s,x1)
      x1 = -s * H(1,2) +c * H(2,2)
      call d_compute_ct(x1,x2,c2,s2,r)
      ! Form Q
      call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c,s)
      call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c2,s2)
      ! Apply to (A,B) and compute Z to make B upper triangular
      ! Q1
      call d_apply_ct_l(B(i,j:j+2),B(i+1,j:j+2),c,s)
      call d_apply_ct_l(A(i,j:j+2),A(i+1,j:j+2),c,s)
      ! Z1
      call d_compute_ct(B(i+1,j+1),B(i+1,j),c,s,r)
      call d_apply_ct_r(B(i:i+2,j),B(i:i+2,j+1),c,s)
      call d_apply_ct_r(A(i:i+2,j),A(i:i+2,j+1),c,s)
      B(i+1,j) = dzero
      ! Q2
      call d_apply_ct_l(B(i+1,j:j+2),B(i+2,j:j+2),c2,s2)
      call d_apply_ct_l(A(i+1,j:j+2),A(i+2,j:j+2),c2,s2)
      ! Z2
      call d_compute_ct(B(i+2,j+2),B(i+2,j+1),c2,s2,r)
      call d_apply_ct_r(B(i:i+2,j+1),B(i:i+2,j+2),c2,s2)
      call d_apply_ct_r(A(i:i+2,j+1),A(i:i+2,j+2),c2,s2)
      B(i+2,j+1) = dzero
      ! Form Z
      call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c,s)
      call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c2,s2)
      A(i+2,j:j+1) = dzero
    else
      ! Make A upper triangular and reverse roles of A and B
      call d_compute_ct(A(i+1,j+1),A(i+2,j+1),c,s,r)
      call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c,s)
      call d_apply_ct_l(A(i+1,j:j+2),A(i+2,j:j+2),c,s)
      call d_apply_ct_l(B(i+1,j:j+2),B(i+2,j:j+2),c,s)
      A(i+2,j+1) = dzero
      ! Continue
      H = A(i,j) * B(i:i+2,j+1:j+2) - B(i,j) * A(i:i+2,j+1:j+2)
      call d_compute_ct(H(3,2),H(3,1),c,s,x2)
      call d_apply_ct_r(H(1:2,1),H(1:2,2),c,s)
      ! Compute Q
      call d_compute_ct(H(1,1),H(2,1),c,s,x1)
      x1 = -s * H(1,2) + c * H(2,2)
      call d_compute_ct(x1,x2,c2,s2,r)
      ! Form Q
      call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c,s)
      call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c2,s2)
      ! Apply to (A,B) and compute Z to make A upper triangular
      ! Q1
      call d_apply_ct_l(B(i,j:j+2),B(i+1,j:j+2),c,s)
      call d_apply_ct_l(A(i,j:j+2),A(i+1,j:j+2),c,s)
      ! Z1
      call d_compute_ct(A(i+1,j+1),A(i+1,j),c,s,r)
      call d_apply_ct_r(B(i:i+2,j),B(i:i+2,j+1),c,s)
      call d_apply_ct_r(A(i:i+2,j),A(i:i+2,j+1),c,s)
      ! Q2
      call d_apply_ct_l(B(i+1,j:j+2),B(i+2,j:j+2),c2,s2)
      call d_apply_ct_l(A(i+1,j:j+2),A(i+2,j:j+2),c2,s2)
      ! Z2
      call d_compute_ct(A(i+2,j+2),A(i+2,j+1),c2,s2,r)
      call d_apply_ct_r(B(i:i+2,j+1),B(i:i+2,j+2),c2,s2)
      call d_apply_ct_r(A(i:i+2,j+1),A(i:i+2,j+2),c2,s2)
      ! Form Z
      call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c,s)
      call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c2,s2)
      ! Finish off by making B triangular
      call d_compute_ct(B(i,j),B(i+1,j),c,s,r)
      call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c,s)
      call d_apply_ct_l(A(i,j:j+2),A(i+1,j:j+2),c,s)
      call d_apply_ct_l(B(i,j:j+2),B(i+1,j:j+2),c,s)
      B(i+1:i+2,j) = dzero
      B(i+2,j+1) = dzero
      A(i+2,j:j+1) = dzero
    end if

    if (compAB) then
      ! Row update
      if (j+3 .le. tstp) then
        call DGEMM('N','N',3,tstp-j-2,3,done,Q4(1,1),4,A(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        A(i:i+2,j+3:tstp) = Ru(1:3,1:tstp-j-2)
        call DGEMM('N','N',3,tstp-j-2,3,done,Q4(1,1),4,B(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        B(i:i+2,j+3:tstp) = Ru(1:3,1:tstp-j-2)
      end if
      ! Column update
      if (i-1 .ge. tstrt+1) then
        call DGEMM('N','N',i-tstrt-1,3,3,done,A(tstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(tstrt+1:i-1,j:j+2) = Cu(1:i-tstrt-1,1:3)
        call DGEMM('N','N',i-tstrt-1,3,3,done,B(tstrt+1,j),LDAB,&
               Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(tstrt+1:i-1,j:j+2) = Cu(1:i-tstrt-1,1:3)
      end if
    else
      ! Row update
      if (j+3 .le. cstp) then
        call DGEMM('N','N',3,cstp-j-2,3,done,Q4(1,1),4,A(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        A(i:i+2,j+3:cstp) = Ru(1:3,1:cstp-j-2)
        call DGEMM('N','N',3,cstp-j-2,3,done,Q4(1,1),4,B(i,j+3),LDAB,&
               dzero,Ru(1,1),LDR)
        B(i:i+2,j+3:cstp) = Ru(1:3,1:cstp-j-2)
      end if
      ! Column update
      if (i-1 .ge. cstrt+1) then
        call DGEMM('N','N',i-cstrt-1,3,3,done,A(cstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(cstrt+1:i-1,j:j+2) = Cu(1:i-cstrt-1,1:3)
        call DGEMM('N','N',i-cstrt-1,3,3,done,B(cstrt+1,j),LDAB,&
              Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(cstrt+1:i-1,j:j+2) = Cu(1:i-cstrt-1,1:3)
      end if
    end if

    if (compQ) then
      call DGEMM('N','T',LDQ,3,3,done,Q(1,qcol),LDQ,Q4(1,1),4,dzero,Cu(1,1),LDC)
      Q(:,qcol:qcol+2) = Cu(1:LDQ,1:3)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,3,3,done,Z(1,zcol),LDZ,Z4(1,1),4,dzero,Cu(1,1),LDC)
      Z(:,zcol:zcol+2) = Cu(1:LDZ,1:3)
    end if
  end subroutine

end module
