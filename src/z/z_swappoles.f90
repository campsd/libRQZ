! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Swap poles in complex Hessenberg pairs
! ___________________________________________________________________
module z_swappoles

use u_parameters
use u_activeparts
use z_memorymgmt
use z_ctransformations

implicit none
private
public z_swap, z_swapmk

contains

  subroutine z_swap(i, j, A, B, compAB, apC, apT, Q, compQ, qcol, &
                    Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! swaps poles (HH) or eigenvalues (TT) i and i+1 in (A,B)
  ! The pencil (A,B) and Q and Z are updated as requested
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
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
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
  !           holds when calling z_swap. c = current, t = total
  ! Q       complex double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: July 24, 2018
    complex(kind=dp), intent(inout) :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)             :: compAB, compQ, compZ
    integer, intent(in)             :: i, j, qcol, zcol
    type(tAp), pointer, intent(in)  :: apC, apT


    complex(kind=dp)                :: s, r
    real(kind=dp)                   :: c
    real(kind=dp)                   :: an, bn
    integer,pointer                 :: cstrt, cstp, tstrt, tstp

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    an = abs(A(i+1,j+1))
    bn = abs(B(i+1,j+1))

    ! Compute Z
    call z_compute_ct(A(i+1,j+1)*B(i,j+1) - B(i+1,j+1)*A(i,j+1), &
                      A(i+1,j+1)*B(i,j) - B(i+1,j+1)*A(i,j), &
                     c,s,r)
    !c = dconjg(c) ! Own ct

    if (compAB) then
     call z_apply_ct_r(A(tstrt+1:i+1,j),A(tstrt+1:i+1,j+1),c,s)
     call z_apply_ct_r(B(tstrt+1:i+1,j),B(tstrt+1:i+1,j+1),c,s)
    else
     call z_apply_ct_r(A(cstrt+1:i+1,j),A(cstrt+1:i+1,j+1),c,s)
     call z_apply_ct_r(B(cstrt+1:i+1,j),B(cstrt+1:i+1,j+1),c,s)
    endif

    if (compZ) then
      call z_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c,s)
    end if

    ! Compute Q
    !if (abs(B(i+2,i+1)) .ge. abs(A(i+2,i+1))) then
    if (bn .ge. an) then
    !if (abs(B(i+2,i)) .ge. abs(A(i+2,i))) then
      call z_compute_ct(B(i,j),B(i+1,j),c,s,r)
    else
      call z_compute_ct(A(i,j),A(i+1,j),c,s,r)
    end if

    if (compAB) then
     call z_apply_ct_l(A(i,j:tstp),A(i+1,j:tstp),c,s)
     call z_apply_ct_l(B(i,j:tstp),B(i+1,j:tstp),c,s)
    else
     call z_apply_ct_l(A(i,j:cstp),A(i+1,j:cstp),c,s)
     call z_apply_ct_l(B(i,j:cstp),B(i+1,j:cstp),c,s)
    end if
    A(i+1,j) = czero
    B(i+1,j) = czero
    if (compQ) then
     !call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),dconjg(c),-s) ! Own ct
     call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),c,-s)
    end if
  end subroutine

  subroutine z_swapmk(i, m, k, A, B, compAB, apC, apT, Q, compQ, qcol, &
                      Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! swaps a batch of m poles starting at pole i up to i+m-1 with poles
  ! i+m up to i+m+k-1.
  ! This moves the entire batch k positions down along the subdiagonal
  ! The pencil and Q and Z are updated as requested via BLAS operations
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (column in A,B) of the first pole in the batch
  ! m       integer [IN]
  !           Batch size of poles
  ! k       integer [IN]
  !           Number of swaps to be made
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [INOUT]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  !           Normally, the following order,
  !              tstrt <= cstrt < i < i+m < cstp <= tstp,
  !           holds when calling z_swapm1. c = current, t = total
  ! Q       complex double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, ..., qcol+m
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, ..., zcol+m
  !___________________________________________________________________________
  ! last edit: August 28, 2018
    complex(kind=dp), intent(inout)   :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    integer, intent(in)               :: i, m, k, qcol, zcol
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT

    integer                           :: j, jj, LDAB, LDQ, LDQM, LDZ, &
                                         LDZM, LDC, LDR
    integer, pointer                  :: cstrt, cstp, tstrt, tstp

    type(tAp), pointer                :: apCI, apTI

    external ZGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call tApInit(apCI,i, i+m+k-1)
    call tApInit(apTI,i, i+m+k-1)

    call z_getLeadingDim(LDAB,A)
    call z_getLeadingDim(LDQ,Q)
    call z_getLeadingDim(LDQM,QM)
    call z_getLeadingDim(LDZ,Z)
    call z_getLeadingDim(LDZM,ZM)
    call z_getLeadingDim(LDC,CuM)
    call z_getLeadingDim(LDR,RuM)

    ! Initialize variables
    ! QM, ZM, RuM, and CuM are provided by z_memorymgmt
    ! (should be allocated first of size m+k)
    QM = czero
    do j=1, m+k
      QM(j,j) = cone
    end do
    ZM = czero
    do j=1, m+k
      ZM(j,j) = cone
    end do

    do j = m, 1, -1
      do jj = 0, k-1
        ! Swap the next pole through the batch of m poles
        call z_swap(i+j+jj,i+j+jj-1, A, B, .false., apCI, apTI, &
                    QM, .true., j+jj, ZM, .true.,j+jj)
      end do
    end do

    ! Update the pencil as requested via BLAS
    if (compAB) then
      ! A row update
      call ZGEMM('C','N',m+k,tstp-i-m-k+1,m+k,cone,QM(1,1),LDQM,&
                 A(i+1,i+m+k),LDAB,czero,RuM(1,1),LDR)
      A(i+1:i+m+k,i+m+k:tstp) = RuM(1:m+k,1:tstp-i-m-k+1)
      ! B row update
      call ZGEMM('C','N',m+k,tstp-i-m-k+1,m+k,cone,QM(1,1),LDQM,&
                 B(i+1,i+m+k),LDAB,czero,RuM(1,1),LDR)
      B(i+1:i+m+k,i+m+k:tstp) = RuM(1:m+k,1:tstp-i-m-k+1)
      ! A column update
      call ZGEMM('N','N',i-tstrt,m+k,m+k,cone,A(tstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      A(tstrt+1:i,i:i+m+k-1) = CuM(1:i-tstrt,1:m+k)
      ! B column update
      call ZGEMM('N','N',i-tstrt,m+k,m+k,cone,B(tstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      B(tstrt+1:i,i:i+m+k-1) = CuM(1:i-tstrt,1:m+k)
    else
      ! A row update
      call ZGEMM('C','N',m+k,cstp-i-m-k+1,m+k,cone,QM(1,1),LDQM,&
                A(i+1,i+m+k),LDAB,czero,RuM(1,1),LDR)
      A(i+1:i+m+k,i+m+k:cstp) = RuM(:,1:cstp-i-m-k+1)
      ! B row update
      call ZGEMM('C','N',m+k,cstp-i-m-k+1,m+k,cone,QM(1,1),LDQM,&
                 B(i+1,i+m+k),LDAB,czero,RuM(1,1),LDR)
      B(i+1:i+m+k,i+m+k:cstp) = RuM(:,1:cstp-i-m-k+1)
      ! A column update
      call ZGEMM('N','N',i-cstrt,m+k,m+k,cone,A(cstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      A(cstrt+1:i,i:i+m+k-1) = CuM(1:i-cstrt,1:m+k)
      ! B column update
      call ZGEMM('N','N',i-cstrt,m+k,m+k,cone,B(cstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      B(cstrt+1:i,i:i+m+k-1) = CuM(1:i-cstrt,1:m+k)
    end if

    if (compQ) then
      call ZGEMM('N','N',LDQ,m+k,m+k,cone,Q(1,qcol),LDQ,QM(1,1),LDQM,&
                czero,CuM,LDC)
      Q(1:LDQ,qcol:qcol+m+k-1) = CuM(1:LDQ,1:m+k)
    end if

    if (compZ) then
      call ZGEMM('N','N',LDZ,m+k,m+k,cone,Z(1,zcol),LDZ,ZM(1,1),LDZM,&
           czero,CuM(1,1),LDC)
      Z(1:LDZ,zcol:zcol+m+k-1) = CuM(1:LDZ,1:m+k)
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)
  end subroutine
end module
