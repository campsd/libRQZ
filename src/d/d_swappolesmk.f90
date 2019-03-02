! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Swap poles in real block-Hessenberg pairs:
!     m poles with next k poles
! ___________________________________________________________________
module d_swappolesmk

use u_parameters
use u_activeparts
use d_ctransformations
use d_memorymgmt
use d_computepoles
use d_swappoles12
use d_swappoles22

implicit none
private
public d_swapmk

contains

  subroutine d_swapmk(i, m, k, A, B, compAB, apC, apT, Q, compQ, qcol, &
                      Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! swaps a batch of m poles starting at pole i up to i+m-1 with poles
  ! i+m up to i+m+k(-1).
  ! This moves the entire batch k(+1) positions down along the subdiagonal
  ! The pencil and Q and Z are updated as requested via BLAS operations
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           Index (column in A,B) of the first pole in the batch
  ! m       integer [IN]
  !           Batch size of poles
  ! k       integer [INOUT]
  !           Number of swaps to be made
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
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
  !           holds when calling d_swapm1. c = current, t = total
  ! Q       double array [INOUT]
  !           Orthonormal equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, ..., qcol+m
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, ..., zcol+m
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    integer, intent(in)               :: i, m, qcol, zcol
    integer, intent(inout)            :: k
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT

    integer                           :: j, jj, LDAB, LDQ, LDQM, LDZ, &
                                         LDZM, LDC, LDR
    integer, pointer                  :: cstrt, cstp, tstrt, tstp

    type(tAp), pointer                :: apCI, apTI
    logical                           :: h1, h2

    external DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)


    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDQM,QM)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDZM,ZM)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)

    if (i+m+k+1 .le. LDAB) then
      if  (abs(A(i+m+k+1,i+m+k-1)) .gt. dzero) then
        k = k + 1
      end if
    end if

    call tApInit(apCI,i, i+m+k-1)
    call tApInit(apTI,i, i+m+k-1)

    QM = dzero
    do j=1, m+k
      QM(j,j) = done
    end do
    ZM = dzero
    do j=1, m+k
      ZM(j,j) = done
    end do


    j = m
    do while (j .ge. 1)
      jj = 0
      h1 = .false.
      if (j .ge. 2) then
        if (abs(A(i+j,i+j-2)) .gt. dzero) then
          h1 = .true.
          ! this pole is a pair of C.C. shifts
          do while (jj .lt. k)
            h2 = .false.
            if (i+j+jj .lt. LDAB-1) then
              if (abs(A(i+j+jj+2,i+j+jj)) .gt. dzero) then
                h2 = .true.
                ! next pole is C.C.
                call d_swap22(i+j+jj-1, i+j+jj-2, A, B, .false., apCI, apTI,&
                       QM, .true., j+jj-1, ZM, .true., j+jj-1)
                jj = jj + 2
              end if
            end if
            if (.not. h2) then
              ! next pole is real
              call d_swap21(i+j+jj-1, i+j+jj-2, A, B, .false., apCI, apTI,&
                     QM, .true., j+jj-1, ZM, .true., j+jj-1)
              jj = jj + 1
            end if
          end do
          j = j - 2
        end if
      end if
      if (.not. h1) then
        ! this is a real pole
        do while (jj .lt. k)
          h2 = .false.
          if (i+j+jj .lt. LDAB-1) then
            if (abs(A(i+j+jj+2,i+j+jj)) .gt. dzero) then
              h2 = .true.
              ! next pole is C.C.
              call d_swap12(i+j+jj, i+j+jj-1, A, B, .false., apCI, apTI,&
                     QM, .true., j+jj,ZM, .true., j+jj)
              jj = jj + 2
            end if
          end if
          if (.not. h2) then
            call d_swap11(i+j+jj, i+j+jj-1, A, B, .false., apCI, apTI,&
                   QM, .true., j+jj,ZM, .true., j+jj)
            jj = jj + 1
          end if
        end do
        j = j - 1
      end if
    end do

    ! Update the pencil as requested via BLAS
    if (compAB) then
      ! A row update
      call DGEMM('T','N',m+k,tstp-i-m-k+1,m+k,done,QM(1,1),LDQM,&
                 A(i+1,i+m+k),LDAB,dzero,Ru(1,1),LDR)
      A(i+1:i+m+k,i+m+k:tstp) = Ru(1:m+k,1:tstp-i-m-k+1)
      ! B row update
      call DGEMM('T','N',m+k,tstp-i-m-k+1,m+k,done,QM(1,1),LDQM,&
                 B(i+1,i+m+k),LDAB,dzero,Ru(1,1),LDR)
      B(i+1:i+m+k,i+m+k:tstp) = Ru(1:m+k,1:tstp-i-m-k+1)
      ! A column update
      call DGEMM('N','N',i-tstrt,m+k,m+k,done,A(tstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 dzero,Cu(1,1),LDC)
      A(tstrt+1:i,i:i+m+k-1) = Cu(1:i-tstrt,1:m+k)
      ! B column update
      call DGEMM('N','N',i-tstrt,m+k,m+k,done,B(tstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 dzero,Cu(1,1),LDC)
      B(tstrt+1:i,i:i+m+k-1) = Cu(1:i-tstrt,1:m+k)
    else
      ! A row update
      call DGEMM('T','N',m+k,cstp-i-m-k+1,m+k,done,QM(1,1),LDQM,&
                A(i+1,i+m+k),LDAB,dzero,Ru(1,1),LDR)
      A(i+1:i+m+k,i+m+k:cstp) = Ru(:,1:cstp-i-m-k+1)
      ! B row update
      call DGEMM('T','N',m+k,cstp-i-m-k+1,m+k,done,QM(1,1),LDQM,&
                 B(i+1,i+m+k),LDAB,dzero,Ru(1,1),LDR)
      B(i+1:i+m+k,i+m+k:cstp) = Ru(:,1:cstp-i-m-k+1)
      ! A column update
      call DGEMM('N','N',i-cstrt,m+k,m+k,done,A(cstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 dzero,Cu(1,1),LDC)
      A(cstrt+1:i,i:i+m+k-1) = Cu(1:i-cstrt,1:m+k)
      ! B column update
      call DGEMM('N','N',i-cstrt,m+k,m+k,done,B(cstrt+1,i),LDAB,ZM(1,1),LDZM,&
                 dzero,Cu(1,1),LDC)
      B(cstrt+1:i,i:i+m+k-1) = Cu(1:i-cstrt,1:m+k)
    end if

    if (compQ) then
      call DGEMM('N','N',LDQ,m+k,m+k,done,Q(1,qcol),LDQ,QM(1,1),LDQM,&
                dzero,Cu,LDC)
      Q(1:LDQ,qcol:qcol+m+k-1) = Cu(1:LDQ,1:m+k)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,m+k,m+k,done,Z(1,zcol),LDZ,ZM(1,1),LDZM,&
           dzero,Cu(1,1),LDC)
      Z(1:LDZ,zcol:zcol+m+k-1) = Cu(1:LDZ,1:m+k)
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)
  end subroutine
end module
