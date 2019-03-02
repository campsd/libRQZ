! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Change poles in complex Hessenberg pairs
! ___________________________________________________________________
module z_setpoles

use u_parameters
use u_activeparts
use z_memorymgmt
use z_ctransformations
use z_swappoles

implicit none
private
public z_set_first_pole, z_set_last_pole, &
       z_set_first_m_poles, z_set_last_m_poles

contains

  subroutine z_set_first_pole(alpha, beta, A, B, compAB, apC, apT, &
                              Q, compQ, qcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the pole of (A,B) in the first position of the current
  ! active part apC to alpha / beta.
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC of the pencil is updated
  ! If compQ is true, the transformation is also applied to
  ! the relevant two columns of Q specified by qcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   complex double [IN]
  !           Defines the new first pole: alpha/beta
  ! beta    complex double [IN]
  !           Defines the new first pole: alpha/beta
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the current active part
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       complex double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+1
  !___________________________________________________________________________
  ! last edit: July 24, 2018
    integer, intent(in)               ::  qcol
    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:)
    complex(kind=dp), intent(in)      ::  alpha, beta
    logical, intent(in)               ::  compQ, compAB
    complex(kind=dp), intent(inout)   ::  Q(:,:)
    type(tAp), pointer, intent(in)    ::  apC, apT

    complex(kind=dp)                  :: v(2), s, r
    real(kind=dp)                     :: c
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    complex(kind=dp), external        :: zdotc

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    ! Compute and normalize the vector of interest
    v = A(cstrt+1:cstrt+2,cstrt+1) * beta - B(cstrt+1:cstrt+2,cstrt+1) * alpha
    s = zdotc(2,v,1,v,1)
    v = v / sqrt(abs(s))

    ! Update pair and equivalence
    if (abs(v(2)) .gt. mp) then
      call z_compute_ct(v(1),v(2),c,s,r)
      if (compAB) then
        call z_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c,s)
        call z_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c,s)
      else
        call z_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c,s)
        call z_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c,s)
      endif

      if (compQ) then
        !call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),dconjg(c),-s) ! Own ct
        call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),c,-s)
      endif
    end if
  end subroutine

  subroutine z_set_last_pole(alpha, beta, A, B, compAB, apC, apT, &
                             Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the pole of (A,B) in the last position of apC to
  ! alpha / beta
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC is updated
  ! If compZ is true, the ctransformation is also applied to
  ! the columns of Z specified by zcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   complex double [IN]
  !           Defines the new last pole: alpha/beta
  ! beta    complex double [IN]
  !           Defines the new last pole: alpha/beta
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
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: July 24, 2018
    integer, intent(in)               ::  zcol
    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:)
    complex(kind=dp), intent(in)      ::  alpha, beta
    logical, intent(in)               ::  compAB, compZ
    complex(kind=dp), intent(inout)   ::  Z(:,:)
    type(tAp), pointer, intent(in)    ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    complex(kind=dp)                  ::  v(2), s, r
    real(kind=dp)                     ::  c
    integer,pointer                   ::  cstrt, cstp, tstrt, tstp
    complex(kind=dp), external        ::  zdotc

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    ! Compute and normalize the vector of interest
    v = A(cstp,cstp-1:cstp) * beta - B(cstp,cstp-1:cstp) * alpha
    s = zdotc(2,v,1,v,1)
    v = v / sqrt(abs(s))

    if (abs(v(2)) > mp) then
      call z_compute_ct(v(2),v(1),c,s,r)
      !c = dconjg(c) ! Own ct
      if (compAB) then
          call z_apply_ct_r(A(tstrt+1:cstp,cstp-1),A(tstrt+1:cstp,cstp),c,s)
          call z_apply_ct_r(B(tstrt+1:cstp,cstp-1),B(tstrt+1:cstp,cstp),c,s)
      else
          call z_apply_ct_r(A(cstrt+1:cstp,cstp-1),A(cstrt+1:cstp,cstp),c,s)
          call z_apply_ct_r(B(cstrt+1:cstp,cstp-1),B(cstrt+1:cstp,cstp),c,s)
      endif

      if (compZ) then
        call z_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c,s)
      end if
    end if
  end subroutine

  subroutine z_set_first_m_poles(alpha,beta,A,B,compAB,apC, apT,&
                                 Q,compQ,qcol,Z,compZ,zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  !
  ! Sets the first m poles of (A,B) starting at the first position of
  ! the current active part apC to the m poles as specified by the arrays:
  !                    alpha./beta
  ! If compAB is true, the complete pencil is updated based on apT, if not,
  ! only the active part apC of the pencil is updated
  ! If compQ is true, the matrix Q is updated accordingly. If compZ is true,
  ! the matrix Z is updated accordingly
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   complex double array [IN]
  !           Defines the new first pole: alpha/beta
  ! beta    complex double array [IN]
  !           Defines the new first pole: alpha/beta
  !           alpha and beta are expected to be of the same
  !           length and this length m are the number of poles
  !           to be introduced
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       complex double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+m+1
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+m
  !___________________________________________________________________________
  ! last edit: July 24, 2019

    complex(kind=dp), intent(inout)   :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    complex(kind=dp), intent(in)      :: alpha(:), beta(:)
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT
    integer, intent(in)               :: qcol, zcol

    integer                           :: m, i, j, LDAB, LDQ, LDQM, LDZ, &
                                         LDZM, LDC, LDR
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external ZGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    m = size(alpha) ! TODO`

    call z_getLeadingDim(LDAB,A)
    call z_getLeadingDim(LDQ,Q)
    call z_getLeadingDim(LDQM,QM)
    call z_getLeadingDim(LDZ,Z)
    call z_getLeadingDim(LDZM,ZM)
    call z_getLeadingDim(LDC,CuM)
    call z_getLeadingDim(LDR,RuM)

    ! Initialize internal variables
    call tApInit(apCI,cstrt, cstrt+m)
    call tApInit(apTI,cstrt, cstrt+m)
    QM = czero
    do i=1, m+1
        QM(i,i) = cone
    end do
    ZM = czero
    do i=1, m
        ZM(i,i) = cone
    end do

    do i = m, 1, -1
        ! Change the first pole
        call z_set_first_pole(alpha(i), beta(i), A, B, .true., apCI,apTI, &
                              QM, .true., 1)
        ! Swap it to the correct position
        do j=1,i-1
          call z_swap(cstrt+j+1,cstrt+j, A, B, .true., apCI, apTI, &
                      QM, .true., j+1, ZM, .true.,j)
        end do
    end do

    ! Update the rest of the pencil and Q,Z matrices via BLAS
    if (compAB) then
      ! A row update
      call ZGEMM('C','N',m+1,tstp-cstrt-m,m+1,cone,QM(1,1),LDQM,&
                 A(cstrt+1,cstrt+m+1),LDAB,czero,RuM(1,1),LDR)
      A(cstrt+1:cstrt+m+1,cstrt+m+1:tstp) = RuM(1:m+1,1:tstp-cstrt-m)
      ! B row update
      call ZGEMM('C','N',m+1,tstp-cstrt-m,m+1,cone,QM(1,1),LDQM,&
                 B(cstrt+1,cstrt+m+1),LDAB,czero,RuM(1,1),LDR)
      B(cstrt+1:cstrt+m+1,cstrt+m+1:tstp) = RuM(1:m+1,1:tstp-cstrt-m)
      if (cstrt .gt. 0) then
        ! A column update
        call ZGEMM('N','N',cstrt-tstrt,m,m,cone,A(tstrt+1,cstrt+1),LDAB,&
                   ZM(1,1),LDZM,czero,CuM(1,1),LDC)
        A(tstrt+1:cstrt,cstrt+1:cstrt+m) = CuM(1:cstrt-tstrt,1:m)
        ! B column update
        call ZGEMM('N','N',cstrt-tstrt,m,m,cone,B(tstrt+1,cstrt+1),LDAB,&
                   ZM(1,1),LDZM,czero,CuM(1,1),LDC)
        B(tstrt+1:cstrt,cstrt+1:cstrt+m) = CuM(1:cstrt-tstrt,1:m)
      end if
    else
      ! A row update
      call ZGEMM('C','N',m+1,cstp-cstrt-m,m+1,cone,QM(1,1),LDQM,&
                 A(cstrt+1,cstrt+m+1),LDAB,czero,RuM(1,1),LDR)
      A(cstrt+1:cstrt+m+1,cstrt+m+1:cstp) = RuM(1:m+1,1:cstp-cstrt-m)
      ! B row update
      call ZGEMM('C','N',m+1,cstp-cstrt-m,m+1,cone,QM(1,1),LDQM,&
                 B(cstrt+1,cstrt+m+1),LDAB,czero,RuM(1,1),LDR)
      B(cstrt+1:cstrt+m+1,cstrt+m+1:cstp) = RuM(1:m+1,1:cstp-cstrt-m)
    end if

    if (compQ) then
      call ZGEMM('N','N',LDQ,m+1,m+1,cone,Q(1,qcol),LDQ,&
                 QM(1,1),LDQM,czero,CuM,LDC)
      Q(1:LDQ,qcol:qcol+m) = CuM(1:LDQ,1:m+1)
    end if

    if (compZ) then
      call ZGEMM('N','N',LDZ,m,m,cone,Z(1,zcol),LDZ,ZM(1,1),LDZM,czero,CuM,LDC)
      Z(1:LDZ,zcol:zcol+m-1) = CuM(1:LDZ,1:m)
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)

  end subroutine

  subroutine z_set_last_m_poles(alpha,beta,A,B,compAB,apC,apT,&
                                Q,compQ,qcol,Z,compZ,zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  !
  ! Sets the last m poles of (A,B) ending at position (stpidx,stpidx-1) to
  ! the m poles as specified by the arrays alpha./beta
  ! If compAB is true, the complete pencil is updated, if not, only the active
  ! part of the pencil is updated
  ! If compQ is true, the matrix Q is updated accordingly. If compZ is true,
  ! the matrix Z is updated accordingly
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   complex double array [IN]
  !           Defines the new last pole: alpha/beta
  ! beta    complex double array [IN]
  !           Defines the new last pole: alpha/beta
  !           alpha and beta are expected to be of the same
  !           length and this length m are the number of poles
  !           to be introduced
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       complex double array [INOUT]
  !           Unitary equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+m
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+m-1
  !___________________________________________________________________________
  ! last edit: July 24, 2019
    complex(kind=dp), intent(inout)   :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    complex(kind=dp), intent(in)      :: alpha(:), beta(:)
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT
    integer, intent(in)               :: qcol, zcol

    integer                           :: m, i, j, LDAB, LDQ, LDQM, LDZ,&
                                         LDZM, LDC, LDR
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external ZGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    m = size(alpha) ! TODO
    call z_getLeadingDim(LDAB,A)
    call z_getLeadingDim(LDQ,Q)
    call z_getLeadingDim(LDQM,QM)
    call z_getLeadingDim(LDZ,Z)
    call z_getLeadingDim(LDZM,ZM)
    call z_getLeadingDim(LDC,CuM)
    call z_getLeadingDim(LDR,RuM)

    call tApInit(apCI,cstp-m, cstp)
    call tApInit(apTI,cstp-m, cstp)
    ! Initialize variables
    QM = czero
    do i=1, m
      QM(i,i) = cone
    end do
    ZM = czero
    do i=1, m+1
      ZM(i,i) = cone
    end do

    do i = m, 1, -1
      ! Change the last pole
      call z_set_last_pole(alpha(i), beta(i), A, B, .false., apCI, apTI,&
                           ZM, .true., m)
      ! Swap it to the correct position
      do j=1,i-1
          call z_swap(cstp-j,cstp-j-1, A, B, .false., apCI, apTI, &
                      QM, .true., m-j, ZM, .true.,m-j)
      end do
    end do

    ! Update the rest of the pencil and Q,Z matrices via BLAS
    if (compAB) then
      ! A column update
      call ZGEMM('N','N',cstp-m-tstrt,m+1,m+1,cone,A(tstrt+1,cstp-m),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      A(tstrt+1:cstp-m,cstp-m:cstp) = CuM(1:cstp-m-tstrt,1:m+1)
      ! B column update
      call ZGEMM('N','N',cstp-m-tstrt,m+1,m+1,cone,B(tstrt+1,cstp-m),LDAB,ZM(1,1),LDZM,&
                 czero,CuM(1,1),LDC)
      B(tstrt+1:cstp-m,cstp-m:cstp) = CuM(1:cstp-m-tstrt,1:m+1)
      if (cstp .lt. tstp) then
        ! A row update
        call ZGEMM('C','N',m,tstp-cstp,m,cone,QM(1,1),LDQM,A(cstp-m+1,cstp+1),&
                   LDAB,czero,RuM(1,1),LDR)
        A(cstp-m+1:cstp,cstp+1:tstp) = RuM(1:m,1:tstp-cstp)
        ! B row update
        call ZGEMM('C','N',m,tstp-cstp,m,cone,QM(1,1),LDQM,B(cstp-m+1,cstp+1),&
                   LDAB,czero,RuM(1,1),LDR)
        B(cstp-m+1:cstp,cstp+1:tstp) = RuM(1:m,1:tstp-cstp)
      end if
    else
      ! A column update
      call ZGEMM('N','N',cstp-cstrt-m,m+1,m+1,cone,A(cstrt+1,cstp-m),&
                 LDAB,ZM(1,1),LDZM,czero,CuM(1,1),LDC)
      A(cstrt+1:cstp-m,cstp-m:cstp) = CuM(1:cstp-cstrt-m,1:m+1)
      ! B column update
      call ZGEMM('N','N',cstp-cstrt-m,m+1,m+1,cone,B(cstrt+1,cstp-m),&
                 LDAB,ZM(1,1),LDZM,czero,CuM(1,1),LDC)
      B(cstrt+1:cstp-m,cstp-m:cstp) = CuM(1:cstp-cstrt-m,1:m+1)
    end if

    if (compQ) then
      call ZGEMM('N','N',LDQ,m,m,cone,Q(1,qcol),LDQ,QM(1,1),LDQM,czero,CuM,LDC)
      Q(1:LDQ,qcol:qcol+m-1) = CuM(1:LDQ,1:m)
    end if

    if (compZ) then
      call ZGEMM('N','N',LDZ,m+1,m+1,cone,Z(1,zcol),LDZ,ZM(1,1),LDZM,&
                 czero,CuM,LDC)
      Z(1:LDZ,zcol:zcol+m) = CuM(1:LDZ,1:m+1)
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)
  end subroutine
end module
