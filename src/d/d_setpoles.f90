! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Change poles in real block-Hessenberg pairs
! ___________________________________________________________________
module d_setpoles

use u_parameters
use u_activeparts
use d_ctransformations
use d_computepoles
use d_memorymgmt
use d_swappoles12
use d_swappoles22
use d_deflations

implicit none
private
public d_first_single_pole, d_first_double_pole, &
       d_last_single_pole, d_last_double_pole, &
       d_first_double_to_inf, d_last_double_to_inf, &
       d_first_m_poles, d_last_m_poles

contains

  subroutine d_first_single_pole(alpha, beta, A, B, compAB, apC, apT, &
                                 Q, compQ, qcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the pole of (A,B) in the first position of the current
  ! active part apC to alpha / beta.
  ! It is assumed that the first pole is currently a single pole.
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC of the pencil is updated
  ! If compQ is true, the transformation is also applied to
  ! the relevant two columns of Q specified by qcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   double [IN]
  !           Defines the new first pole: alpha/beta
  ! beta    double [IN]
  !           Defines the new first pole: alpha/beta
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
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
  ! last edit: September 5, 2018
    integer, intent(in)               ::  qcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    real(kind=dp), intent(in)         ::  alpha, beta
    logical, intent(in)               ::  compQ, compAB
    real(kind=dp), intent(inout)      ::  Q(:,:)
    type(tAp), pointer, intent(in)    ::  apC, apT

    real(kind=dp)                     :: v(2), c, s, r
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: LDQ
    real(kind=dp), external           :: ddot

    call d_getLeadingDim(LDQ,Q)
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    ! Compute and normalize the vector of interest
    v = A(cstrt+1:cstrt+2,cstrt+1) * beta - B(cstrt+1:cstrt+2,cstrt+1) * alpha
    s = ddot(2,v,1,v,1)
    v = v / sqrt(abs(s))

    ! Update pair and equivalence
    if (abs(v(2)) .gt. mp) then
      call d_compute_ct(v(1),v(2),c,s,r)
      if (compAB) then
        call d_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c,s)
        call d_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c,s)
      else
        call d_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c,s)
        call d_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c,s)
      endif

      if (compQ) then
        call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),c,-s)
      endif
    end if
  end subroutine

  subroutine d_first_double_pole(mu, A, B, compAB, apC, apT, Q, compQ, qcol, &
                Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the pole of (A,B) in the first and second position of the current
  ! active part apC to mu and conjg(mu).
  ! It is assumed that the first and second pole are currently
  ! either two consecutive single poles or a double pole.
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC of the pencil is updated
  ! If compQ is true, the transformation is also applied to
  ! the relevant two columns of Q specified by qcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! mu      complex double [IN]
  !           Defines the new first and second pole: mu and conjg(mu)
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
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
  ! last edit: January 30, 2019
    integer, intent(in)               ::  qcol, zcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    complex(kind=dp), intent(in)      ::  mu
    logical, intent(in)               ::  compQ, compZ, compAB
    real(kind=dp), intent(inout)      ::  Q(:,:), Z(:,:)
    type(tAp), pointer, intent(in)    ::  apC, apT

    real(kind=dp)                     :: c1, s1, c2, s2, c3, s3, r, l1, l2, &
                                         RS(3)
    complex(kind=dp)                  :: V(3), V2(3), M(3,2), a1, a2, &
                                         b1, b2, WORK(6)
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: INFO, LDQ
    logical                           :: isreal
    real(kind=dp), external           :: DZNRM2
    external                          :: ZGEMV, ZGELS

    call d_getLeadingDim(LDQ,Q)
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    if (A(cstrt+3,cstrt+1) .gt. dzero) then
      call d_eigenvalues2x2(A,B,cstrt+2,cstrt+1,l1,l2,isreal)
      if (isreal) then
        call d_make_hess(l1,cstrt+2,cstrt+1,A,B,compAB,apC,apT,Q,compQ,qcol+1,&
                Z,compZ,zcol)
        a1 = cmplx(A(cstrt+2,cstrt+1),dzero,kind=dp)
        a2 = cmplx(A(cstrt+3,cstrt+2),dzero,kind=dp)
        b1 = cmplx(B(cstrt+2,cstrt+1),dzero,kind=dp)
        b2 = cmplx(B(cstrt+3,cstrt+2),dzero,kind=dp)
      else
        a1 = cmplx(l1,l2,kind=dp)
        a2 = cmplx(l1,-l2,kind=dp)
        b1 = cone
        b2 = cone
      end if
    else
      a1 = cmplx(A(cstrt+2,cstrt+1),dzero,kind=dp)
      a2 = cmplx(A(cstrt+3,cstrt+2),dzero,kind=dp)
      b1 = cmplx(B(cstrt+2,cstrt+1),dzero,kind=dp)
      b2 = cmplx(B(cstrt+3,cstrt+2),dzero,kind=dp)
      isreal = .true.
    end if

    ! Compute the vector
    V = czero
    V(1) = cone

    ! First pole
    M = b1 * A(cstrt+1:cstrt+3,cstrt+1:cstrt+2) - &
        a1 * B(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    call ZGELS('N',3,2,1,M(1,1),3,V,3,WORK,6,INFO)
    if (INFO .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
      '[0m solving first LS system for first double pole failed.'
      write(*,*) 'Matrix: ', M
      stop
    end if
    M = A(cstrt+1:cstrt+3,cstrt+1:cstrt+2) - &
        mu * B(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    call ZGEMV('N',3,2,done,M,3,V,1,czero,V2,1)
    s1 = DZNRM2(3,V2,1)
    V = V2 / s1
    ! Second pole
    M = b2 * A(cstrt+1:cstrt+3,cstrt+1:cstrt+2) - &
        a2 * B(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    call ZGELS('N',3,2,1,M(1,1),3,V,3,WORK,6,INFO)
    if (INFO .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
      '[0m solving second system for first double pole failed.'
      write(*,*) 'Matrix: ', M
      stop
    end if
    M = A(cstrt+1:cstrt+3,cstrt+1:cstrt+2) - &
        conjg(mu) * B(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    call ZGEMV('N',3,2,done,M,3,V,1,czero,V2,1)
    s1 = DZNRM2(3,V2,1)
    V = V2 / s1

    ! Compute the transformation
    call d_compute_ct(real(v(2)),real(v(3)),c1,s1,r)
    call d_compute_ct(real(v(1)),r,c2,s2,l1)

    ! Make B block upper triangular
    RS = B(cstrt+1:cstrt+3,cstrt+1)
    l1 = c1 * RS(2) + s1 * RS(3)
    RS(3) = -s1 * RS(2) + c1 * RS(3)
    RS(2) = -s2 * RS(1) + c2 * l1
    call d_compute_ct(RS(2),RS(3),c3,s3,r)

    ! Apply the transformation
    if (compAB) then
      call d_apply_ct_l(A(cstrt+2,cstrt+1:tstp),A(cstrt+3,cstrt+1:tstp),c1,s1)
      call d_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c2,s2)
      call d_apply_ct_l(A(cstrt+2,cstrt+1:tstp),A(cstrt+3,cstrt+1:tstp),c3,s3)
      call d_apply_ct_l(B(cstrt+2,cstrt+1:tstp),B(cstrt+3,cstrt+1:tstp),c1,s1)
      call d_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c2,s2)
      call d_apply_ct_l(B(cstrt+2,cstrt+1:tstp),B(cstrt+3,cstrt+1:tstp),c3,s3)
    else
      call d_apply_ct_l(A(cstrt+2,cstrt+1:cstp),A(cstrt+3,cstrt+1:cstp),c1,s1)
      call d_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c2,s2)
      call d_apply_ct_l(A(cstrt+2,cstrt+1:cstp),A(cstrt+3,cstrt+1:cstp),c3,s3)
      call d_apply_ct_l(B(cstrt+2,cstrt+1:cstp),B(cstrt+3,cstrt+1:cstp),c1,s1)
      call d_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c2,s2)
      call d_apply_ct_l(B(cstrt+2,cstrt+1:cstp),B(cstrt+3,cstrt+1:cstp),c3,s3)
    end if
    B(cstrt+3,cstrt+1) = dzero

    if (compQ) then
      call d_apply_ct_r(Q(1:LDQ,qcol+1),Q(1:LDQ,qcol+2),c1,-s1)
      call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),c2,-s2)
      call d_apply_ct_r(Q(1:LDQ,qcol+1),Q(1:LDQ,qcol+2),c3,-s3)
    end if
  end subroutine

  subroutine d_first_double_to_inf(A, B, compAB, apC, apT, Q, compQ, qcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Shifts the first pole at the start of the current active part to infinity.
  ! This first pole is assumed to be a double pole.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
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
  !              qcol, qcol+1, qcol+2
  !___________________________________________________________________________
  ! last edit: September 13, 2018
    integer, intent(in)               ::  qcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    logical, intent(in)               ::  compQ, compAB
    real(kind=dp), intent(inout)      ::  Q(:,:)
    type(tAp), pointer, intent(in)    ::  apC, apT

    real(kind=dp)                     :: c1, s1, c2, s2, r
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: LDQ

    call d_getLeadingDim(LDQ,Q)
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_compute_ct(B(cstrt+1,cstrt+1),B(cstrt+2,cstrt+1),c1,s1,r)

    if (compAB) then
      call d_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c1,s1)
      call d_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c1,s1)
    else
      call d_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c1,s1)
      call d_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c1,s1)
    end if
    B(cstrt+2,cstrt+1) = dzero

    call d_compute_ct(A(cstrt+2,cstrt+1),A(cstrt+3,cstrt+1),c2,s2,r)

    if (compAB) then
      call d_apply_ct_l(A(cstrt+2,cstrt+1:tstp),A(cstrt+3,cstrt+1:tstp),c2,s2)
      call d_apply_ct_l(B(cstrt+2,cstrt+2:tstp),B(cstrt+3,cstrt+2:tstp),c2,s2)
    else
      call d_apply_ct_l(A(cstrt+2,cstrt+1:cstp),A(cstrt+3,cstrt+1:cstp),c2,s2)
      call d_apply_ct_l(B(cstrt+2,cstrt+2:cstp),B(cstrt+3,cstrt+2:cstp),c2,s2)
    end if
    A(cstrt+3,cstrt+1) = dzero

    if (compQ) then
      call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),c1,-s1)
      call d_apply_ct_r(Q(1:LDQ,qcol+1),Q(1:LDQ,qcol+2),c2,-s2)
    end if

  end subroutine

  subroutine d_last_single_pole(alpha, beta, A, B, compAB, apC, apT, &
                                Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the pole of (A,B) in the last position of apC to
  ! alpha / beta
  ! It is assumed that this last pole is currently a single
  ! pole.
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC is updated
  ! If compZ is true, the ctransformation is also applied to
  ! the columns of Z specified by zcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! alpha   double [IN]
  !           Defines the new last pole: alpha/beta
  ! beta    double [IN]
  !           Defines the new last pole: alpha/beta
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
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    integer, intent(in)               ::  zcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    real(kind=dp), intent(in)         ::  alpha, beta
    logical, intent(in)               ::  compAB, compZ
    real(kind=dp), intent(inout)      ::  Z(:,:)
    type(tAp), pointer, intent(in)    ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     ::  v(2), c, s, r
    integer,pointer                   ::  cstrt, cstp, tstrt, tstp
    real(kind=dp), external           ::  ddot

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    ! Compute and normalize the vector of interest
    v = A(cstp,cstp-1:cstp) * beta - B(cstp,cstp-1:cstp) * alpha
    s = ddot(2,v,1,v,1)
    v = v / sqrt(abs(s))

    if (abs(v(2)) > mp) then
      call d_compute_ct(v(2),v(1),c,s,r)
      if (compAB) then
          call d_apply_ct_r(A(tstrt+1:cstp,cstp-1),A(tstrt+1:cstp,cstp),c,s)
          call d_apply_ct_r(B(tstrt+1:cstp,cstp-1),B(tstrt+1:cstp,cstp),c,s)
      else
          call d_apply_ct_r(A(cstrt+1:cstp,cstp-1),A(cstrt+1:cstp,cstp),c,s)
          call d_apply_ct_r(B(cstrt+1:cstp,cstp-1),B(cstrt+1:cstp,cstp),c,s)
      endif

      if (compZ) then
        call d_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c,s)
      end if
    end if
  end subroutine

  subroutine d_last_double_pole(mu, A, B, compAB, apC, apT, Q, compQ, qcol, &
               Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Sets the poles of (A,B) in the last and second-to-last position
  ! of apC to mu and conjg(mu).
  ! It is assumed that these poles are currently either consecutive
  ! single poles or a double pole.
  ! If compAB is true, the complete pencil is updated based on apT,
  ! if not only the active part apC is updated
  ! If compZ is true, the ctransformation is also applied to
  ! the columns of Z specified by zcol
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! mu      complex double [IN]
  !           Defines the new last poles: mu and conjg(mu)
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
  ! Z       double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: January 30, 2019
    integer, intent(in)               ::  zcol, qcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    complex(kind=dp), intent(in)      ::  mu
    logical, intent(in)               ::  compAB, compQ, compZ
    real(kind=dp), intent(inout)      ::  Z(:,:), Q(:,:)
    type(tAp), pointer, intent(in)    ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     :: c1, s1, c2, s2, c3, s3, r, l1, l2, &
                                         RS(3)
    complex(kind=dp)                  :: V(3), V2(3), M(2,3), a1, a2, &
                                         b1, b2, WORK(6)
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: INFO
    logical                           :: isreal
    real(kind=dp), external           :: DZNRM2
    external                          :: ZGEMV, ZGELS

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    if (abs(A(cstp,cstp-2)) .gt. dzero) then
      call d_eigenvalues2x2(A,B,cstp-1,cstp-2,l1,l2,isreal)
      if (isreal) then
        call d_make_hess(l1,cstp-1,cstp-2,A,B,compAB,apC,apT,Q,compQ,qcol,&
               Z,compZ,zcol)
        a1 = cmplx(A(cstp,cstp-1),dzero,kind=dp)
        a2 = cmplx(A(cstp-1,cstp-2),dzero,kind=dp)
        b1 = cmplx(B(cstp,cstp-1),dzero,kind=dp)
        b2 = cmplx(B(cstp-1,cstp-2),dzero,kind=dp)
      else
        a1 = cmplx(l1,l2,kind=dp)
        a2 = cmplx(l1,-l2,kind=dp)
        b1 = cone
        b2 = cone
      end if
    else
      a1 = cmplx(A(cstp,cstp-1),dzero,kind=dp)
      a2 = cmplx(A(cstp-1,cstp-2),dzero,kind=dp)
      b1 = cmplx(B(cstp,cstp-1),dzero,kind=dp)
      b2 = cmplx(B(cstp-1,cstp-2),dzero,kind=dp)
      isreal = .true.
    end if

    ! Compute the vector
    V = czero
    V(3) = cone

    ! First pole
    M = b1 * A(cstp-1:cstp,cstp-2:cstp) - &
        a1 * B(cstp-1:cstp,cstp-2:cstp)
    call ZGELS('C',2,3,1,M(1,1),2,V,3,WORK,6,INFO)
    if (INFO .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
      '[0m solving LS system for last double pole failed.'
      write(*,*) 'Matrix: ', M
      stop
    end if
    M = A(cstp-1:cstp,cstp-2:cstp) - &
        mu * B(cstp-1:cstp,cstp-2:cstp)
    call ZGEMV('C',2,3,done,M(1,1),2,V(1),1,czero,V2(1),1)
    s1 = DZNRM2(3,V2,1)
    V = V2 / s1
    ! Second pole
    M = b2 * A(cstp-1:cstp,cstp-2:cstp) - &
        a2 * B(cstp-1:cstp,cstp-2:cstp)
    call ZGELS('C',2,3,1,M(1,1),2,V,3,WORK,6,INFO)
    if (INFO .ne. 0) then
      write(*,*) ' '//achar(27)//'[93m warning '//achar(27)//&
      '[0m solving LS system for last double pole failed.'
      write(*,*) 'Matrix: ', M
      stop
    end if
    M = A(cstp-1:cstp,cstp-2:cstp) - &
        conjg(mu) * B(cstp-1:cstp,cstp-2:cstp)
    call ZGEMV('C',2,3,done,M(1,1),2,V(1),1,czero,V2(1),1)
    s1 = DZNRM2(3,V2,1)
    V = V2 / s1
    ! Compute the transformation
    call d_compute_ct(real(V(2)),real(V(1)),c1,s1,r)
    call d_compute_ct(real(V(3)),r,c2,s2,l1)

    ! Make B block upper triangular
    RS = B(cstp,cstp-2:cstp)
    l1 = s1 * RS(1) + c1 * RS(2)
    RS(1) = c1 * RS(1) - s1 * RS(2)
    RS(2) = c2 * l1 - s2 * RS(3)
    call d_compute_ct(RS(2),RS(1),c3,s3,r)

    ! Apply the transformation
    if (compAB) then
      call d_apply_ct_r(A(tstrt+1:cstp,cstp-2),A(tstrt+1:cstp,cstp-1),c1,s1)
      call d_apply_ct_r(A(tstrt+1:cstp,cstp-1),A(tstrt+1:cstp,cstp),c2,s2)
      call d_apply_ct_r(A(tstrt+1:cstp,cstp-2),A(tstrt+1:cstp,cstp-1),c3,s3)
      call d_apply_ct_r(B(tstrt+1:cstp,cstp-2),B(tstrt+1:cstp,cstp-1),c1,s1)
      call d_apply_ct_r(B(tstrt+1:cstp,cstp-1),B(tstrt+1:cstp,cstp),c2,s2)
      call d_apply_ct_r(B(tstrt+1:cstp,cstp-2),B(tstrt+1:cstp,cstp-1),c3,s3)
    else
      call d_apply_ct_r(A(cstrt+1:cstp,cstp-2),A(cstrt+1:cstp,cstp-1),c1,s1)
      call d_apply_ct_r(A(cstrt+1:cstp,cstp-1),A(cstrt+1:cstp,cstp),c2,s2)
      call d_apply_ct_r(A(cstrt+1:cstp,cstp-2),A(cstrt+1:cstp,cstp-1),c3,s3)
      call d_apply_ct_r(B(cstrt+1:cstp,cstp-2),B(cstrt+1:cstp,cstp-1),c1,s1)
      call d_apply_ct_r(B(cstrt+1:cstp,cstp-1),B(cstrt+1:cstp,cstp),c2,s2)
      call d_apply_ct_r(B(cstrt+1:cstp,cstp-2),B(cstrt+1:cstp,cstp-1),c3,s3)
    end if
    B(cstp,cstp-2) = dzero

    if (compZ) then
      call d_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c1,s1)
      call d_apply_ct_r(Z(:,zcol+1),Z(:,zcol+2),c2,s2)
      call d_apply_ct_r(Z(:,zcol),Z(:,zcol+1),c3,s3)
    end if

  end subroutine

  subroutine d_last_double_to_inf( A, B, compAB, apC, apT, Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Shifts the last pole at the end of the current active part to infinity.
  ! This last pole is assumed to be a double pole.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the current active part
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zqcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1, zcol+2
  !___________________________________________________________________________
  ! last edit: September 13, 2018
    integer, intent(in)               ::  zcol
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:)
    logical, intent(in)               ::  compZ, compAB
    real(kind=dp), intent(inout)      ::  Z(:,:)
    type(tAp), pointer, intent(in)    ::  apC, apT

    real(kind=dp)                     :: c1, s1, c2, s2, r
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: LDZ

    call d_getLeadingDim(LDZ,Z)
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_compute_ct(B(cstp,cstp),B(cstp,cstp-1),c1,s1,r)

    if (compAB) then
      call d_apply_ct_r(A(tstrt+1:cstp,cstp-1),A(tstrt+1:cstp,cstp),c1,s1)
      call d_apply_ct_r(B(tstrt+1:cstp,cstp-1),B(tstrt+1:cstp,cstp),c1,s1)
    else
      call d_apply_ct_r(A(cstrt+1:cstp,cstp-1),A(cstrt+1:cstp,cstp),c1,s1)
      call d_apply_ct_r(B(cstrt+1:cstp,cstp-1),B(cstrt+1:cstp,cstp),c1,s1)
    end if
    B(cstp,cstp-1) = dzero

    call d_compute_ct(A(cstp,cstp-1),A(cstp,cstp-2),c2,s2,r)

    if (compAB) then
      call d_apply_ct_r(A(tstrt+1:cstp,cstp-2),A(tstrt+1:cstp,cstp-1),c2,s2)
      call d_apply_ct_r(B(tstrt+1:cstp,cstp-2),B(tstrt+1:cstp,cstp-1),c2,s2)
    else
      call d_apply_ct_r(A(cstrt+1:cstp,cstp-2),A(cstrt+1:cstp,cstp-1),c2,s2)
      call d_apply_ct_r(B(cstrt+1:cstp,cstp-2),B(cstrt+1:cstp,cstp-1),c2,s2)
    end if
    A(cstp,cstp-2) = dzero

    if (compZ) then
      call d_apply_ct_r(Z(1:LDZ,zcol+1),Z(1:LDZ,zcol+2),c1,s1)
      call d_apply_ct_r(Z(1:LDZ,zcol),Z(1:LDZ,zcol+1),c2,s2)
    end if
  end subroutine

  subroutine d_first_m_poles(m,alpha,beta,k,mur,mui,l,A,B,compAB,apC, apT,&
                Q, compQ, qcol, Z, compZ, zcol )
  ! DESCRIPTION
  !___________________________________________________________________________
  !
  ! Sets the first k+2*l(+1) poles of (A,B) starting at the first position of
  ! the current active part apC to the poles as specified by the arrays:
  !                    k: alpha./beta, and, 2*l: mur +-1i* mui
  ! If compAB is true, the complete pencil is updated based on apT, if not,
  ! only the active part apC of the pencil is updated
  ! If compQ is true, the matrix Q is updated accordingly. If compZ is true,
  ! the matrix Z is updated accordingly
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! m       integer [OUT]
  !           Total number of shifts that was introduced: m = k + 2*l (+1)
  !           Can be 1 larger than requested IF it did not match
  ! alpha   double array [IN]
  !           Defines k real shifts: alpha/beta
  !           Array with at least k elements
  ! beta    double array [IN]
  !           Defines k real shifts: alpha/beta
  !           Array with at least k elements
  ! k       integer [IN]
  !           Number of real shifts
  ! mur     double array [IN]
  !           Real parts of l pairs of C.C. shifts: mur +- 1i * muc
  !           Array with at least l elements
  ! mui     double array [IN]
  !           Imaginary parts of l pairs of C.C. shifts: mur +- 1i * muc
  !           Array with at least l elements
  ! l       integer [IN]
  !           Number of C.C. pairs to be introduced
  !           (total of 2*l complex shifts)
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       double array [INOUT]
  !           Orthonormal equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+m+1
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+m
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    real(kind=dp), intent(in)         :: alpha(:), beta(:),mur(:),mui(:)
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT
    integer, intent(in)               :: qcol, zcol
    integer, intent(inout)            :: m, k, l

    integer                           :: i, j, LDAB, LDQ, LDQM, LDZ, &
                                         LDZM, LDC, LDR
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    m = k + 2*l
    if (abs(A(cstrt+m+2,cstrt+m)) .gt. dzero) then
      m = m + 1
    end if

    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDQM,QM)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDZM,ZM)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)

    ! Initialize internal variables
    call tApInit(apCI,cstrt, cstrt+m)
    call tApInit(apTI,cstrt, cstrt+m)

    QM = dzero
    do i=1, m+1
        QM(i,i) = done
    end do

    ZM = dzero
    do i=1, m
        ZM(i,i) = done
    end do

    do while (k+l .gt. 0)
      ! Check what the current first pole is
      if (abs(A(cstrt+3,cstrt+1)) .gt. dzero) then
        ! Current first is c.c.
        if (l .gt. 0) then
          ! We can introduce c.c
          call d_first_double_pole(cmplx(mur(l),mui(l),kind=dp), A, B, .true., &
                 apCI, apTI, QM, .true., 1, ZM, .true., 1)
          l = l - 1
          j = 1
          do while (j .le. k+2*l)
            if (abs(A(cstrt+j+4,cstrt+j+2)) .gt. dzero) then
              call d_swap22(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 2
            else
              call d_swap21(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 1
            end if
          end do
        else
          ! Shift to Inf
          call d_first_double_to_inf(A, B, .true., apCI, apTI, QM, .true., 1)
          call d_first_single_pole(alpha(k), beta(k), A, B, .true., apCI, apTI,&
                              QM, .true., 1)
          call d_swap11(cstrt+2, cstrt+1, A, B, .true., apCI, apTI, &
                   QM, .true., 2, ZM, .true., 1)
          if (k .gt. 1) then
            ! introduce two different reals
            call d_first_single_pole(alpha(k-1), beta(k-1), A, B, .true., &
                   apCI, apTI, QM, .true., 1)
            k = k - 2
          else
            ! introduce twice the same real (should be the last)
            call d_first_single_pole(alpha(k), beta(k), A, B, .true., &
                    apCI, apTI, QM, .true., 1)
            k = k - 1
          end if
          j = 1
          do while (j .le. k+2*l)
            if (abs(A(cstrt+j+4,cstrt+j+2)) .gt. dzero) then
              call d_swap12(cstrt+j+2, cstrt+j+1, A, B, .true., apCI, apTI, &
                     QM, .true., j+2, ZM, .true., j+1)
              call d_swap12(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 2
            else
              call d_swap11(cstrt+j+2, cstrt+j+1, A, B, .true., apCI, apTI, &
                     QM, .true., j+2, ZM, .true., j+1)
              call d_swap11(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 1
            end if
          end do
        end if
      else
        ! Current first is real
        if (k .gt. 0) then
          ! Set real
          call d_first_single_pole(alpha(k), beta(k), A, B, .true., apCI, apTI,&
                                   QM, .true., 1)
          ! Swap to the back
          j = 1
          do while (j .lt. k+2*l)
            if (abs(A(cstrt+j+3,cstrt+j+1)) .gt. dzero) then
              call d_swap12(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 2
            else
              call d_swap11(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 1
            end if
          end do
          k = k - 1 !Update counter
        else
          ! l is greater than 0 (otherwise this is never reached)
          if (abs(A(cstrt+4,cstrt+2)) .gt. dzero) then
            call d_swap12(cstrt+2, cstrt+1, A, B, .true., apCI, apTI, &
                   QM, .true.,2, ZM, .true., 1)
          end if
          ! Set cc
          call d_first_double_pole(cmplx(mur(l),mui(l),kind=dp),A, B, .true., &
                 apCI,apTI,QM, .true. ,1,ZM,.true.,1)
          j = 1

          do while (j .lt. 2*l-1)
            if (abs(A(cstrt+j+4,cstrt+j+2)) .gt. dzero) then
              call d_swap22(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 2
            else
              call d_swap21(cstrt+j+1, cstrt+j, A, B, .true., apCI, apTI, &
                     QM, .true., j+1, ZM, .true., j)
              j = j + 1
            end if
          end do
          l = l - 1 !Update counter
        end if
      end if
    end do

    ! Update the remainder of the pencil
    if (compAB) then
      ! A row update
      call DGEMM('T','N',m+1,tstp-cstrt-m,m+1,done,QM(1,1),LDQM,&
                 A(cstrt+1,cstrt+m+1),LDAB,dzero,Ru(1,1),LDR)
      A(cstrt+1:cstrt+m+1,cstrt+m+1:tstp) = Ru(1:m+1,1:tstp-cstrt-m)
      ! B row update
      call DGEMM('T','N',m+1,tstp-cstrt-m,m+1,done,QM(1,1),LDQM,&
                 B(cstrt+1,cstrt+m+1),LDAB,dzero,Ru(1,1),LDR)
      B(cstrt+1:cstrt+m+1,cstrt+m+1:tstp) = Ru(1:m+1,1:tstp-cstrt-m)
      if (cstrt .gt. 0) then
        ! A column update
        call DGEMM('N','N',cstrt-tstrt,m,m,done,A(tstrt+1,cstrt+1),LDAB,&
                   ZM(1,1),LDZM,dzero,Cu(1,1),LDC)
        A(tstrt+1:cstrt,cstrt+1:cstrt+m) = Cu(1:cstrt-tstrt,1:m)
        ! B column update
        call DGEMM('N','N',cstrt-tstrt,m,m,done,B(tstrt+1,cstrt+1),LDAB,&
                   ZM(1,1),LDZM,dzero,Cu(1,1),LDC)
        B(tstrt+1:cstrt,cstrt+1:cstrt+m) = Cu(1:cstrt-tstrt,1:m)
      end if
    else
      ! A row update
      call DGEMM('N','N',m+1,cstp-cstrt-m,m+1,done,QM(1,1),LDQM,&
                 A(cstrt+1,cstrt+m+1),LDAB,dzero,Ru(1,1),LDR)
      A(cstrt+1:cstrt+m+1,cstrt+m+1:cstp) = Ru(1:m+1,1:cstp-cstrt-m)
      ! B row update
      call DGEMM('N','N',m+1,cstp-cstrt-m,m+1,done,QM(1,1),LDQM,&
                 B(cstrt+1,cstrt+m+1),LDAB,dzero,Ru(1,1),LDR)
      B(cstrt+1:cstrt+m+1,cstrt+m+1:cstp) = Ru(1:m+1,1:cstp-cstrt-m)
    end if

    if (compQ) then
      call DGEMM('N','N',LDQ,m+1,m+1,done,Q(1,qcol),LDQ,&
                 QM(1,1),LDQM,dzero,Cu,LDC)
      Q(1:LDQ,qcol:qcol+m) = Cu(1:LDQ,1:m+1)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,m,m,done,Z(1,zcol),LDZ,ZM(1,1),LDZM,dzero,Cu,LDC)
      Z(1:LDZ,zcol:zcol+m-1) = Cu(1:LDZ,1:m)
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)
  end subroutine

  subroutine d_last_m_poles(m,bf,alpha,beta,k,mur,mui,l,A,B,compAB,apC, apT, &
               Q, compQ, qcol, Z, compZ, zcol )
  ! DESCRIPTION
  !___________________________________________________________________________
  !
  ! Sets the last m=k+2*l poles of (A,B) starting at the last position of
  ! the current active part apC to the poles as specified by the arrays:
  !                    k: alpha./beta, and, 2*l: mur +-1i* mui
  ! If compAB is true, the complete pencil is updated, if not, only the active
  ! part of the pencil is updated
  ! If compQ is true, the matrix Q is updated accordingly. If compZ is true,
  ! the matrix Z is updated accordingly
  ! Deflations are monitored during the pole introduction. Active part can
  ! shift.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! m       integer [OUT]
  !           Total number of shifts that was introduced: m = k + 2*l (+1)
  !           Can be 1 larger than requested IF it did not match
  ! bf      integer [IN]
  !           Buffer size to compensate for deflations that might occur during
  !           pole introductions. QM, ZM should be of dimension m+bf+1
  ! alpha   double array [IN]
  !           Defines k real shifts: alpha/beta
  !           Array with at least k elements
  ! beta    double array [IN]
  !           Defines k real shifts: alpha/beta
  !           Array with at least k elements
  ! k       integer [IN]
  !           Number of real shifts
  ! mur     double array [IN]
  !           Real parts of l pairs of C.C. shifts: mur +- 1i * muc
  !           Array with at least l elements
  ! mui     double array [IN]
  !           Imaginary parts of l pairs of C.C. shifts: mur +- 1i * muc
  !           Array with at least l elements
  ! l       integer [IN]
  !           Number of C.C. pairs to be introduced
  !           (total of 2*l complex shifts)
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [IN]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       double array [INOUT]
  !           Orthonormal equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN] ! TODO qcol should point at the last one
  !           Start index of the columns of Q that are updated:
  !              qcol, qcol+m+bf-1
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+m+bf
  !___________________________________________________________________________
  ! last edit: February 20, 2019
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    real(kind=dp), intent(in)         :: alpha(:), beta(:),mur(:),mui(:)
    logical, intent(in)               :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)    :: apC, apT
    integer, intent(in)               :: qcol, zcol
    integer, intent(inout)            :: m, k, l, bf

    logical                           :: defl, defli
    integer                           :: i, LDAB, LDQ, LDQM, LDZ, &
                                         LDZM, LDC, LDR, mi, &
                                         cstpe, cendm, fendm, &
                                         ndef, nbuf, nuc, nc, sz, cdif
    integer,pointer                   :: cstrt, cstp, tstrt, tstp, &
                                         cstrti, cstpi
    type(tAp), pointer                :: apCI, apTI

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

    m = k + 2*l
    ! Initialize internal variables
    call tApInit(apCI,cstp-m-bf, cstp)
    call tApInit(apTI,cstp-m-bf, cstp)
    call tApGet(apCI,cstrti,cstpi)
    cstpe = cstpi

    ! Set entire matrices to identity
    QM = dzero
    do i=1, LDQM
        QM(i,i) = done
    end do
    ZM = dzero
    do i=1, LDZM
        ZM(i,i) = done
    end do

    ! Ensure m covers integer number of blocks
    mi = m
    if (abs(A(cstp-m+1,cstp-m-1)) .gt. dzero) then
      m = m + 1
      bf = bf - 1
    end if
    ! Set all counters
    defl = .false.
    ndef = 0
    nbuf = bf
    nc = 0
    nuc = m
    cendm = LDZM ! current end index in QM,ZM corresponding with cstpi
    fendm = LDZM

    ! Phase 1: Introduce all l double poles
    do while (l .gt. 0)
      sz = 2
      call moveLucSzToEnd(A,B,apCI,apTI,sz,cendm,nc,nuc,nbuf)
      if (sz .gt. 0) then
        ! The two last poles are now unchanged
        call d_stop_deflation_double(defli, A, B, .true., &
             apCI, apTI, QM, .true., cendm-1, ZM, .true., cendm-2)
        if (defli) then
          defl = .true.
          cdif = cstpe-cstpi
          ndef = ndef + cdif
          nbuf = nbuf - cdif
          nuc = nuc + cdif
          cendm = cendm - cdif
          cstpe = cstpi
          if (nbuf .lt. 3) exit
          if (cdif .eq. 1) then
            nc = nc + 1
          end if
        else
          nuc = nuc - sz
          nc = nc + sz
          call d_last_double_pole(cmplx(mur(l),mui(l),kind=dp), A, B, &
                 .true., apCI, apTI, QM, .true., cendm-1, ZM, .true., cendm-2)
          l = l - 1
        end if
      else
        ! sz = - 1
        exit
      end if
    ! End of phase 1
    end do

    ! Phase 2: Introduce all k single poles
    do while (k .gt. 0)
      sz = 1
      call moveLucSzToEnd(A,B,apCI,apTI,sz,cendm,nc,nuc,nbuf)
      if (sz .gt. 0) then
        if (sz .eq. 1) then
          call d_stop_deflation_single(defli, A, B, .true., &
                  apCI, apTI, ZM, .true., cendm-1)
          if (defli) then
            defl = .true.
            cdif = cstpe-cstpi
            ndef = ndef + cdif
            nbuf = nbuf - cdif
            nuc = nuc + cdif
            cendm = cendm - cdif
            cstpe = cstpi
            if (nbuf .lt. 2) exit
            if (cdif .eq. 2) then
              nc = nc - 1
            end if
          else
            nuc = nuc - 1
            nc = nc + 1
            call d_last_single_pole(alpha(k), beta(k), A, B, .true., &
                 apCI, apTI, ZM, .true., cendm-1)
            k = k - 1
          end if
        else
          ! sz = 2
          call d_stop_deflation_double(defli, A, B, .true., &
                apCI, apTI, QM, .true., cendm-1, ZM, .true., cendm-2)
          if (defli) then
            defl = .true.
            cdif = cstpe-cstpi
            ndef = ndef + cdif
            nbuf = nbuf - cdif
            nuc = nuc + cdif
            cendm = cendm - cdif
            cstpe = cstpi
            if (nbuf .lt. 2) exit
            if (cdif .eq. 1) then
              nc = nc + 1
            end if
          else
            nuc = nuc - 2
            nc = nc + 2
            ! Can be a double pole
            if (abs(A(cstpi,cstpi-2)) .gt. dzero) then
              call d_last_double_to_inf( A, B, .true., apCI, apTI, &
                   ZM, .true., cendm-2)
            end if
            call d_last_single_pole(alpha(k), beta(k), A, B, .true., &
               apCI, apTI, ZM, .true., cendm-1)
            k = k - 1
            call d_swap11(cstpi-1, cstpi-2, A, B, .true., apCI, apTI, &
               QM, .true., cendm-1, ZM, .true., cendm-2)
            if (k .eq. 0) k = 1
            call d_last_single_pole(alpha(k), beta(k), A, B, .true., &
               apCI, apTI, ZM, .true., cendm-1)
            k = k - 1
          end if
        end if
      else
        ! sz = -1
        exit
      end if
    end do

    ! Update the rest of the pencil and Q,Z matrices via BLAS
    if (compAB) then
      ! A column update
      call DGEMM('N','N',cstp-m-bf-tstrt,m+ndef+1,m+ndef+1,done,A(tstrt+1,cstp-m-ndef),LDAB,&
            ZM(LDZM-m-ndef,LDZM-m-ndef), LDZM, dzero, Cu(1,1), LDC)
      A(tstrt+1:cstp-m-bf,cstp-m-ndef:cstp) = Cu(1:cstp-m-bf-tstrt,1:m+ndef+1)
      ! B column update
      call DGEMM('N','N',cstp-m-bf-tstrt,m+ndef+1,m+ndef+1,done,B(tstrt+1,cstp-m-ndef),LDAB,&
            ZM(LDZM-m-ndef,LDZM-m-ndef), LDZM, dzero, Cu(1,1), LDC)
      B(tstrt+1:cstp-m-bf,cstp-m-ndef:cstp) = Cu(1:cstp-m-bf-tstrt,1:m+ndef+1)
      if (cstp .lt. tstp) then
        ! A row update
        call DGEMM('T','N', m+ndef, tstp-cstp, m+ndef, done, QM(LDQM-m-ndef+1,LDQM-m-ndef+1), LDQM,&
              A(cstp-m-ndef+1,cstp+1), LDAB, dzero, Ru(1,1), LDR)
        A(cstp-m-ndef+1:cstp,cstp+1:tstp) = Ru(1:m+ndef,1:tstp-cstp)
        ! B row update
        call DGEMM('T','N',m+ndef,tstp-cstp,m+ndef,done,QM(LDQM-m-ndef+1,LDQM-m-ndef+1),LDQM,&
              B(cstp-m-ndef+1,cstp+1), LDAB, dzero, Ru(1,1), LDR)
        B(cstp-m-ndef+1:cstp,cstp+1:tstp) = Ru(1:m+ndef,1:tstp-cstp)
      end if
    else
      ! A column update
      call DGEMM('N','N',cstp-m-bf-cstrt,m+ndef+1,m+ndef+1,done,A(cstrt+1,cstp-m-ndef),LDAB,&
            ZM(LDZM-m-ndef,LDZM-m-ndef), LDZM, dzero, Cu(1,1), LDC)
      A(cstrt+1:cstp-m-bf,cstp-m-ndef:cstp) = Cu(1:cstp-m-bf-cstrt,1:m+ndef+1)
      ! B column update
      call DGEMM('N','N',cstp-m-bf-cstrt,m+ndef+1,m+ndef+1,done,B(cstrt+1,cstp-m-ndef),LDAB,&
            ZM(LDZM-m-ndef,LDZM-m-ndef), LDZM, dzero, Cu(1,1), LDC)
      B(cstrt+1:cstp-m-bf,cstp-m-ndef:cstp) = Cu(1:cstp-m-bf-cstrt,1:m+ndef+1)
    end if

    if (compQ) then
      call DGEMM('N','N',LDQ,m+ndef,m+ndef,done,Q(1,qcol-m-ndef+1),LDQ,QM(LDQM-m-ndef+1,LDQM-m-ndef+1),LDQM,dzero,Cu,LDC)
      Q(1:LDQ,qcol-m-ndef+1:qcol) = Cu(1:LDQ,1:m+ndef)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,m+ndef+1,m+ndef+1,done,Z(1,zcol-m-ndef),LDZ,ZM(LDZM-m-ndef,LDZM-m-ndef),LDZM,&
                 dzero,Cu,LDC)
      Z(1:LDZ,zcol-m-ndef:zcol) = Cu(1:LDZ,1:m+ndef+1)
    end if

    if (defl) then
      cstp = cstpi
    end if

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)

  end subroutine

  recursive subroutine moveLucSzToEnd(A,B,apCI,apTI,sz,cendm,nc,nuc,nbuf)
    ! Move last unchanged pole of the desired size
    ! to the end of the range
    !
    ! (A,B) : pencil
    ! apCI, apTI : active range (contains cstpi)
    ! sz:  specifies if we want a single or double pole
    ! cendm: current end in QM,ZM
    ! nc: number of changed poles
    ! nuc: number of unchanged poles
    ! nbuf: remaining buffer size

    real(kind=dp), intent(inout)    :: A(:,:), B(:,:)
    type(tAp), pointer, intent(in)  :: apCI, apTI
    integer, intent(in)             :: cendm, nc, nuc, nbuf
    integer, intent(inout)          :: sz

    integer :: j, jm ! index in A,B and QM, ZM
    integer :: jmin ! minimal column index in A,B
    integer, pointer :: cstrti, cstpi
    logical :: ds

    call tApGet(apCI,cstrti,cstpi)
    j = cstpi - nc - 1
    jm = cendm - nc - 1
    jmin = cstpi - nc - nuc - 1

    if ((sz .eq. 2) .and. (nuc .gt. 1)) then
      ! Look for 2 consecutive 1x1 or a 2x2 (can be found in last 3 poles)
      j = j - 1
      jm = jm - 1
      if (abs(A(j+2,j)) .gt. dzero) then
        ! At column j starts a double pole
        ds = .false.
      else
        ! At column j+1 is a single pole
        if ((nuc .gt. 2) .and. (abs(A(j+1,j-1)) .gt. dzero)) then
          ! At columns j-1:j is a double pole, so we shift one more
          j = j - 1
          jm = jm - 1
          ds = .false.
        else
          ! At column j is also a single pole
          ds = .true.
        end if
      end if
      ! Swap the identified one to the front
      do while (jm .lt. cendm-2)
        if (ds) then
          ! Swapping 2 1x1
          if (jm .lt. cendm-3) then
            if (abs(A(j+4,j+2)) .gt. dzero) then
              call d_swap12(j+2, j+1, A, B, .true., apCI, apTI, &
                   QM, .true., jm+2, ZM, .true., jm+1)
              call d_swap12(j+1, j, A, B, .true., apCI, apTI, &
                   QM, .true., jm+1, ZM, .true., jm)
              j = j + 2
              jm = jm + 2
            else
              call d_swap11(j+2, j+1, A, B, .true., apCI, apTI, &
                   QM, .true., jm+2, ZM, .true., jm+1)
              call d_swap11(j+1, j, A, B, .true., apCI, apTI, &
                   QM, .true., jm+1, ZM, .true., jm)
              j = j + 1
              jm = jm + 1
            end if
          else
            call d_swap11(j+2, j+1, A, B, .true., apCI, apTI, &
                 QM, .true., jm+2, ZM, .true., jm+1)
            call d_swap11(j+1, j, A, B, .true., apCI, apTI, &
                 QM, .true., jm+1, ZM, .true., jm)
            j = j + 1
            jm = jm + 1
          end if
        else
          ! Swapping 2x2 block
          if (jm .lt. cendm-3) then
            if (abs(A(j+4,j+2)) .gt. dzero) then
              call d_swap22(j+1, j, A, B, .true., apCI, apTI, &
                   QM, .true., jm+1, ZM, .true., jm)
              j = j + 2
              jm = jm + 2
            else
              call d_swap21(j+1, j, A, B, .true., apCI, apTI, &
                   QM, .true., jm+1, ZM, .true., jm)
              j = j + 1
              jm = jm + 1
            end if
          else
            call d_swap21(j+1, j, A, B, .true., apCI, apTI, &
                 QM, .true., jm+1, ZM, .true., jm)
            j = j + 1
            jm = jm + 1
          end if
        end if
      end do
    elseif ((sz .eq. 1) .and. (nuc .gt. 0)) then
      ! Look for 1x1 if not available a 2x2
      do while (jm .gt. jmin)
        if (abs(A(j+1,j-1)) .gt. dzero) then
          j = j - 2
          jm = jm - 2
        else
          exit
        end if
      end do
      if (abs(A(j+1,j-1)) .gt. dzero) then
        ! We acctually did not find a 1x1
        sz = 2
        call moveLucSzToEnd(A,B,apCI,apTI,sz,cendm,nc,nuc,nbuf)
        return
      end if
      ! We swap the 1x1 to the end
      do while (jm .lt. cendm-1)
        if (jm .lt. cendm-2) then
          if (abs(A(j+3,j+1)) .gt. dzero) then
            call d_swap12(j+1, j, A, B, .true., apCI, apTI, &
                  QM, .true., jm+1, ZM, .true., jm)
            j = j + 2
            jm = jm + 2
          else
            call d_swap11(j+1, j, A, B, .true., apCI, apTI, &
                  QM, .true., jm+1, ZM, .true., jm)
            j = j + 1
            jm = jm + 1
          end if
        else
          call d_swap11(j+1, j, A, B, .true., apCI, apTI, &
                QM, .true., jm+1, ZM, .true., jm)
          j = j + 1
          jm = jm + 1
        end if
      end do
    else
      sz = -1 ! We cannot find a desired option
    end if
  end subroutine

end module
