! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Deflation monitoring
! ___________________________________________________________________
module d_deflations

use u_activeparts
use u_parameters
use d_ctransformations
use d_memorymgmt
use d_computepoles

implicit none
private
public d_check_interior_deflations, d_start_deflation_single, &
       d_start_deflation_double, d_stop_deflation_single, &
       d_stop_deflation_double
contains

  subroutine d_check_interior_deflations(defl,A, B, compAB, apC, apT, &
              Q, compQ, qcol, Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Check for interior deflations in active part defined by apC.
  ! When deflated poles are encountered, the pencil is updated,
  ! i.e., these poles are explicitly set to zero, and the active
  ! parts list apC is updated.The boolean defl returns true when
  ! at least 1 deflation is found.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [IN]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [INOUT]
  !           Active parts linked list.
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
  ! last edit: January 11, 2018
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    type(tAp), pointer, intent(inout) :: apC
    type(tAp), pointer, intent(in)    :: apT
    logical, intent(out)              :: defl
    logical, intent(in)               :: compAB, compQ, compZ
    integer, intent(in)               :: qcol, zcol

    integer                           :: i, j, cpy, cpystrt
    integer, pointer                  :: cstrt, cstp, tstrt, tstp
    real(kind=dp)                     :: l1, l2
    logical                           :: isreal, defli
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    defl = .false.
    i = cstrt+1
    cpystrt = cstrt
    do while (i .le. cstp-2)
      if (abs(A(i+2,i)) .gt. dzero) then
        ! double pole
        call d_eigenvalues2x2(A,B,i+1,i,l1,l2,isreal)
        if (isreal) then
          ! Not C.C.
          call d_make_hess(l1,i+1,i,A,B,compAB,apC,apT,Q,compQ,qcol+i-cpystrt,&
                 Z,compZ,zcol+i-cpystrt-1)
          call d_interior_deflation_sng(i,defli,A,B)
          if (defli) then
            defl = .true.
            ! update active parts
            cpy = cstp
            call tApPut(apC,strt=cstrt, stp=i)
            call tApInsertAfter(apC,strt=i,stp=cpy)
            apC => tApNext(apC)
            call tApGet(apC,cstrt,cstp) ! this does update cstrt,cstp
                                        !(but does not affect the loop counter)
          end if
          i = i + 1
        else
          j = i
          call d_interior_deflation_dbl(j,defli,A,B)
          if (defli) then
            defl = .true.
            ! update active parts
            cpy = cstp
            call tApPut(apC,strt=cstrt, stp=i+j)
            call tApInsertAfter(apC,strt=i+j,stp=cpy)
            apC => tApNext(apC)
            call tApGet(apC,cstrt,cstp) ! this does update cstrt,cstp
                                        !(but does not affect the loop counter)
          end if
          i = i + 2
        end if
      else
        call d_interior_deflation_sng(i,defli,A,B)
        if (defli) then
          defl = .true.
          ! update active parts
          cpy = cstp
          call tApPut(apC,strt=cstrt, stp=i)
          call tApInsertAfter(apC,strt=i,stp=cpy)
          apC => tApNext(apC)
          call tApGet(apC,cstrt,cstp) ! this does update cstrt,cstp
                                      !(but does not affect the loop counter)
        end if
        i = i + 1
      end if
    end do
  end subroutine d_check_interior_deflations

  subroutine d_interior_deflation_sng(i,defl,A,B)
    ! DESCRIPTION
    !___________________________________________________________________________
    ! Test if the single interior pole at position (i+1,i) in (A,B)
    ! is deflated
    !
    ! ARGUMENTS
    !___________________________________________________________________________
    ! i       integer [IN]
    !           column index of single pole
    ! defl    boolean [INOUT]
    !           returns true of the ith pole can be deflated.
    !           If so, the respective entries in A,B are set to zero.
    ! A       double array [INOUT]
    !           block-Hessenberg matrix
    ! B       double array [INOUT]
    !           Hessenberg matrix
    !
    ! NOTA
    !___________________________________________________________________________
    ! This is a private subroutine called from  d_check_interior_deflations
    !___________________________________________________________________________
    ! last edit: January 11, 2019
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:)
    integer, intent(in)               :: i
    logical, intent(inout)            :: defl

    defl = .false.
    if ((abs(A(i+1,i)) .lt. mp * (abs(A(i,i)) + abs(A(i+1,i+1)))) .and. &
      (abs(B(i+1,i)) .lt. mp * (abs(B(i,i)) + abs(B(i+1,i+1))))) then
      ! Set to zero
      defl = .true.
      A(i+1,i) = dzero
      B(i+1,i) = dzero
    end if
  end subroutine

  subroutine d_interior_deflation_dbl(i,defl,A,B)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Test if the double interior pole at position (i+1:i+2,i:i+1)
  ! in (A,B) is deflated
  ! ARGUMENTS
  !___________________________________________________________________________
  ! i       integer [IN]
  !           column index of double pole
  ! defl    boolean [INOUT]
  !           returns true of the ith pole can be deflated.
  !           If so, the respective entries in A,B are set to zero.
  ! A       double array [INOUT]
  !           block-Hessenberg matrix
  ! B       double array [INOUT]
  !           Hessenberg matrix
  !
  ! NOTA
  !___________________________________________________________________________
  ! This is a private subroutine called from  d_check_interior_deflations
  !___________________________________________________________________________
  ! last edit: February 7, 2019
    real(kind=dp), intent(inout)      :: A(:,:), B(:,:)
    logical, intent(inout)            :: defl
    integer, intent(inout)            :: i

    defl = .false.
    if (((abs(A(i+2,i)) + abs(A(i+2,i+1))) .lt. mp * &
        (abs(A(i+1,i+1)) + abs(A(i+2,i+2)))) .and. &
        (abs(B(i+2,i+1)) .lt. mp * &
        (abs(B(i+1,i+1)) + abs(B(i+2,i+2))))) then
        defl = .true.
        B(i+2,i:i+1) = dzero
        A(i+2,i:i+1) = dzero
        i = 1 ! Ensures that active parts are handled in calling function
    elseif (((abs(A(i+1,i)) + abs(A(i+2,i))) .lt. mp * &
        (abs(A(i,i)) + abs(A(i+1,i+1)))) .and. &
        (abs(B(i+1,i)) .lt. mp * &
        (abs(B(i,i)) + abs(B(i+1,i+1))))) then
        defl = .true.
        B(i+1:i+2,i) = dzero
        A(i+1:i+2,i) = dzero
        i = 0 ! Ensures that active parts are handled in calling function
      end if
  end subroutine

  subroutine d_start_deflation_single(defl, A, B, compAB, apC, apT, &
                                     Q, compQ, qcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Check for deflation at the start of the active area
  ! The pencil (A,B) and the Q matrix is updated as requested if a deflation
  ! is found.
  ! The current active part apC is also updated, which also updates
  ! it in the caller function.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator (true if deflated, else false)
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B) as defined by
  !           apT. If not, then only to the active part apC
  ! apC     type(tAp) [INOUT]
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
  ! last edit: September 12, 2018
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:)
    logical, intent(in)               ::  compAB, compQ
    integer, intent(in)               ::  qcol
    logical, intent(out)              ::  defl
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     :: c,s,r,t, tt, u(2)
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: LDQ

    call d_getLeadingDim(LDQ,Q)

    ! Get current and total ranges
    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    defl = .false.
    ! first check for direct deflation
    if ((abs(A(cstrt+2,cstrt+1)) .lt. &
         mp*(abs(A(cstrt+1,cstrt+1)) + abs(A(cstrt+2,cstrt+2)))) .and. &
        (abs(B(cstrt+2,cstrt+1)) .lt. &
         mp*(abs(B(cstrt+1,cstrt+1)) + abs(B(cstrt+2,cstrt+2))))) then
        defl = .true.
        A(cstrt+2,cstrt+1) = dzero
        B(cstrt+2,cstrt+1) = dzero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt+1,stp=cstp)
    elseif (abs(A(cstrt+1,cstrt+1) * B(cstrt+2,cstrt+1) - &
                A(cstrt+2,cstrt+1) * B(cstrt+1,cstrt+1)) .lt. &
           mp*(abs(A(cstrt+2,cstrt+2)) + abs(B(cstrt+2,cstrt+2)))) then
      if (abs(A(cstrt+2,cstrt+1)) .ge. abs(B(cstrt+2,cstrt+1))) then
        ! Compute rotation based on A
        call d_compute_ct(A(cstrt+1,cstrt+1),A(cstrt+2,cstrt+1),c,s,r)
        ! Apply to B
        r = B(cstrt+1,cstrt+1)
        t = B(cstrt+2,cstrt+1)
        u = B(cstrt+1:cstrt+2,cstrt+2)
      else
        ! Compute rotation based on B
        call d_compute_ct(B(cstrt+1,cstrt+1),B(cstrt+2,cstrt+1),c,s,r)
        ! Apply to A
        r = A(cstrt+1,cstrt+1)
        t = A(cstrt+2,cstrt+1)
        u = A(cstrt+1:cstrt+2,cstrt+2)
      end if

      tt = -s * r + c * t
      r = c * r + s * t
      t = tt
      u(2) = -s * u(1) + c * u(2)

      if (abs(t) .lt. mp*(abs(r) + abs(u(2)))) then
        defl = .true.
        if (compAB) then
          call d_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c,s)
          call d_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c,s)
        else
          call d_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c,s)
          call d_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c,s)
        end if
        if (compQ) then
          call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),c,-s)
        end if
        A(cstrt+2,cstrt+1) = dzero
        B(cstrt+2,cstrt+1) = dzero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt+1,stp=cstp)
      end if
    end if

    ! Double check if a 2x2 block has separated
    if (cstp-cstrt .gt. 3) then
      if (.not. (defl .or. (abs(A(cstrt+4,cstrt+2)) .gt. dzero))) then
      call d_interior_deflation_sng(cstrt+2,defl,A,B)
      if (defl) then
        ! update active parts
        call tApInsertAfter(apC,strt=cstrt+2,stp=cstp)
        call tApPut(apC,strt=cstrt, stp=cstrt+2)
      end if
     end if
    end if
  end subroutine d_start_deflation_single

  subroutine d_start_deflation_double(defl, A, B, compAB, apC, apT, &
                                     Q, compQ, qcol, Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Check for deflation at the start of the active area. It is assumed
  ! that there is a double pole at the start. It is also checked if the double
  ! pole actually is a complex conjugate pair. If not, the Hessenberg form
  ! is restored and a single deflation is checked.
  ! The pencil (A,B) and the Q and Z matrices are updated as requested
  ! if a deflation is found.
  ! The current active part apC is also updated, which also updates
  ! it in the caller function.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator (true if deflated, else false)
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B) as defined by
  !           apT. If not, then only to the active part apC
  ! apC     type(tAp) [INOUT]
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
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: September 12, 2018
  real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
  logical, intent(in)               ::  compAB, compQ, compZ
  integer, intent(in)               ::  qcol, zcol
  logical, intent(out)              ::  defl
  type(tAp), pointer, intent(inout) ::  apC
  type(tAp), pointer, intent(in)    ::  apT

  real(kind=dp)                     :: c,c2,s,s2,r,t
  logical                           :: isreal
  integer,pointer                   :: cstrt, cstp, tstrt, tstp
  integer                           :: LDAB, LDR, LDC, LDQ
  external                          :: DGEMM

  call d_getLeadingDim(LDAB,A)
  call d_getLeadingDim(LDR,Ru)
  call d_getLeadingDim(LDC,Cu)
  call d_getLeadingDim(LDQ,Q)

  ! Get current and total ranges
  call tApGet(apC,cstrt,cstp)
  call tApGet(apT,tstrt,tstp)

  defl = .false.

  ! Check first if it actually a cc pole
  call d_eigenvalues2x2(A,B,cstrt+2,cstrt+1,c,c2,isreal)
  if (isreal) then
    call d_make_hess(c,cstrt+2,cstrt+1,A,B,compAB,apC,apT,Q,compQ,qcol+1,&
                       Z,compZ,zcol)
    call d_start_deflation_single(defl,A,B,compAB,apC,apT,Q,compQ,qcol)
  else
    ! Compute the QR factorization
    ! We create a copy of the leading block for the deflation test
    Q4(1:3,1:2) = A(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    Z4(1:3,1:2) = B(cstrt+1:cstrt+3,cstrt+1:cstrt+2)
    ! Create zero in B(cstrt+2,cstrt+1)
    call d_compute_ct(Z4(1,1),Z4(2,1),c,s,r)
    ! Apply transformation
    Z4(1,1) = r
    Z4(2,1) = dzero
    t = -s * Z4(1,2) + c * Z4(2,2)
    Z4(1,2) = c * Z4(1,2) + s * Z4(2,2)
    Z4(2,2) = t
    call d_apply_ct_l(Q4(1,1:2),Q4(2,1:2),c,s)
    ! Try to zero row cstrt+3 all at once
    if ((abs(Z4(3,2)) .gt. abs(Q4(3,1))) .and. &
       (abs(Z4(3,2)) .gt. abs(Q4(3,2)))) then
      call d_compute_ct(Z4(2,2),Z4(3,2),c2,s2,r)
      Z4(2,2) = r
      Z4(3,2) = dzero
    else
      if (abs(Q4(3,1)) .gt. abs(Q4(3,2))) then
        call d_compute_ct(Q4(2,1),Q4(3,1),c2,s2,r)
      else
        call d_compute_ct(Q4(2,2),Q4(3,2),c2,s2,r)
      end if
      t = -s2 * Z4(2,2) + c2 * Z4(3,2)
      Z4(2,2) = c2 * Z4(2,2) + s2 * Z4(3,2)
      Z4(3,2) = t
    end if
    call d_apply_ct_l(Q4(2,1:2),Q4(3,1:2),c2,s2)
    ! Test if deflated
    if (((abs(Q4(3,1)) + abs(Q4(3,2))) .lt. &
          mp * (abs(Q4(1,1)) + abs(Q4(2,2)) + abs(Q4(1,2)) &
          + abs(Q4(2,1)))) .and. (abs(Z4(3,2)) .lt. &
           mp * (abs(Z4(1,1)) + abs(Z4(2,2)) + abs(Z4(1,2))))) then
        defl = .true.
        ! Construct Q
        Q4 = dzero
        Q4(1,1) = done
        Q4(2,2) = done
        Q4(3,3) = done
        call d_apply_ct_l(Q4(1,1:3),Q4(2,1:3),c,s)
        call d_apply_ct_l(Q4(2,1:3),Q4(3,1:3),c2,s2)
    end if

    if (defl) then
      if (compAB) then
        call DGEMM('N','N',3,tstp-cstrt,3,done,Q4(1,1),4,A(cstrt+1,cstrt+1),&
                       LDAB,dzero,Ru(1,1),LDR)
        A(cstrt+1:cstrt+3,cstrt+1:tstp) = Ru(1:3,1:tstp-cstrt)
        call DGEMM('N','N',3,tstp-cstrt,3,done,Q4(1,1),4,B(cstrt+1,cstrt+1),&
                       LDAB,dzero,Ru(1,1),LDR)
        B(cstrt+1:cstrt+3,cstrt+1:tstp) = Ru(1:3,1:tstp-cstrt)
      else
        call DGEMM('N','N',3,cstp-cstrt,3,done,Q4(1,1),4,A(cstrt+1,cstrt+1),&
                       LDAB,dzero,Ru(1,1),LDR)
        A(cstrt+1:cstrt+3,cstrt+1:cstp) = Ru(1:3,1:cstp-cstrt)
        call DGEMM('N','N',3,cstp-cstrt,3,done,Q4(1,1),4,B(cstrt+1,cstrt+1),&
                       LDAB,dzero,Ru(1,1),LDR)
        B(cstrt+1:cstrt+3,cstrt+1:cstp) = Ru(1:3,1:cstp-cstrt)
      end if
      ! Put elements to zero
      A(cstrt+3,cstrt+1:cstrt+2) = dzero
      B(cstrt+3,cstrt+1:cstrt+2) = dzero
      B(cstrt+2,cstrt+1) = dzero
      if (compQ) then
        call DGEMM('N','T',LDQ,3,3,cone,Q(1,qcol),LDQ,Q4(1,1),4,&
                      czero,Cu(1,1),LDC)
        Q(1:LDQ,qcol:qcol+2) = Cu(1:LDQ,1:3)
      end if
      ! This also updates it in the caller function
      call tApPut(apC,strt=cstrt+2,stp=cstp)
    end if
  end if
  end subroutine d_start_deflation_double

  subroutine d_stop_deflation_single(defl, A, B, compAB, apC, apT,&
                                    Z, compZ, zcol )
  !___________________________________________________________________________
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Check for deflation at the end of the active area
  ! The pencil (A,B) and the Z matrix is updated as requested if a deflation
  ! is found.
  ! The current active part apC is also updated, which also updates
  ! it in the caller function.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator (true if deflated, else false)
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the current active part
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
  ! last edit: September 12, 2018
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compZ
    integer, intent(in)               ::  zcol
    logical, intent(out)              ::  defl
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     :: c, s, r,t, tt, u(2)
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    integer                           :: LDZ

    call d_getLeadingDim(LDZ,Z)

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    defl = .false.
    ! first check for direct deflation
    if ((abs(A(cstp,cstp-1)) .lt. &
         mp*(abs(A(cstp,cstp)) + abs(A(cstp-1,cstp-1)))) .and. &
        (abs(B(cstp,cstp-1)) .lt. &
         mp*(abs(B(cstp,cstp)) + abs(B(cstp-1,cstp-1))))) then
      defl = .true.
      A(cstp,cstp-1) = dzero
      B(cstp,cstp-1) = dzero
      ! This also updates it in the caller function
      call tApPut(apC,strt=cstrt,stp=cstp-1)
    elseif (abs(A(cstp,cstp) * B(cstp,cstp-1) - &
                A(cstp,cstp-1) * B(cstp,cstp)) .lt. &
           mp*(abs(A(cstp-1,cstp-1)) + abs(B(cstp-1,cstp-1)))) then
      if (abs(A(cstp,cstp-1)) .ge. abs(B(cstp,cstp-1))) then
        ! Compute rotation based on A
        call d_compute_ct(A(cstp,cstp),A(cstp,cstp-1),c,s,r)
        ! Apply to B
        r = B(cstp,cstp)
        t = B(cstp,cstp-1)
        u = B(cstp-1,cstp-1:cstp)
      else
        ! Compute rotation based on B
        call d_compute_ct(B(cstp,cstp),B(cstp,cstp-1),c,s,r)
        ! Apply to A
        r = A(cstp,cstp)
        t = A(cstp,cstp-1)
        u = A(cstp-1,cstp-1:cstp)
      end if

      tt = -s * r + c * t
      r = c * r + s * t
      t = tt
      u(2) =  -s * u(1) + c * u(2)

      if (abs(t) .lt. mp*(abs(r) + abs(u(2)))) then
        defl = .true.
        if (compAB) then
          call d_apply_ct_r(A(tstrt+1:cstp,cstp-1),A(tstrt+1:cstp,cstp),c,s)
          call d_apply_ct_r(B(tstrt+1:cstp,cstp-1),B(tstrt+1:cstp,cstp),c,s)
        else
          call d_apply_ct_r(A(cstrt+1:cstp,cstp-1),A(cstrt+1:cstp,cstp),c,s)
          call d_apply_ct_r(B(cstrt+1:cstp,cstp-1),B(cstrt+1:cstp,cstp),c,s)
        endif
        if (compZ) then
          call d_apply_ct_r(Z(1:LDZ,zcol),Z(1:LDZ,zcol+1),c,s)
        end if
        A(cstp,cstp-1) = dzero
        B(cstp,cstp-1) = dzero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt,stp=cstp-1)
      end if
    end if

    ! Double check if a 2x2 block has separated
    if (cstp-cstrt .gt. 3) then
      if (.not. (defl .or. (abs(A(cstp-1,cstp-3)) .gt. dzero))) then
      call d_interior_deflation_sng(cstp-2,defl,A,B)
      if (defl) then
        ! update active parts
        call tApInsertAfter(apC,strt=cstp-2,stp=cstp)
        call tApPut(apC,strt=cstrt, stp=cstp-2)
      end if
     end if
    end if
  end subroutine d_stop_deflation_single

  subroutine d_stop_deflation_double(defl, A, B, compAB, apC, apT,&
                                    Q, compQ, qcol, Z, compZ, zcol )
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Check for deflation at the end of the active area. It is assumed
  ! that there is a double pole at the end. It is also checked if the double
  ! pole actually is a complex conjugate pair. If not, the Hessenberg form
  ! is restored and a single deflation is checked.
  ! The pencil (A,B) and the Q and Z matrices are updated as requested
  ! if a deflation is found.
  ! The current active part apC is also updated, which also updates
  ! it in the caller function.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator (true if deflated, else false)
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B) as defined by
  !           apT. If not, then only to the active part apC
  ! apC     type(tAp) [INOUT]
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
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1, zcol+2
  !___________________________________________________________________________
  ! last edit: January 22, 2019
  real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:),Z(:,:)
  logical, intent(in)               ::  compAB, compQ, compZ
  integer, intent(in)               ::  qcol, zcol
  logical, intent(out)              ::  defl
  type(tAp), pointer, intent(inout) ::  apC
  type(tAp), pointer, intent(in)    ::  apT

  real(kind=dp)                     :: c, c2, s, s2, r,t
  logical                           :: isreal
  integer,pointer                   :: cstrt, cstp, tstrt, tstp
  integer                           :: LDAB, LDC, LDZ

  external                          :: DGEMM

  call d_getLeadingDim(LDAB,A)
  call d_getLeadingDim(LDC,Cu)
  call d_getLeadingDim(LDZ,Z)

  call tApGet(apC,cstrt,cstp)
  call tApGet(apT,tstrt,tstp)

  defl = .false.

  ! Check first if it is actually a cc pole
  call d_eigenvalues2x2(A,B,cstp-1,cstp-2,c,c2,isreal)
  if (isreal) then
    call d_make_hess(c,cstp-1,cstp-2,A,B,compAB,apC,apT,Q,compQ,qcol,&
                       Z,compZ,zcol)
    call d_stop_deflation_single(defl,A,B,compAB,apC,apT,Z,compZ,zcol)
  else
    ! Compute the RQ factorization
    ! We create a copy of the trailing block for the deflation test
    Q4(1:2,1:3) = A(cstp-1:cstp,cstp-2:cstp)
    Z4(1:2,1:3) = B(cstp-1:cstp,cstp-2:cstp)
    ! Create zero in B(cstp,cstp-1)
    call d_compute_ct(Z4(2,3),Z4(2,2),c,s,r)
    ! Apply transformation
    Z4(2,3) = r
    Z4(2,2) = dzero
    t = s * Z4(1,2) + c * Z4(1,3)
    Z4(1,2) = c * Z4(1,2) - s * Z4(1,3)
    Z4(1,3) = t
    call d_apply_ct_r(Q4(1:2,2),Q4(1:2,3),c,s)
    ! Try to zero column cstp-2 all at once
    if ((abs(Z4(1,1)) .gt. abs(Q4(1,1))) .and. &
       (abs(Z4(1,1)) .gt. abs(Q4(2,1)))) then
      call d_compute_ct(Z4(1,2),Z4(1,1),c2,s2,r)
      Z4(1,1) = dzero
      Z4(1,2) = r
    else
      if (abs(Q4(1,1)) .gt. abs(Q4(2,1))) then
        call d_compute_ct(Q4(1,2),Q4(1,1),c2,s2,r)
      else
        call d_compute_ct(Q4(2,2),Q4(2,1),c2,s2,r)
      end if
      t = s2 * Z4(1,1) + c2 * Z4(1,2)
      Z4(1,1) = c2 * Z4(1,1) - s2 * Z4(1,2)
      Z4(1,2) = t
    end if
    call d_apply_ct_r(Q4(1:2,1),Q4(1:2,2),c2,s2)
    ! Test if deflated
    if (((abs(Q4(1,1)) + abs(Q4(2,1))) .lt. &
          mp * (abs(Q4(1,2)) + abs(Q4(2,3)) + abs(Q4(1,3)) &
          + abs(Q4(2,2)))) .and. (abs(Z4(1,1)) .lt. &
           mp * (abs(Z4(1,2)) + abs(Z4(2,3))+ abs(Z4(1,3))))) then
        defl = .true.
        ! Construct Z
        Z4 = dzero
        Z4(1,1) = done
        Z4(2,2) = done
        Z4(3,3) = done
        call d_apply_ct_r(Z4(1:3,2),Z4(1:3,3),c,s)
        call d_apply_ct_r(Z4(1:3,1),Z4(1:3,2),c2,s2)
    end if

    if (defl) then
      if (compAB) then
        call DGEMM('N','N',cstp-tstrt,3,3,done,A(tstrt+1,cstp-2),LDAB,&
                  Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(tstrt+1:cstp,cstp-2:cstp) = Cu(1:cstp-tstrt,1:3)
        call DGEMM('N','N',cstp-tstrt,3,3,done,B(tstrt+1,cstp-2),LDAB,&
                  Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(tstrt+1:cstp,cstp-2:cstp) = Cu(1:cstp-tstrt,1:3)
      else
        call DGEMM('N','N',cstp-cstrt,3,3,done,A(cstrt+1,cstp-2),LDAB,&
                 Z4(1,1),4,dzero,Cu(1,1),LDC)
        A(cstrt+1:cstp,cstp-2:cstp) = Cu(1:cstp-cstrt,1:3)
        call DGEMM('N','N',cstp-cstrt,3,3,done,B(cstrt+1,cstp-2),LDAB,&
                 Z4(1,1),4,dzero,Cu(1,1),LDC)
        B(cstrt+1:cstp,cstp-2:cstp) = Cu(1:cstp-cstrt,1:3)
      end if
      ! Put elements to zero
      A(cstp-1:cstp,cstp-2) = dzero
      B(cstp-1:cstp,cstp-2) = dzero
      B(cstp,cstp-1) = dzero
      if (compZ) then
        call DGEMM('N','N',LDZ,3,3,done,Z(1,zcol),LDZ,Z4(1,1),4,&
                   dzero,Cu(1,1),LDC)
        Z(1:LDZ,zcol:zcol+2) = Cu(1:LDZ,1:3)
      end if
      ! This also updates it in the caller function
      call tApPut(apC,strt=cstrt,stp=cstp-2)
    end if
  end if
  end subroutine d_stop_deflation_double
end module
