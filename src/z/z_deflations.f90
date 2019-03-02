! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   Deflation monitoring
! ___________________________________________________________________
module z_deflations

use u_activeparts
use u_parameters
use z_ctransformations

implicit none
private
public z_check_interior_deflations, z_check_start_deflation, &
       z_check_stop_deflation

contains

  subroutine z_check_interior_deflations(defl,A, B, apC)
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
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! apC     type(tAp) [INOUT]
  !           Active parts linked list.
  !___________________________________________________________________________
  ! last edit: July 30, 2018
    complex(kind=dp), intent(inout)   :: A(:,:), B(:,:)
    type(tAp), pointer, intent(inout) :: apC
    logical, intent(out)              :: defl

    integer                           :: i, cpy
    integer, pointer                  :: cstrt, cstp

    call tApGet(apC,cstrt,cstp)
    defl = .false.
    do i = cstrt+1, cstp-1
      if ((abs(A(i+1,i)) .lt. mp * (abs(A(i,i)) + abs(A(i+1,i+1)))) .and. &
        (abs(B(i+1,i)) .lt. mp * (abs(B(i,i)) + abs(B(i+1,i+1))))) then
        ! Set to zero
        defl = .true.
        A(i+1,i) = czero
        B(i+1,i) = czero
        ! update active parts
        cpy = cstp
        call tApPut(apC,strt=cstrt, stp=i)
        call tApInsertAfter(apC,strt=i,stp=cpy)
        !apC => apC%nextElem
        apC => tApNext(apC)
        call tApGet(apC,cstrt,cstp) ! this does update cstrt,cstp
                                    !(but does not affect the loop counter)
       end if
    end do
  end subroutine z_check_interior_deflations

  subroutine z_check_start_deflation(defl, A, B, compAB, apC, apT, &
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
  ! last edit: July 24, 2018
    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Q(:,:)
    logical, intent(in)               ::  compAB, compQ
    integer, intent(in)               ::  qcol
    logical, intent(out)              ::  defl
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    complex(kind=dp)                  :: s,r,t, tt, u(2)
    real(kind=dp)                     :: c
    integer,pointer                   :: cstrt, cstp, tstrt, tstp

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
        A(cstrt+2,cstrt+1) = czero
        B(cstrt+2,cstrt+1) = czero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt+1,stp=cstp)
    elseif (abs(A(cstrt+1,cstrt+1) * B(cstrt+2,cstrt+1) - &
                A(cstrt+2,cstrt+1) * B(cstrt+1,cstrt+1)) .lt. &
           mp*(abs(A(cstrt+2,cstrt+2)) + abs(B(cstrt+2,cstrt+2)))) then
      if (abs(A(cstrt+2,cstrt+1)) .ge. abs(B(cstrt+2,cstrt+1))) then
        ! Compute rotation based on A
        call z_compute_ct(A(cstrt+1,cstrt+1),A(cstrt+2,cstrt+1),c,s,r)
        ! Apply to B
        r = B(cstrt+1,cstrt+1)
        t = B(cstrt+2,cstrt+1)
        u = B(cstrt+1:cstrt+2,cstrt+2)
      else
        ! Compute rotation based on B
        call z_compute_ct(B(cstrt+1,cstrt+1),B(cstrt+2,cstrt+1),c,s,r)
        ! Apply to A
        r = A(cstrt+1,cstrt+1)
        t = A(cstrt+2,cstrt+1)
        u = A(cstrt+1:cstrt+2,cstrt+2)
      end if

      !tt = -dconjg(s) * r + dconjg(c) * t ! For use with own ct
      tt = -dconjg(s) * r + c * t
      r = c * r + s * t
      t = tt

      !u(2) = -dconjg(s) * u(1) + dconjg(c) * u(2) ! For use with own ct
      u(2) = -dconjg(s) * u(1) + c * u(2)

      if (abs(t) .lt. mp*(abs(r) + abs(u(2)))) then
        defl = .true.
        if (compAB) then
          call z_apply_ct_l(A(cstrt+1,cstrt+1:tstp),A(cstrt+2,cstrt+1:tstp),c,s)
          call z_apply_ct_l(B(cstrt+1,cstrt+1:tstp),B(cstrt+2,cstrt+1:tstp),c,s)
        else
          call z_apply_ct_l(A(cstrt+1,cstrt+1:cstp),A(cstrt+2,cstrt+1:cstp),c,s)
          call z_apply_ct_l(B(cstrt+1,cstrt+1:cstp),B(cstrt+2,cstrt+1:cstp),c,s)
        end if
        if (compQ) then
          !call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),conjg(c),-s) ! With own ct
          call z_apply_ct_r(Q(:,qcol),Q(:,qcol+1),c,-s)
        end if
        A(cstrt+2,cstrt+1) = czero
        B(cstrt+2,cstrt+1) = czero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt+1,stp=cstp)
      end if
    end if
  end subroutine z_check_start_deflation

  subroutine z_check_stop_deflation(defl, A, B, compAB, apC, apT,&
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
  ! last edit: June 7, 2018
    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compZ
    integer, intent(in)               ::  zcol
    logical, intent(out)              ::  defl
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    complex(kind=dp)                  :: s,r,t, tt, u(2)
    real(kind=dp)                     :: c
    integer,pointer                   :: cstrt, cstp, tstrt, tstp

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    defl = .false.
    ! first check for direct deflation
    if ((abs(A(cstp,cstp-1)) .lt. &
         mp*(abs(A(cstp,cstp)) + abs(A(cstp-1,cstp-1)))) .and. &
        (abs(B(cstp,cstp-1)) .lt. &
         mp*(abs(B(cstp,cstp)) + abs(B(cstp-1,cstp-1))))) then
      defl = .true.
      A(cstp,cstp-1) = czero
      B(cstp,cstp-1) = czero
      ! This also updates it in the caller function
      call tApPut(apC,strt=cstrt,stp=cstp-1)
    elseif (abs(A(cstp,cstp) * B(cstp,cstp-1) - &
                A(cstp,cstp-1) * B(cstp,cstp)) .lt. &
           mp*(abs(A(cstp-1,cstp-1)) + abs(B(cstp-1,cstp-1)))) then
      if (abs(A(cstp,cstp-1)) .ge. abs(B(cstp,cstp-1))) then
        ! Compute rotation based on A
        call z_compute_ct(A(cstp,cstp),A(cstp,cstp-1),c,s,r)
        !c = dconjg(c) ! Own ct
        ! Apply to B
        r = B(cstp,cstp)
        t = B(cstp,cstp-1)
        u = B(cstp-1,cstp-1:cstp)
      else
        ! Compute rotation based on B
        call z_compute_ct(B(cstp,cstp),B(cstp,cstp-1),c,s,r)
        !c = dconjg(c) ! Own ct
        ! Apply to A
        r = A(cstp,cstp)
        t = A(cstp,cstp-1)
        u = A(cstp-1,cstp-1:cstp)
      end if

      tt = -dconjg(s) * r + c * t
      !r = dconjg(c) * r + s * t ! Own ct
      r = c * r + s * t
      t = tt
      u(2) =  -dconjg(s) * u(1) + c * u(2)

      if (abs(t) .lt. mp*(abs(r) + abs(u(2)))) then
        defl = .true.
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
        A(cstp,cstp-1) = czero
        B(cstp,cstp-1) = czero
        ! This also updates it in the caller function
        call tApPut(apC,strt=cstrt,stp=cstp-1)
      end if
    end if
  end subroutine z_check_stop_deflation

end module z_deflations
