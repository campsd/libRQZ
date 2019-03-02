! libRQZ
! author: daan.camps@cs.kuleuven.be
! Description:
!   Implements the routines for selecting the shifts and poles during
!   the iteration.
! ___________________________________________________________________
module d_computepoles
  use u_parameters
  use u_activeparts
  use d_ctransformations
  use d_memorymgmt

  implicit none
  private
  public d_get_shift, d_get_pole, d_eigenvalues2x2, d_make_hess, &
         d_ccstdform, d_normalizeSchur

contains
  subroutine d_get_shift(A,B,idx,alpha,beta,isreal)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes the eigenvalues of the trailing two by two block.
  ! If these are real,
  !  alpha/beta is set equal to the one closest to the last
  !  diagonal element and isreal is set to true
  ! If this is a complex conjugate pair,
  !  alpha contains the real part and beta the imaginary part,
  !  isreal is set to false
  !
  ! The user can overwrite this function to use a different shift
  ! selection procedure.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [IN]
  !           Upper block-Hessenberg matrix A
  ! B       double array [IN]
  !           Upper block-Hessenberg matrix B
  ! idx     integer [IN]
  !           Index indicating the two-by-two block,
  !             idx-1:idx,idx-1:idx
  ! alpha   double [OUT]
  !           Defines the shift. If isreal is false: real
  !           part of the shift, else numerator
  ! beta    double [OUT]
  !           Defines the shift. If isreal is false: imaginary
  !           part of the shift, else denominator
  ! isreal  boolean [OUT]
  !           Defines if shift should be treated as a single
  !           real-valued shift or a pair of complex conjugate
  !           shifts.
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp),intent(in)     :: A(:,:), B(:,:)
    integer, intent(in)          :: idx
    real(kind=dp),intent(out)    :: alpha, beta
    logical, intent(out)         :: isreal

    call d_eigenvalues2x2(A,B,idx-1,idx-1,alpha,beta,isreal)

    ! Select best one if real
    if (isreal) then
      if (abs(B(idx,idx)*alpha-A(idx,idx)) .lt. &
          abs(B(idx,idx)*beta-A(idx,idx))) then
       beta = done
      else
       alpha = beta
       beta = done
      endif
    end if
  end subroutine

  subroutine d_get_pole(A,B,idx,alpha,beta,isreal)
  !
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes the eigenvalues of the leading two by two block.
  ! If these are real,
  !  alpha/beta is set equal to the one closest to the first
  !  diagonal element and isreal is set to true
  ! If this is a complex conjugate pair,
  !  alpha contains the real part and beta the imaginary part,
  !  isreal is set to false
  !
  ! The user can overwrite this function to use a different pole
  ! selection procedure.
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [IN]
  !           Upper block-Hessenberg matrix A
  ! B       double array [IN]
  !           Upper block-Hessenberg matrix B
  ! idx     integer [IN]
  !           Index indicating the two-by-two block,
  !             idx:idx+1,idx:idx+1
  ! alpha   double [OUT]
  !           Defines the shift. If isreal is false: real
  !           part of the shift, else numerator
  ! beta    double [OUT]
  !           Defines the shift. If isreal is false: imaginary
  !           part of the shift, else denominator
  ! isreal  boolean [OUT]
  !           Defines if shift should be treated as a single
  !           real-valued shift or a pair of complex conjugate
  !           shifts.
  !___________________________________________________________________________
  ! last edit: September 5, 2018
    real(kind=dp),intent(in)     :: A(:,:), B(:,:)
    integer, intent(in)          :: idx
    real(kind=dp),intent(out)    :: alpha, beta
    logical, intent(out)         :: isreal

    call d_eigenvalues2x2(A,B,idx,idx,alpha,beta,isreal)

    ! Select best one if real
    if (isreal) then
      if (abs(B(idx,idx)*alpha-A(idx,idx)) .lt. &
          abs(B(idx,idx)*beta-A(idx,idx))) then
       beta = done
      else
       alpha = beta
       beta = done
      endif
    end if
  end subroutine

  subroutine d_eigenvalues2x2(A,B,i,j,l1,l2,isreal)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Computes the eigenvalues of a 2-by-2 block
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [IN]
  !           array block-Hessenberg
  ! B       double array [IN]
  !           array block-Hessenberg
  ! i       integer [IN]
  !           eigenvalues of range i:i+1, j:j+1 are computed
  ! j       integer [IN]
  !           eigenvalues of range i:i+1, j:j+1 are computed
  ! l1      double [OUT]
  !          If isreal is true, real eigenvalue of (A,B)
  !          If isreal is false, real part of eigenvalue
  ! l2      double [OUT]
  !          If isreal is true, real eigenvalue of (A,B)
  !          If isreal is false, imaginary part of eigenvalue
  ! isreal  boolean [OUT]
  !           Indicates if it are two real eigenvalues or
  !           a pair of complex conjugate eigenvalues
  !
  !___________________________________________________________________________
  ! last edit: September 7, 2018
  ! WRAPPER FOR LAPACK DLAG2
    real(kind=dp), intent(in)   :: A(:,:), B(:,:)
    integer, intent(in)         :: i,j
    real(kind=dp), intent(out)  :: l1, l2
    logical, intent(out)        :: isreal

    integer                     :: LDAB
    real(kind=dp)               :: S1, S2, WI, AC(2,2), BC(2,2), c, s, r
    external                    :: DLAG2

    call d_getLeadingDim(LDAB,A)

    if (abs(B(i+1,j)) .gt. dzero ) then
      AC = A(i:i+1,j:j+1)
      BC = B(i:i+1,j:j+1)
      call d_compute_ct(BC(1,1),BC(2,1),c,s,r)
      call d_apply_ct_l(AC(1,1:2),AC(2,1:2),c,s)
      call d_apply_ct_l(BC(1,1:2),BC(2,1:2),c,s)
      BC(2,1) = dzero
      call DLAG2(AC, 2, BC, 2,sfmin, S1,S2, l1, l2, WI )
    else
      call DLAG2(A(i,j), LDAB, B(i,j), LDAB,sfmin, S1,S2, l1, l2, WI )
    end if

    if (WI .gt. dzero) then
      isreal = .false.
      l1 = l1 /S1
      l2 = WI /S1
    else
      isreal = .true.
      l1 = l1 /S1
      l2 = l2 /S2
    endif

  end subroutine

  subroutine d_make_hess(l1,i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Restores the Hessenberg shape of a 2x2 block that has real eigenvalues,
  ! one of which is equal to l1
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  !
  ! l1         double [IN]
  !             One of two eigenvalues of the block i:i+1,j:j+1
  ! i       integer [IN]
  !           Row index of 2x2 block: i:i+1
  ! j       integer [IN]
  !           Column index of 2x2 block: j:j+1
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [INOUT]
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
  !              qcol, qcol+1
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: January 11, 2019

    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    real(kind=dp), intent(in)       :: l1 ! eigenvalue
    integer, intent(in)             :: i,j, qcol, zcol
    logical, intent(in)             :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)  :: apC, apT

    real(kind=dp)                   :: c,s,r
    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDQ, LDZ

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)

    Z4(1:2,1:2) = A(i:i+1,j:j+1) - l1 * B(i:i+1,j:j+1)

    if (abs(Z4(1,2)) .gt. abs(Z4(2,2))) then
      call d_compute_ct(Z4(1,2),Z4(1,1),c,s,r)
    else
      call d_compute_ct(Z4(2,2),Z4(2,1),c,s,r)
    end if

    if (compAB) then
      call d_apply_ct_r(A(tstrt+1:i+1,j),A(tstrt+1:i+1,j+1),c,s)
      call d_apply_ct_r(B(tstrt+1:i+1,j),B(tstrt+1:i+1,j+1),c,s)
    else
      call d_apply_ct_r(A(cstrt+1:i+1,j),A(cstrt+1:i+1,j+1),c,s)
      call d_apply_ct_r(B(cstrt+1:i+1,j),B(cstrt+1:i+1,j+1),c,s)
    end if

    if (compZ) then
      call d_apply_ct_r(Z(1:LDZ,zcol),Z(1:LDZ,zcol+1),c,s)
    end if

    if (abs(B(i+1,j)) .ge. abs(A(i+1,j))) then
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
      call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),c,-s)
    end if

  end subroutine

  subroutine d_ccstdform(i,j,A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Restores the standard form of a 2x2 block
  ! containing complex conjugate eigenvalues:
  !     B will be made diagonal with positive elements
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  !
  ! i       integer [IN]
  !           Row index of 2x2 block: i:i+1
  ! j       integer [IN]
  !           Column index of 2x2 block: j:j+1
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper block-Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [INOUT]
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
  !              qcol, qcol+1
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1
  !___________________________________________________________________________
  ! last edit: January 15, 2019

    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    integer, intent(in)             :: i,j, qcol, zcol
    logical, intent(in)             :: compAB, compQ, compZ
    type(tAp), pointer, intent(in)  :: apC, apT

    real(kind=dp)                   :: c,s,r,cl,sl,cr,sr,g,h,s1,s2
    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDQ, LDZ

    external  ::  DLASV2

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)

    ! First rotate B if necessary
    if (abs(B(i+1,j)) .gt. dzero) then
      call d_compute_ct(B(i,j),B(i+1,j),c,s,r)
      g = c * B(i,j+1) + s * B(i+1,j+1)
      h = -s * B(i,j+1) + c *B(i+1,j+1)
      ! Then make it diagonal
      call DLASV2(r,g,h,s2,s1,sr,cr,sl,cl)
      sr = -sr
      ! Combine left rotations
      r = cl*c - sl*s
      sl = cl*s + sl*c
      cl = r
    else
      call DLASV2(B(i,j),B(i,j+1),B(i+1,j+1),s2,s1,sr,cr,sl,cl)
      sr = -sr
    end if

    ! Set B block
    B(i,j) = s1
    B(i+1,j) = dzero
    B(i,j+1) = dzero
    B(i+1,j+1) = s2

    ! Update everything
    if (compAB) then
      call d_apply_ct_l(A(i,j:tstp),A(i+1,j:tstp),cl,sl)
      call d_apply_ct_l(B(i,j+2:tstp),B(i+1,j+2:tstp),cl,sl)
      call d_apply_ct_r(A(tstrt+1:i+1,j),A(tstrt+1:i+1,j+1),cr,sr)
      call d_apply_ct_r(B(tstrt+1:i-1,j),B(tstrt+1:i-1,j+1),cr,sr)
    else
      call d_apply_ct_l(A(i,j:cstp),A(i+1,j:cstp),cl,sl)
      call d_apply_ct_l(B(i,j+2:cstp),B(i+1,j+2:cstp),cl,sl)
      call d_apply_ct_r(A(cstrt+1:i+1,j),A(cstrt+1:i+1,j+1),cr,sr)
      call d_apply_ct_r(B(cstrt+1:i-1,j),B(cstrt+1:i-1,j+1),cr,sr)
    end if

    if (compQ) then
      call d_apply_ct_r(Q(1:LDQ,qcol),Q(1:LDQ,qcol+1),cl,-sl)
    end if

    if (compZ) then
      call d_apply_ct_r(Z(1:LDZ,zcol),Z(1:LDZ,zcol+1),cr,sr)
    end if

  end subroutine



  subroutine d_normalizeSchur(A,B,compAB,apC,apT,Q,compQ,qcol,Z,compZ,zcol,sa)
  !
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Normalize all 2x2 blocks in the generalized real Schur form
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  !
  ! A       double array [INOUT]
  !           Upper block-Hessenberg matrix A
  ! B       double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
  ! apC     type(tAp) [INOUT]
  !           Current active part
  ! apT     type(tAp) [IN]
  !           Total part of the pencil the method is acting on
  !           (used when compAB is true)
  ! Q       double array [INOUT]
  !           Orthonormal equivalence matrix Q
  ! compQ   boolean [IN]
  !           If true, the left Schur vectors are updated
  ! qcol    integer [IN]
  !           Start index of the columns of Q that are updated
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated
  ! sa      boolean [IN]
  !           Standalone usage: if true, all cc blocks are normalized,
  !           if false, all 2x2 blocks are tested for cc and corrected.
  !___________________________________________________________________________
  ! last edit: January 16, 2019
    real(kind=dp), intent(inout)    :: A(:,:), B(:,:), Q(:,:), Z(:,:)
    integer, intent(in)             :: qcol, zcol
    logical, intent(in)             :: compAB, compQ, compZ, sa
    type(tAp), pointer, intent(in)  :: apC, apT

    integer,pointer                 :: cstrt, cstp, tstrt, tstp
    integer                         :: LDQ, LDZ
    real(kind=dp)                   :: mu, nu
    integer                         :: i
    logical                         :: isreal

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)

    i = tstrt+1
    do while (i .lt. tstp)
      if (abs(A(i+1,i)) .gt. dzero) then
        call d_eigenvalues2x2(A,B,i,i,mu,nu,isreal)
        if (isreal) then
          call d_make_hess(mu,i,i,A,B,compAB,apC,apT,&
                Q,compQ,qcol+i-1,Z,compZ,zcol+i-1)
        elseif (sa .and. (abs(B(i+1,i)) .gt. dzero) .or. (abs(B(i,i+1)) .gt. dzero)) then
          call d_ccstdform(i,i,A,B,compAB,apC,apT,&
                Q,compQ,qcol+i-1,Z,compZ,zcol+i-1)
        end if
        i = i + 2
      else
        i = i + 1
      end if
    end do
  end subroutine
end module
