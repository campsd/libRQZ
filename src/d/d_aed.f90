! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   aggressive early deflation
! ___________________________________________________________________
module d_aed

  use u_activeparts
  use u_parameters
  use d_memorymgmt
  use d_ctransformations
  use d_swappoles12
  use d_swappoles22
  use d_RQZ

  implicit none
  private
  public d_stop_aed, d_start_aed

contains

  subroutine d_start_aed(defl, w, A, B, compAB, apC, apT, Q, compQ, qcol, &
                         Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Perform aggressive early deflation at the front of the active part with
  ! window length w
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator
  ! w       integer [IN]
  !           Deflation window size
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
  !              qcol, qcol+w
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+w-1
  !___________________________________________________________________________
  ! last edit: January 31, 2019
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compQ, compZ
    logical, intent(out)              ::  defl
    integer, intent(in)               ::  w, qcol, zcol
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     :: amax, ssrc, c, s, r
    integer                           :: i, j, nd, nnd, LDAB, LDQ, LDQW, LDZ, &
                                         LDZW, LDC, LDR, nndth, nnt, wi
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external                          :: DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    wi = w
    if (abs(A(cstrt+w+1,cstrt+w-1)) .gt. dzero) then
      wi = wi + 1 ! we take the entire block
    end if

    call tApInit(apCI,cstrt,cstrt+wi)
    call tApInit(apTI,cstrt,cstrt+wi)

    defl = .false.
    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDQW,QW)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDZW,ZW)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)
    nndth = (wi/3) + 1 ! threshold

    QW = dzero
    do i=1, wi
      QW(i,i) = done
    end do

    ZW = dzero
    do i=1, wi
      ZW(i,i) = done
    end do

    ! Deflation logic
    ! Check if block after deflation window is cc or not
    if (abs(A(cstrt+wi+2,cstrt+wi)) .gt. dzero) then
      amax = max(abs(A(cstrt+wi+1,cstrt+wi)),abs(A(cstrt+wi+2,cstrt+wi)))
    else
      amax = abs(A(cstrt+wi+1,cstrt+wi))
    end if

    ! Reduce to Schur form
    call d_RQZ2(A, B, .true., cstrt, cstrt+wi, QW, .true., 1, &
                 ZW, .true., 1, .false.)

    ! Now we'll test eigenvalues in the window for deflations until too
    ! many non deflatables are detected
    nd = 0 ! number of deflations detected
    nnd = 0 ! number of non-deflatables
    nnt = 0 ! number of not tested

    do while ((nd + nnd) .lt. wi-1)
      ! Compute spikes
      Ru(1,1:wi) = ZW(wi,1:wi) * amax ! spike A
      Ru(2,1:wi) = ZW(wi,1:wi) * B(cstrt+wi+1,cstrt+wi) ! spike B

      ! Check if the eigenvalue at the end of the spike is cc or real
      if (abs(A(cstrt+nd+2,cstrt+nd+1)) .gt. dzero ) then
        ! It is cc, so two consecutive elements need to be small
        if ( (abs(Ru(1,nd+1)) + abs(Ru(1,nd+2))) .lt. mp * &
           (abs(A(cstrt+nd+1,cstrt+nd+1)) + abs(A(cstrt+nd+2,cstrt+nd+2))) .and. &
             (abs(Ru(2,nd+1)) + abs(Ru(2,nd+2))) .lt. mp * &
          (abs(B(cstrt+nd+1,cstrt+nd+1)) + abs(B(cstrt+nd+2,cstrt+nd+2)))) then
          ! This one can be deflated
          defl = .true.
          nd = nd + 2
        else
          ! This one cannot be deflated
          nnd = nnd + 2
        end if
      else
        ! It is real so only one element needs to be small
        if ( (abs(Ru(1,nd+1)) .lt. mp * abs(A(cstrt+nd+1,cstrt+nd+1))) .and. &
             (abs(Ru(2,nd+1)) .lt. mp * abs(B(cstrt+nd+1,cstrt+nd+1)))) then
          ! This one can be deflated
          defl = .true.
          nd = nd + 1
       else
         ! This one cannot be deflated
         nnd = nnd + 1
       end if
      end if

      ! Exit if too many tests failed
      if (nnd .gt. nndth) then
        nnt = wi - nd - nnd
        exit
      end if

      if (nd + nnd .lt. wi-1) then
        ! We now move the eigenvalue at position
        ! cstrt + nd + nnd + 1 to position cstrt + nd + 1
        j = cstrt + nd + nnd + 1 ! Index in large matrix
        i = nd + nnd + 1  ! Index in small matrix
        if (abs(A(j+1,j)) .gt. dzero) then
          ! This is a 2x2 block
          do while (i .gt. nd + 1)
            if (i-2 .gt. nd) then
              if (abs(A(j-1,j-2)) .gt. dzero) then
                call d_swap22(j-2, j-2, A, B, .true., apCI, apTI, &
                      QW, .true., i-2, ZW, .true., i-2)
                j = j - 2
                i = i - 2
              else
                call d_swap12(j-1, j-1, A, B, .true., apCI, apTI, &
                      QW, .true., i-1, ZW, .true., i-1)
                j = j - 1
                i = i - 1
              end if
           else
             call d_swap12(j-1, j-1, A, B, .true., apCI, apTI, &
                   QW, .true., i-1, ZW, .true., i-1)
             j = j - 1
             i = i - 1
           end if
          end do
        else
          ! This is a 1x1 block
          do while (i .gt. nd + 1)
            if (i-2 .gt. nd) then
              if (abs(A(j-1,j-2)) .gt. dzero) then
                call d_swap21(j-2, j-2, A, B, .true., apCI, apTI, &
                      QW, .true., i-2, ZW, .true., i-2)
                j = j - 2
                i = i - 2
              else
                call d_swap11(j-1, j-1, A, B, .true., apCI, apTI, &
                      QW, .true., i-1, ZW, .true., i-1)
                j = j - 1
                i = i - 1
              end if
            else
              call d_swap11(j-1, j-1, A, B, .true., apCI, apTI, &
                    QW, .true., i-1, ZW, .true., i-1)
              j = j - 1
              i = i - 1
            end if
          end do
        end if
      end if
    end do

    if (nd + nnd .eq. wi-1) then
      ! Test the last one, which is a single one
      if ( (abs(Ru(1,nd+1)) .lt. mp * abs(A(cstrt+nd+1,cstrt+nd+1))) .and. &
           (abs(Ru(2,nd+1)) .lt. mp * abs(B(cstrt+nd+1,cstrt+nd+1)))) then
        ! This one can be deflated
        defl = .true.
        nd = nd + 1
      else
       ! This one cannot be deflated
       nnd = nnd + 1
      end if
    end if

    ! Restore the block-Hessenberg shape
    if (amax .gt. abs(B(cstrt+wi+1,cstrt+wi))) then
      if (abs(A(cstrt+wi+2,cstrt+wi)) .gt. abs(A(cstrt+wi+1,cstrt+wi))) then
        ssrc = A(cstrt+wi+2,cstrt+wi)
      else
        ssrc = A(cstrt+wi+1,cstrt+wi)
      end if
    else
      ssrc = B(cstrt+wi+1,cstrt+wi)
    endif

    Ru(3,1:wi) = ZW(wi,1:wi) * ssrc
    do j = nd+1, wi-2
      ! Compute the rotation
      call d_compute_ct(Ru(3,j+1),Ru(3,j),c,s,r)
      ! Apply to (A,B)
      call d_apply_ct_r(A(cstrt+1:cstrt+j+2,cstrt+j), &
                        A(cstrt+1:cstrt+j+2,cstrt+j+1),c,s)
      call d_apply_ct_r(B(cstrt+1:cstrt+j+2,cstrt+j), &
                        B(cstrt+1:cstrt+j+2,cstrt+j+1),c,s)
      ! Apply to ZW
      call d_apply_ct_r(ZW(:,j),ZW(:,j+1),c,s)
      Ru(3,j) = dzero
      Ru(3,j+1) = r
    end do

    if (nd+1 .lt. wi) then
      call d_compute_ct(Ru(3,wi),Ru(3,wi-1),c,s,r)
      ! Apply to (A,B)
      call d_apply_ct_r(A(cstrt+1:cstrt+wi,cstrt+wi-1), &
                        A(cstrt+1:cstrt+wi,cstrt+wi),c,s)
      call d_apply_ct_r(B(cstrt+1:cstrt+wi,cstrt+wi-1), &
                        B(cstrt+1:cstrt+wi,cstrt+wi),c,s)
      ! Apply to ZW
      call d_apply_ct_r(ZW(:,wi-1),ZW(:,wi),c,s)

      if (r/ssrc .lt. dzero) then
        A(cstrt+wi+1:cstrt+wi+2,cstrt+wi) = - A(cstrt+wi+1:cstrt+wi+2,cstrt+wi)
        B(cstrt+wi+1,cstrt+wi) = - B(cstrt+wi+1,cstrt+wi)
      end if
    end if

    ! Update the remainder of the pencil
    if (compAB) then
      if (cstrt .gt. 0) then
        ! A column update
        call DGEMM('N','N',cstrt-tstrt,wi,wi,done,A(tstrt+1,cstrt+1),LDAB,&
                  ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
        A(tstrt+1:cstrt,cstrt+1:cstrt+wi) = Cu(1:cstrt-tstrt,1:wi)
        ! B column update
        call DGEMM('N','N',cstrt-tstrt,wi,wi,done,B(tstrt+1,cstrt+1),LDAB,&
                  ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
        B(tstrt+1:cstrt,cstrt+1:cstrt+wi) = Cu(1:cstrt-tstrt,1:wi)
      end if
      ! A row update
      call DGEMM('T','N',wi,tstp-cstrt-wi,wi,done,QW(1,1),LDQW,&
                A(cstrt+1,cstrt+wi+1),LDAB,dzero,Ru(1,1),LDR)
      A(cstrt+1:cstrt+wi,cstrt+wi+1:tstp) = Ru(1:wi,1:tstp-cstrt-wi)
      ! B row update
      call DGEMM('T','N',wi,tstp-cstrt-wi,wi,done,QW(1,1),LDQW,&
                B(cstrt+1,cstrt+wi+1),LDAB,dzero,Ru(1,1),LDR)
      B(cstrt+1:cstrt+wi,cstrt+wi+1:tstp) = Ru(1:wi,1:tstp-cstrt-wi)
    else
      ! A row update
      call DGEMM('T','N',wi,cstp-cstrt-wi,wi,done,QW(1,1),LDQW,&
                A(cstrt+1,cstrt+wi+1),LDAB,dzero,Ru(1,1),LDR)
      A(cstrt+1:cstrt+wi,cstrt+wi+1:cstp) = Ru(1:wi,1:cstp-cstrt-wi)
      ! B row update
      call DGEMM('T','N',wi,cstp-cstrt-wi,wi,done,QW(1,1),LDQW,&
                B(cstrt+1,cstrt+wi+1),LDAB,dzero,Ru(1,1),LDR)
      B(cstrt+1:cstrt+wi,cstrt+wi+1:cstp) = Ru(1:wi,1:cstp-cstrt-wi)
    end if

    if (compQ) then
      call DGEMM('N','N',LDQ,wi,wi,done,Q(1,qcol),LDQ,QW(1,1),LDQW,dzero,Cu,LDC)
      Q(1:LDQ,qcol:qcol+wi-1) = Cu(1:LDQ,1:wi)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,wi,wi,done,Z(1,zcol),LDZ,ZW(1,1),LDZW,dzero,Cu,LDC)
      Z(1:LDZ,zcol:zcol+wi-1) = Cu(1:LDZ,1:wi)
    end if

    ! Update indices (also updates it in main)
    call tApPut(apC,strt=cstrt+nd,stp=cstp)

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)

  end subroutine

  subroutine d_stop_aed(defl, w, A, B, compAB, apC, apT, Q, compQ, qcol,&
                        Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Perform aggressive early deflation at the end of the active part with
  ! window length w
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! defl    boolean [OUT]
  !           Deflation indicator
  ! w       integer [IN]
  !           Deflation window size
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
  !              qcol, qcol+1, ...,qcol+w-1
  ! Z       double array [INOUT]
  !           Orthonormal equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+1, ..., zcol+w-1
  !___________________________________________________________________________
  ! last edit: January 31, 2019
    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compQ, compZ
    logical, intent(out)              ::  defl
    integer, intent(in)               ::  w, qcol, zcol
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    real(kind=dp)                     :: amax, c, s, r, ssrc
    integer                           :: i, j, nd, nnd, LDAB, LDQ, LDQW, &
                                         LDZ, LDZW, LDC, LDR, nndth, nnt, wi,&
                                         qci,zci
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external                          :: DGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)
    wi = w
    qci = qcol
    zci = zcol

    if (abs(A(cstp-w+2,cstp-w)) .gt. dzero) then
      wi = wi + 1 ! we take the entire block
      qci = qci - 1
      zci = zci - 1
    end if

    call tApInit(apCI,cstp-wi,cstp)
    call tApInit(apTI,cstp-wi,cstp)

    defl = .false.
    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDQW,QW)
    call d_getLeadingDim(LDZ,Z)
    call d_getLeadingDim(LDZW,ZW)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)
    nndth = (wi/3) + 1 ! threshold

    QW = dzero
    do i=1, wi
      QW(i,i) = done
    end do

    ZW = dzero
    do i=1, wi
      ZW(i,i) = done
    end do

    ! Deflation logic
    ! Check if block before deflation window is cc pole or not
    if (abs(A(cstp-wi+1,cstp-wi-1)) .gt. dzero) then
      amax = max(abs(A(cstp-wi+1,cstp-wi-1)),abs(A(cstp-wi+1,cstp-wi)))
    else
      amax = abs(A(cstp-wi+1,cstp-wi))
    end if

    ! Reduce to Schur form
    call d_RQZ2(A, B, .true., cstp-wi, cstp, QW, .true., 1, &
                 ZW, .true., 1, .false.)

    ! Now we'll test eigenvalues in the window for deflations until too
    ! many non deflatables are detected
    nd = 0 ! number of deflations detected
    nnd = 0 ! number of non-deflatables
    nnt = 0 ! number of not tested

    do while ((nd + nnd) .lt. wi-1)
      ! Compute spikes and store then in Ru array
      Ru(1:wi,1) = QW(1,1:wi) * amax ! spike A
      Ru(1:wi,2) = QW(1,1:wi) * B(cstp-wi+1,cstp-wi) ! spike B

      ! Check if the eigenvalue at the end of the spike is cc or real
      if (abs(A(cstp-nd,cstp-nd-1)) .gt. dzero) then
        ! It is cc, so two consecutive elements need to be small
        if( (abs(Ru(wi-nd,1)) + abs(Ru(wi-nd-1,1))) .lt. mp * &
              (abs(A(cstp-nd,cstp-nd)) + abs(A(cstp-nd-1,cstp-nd-1))) .and. &
              (abs(Ru(wi-nd,2)) + abs(Ru(wi-nd-1,2))) .lt. mp * &
              (abs(B(cstp-nd,cstp-nd)) + abs(B(cstp-nd-1,cstp-nd-1)))) then
          ! This one can be deflated
          defl = .true.
          nd = nd + 2
        else
          ! This one cannot be deflated
          nnd = nnd + 2
        end if
      else
        ! It is real, so only one element needs to be small
        if ( (abs(Ru(wi-nd,1)) .lt. mp * abs(A(cstp-nd,cstp-nd))) .and. &
           (abs(Ru(wi-nd,2)) .lt. mp * abs(B(cstp-nd,cstp-nd))) ) then
          ! This one can be deflated
          defl = .true.
          nd = nd + 1
       else
         ! This one cannot be deflated
         nnd = nnd + 1
       end if
      end if

      ! Exit if too many deflation tests failed
      if (nnd .gt. nndth) then
        nnt = wi - nd - nnd
        exit
      end if

      if (nd + nnd .lt. wi-1) then
        ! We now move the eigenvalue at position
        ! cstp-nd-nnd to position cstp-nd
        j = cstp-nd-nnd ! Index in large matrix
        i = wi-nd-nnd ! Index in small matrix
        if (abs(A(j,j-1)) .gt. dzero) then
          ! This is a 2x2 block
          j = j - 1
          i = i - 1
          do while (i .lt. wi-nd-1)

            if (i+2 .lt. wi-nd) then
              if (abs(A(j+3,j+2)) .gt. dzero) then
                call d_swap22(j, j, A, B, .true., apCI, apTI, &
                      QW, .true., i, ZW, .true., i)
                j = j + 2
                i = i + 2
              else
                call d_swap21(j, j, A, B, .true., apCI, apTI, &
                       QW, .true., i, ZW, .true., i)
                j = j + 1
                i = i + 1
              end if
            else
              call d_swap21(j, j, A, B, .true., apCI, apTI, &
                     QW, .true., i, ZW, .true., i)
              j = j + 1
              i = i + 1
            end if
          end do
        else
          ! This is a 1x1 block
          do while (i .lt. wi-nd)
            if (i+1 .lt. wi-nd) then
              if (abs(A(j+2,j+1)) .gt. dzero) then
                call d_swap12(j, j, A, B, .true., apCI, apTI, &
                      QW, .true., i, ZW, .true., i)
                j = j + 2
                i = i + 2
              else
                call d_swap11(j, j, A, B, .true., apCI, apTI, &
                      QW, .true., i, ZW, .true., i)
                j = j + 1
                i = i + 1
              end if
            else
              call d_swap11(j, j, A, B, .true., apCI, apTI, &
                    QW, .true., i, ZW, .true., i)
              j = j + 1
              i = i + 1
            end if
          end do
        end if
      end if
    end do

    if (nd + nnd .eq. wi-1) then
      ! test the last one, which is a single one
      if ( (abs(Ru(wi-nd,1)) .lt. mp * abs(A(cstp-nd,cstp-nd))) .and. &
         (abs(Ru(wi-nd,2)) .lt. mp * abs(B(cstp-nd,cstp-nd))) ) then
        ! This one can be deflated
        defl = .true.
        nd = nd + 1
     else
       ! This one cannot be deflated
       nnd = nnd + 1
     end if
    end if

    ! Restore the (block)-Hessenberg shape
    if (amax .gt. abs(B(cstp-wi+1,cstp-wi))) then
      if (abs(A(cstp-wi+1,cstp-wi)) .gt. abs(A(cstp-wi+1,cstp-wi-1))) then
        ssrc = A(cstp-wi+1,cstp-wi)
      else
        ssrc = A(cstp-wi+1,cstp-wi-1)
      end if
    else
      ssrc = B(cstp-wi+1,cstp-wi)
    end if

    Ru(1:wi,3) = QW(1,1:wi) * ssrc
    do j = wi-nd-1,2,-1
      ! Compute the rotation
      call d_compute_ct(Ru(j,3),Ru(j+1,3),c,s,r)
      ! Apply to (A,B)
      call d_apply_ct_l(A(cstp+j-wi,cstp+j-wi-1:cstp), &
                        A(cstp+j-wi+1,cstp+j-wi-1:cstp),c,s)
      call d_apply_ct_l(B(cstp+j-wi,cstp+j-wi:cstp), &
                        B(cstp+j-wi+1,cstp+j-wi:cstp),c,s)
      ! Appl to QW
      call d_apply_ct_r(QW(1:wi+1,j),QW(1:wi+1,j+1),c,-s)
      Ru(j,3) = r
      Ru(j+1,3) = dzero
    end do

    if (wi-nd-1 .gt. 0) then
      call d_compute_ct(Ru(1,3),Ru(2,3),c,s,r)

      ! Apply to (A,B)
      call d_apply_ct_l(A(cstp-wi+1,cstp-wi+1:cstp), &
                      A(cstp-wi+2,cstp-wi+1:cstp),c,s)
      call d_apply_ct_l(B(cstp-wi+1,cstp-wi+1:cstp), &
                      B(cstp-wi+2,cstp-wi+1:cstp),c,s)
      ! Apply to QW
      call d_apply_ct_r(QW(1:wi+1,1),QW(1:wi+1,2),c,-s)

      if (r/ssrc .lt. dzero) then
        A(cstp-wi+1,cstp-wi-1:cstp-wi) = -A(cstp-wi+1,cstp-wi-1:cstp-wi)
        B(cstp-wi+1,cstp-wi) = -B(cstp-wi+1,cstp-wi)
      end if
    end if

    ! Update the remainder of the pencil
    if (compAB) then
      ! A column update
      call DGEMM('N','N',cstp-wi-tstrt,wi,wi,done,A(tstrt+1,cstp-wi+1),LDAB,&
                 ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
      A(tstrt+1:cstp-wi,cstp-wi+1:cstp) = Cu(1:cstp-wi-tstrt,1:wi)
      ! B column update
      call DGEMM('N','N',cstp-wi-tstrt,wi,wi,done,B(tstrt+1,cstp-wi+1),LDAB,&
                 ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
      B(tstrt+1:cstp-wi,cstp-wi+1:cstp) = Cu(1:cstp-wi-tstrt,1:wi)
      if (cstp .lt. tstp) then
        ! A row update
        call DGEMM('T','N',wi,tstp-cstp,wi,done,QW(1,1),LDQW,&
                   A(cstp-wi+1,cstp+1),LDAB,dzero,Ru(1,1),LDR)
        A(cstp-wi+1:cstp,cstp+1:tstp) = Ru(1:wi,1:tstp-cstp)
        ! B row update
        call DGEMM('T','N',wi,tstp-cstp,wi,done,QW(1,1),LDQW,&
                   B(cstp-wi+1,cstp+1),LDAB,dzero,Ru(1,1),LDR)
        B(cstp-wi+1:cstp,cstp+1:tstp) = Ru(1:wi,1:tstp-cstp)
      end if
    else
      ! A column update
      call DGEMM('N','N',cstp-wi-cstrt,wi,wi,done,A(cstrt+1,cstp-wi+1),&
                LDAB,ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
      A(cstrt+1:cstp-wi,cstp-wi+1:cstp) = Cu(1:cstp-wi-cstrt,1:wi)
      ! B column update
      call DGEMM('N','N',cstp-wi-cstrt,wi,wi,done,A(cstrt+1,cstp-wi+1),&
                LDAB,ZW(1,1),LDZW,dzero,Cu(1,1),LDC)
      A(cstrt+1:cstp-wi,cstp-wi+1:cstp) = Cu(1:cstp-wi-cstrt,1:wi)
    end if

    if (compQ) then
      call DGEMM('N','N',LDQ,wi,wi,done,Q(1,qci),LDQ,QW(1,1),LDQW,dzero,Cu,LDC)
      Q(1:LDQ,qci:qci+wi-1) = Cu(1:LDQ,1:wi)
    end if

    if (compZ) then
      call DGEMM('N','N',LDZ,wi,wi,done,Z(1,zci),LDZ,&
                ZW(1,1),LDZW,dzero,Cu,LDC)
      Z(1:LDZ,zci:zci+wi-1) = Cu(1:LDZ,1:wi)
    end if

    ! Update indices (also updates it in main)
    call tApPut(apC,strt=cstrt,stp=cstp-nd)

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)

  end subroutine
end module
