! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   aggressive early deflation
! ___________________________________________________________________
module z_aed

  use u_activeparts
  use u_parameters
  use z_memorymgmt
  use z_swappoles
  use z_RQZ

  implicit none
  private
  public z_start_aed, z_stop_aed

contains

  subroutine z_start_aed(defl, w, A, B, compAB, apC, apT, Q, compQ, qcol, &
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
  ! A       complex double array [INOUT]
  !           Upper Hessenberg matrix A
  ! B       complex double array [INOUT]
  !           Upper Hessenberg matrix B
  ! compAB  boolean [INOUT]
  !           If true, the equivalences are applied to (A,B). If not,
  !           then only to the active parts
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
  !              qcol, qcol+w
  ! Z       complex double array [INOUT]
  !           Unitary equivalence matrix Z
  ! compZ   boolean [IN]
  !           If true, the right Schur vectors are updated
  ! zcol    integer [IN]
  !           Start index of the columns of Z that are updated:
  !              zcol, zcol+w-1
  !___________________________________________________________________________
  ! last edit: September 4, 2018
  complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
  logical, intent(in)               ::  compAB, compQ, compZ
  logical, intent(out)              ::  defl
  integer, intent(in)               ::  w, qcol, zcol
  type(tAp), pointer, intent(inout) ::  apC
  type(tAp), pointer, intent(in)    ::  apT

  integer                           :: i, j, nd, nnd, LDAB, LDQ, LDQW, LDZ, &
                                       LDZW, LDC, LDR, nndth, nnt
  integer,pointer                   :: cstrt, cstp, tstrt, tstp
  type(tAp), pointer                :: apCI, apTI

  external ZGEMM

  call tApGet(apC,cstrt,cstp)
  call tApGet(apT,tstrt,tstp)

  call tApInit(apCI,cstrt,cstrt+w)
  call tApInit(apTI,cstrt,cstrt+w)

  defl = .false.
  call z_getLeadingDim(LDAB,A)
  call z_getLeadingDim(LDQ,Q)
  call z_getLeadingDim(LDQW,QW)
  call z_getLeadingDim(LDZ,Z)
  call z_getLeadingDim(LDZW,ZW)
  call z_getLeadingDim(LDC,CuW)
  call z_getLeadingDim(LDR,RuW)
  nndth = (w/3) + 1 ! threshold

  QW = czero
  do i=1, w+1
    QW(i,i) = cone
  end do

  ZW = czero
  do i=1, w
    ZW(i,i) = cone
  end do

  ! Transform the part in the deflation window to Schur form with
  ! single shift RQZ
  call z_RQZ1(A, B, .true., cstrt, cstrt+w, QW, .true., 1, ZW, .true., 1, &
              RuW(1:w,1), RuW(1:w,2))

  ! Now we'll test eigenvalues in the window for deflations until too
  ! many non deflatables are detected
  nd = 0 ! number of deflations detected
  nnd = 0 ! number of non-deflatables
  nnt = 0 ! number of not tested

  do i=1,w-1
    ! Compute spikes and store them in RuW array
    RuW(1,1:w) = ZW(w,1:w) * A(cstrt+w+1,cstrt+w) ! spike A
    RuW(2,1:w) = ZW(w,1:w) * B(cstrt+w+1,cstrt+w) ! spike B

    if ( (abs(RuW(1,nd+1)) .lt. mp * abs(A(cstrt+nd+1,cstrt+nd+1))) .and. &
         (abs(RuW(2,nd+1)) .lt. mp * abs(B(cstrt+nd+1,cstrt+nd+1))) ) then
      ! This one can be deflated
      defl = .true.
      nd = nd + 1
    else
      ! This one cannot be deflated
      nnd = nnd + 1
      if (nnd .gt. nndth) then
        nnt = w - nd - nnd
        exit
      end if
      ! we move another one upfront
      do j =  cstrt+nd+nnd,cstrt+nd+1,-1
        call z_swap(j, j, A, B, .true., apCI, apTI, QW, .true., &
                    j-cstrt, ZW, .true., j-cstrt)
      end do
    end if
  end do

  ! Compute spikes and store them in RuW array
  RuW(1,1:w) = ZW(w,1:w) * A(cstrt+w+1,cstrt+w) ! spike A
  RuW(2,1:w) = ZW(w,1:w) * B(cstrt+w+1,cstrt+w) ! spike B

  ! Also test the last one if all are ok
  if (nnt .eq. 0) then
    if ( (abs(RuW(1,nd+1)) .lt. mp * abs(A(cstrt+nd+1,cstrt+nd+1))) .and. &
         (abs(RuW(2,nd+1)) .lt. mp * abs(B(cstrt+nd+1,cstrt+nd+1))) ) then
      ! This one can be deflated
      nd = nd + 1
    end if
  end if

  ! Set first nd elements in the spike to zero and update pencil
  RuW(1,1:nd) = czero
  A(cstrt+w+1,cstrt+1:cstrt+w) = RuW(1,1:w)
  RuW(2,1:nd) = czero
  B(cstrt+w+1,cstrt+1:cstrt+w) = RuW(2,1:w)

  ! Restore Hessenberg shape if required via permutation
  ! Warning: need to change the poles afterwards, this is done in RQZtp !
  if (nd .lt. w) then
    ! QW
    RuW(1:w+1,1:w-nd-1) = QW(1:w+1,nd+2:w)
    QW(1:w+1,nd+2) = QW(1:w+1,w+1)
    QW(1:w+1,nd+3:w+1) = RuW(1:w+1,1:w-nd-1)
    ! A
    RuW(1:w-nd-1,1:w-nd) = A(cstrt+nd+2:cstrt+w,cstrt+nd+1:cstrt+w)
    A(cstrt+nd+2,cstrt+nd+1:cstrt+w) = A(cstrt+w+1,cstrt+nd+1:cstrt+w)
    A(cstrt+nd+3:cstrt+w+1,cstrt+nd+1:cstrt+w) = RuW(1:w-nd-1,1:w-nd)
    ! B
    RuW(1:w-nd-1,1:w-nd) = B(cstrt+nd+2:cstrt+w,cstrt+nd+1:cstrt+w)
    B(cstrt+nd+2,cstrt+nd+1:cstrt+w) = B(cstrt+w+1,cstrt+nd+1:cstrt+w)
    B(cstrt+nd+3:cstrt+w+1,cstrt+nd+1:cstrt+w) = RuW(1:w-nd-1,1:w-nd)
  end if

  ! Update the remainder of the pencil
  if (compAB) then
    if (cstrt .gt. 0) then
      ! A column update
      call ZGEMM('N','N',cstrt-tstrt,w,w,cone,A(tstrt+1,cstrt+1),LDAB,&
                ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      A(tstrt+1:cstrt,cstrt+1:cstrt+w) = CuW(1:cstrt-tstrt,1:w)
      ! B column update
      call ZGEMM('N','N',cstrt-tstrt,w,w,cone,B(tstrt+1,cstrt+1),LDAB,&
                ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      B(tstrt+1:cstrt,cstrt+1:cstrt+w) = CuW(1:cstrt-tstrt,1:w)
    end if
    ! A row update
    call ZGEMM('C','N',w+1,tstp-cstrt-w,w+1,cone,QW(1,1),LDQW,&
              A(cstrt+1,cstrt+w+1),LDAB,czero,RuW(1,1),LDR)
    A(cstrt+1:cstrt+w+1,cstrt+w+1:tstp) = RuW(1:w+1,1:tstp-cstrt-w)
    ! B row update
    call ZGEMM('C','N',w+1,tstp-cstrt-w,w+1,cone,QW(1,1),LDQW,&
              B(cstrt+1,cstrt+w+1),LDAB,czero,RuW(1,1),LDR)
    B(cstrt+1:cstrt+w+1,cstrt+w+1:tstp) = RuW(1:w+1,1:tstp-cstrt-w)
  else
    ! A row update
    call ZGEMM('C','N',w+1,cstp-cstrt-w,w+1,cone,QW(1,1),LDQW,&
              A(cstrt+1,cstrt+w+1),LDAB,czero,RuW(1,1),LDR)
    A(cstrt+1:cstrt+w+1,cstrt+w+1:cstp) = RuW(1:w,1:cstp-cstrt-w)
    ! B row update
    call ZGEMM('C','N',w+1,cstp-cstrt-w,w+1,cone,QW(1,1),LDQW,&
              B(cstrt+1,cstrt+w+1),LDAB,czero,RuW(1,1),LDR)
    B(cstrt+1:cstrt+w+1,cstrt+w+1:cstp) = RuW(1:w,1:cstp-cstrt-w)
  end if

  if (compQ) then
    call ZGEMM('N','N',LDQ,w+1,w+1,cone,Q(1,qcol),LDQ,QW(1,1),LDQW,czero,CuW,LDC)
    Q(1:LDQ,qcol:qcol+w) = CuW(1:LDQ,1:w+1)
  end if

  if (compZ) then
    call ZGEMM('N','N',LDZ,w,w,cone,Z(1,zcol),LDZ,ZW(1,1),LDZW,czero,CuW,LDC)
    Z(1:LDZ,zcol:zcol+w-1) = CuW(1:LDZ,1:w)
  end if

  ! Update indices (also updates it in main)
  call tApPut(apC,strt=cstrt+nd,stp=cstp)

  ! Free everything
  call tApFree(apCI)
  call tApFree(apTI)
  end subroutine z_start_aed

  subroutine z_stop_aed(defl, w, A, B, compAB, apC, apT, Q, compQ, qcol,&
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
    ! A       complex double array [INOUT]
    !           Upper Hessenberg matrix A
    ! B       complex double array [INOUT]
    !           Upper Hessenberg matrix B
    ! compAB  boolean [INOUT]
    !           If true, the equivalences are applied to (A,B). If not,
    !           then only to the active parts
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
    !              qcol, qcol+1, ...,qcol+w-1
    ! Z       complex double array [INOUT]
    !           Unitary equivalence matrix Z
    ! compZ   boolean [IN]
    !           If true, the right Schur vectors are updated
    ! zcol    integer [IN]
    !           Start index of the columns of Z that are updated:
    !              zcol, zcol+1, ..., zcol+w
    !___________________________________________________________________________
    ! last edit: September 4, 2018
    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compQ, compZ
    logical, intent(out)              ::  defl
    integer, intent(in)               ::  w, qcol, zcol
    type(tAp), pointer, intent(inout) ::  apC
    type(tAp), pointer, intent(in)    ::  apT

    integer                           :: i, j, nd, nnd, LDAB, LDQ, LDQW, &
                                         LDZ, LDZW, LDC, LDR, nndth, nnt
    integer,pointer                   :: cstrt, cstp, tstrt, tstp
    type(tAp), pointer                :: apCI, apTI

    external ZGEMM

    call tApGet(apC,cstrt,cstp)
    call tApGet(apT,tstrt,tstp)

    call tApInit(apCI,cstp-w,cstp)
    call tApInit(apTI,cstp-w,cstp)

    defl = .false.
    call z_getLeadingDim(LDAB,A)
    call z_getLeadingDim(LDQ,Q)
    call z_getLeadingDim(LDQW,QW)
    call z_getLeadingDim(LDZ,Z)
    call z_getLeadingDim(LDZW,ZW)
    call z_getLeadingDim(LDC,CuW)
    call z_getLeadingDim(LDR,RuW)
    nndth = (w/2) + 1 ! threshold

    QW = czero
    do i=1, w
      QW(i,i) = cone
    end do

    ZW = czero
    do i=1, w+1
      ZW(i,i) = cone
    end do

    ! Transform the part in the deflation window to Schur form with single shift RQZ
    call z_RQZ1(A, B, .true., cstp-w, cstp, QW, .true., 1, ZW, .true., 2, &
                RuW(1:w,1), RuW(1:w,2))

    ! Now we'll test eigenvalues in the window for deflations until too
    ! many non deflatables are detected
    nd = 0 ! number of deflations detected
    nnd = 0 ! number of non-deflatables
    nnt = 0 ! number of not tested

    do i=w-1,1,-1
      ! Compute spikes and store them in RuW array
      RuW(1:w,1) = dconjg(QW(1,1:w)) * A(cstp-w+1,cstp-w) ! spike A
      RuW(1:w,2) = dconjg(QW(1,1:w)) * B(cstp-w+1,cstp-w) ! spike B

      if ( (abs(RuW(w-nd,1)) .lt. mp * abs(A(cstp-nd,cstp-nd))) .and. &
           (abs(RuW(w-nd,2)) .lt. mp * abs(B(cstp-nd,cstp-nd))) ) then
        ! This one can be deflated
        defl = .true.
        nd = nd + 1
      else
        ! This one cannot be deflated
        nnd = nnd + 1
        if (nnd .gt. nndth) then
          nnt = w - nd - nnd
          exit
        end if
        ! we move another one upfront
        do j =  cstp-nd-nnd, cstp-nd-1
          call z_swap(j, j, A, B, .true., apCI, apTI, QW, .true., j+w-cstp, &
                      ZW, .true., j+w-cstp+1)
        end do
      end if
    end do

    ! Compute spikes and store them in RuW array
    RuW(1:w,1) = dconjg(QW(1,1:w)) * A(cstp-w+1,cstp-w) ! spike A
    RuW(1:w,2) = dconjg(QW(1,1:w)) * B(cstp-w+1,cstp-w) ! spike B

    ! Also test the last one if all are ok
    if (nnt .eq. 0) then
      if ( (abs(RuW(w-nd,1)) .lt. mp * abs(A(cstp-nd,cstp-nd))) .and. &
           (abs(RuW(w-nd,2)) .lt. mp * abs(B(cstp-nd,cstp-nd))) ) then
        ! This one can be deflated
        nd = nd + 1
      end if
    end if

    ! Set last nd elements in the spike to zero and update pencil
    RuW(w-nd+1:w,1) = czero
    A(cstp-w+1:cstp,cstp-w) = RuW(1:w,1)
    RuW(w-nd+1:w,2) = czero
    B(cstp-w+1:cstp,cstp-w) = RuW(1:w,2)

    ! Restore Hessenberg shape if required via permutation
    ! Warning: need to change the poles afterwards (done in RQZtp) !
    if (w-nd .gt. 1) then
      ! ZW
      RuW(1:w+1,1:w-nd-1) = ZW(1:w+1,2:w-nd)
      ZW(1:w+1,w-nd) = ZW(1:w+1,1)
      ZW(1:w+1,1:w-nd-1) = RuW(1:w+1,1:w-nd-1)
      ! A
      RuW(1:w-nd,1:w-nd-1) = A(cstp-w+1:cstp-nd,cstp-w+1:cstp-nd-1)
      A(cstp-w+1:cstp-nd,cstp-nd-1) = A(cstp-w+1:cstp-nd,cstp-w)
      A(cstp-w+1:cstp-nd,cstp-w:cstp-nd-2) = RuW(1:w-nd,1:w-nd-1)
      ! B
      RuW(1:w-nd,1:w-nd-1) = B(cstp-w+1:cstp-nd,cstp-w+1:cstp-nd-1)
      B(cstp-w+1:cstp-nd,cstp-nd-1) = B(cstp-w+1:cstp-nd,cstp-w)
      B(cstp-w+1:cstp-nd,cstp-w:cstp-nd-2) = RuW(1:w-nd,1:w-nd-1)
    end if

    ! Update the remainder of the pencil
    if (compAB) then
      ! A column update
      call ZGEMM('N','N',cstp-w-tstrt,w+1,w+1,cone,A(tstrt+1,cstp-w),LDAB,&
                 ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      A(tstrt+1:cstp-w,cstp-w:cstp) = CuW(1:cstp-w-tstrt,1:w+1)
      ! B column update
      call ZGEMM('N','N',cstp-w-tstrt,w+1,w+1,cone,B(tstrt+1,cstp-w),LDAB,&
                 ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      B(tstrt+1:cstp-w,cstp-w:cstp) = CuW(1:cstp-w-tstrt,1:w+1)
      if (cstp .lt. tstp) then
        ! A row update
        call ZGEMM('C','N',w,tstp-cstp,w,cone,QW(1,1),LDQW,&
                   A(cstp-w+1,cstp+1),LDAB,czero,RuW(1,1),LDR)
        A(cstp-w+1:cstp,cstp+1:tstp) = RuW(1:w,1:tstp-cstp)
        ! B row update
        call ZGEMM('C','N',w,tstp-cstp,w,cone,QW(1,1),LDQW,&
                   B(cstp-w+1,cstp+1),LDAB,czero,RuW(1,1),LDR)
        B(cstp-w+1:cstp,cstp+1:tstp) = RuW(1:w,1:tstp-cstp)
      end if
    else
      ! A column update
      call ZGEMM('N','N',cstp-w-cstrt,w+1,w+1,cone,A(cstrt+1,cstp-w),&
                LDAB,ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      A(cstrt+1:cstp-w,cstp-w:cstp) = CuW(1:cstp-w-cstrt,1:w+1)
      ! B column update
      call ZGEMM('N','N',cstp-w-cstrt,w+1,w+1,cone,A(cstrt+1,cstp-w),&
                LDAB,ZW(1,1),LDZW,czero,CuW(1,1),LDC)
      A(cstrt+1:cstp-w,cstp-w:cstp) = CuW(1:cstp-w-cstrt,1:w+1)
    end if

    if (compQ) then
      call ZGEMM('N','N',LDQ,w,w,cone,Q(1,qcol),LDQ,QW(1,1),LDQW,czero,CuW,LDC)
      Q(1:LDQ,qcol:qcol+w-1) = CuW(1:LDQ,1:w)
    end if

    if (compZ) then
      call ZGEMM('N','N',LDZ,w+1,w+1,cone,Z(1,zcol),LDZ,&
                ZW(1,1),LDZW,czero,CuW,LDC)
      Z(1:LDZ,zcol:zcol+w) = CuW(1:LDZ,1:w+1)
    end if

    ! Update indices (also updates it in main)
    call tApPut(apC,strt=cstrt,stp=cstp-nd)

    ! Free everything
    call tApFree(apCI)
    call tApFree(apTI)

  end subroutine z_stop_aed


end module z_aed
