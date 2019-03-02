! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!  Tightly-packed RQZ iteration for real block-Hessenberg pairs with
!  aggressive early deflation
! ___________________________________________________________________
module d_RQZtp
  use u_parameters
  use u_activeparts
  use d_memorymgmt
  use d_setpoles
  use d_swappoles12
  use d_swappoles22
  use d_swappolesmk
  use d_computepoles
  use d_deflations
  use d_aed
  use d_RQZ

  implicit none
  private
  public d_RQZm

contains
  subroutine d_RQZm(A,B, compAB, strtidx, stpidx, Q, compQ, qcol, &
                     Z, compZ, zcol)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Implements the multishift, multipole RQZ method with aggressive
  ! early deflation
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       double array [INOUT]
  !           Matrix in upper block-Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper blocktriangular form and is equivalent with A
  !           if compAB is true.
  !
  ! B       double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with B
  !           if compAB is true.
  !
  ! compAB  boolean       [IN]
  !           If true, the equivalences are applied to (A,B) such
  !           that the output (block)triangular pair is equivalent to
  !           the input pair. If not, the output is not
  !           equivalent but the eigenvalues of the original pair
  !           are available on the diagonal.
  !
  ! strtidx integer     [IN]
  !           Start index of part of the matrix pencil that needs
  !           to be reduced to Schur form. Provide '0' if start is
  !           the beginning of the pencil.
  !
  ! stpidx  integer     [IN]
  !           Stop index of the part of the matrix pencil that needs
  !           to be reduced to Schur form. Provide N (matrix size) if
  !           stop is the end of the pencil. The strt and stp parameters
  !           need to satisfy:
  !               0 <= strtidx < stpidx <= N
  !
  ! Q       double array [INOUT]
  !           Orthonormal matrix containing the left Schur
  !           vectors. The number of columns of Q should be at least
  !           stpidx-strtidx. The Q array is not initialized in this routine.
  !
  ! compQ   boolean       [IN]
  !           If true, the left Schur vectors are computed.
  !
  ! qcol    integer       [IN]
  !           Column index from whereon the equivalence transformations
  !           are applied. If they should be applied from the beginning
  !           of the Q array, qcol needs to be set to 1. If qcol > 1, then
  !           the Q array should have at least qcol+stpidx-strtidx-1 columns.
  !
  ! Z       double array [INOUT]
  !           Orthonormal matrix containing the right Schur
  !           vectors. The number of columns of Z should be at least
  !           stp-strt. The Z array is not initialized in this routine.
  !
  ! compZ   boolean       [IN]
  !           If true, the right Schur vectors are computed.
  !
  ! zcol    integer       [IN]
  !           Column index from whereon the equivalence transformations
  !           are applied. If they should be applied from the beginning
  !           of the Z array, zcol needs to be set to 1. If zcol > 1, then
  !           the Z array should have at least zcol+stpidx-strtidx-1 columns.
  !
  !
  ! On successful completion of the routine the following equivalences hold:
  !   A(out) = Q^T A(in) Z
  !   B(out) = Q^T B(in) Z
  !___________________________________________________________________________
  ! last edit: January 15, 2019

    real(kind=dp), intent(inout)      ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compQ, compZ
    integer, intent(in)               ::  strtidx, stpidx, qcol, zcol

    ! Internal variables
    type(tAp), pointer                :: apHead, apCurrent, apTotal
    integer                           :: i, j, N, nit, nswaps, qci, zci, mi, &
                                         mc, ki, kc, kco, wei, wec, wsc, wsi, LDAB, &
                                         LDQ,LDQM, LDZ, LDZM, LDC, LDR, &
                                         nreal, ncc, bc, maxbk
    integer,pointer                   :: cstrt, cstp
    logical                           :: problemsolved, defl
    real(kind=dp), allocatable        :: reva(:),revb(:), ccevr(:), ccevc(:)
    external DGEMM

    real(kind=dp) :: DLAMCH

    ! Initialize variables
    sfmin = DLAMCH('S')

    ! apHead stores the active parts throughout the iteration
    call tApInit(apHead,strt=strtidx,stp=stpidx)
    ! apTotal stores the active part at the start and is never modified
    call tApInit(apTotal,strt=strtidx,stp=stpidx)
    ! Q and Z indices
    qci = qcol-strtidx
    zci = zcol-strtidx
    ! Tightly-packed settings and memory
    call d_getSettings(stpidx-strtidx,mi,ki,wei,wsi)
    mc = mi
    kc = ki
    wec = wei
    wsc = wsi
    bc = 6
    maxbk = max(bc,ki)
    ! Leading dimensions
    call d_getLeadingDim(LDAB,A)
    call d_getLeadingDim(LDQ,Q)
    call d_getLeadingDim(LDZ,Z)
    N = max(LDAB,LDQ,LDZ)
    ! Allocate memory for swapping and deflation
    if (wei .gt. 1) then
      call d_allocateAedMem(wei+2)
    end if
    call d_allocateSwap4Mem
    call d_allocateSwapMMem(mi+maxbk+2) ! m can become 1 larger and k as well
    call d_allocateRCMem(mi+maxbk+2,N) ! This should always be the largest
    allocate(reva(mi+1),revb(mi+1),ccevr((mi+1)/2),ccevc((mi+1)/2)) ! For shift intro
    revb = done
    ! Leading dimensions
    call d_getLeadingDim(LDZM,ZM)
    call d_getLeadingDim(LDC,Cu)
    call d_getLeadingDim(LDR,Ru)
    call d_getLeadingDim(LDQM,QM)
    !Performance tracking
    nit = 0
    nswaps = 0
    problemsolved = .false.

    ! Start of main program
    do while (.not. problemsolved)
    ! Main loop is repeated until no active parts are left
    apCurrent => apHead
    do
    ! Loop over all active parts
      call tApGet(apCurrent,cstrt,cstp)
      call d_getSettings(cstp-cstrt,mc,kc,wec,wsc)
      bc = 6
      kco = kc
      if (mc .eq. 1) then
        ! Use the combination shift implementation
        QM = dzero
        do i=1, cstp-cstrt
          QM(i,i) = done
        end do
        ZM = dzero
        do i=1, cstp-cstrt
          ZM(i,i) = done
        end do

        call d_RQZ2(A, B, compAB, cstrt, cstp, QM, .true., 1, &
                     ZM, .true., 1, .false.)

        ! Update the remainder of the pencil
        if (compAB) then
          if (cstp .lt. stpidx) then
              ! Row update A
              call DGEMM('T','N',cstp-cstrt,stpidx-cstp,cstp-cstrt,cone,&
                     QM(1,1),LDQM,A(cstrt+1,cstp+1),LDAB,czero,Ru(1,1),LDR)
              A(cstrt+1:cstp,cstp+1:stpidx) = Ru(1:cstp-cstrt,1:stpidx-cstp)
              ! Row update B
              call DGEMM('T','N',cstp-cstrt,stpidx-cstp,cstp-cstrt,cone,&
                     QM(1,1),LDQM,B(cstrt+1,cstp+1),LDAB,czero,Ru(1,1),LDR)
              B(cstrt+1:cstp,cstp+1:stpidx) = Ru(1:cstp-cstrt,1:stpidx-cstp)
          end if
          if (cstrt .gt. 0) then
            ! Column update A
            call DGEMM('N','N',cstrt-strtidx,cstp-cstrt,cstp-cstrt,cone,&
                      A(strtidx+1,cstrt+1),LDAB,ZM(1,1),LDZM,czero,Cu(1,1),LDC)
            A(strtidx+1:cstrt,cstrt+1:cstp) = Cu(1:cstrt-strtidx,1:cstp-cstrt)
            ! Column update B
            call DGEMM('N','N',cstrt-strtidx,cstp-cstrt,cstp-cstrt,cone,&
                      B(strtidx+1,cstrt+1),LDAB,ZM(1,1),LDZM,czero,Cu(1,1),LDC)
            B(strtidx+1:cstrt,cstrt+1:cstp) = Cu(1:cstrt-strtidx,1:cstp-cstrt)
          end if
        end if

        if (compQ) then
          call DGEMM('N','N',LDQ,cstp-cstrt,cstp-cstrt,cone,&
             Q(1,qci+cstrt),LDQ,QM(1,1),LDQM,czero,Cu,LDC)
          Q(1:LDQ,qci+cstrt:qci+cstp-1) = Cu(1:LDQ,1:cstp-cstrt)
        end if

        if (compZ) then
          call DGEMM('N','N',LDZ,cstp-cstrt,cstp-cstrt,cone,&
             Z(1,zci+cstrt),LDZ,ZM(1,1),LDZM,czero,Cu,LDC)
          Z(1:LDZ,zci+cstrt:zci+cstp-1) = Cu(1:LDZ,1:cstp-cstrt)
        end if

        ! Not an active part anymore, remove from list
        if (.not. (associated(tApPrevious(apCurrent)))) then
          ! apCurrent is the head, update head
          apHead => tApNext(apCurrent) ! can be null(), no issue
        end if
        ! Delete and restart
        call tApDelete(apCurrent)
        apCurrent => apHead
      else
        ! AED at front of pencil
        call d_start_aed(defl, wsc, A, B, compAB, apCurrent, apTotal, &
              Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)
        if (defl) then
          ! Make sure that we're handling the correct range
          call tApGet(apCurrent,cstrt,cstp)
        end if
        ! Regular deflation at front of pencil
        if (abs(A(cstrt+3,cstrt+1)) .gt. dzero) then
          call d_start_deflation_double(defl,A,B,compAB,apCurrent,apTotal,&
                                     Q,compQ,qci+cstrt,Z,compZ,zci+cstrt)
        else
          call d_start_deflation_single(defl,A,B,compAB,apCurrent,apTotal,&
                                     Q,compQ,qci+cstrt)
        end if
        if (defl) then
          ! Make sure that we're handling the correct range
          call tApGet(apCurrent,cstrt,cstp)
        end if
        ! Deflations along subdiagonal
        call d_check_interior_deflations(defl,A, B, compAB, apCurrent, apTotal,&
                Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)
        ! apCurrent is updated; consequently cstrt, and cstp
        if (defl) then
          ! Make sure that we're handling the correct range
          call tApGet(apCurrent,cstrt,cstp)
        end if
        if (cstp - cstrt .gt. mc) then
          ! Compute the shifts as the eigenvalues of the trailing mcxmc block
          ! Create copy (we do not want to change the pencil here)
          QM(1:mc,1:mc) = A(cstp-mc+1:cstp,cstp-mc+1:cstp)
          ZM(1:mc,1:mc) = B(cstp-mc+1:cstp,cstp-mc+1:cstp)
          call d_RQZ2(QM, ZM, .false., 0, mc, Q, .false., 1, &
                       Z, .false., 1, .false.)
          call d_getEigenvalues(QM, ZM, 0, mc, reva,nreal,ccevr,ccevc,ncc)
          ! Introduce the mc shifts in pencil
          call d_first_m_poles(mc,reva,revb,nreal,ccevr,ccevc,ncc,A,B,compAB,&
                  apCurrent, apTotal, Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)

          if (cstp - cstrt .gt. mc) then
            nit = nit + mc
            ! Chase
            ! Take as many steps of size kc as possible
            j = 0
            i = cstrt + 1
            do while (i .le. cstp-mc-kc-1)
              call d_swapmk(i, mc, kc,A, B, compAB, apCurrent, apTotal, &
                Q, compQ, qci+i, Z, compZ, zci+i-1)
              i = i + kc
              nswaps = nswaps + (mc*kc)
              j = j + 1
              ! reset kc
              kc = kco
            end do

            ! Take remaining steps
            kc = cstp-i-mc
            call d_swapmk(i, mc, kc ,A, B, compAB, apCurrent,&
              apTotal, Q, compQ, qci+i,&
              Z, compZ, zci+i-1)

            ! AED at end of pencil
            call d_stop_aed(defl, wec, A, B, compAB, apCurrent, apTotal, &
                  Q, compQ, qci+cstp-wec, Z, compZ, zci+cstp-wec)
            ! Regular deflations at end
            if (.not.(abs(A(cstp,cstp-2)) .gt. dzero)) then
              call d_stop_deflation_single(defl, A, B, compAB,apCurrent, apTotal,&
                                          Z, compZ, zci+cstp-2)
            else
              call d_stop_deflation_double(defl, A, B, compAB,apCurrent, apTotal,&
                                          Q, compQ, qci+cstp-2,Z, compZ, zci+cstp-3)
            end if
            if (cstp - cstrt .gt. mc) then
              ! Compute mc poles as the eigenvalues of the leading mcxmc block
              ! Create copy (we do not want to change the pencil here)
              QM(1:mc,1:mc) = A(cstrt+1:cstrt+mc,cstrt+1:cstrt+mc)
              ZM(1:mc,1:mc) = B(cstrt+1:cstrt+mc,cstrt+1:cstrt+mc)
              call d_RQZ2(QM, ZM, .false., 0, mc, Q, .false., 1, &
                         Z, .false., 1, .false.)
              call d_getEigenvalues(QM, ZM, 0, mc, reva,nreal,ccevr,ccevc,ncc)
              ! Introduce the mc poles (mc should match here)
              call d_last_m_poles(mc,bc,reva,revb,nreal,ccevr,ccevc,ncc,A,B,compAB,&
                apCurrent,apTotal,Q,compQ,cstp,Z,compZ,cstp)
            end if
          end if
        end if
      end if

      ! Check for exit condition
      if (associated(tApNext(apCurrent))) then
        apCurrent => tApNext(apCurrent) ! Set the next one
      else
        exit ! Exit this loop over the active parts
      end if
    ! end of loop over active parts
    end do

    ! Check if entire problem is solved
    call tApGet(apHead,cstrt,cstp)
    if (.not. (associated(cstrt) .or. associated(cstp))) then
      problemsolved = .true.
    end if

    ! end of main loop
    end do

    ! Normalize all 2x2 blocks
    apCurrent => apTotal ! This allows us to run make_hess if needed
    call d_normalizeSchur(A,B,compAB,apCurrent,apTotal,&
             Q,compQ,qci,Z,compZ,zci, .true.)

    ! free everything
    apCurrent => null()
    call tApFree(apHead)
    call tApFree(apTotal)

    call d_deallocateSwapMMem
    call d_deallocateSwap4Mem
    call d_deallocateAedMem
    call d_deallocateRCMem

  end subroutine

  subroutine d_getSettings(N,m,k,we,ws)
    integer, intent(in)             ::  N
    integer, intent(out)            ::  m, k, we, ws

    if (N .lt. ssth) then
      k = 1
      m = 1
      we = 1
      ws = 1
    elseif (N .lt. 150) then
      k = 4
      m = 4
      we = 6
      ws = 4
    elseif (N .lt. 250) then
      k = 8
      m = 8
      we = 10
      ws = 4
    elseif (N .lt. 501) then
      k = 16
      m = 16
      we = 18
      ws = 6
    elseif (N .lt. 1001) then
      k = 32
      m = 32
      we = 34
      ws = 10
    elseif (N .lt. 3000) then
      k = 64
      m = 64
      we = 66
      ws = 16
    elseif (N .lt. 6000) then
      k = 128
      m = 128
      we = 130
      ws = 32
    else
      k = 256
      m = 256
      we = 266
      ws = 48
    end if
  end subroutine
end module
