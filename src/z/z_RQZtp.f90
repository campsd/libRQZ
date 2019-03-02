! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!  Tightly-packed RQZ iteration for complex Hessenberg pairs with
!  aggressive early deflation
! ___________________________________________________________________
module z_RQZtp
  use u_parameters
  use u_activeparts
  use z_memorymgmt
  use z_setpoles
  use z_swappoles
  use z_computepoles
  use z_deflations
  use z_aed
  use z_RQZ

  implicit none
  private
  public z_RQZm

contains
  subroutine z_RQZm(A,B, compAB, strtidx, stpidx, Q, compQ, qcol, &
                     Z, compZ, zcol, alpha, beta)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Implements the tightly-packed shift complex-valued RQZ method
  ! with aggressive early deflation
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! A       complex double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with A
  !           if compAB is true.
  !
  ! B       complex double array [INOUT]
  !           Matrix in upper Hessenberg form on entrance.
  !           On exit, the part defined by strtidx and stpidx is
  !           in upper triangular form and is equivalent with B
  !           if compAB is true.
  !
  ! compAB  boolean       [IN]
  !           If true, the equivalences are applied to (A,B) such
  !           that the output triangular pair is equivalent to
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
  ! Q       complex double array [INOUT]
  !           Unitary matrix containing the left Schur
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
  ! Z       complex double array [INOUT]
  !           Unitary matrix containing the right Schur
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
  ! alpha   complex double array [OUT]
  !           Array of length stpidx-strtidx containing the numerators for
  !           the generalized eigenvalues
  !
  ! beta    complex double array [OUT]
  !           Array of length stpidx-strtidx containing the denominators
  !           for the generalized eigenvalues
  !
  ! On successful completion of the routine the following equivalences hold:
  !   A(out) = Q^* A(in) Z
  !   B(out) = Q^* B(in) Z
  !___________________________________________________________________________
  ! last edit: September 4, 2018

    complex(kind=dp), intent(inout)   ::  A(:,:), B(:,:), Q(:,:), Z(:,:)
    logical, intent(in)               ::  compAB, compQ, compZ
    integer, intent(in)               ::  strtidx, stpidx, qcol, zcol
    complex(kind=dp), intent(out)     ::  alpha(:), beta(:)

    ! Internal variables
    type(tAp), pointer                :: apHead, apCurrent, apTotal
    integer                           :: i, j, N, nit, nswaps, qci, zci, mi, &
                                         mc, ki, kc, wei, wec, wsc, wsi, LDAB, &
                                         LDQ,LDQM, LDZ, LDZM, LDC, LDR
    integer,pointer                   :: cstrt, cstp
    logical                           :: problemsolved, defl

    external ZGEMM

    ! Initialize variables

    ! apHead stores the active parts throughout the iteration
    call tApInit(apHead,strt=strtidx,stp=stpidx)
    ! apTotal stores the active part at the start and is never modified
    call tApInit(apTotal,strt=strtidx,stp=stpidx)
    ! Q and Z indices
    qci = qcol-strtidx
    zci = zcol-strtidx
    ! Tightly-packed settings and memory
    call z_getSettings(stpidx-strtidx,mi,ki,wei,wsi)
    mc = mi
    kc = ki
    wec = wei
    wsc = wsi
    ! Leading dimensions
    call z_getLeadingDim(LDAB,A)
    call z_getLeadingDim(LDQ,Q)
    call z_getLeadingDim(LDZ,Z)
    N = max(LDAB,LDQ,LDZ)
    ! Allocate memory for swapping and deflation
    if (wei .gt. 1) then
      call z_allocateAedMem(wei,N)
    end if
    call z_allocateSwapMem(mi+ki-1,N)
    ! Leading dimensions
    call z_getLeadingDim(LDZM,ZM)
    call z_getLeadingDim(LDC,CuM)
    call z_getLeadingDim(LDR,RuM)
    call z_getLeadingDim(LDQM,QM)
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
      call z_getSettings(cstp-cstrt,mc,kc,wec,wsc)
      if (mc .eq. 1) then
        ! Use the single shift implementation
        QM = czero
        do i=1, cstp-cstrt
          QM(i,i) = cone
        end do
        ZM = czero
        do i=1, cstp-cstrt
          ZM(i,i) = cone
        end do
        call z_RQZ1(A, B, compAB, cstrt, cstp, QM, .true., 1, ZM, .true., 1,&
                    alpha, beta)
        ! Update the remainder of the pencil
        if (compAB) then
          ! Row update A
          call ZGEMM('C','N',cstp-cstrt,stpidx-cstp,cstp-cstrt,cone,&
                     QM(1,1),LDQM,A(cstrt+1,cstp+1),LDAB,czero,RuM(1,1),LDR)
          A(cstrt+1:cstp,cstp+1:stpidx) = RuM(1:cstp-cstrt,1:stpidx-cstp)
          ! Row update B
          call ZGEMM('C','N',cstp-cstrt,stpidx-cstp,cstp-cstrt,cone,&
                     QM(1,1),LDQM,B(cstrt+1,cstp+1),LDAB,czero,RuM(1,1),LDR)
          B(cstrt+1:cstp,cstp+1:stpidx) = RuM(1:cstp-cstrt,1:stpidx-cstp)
          if (cstrt .gt. 0) then
            ! Column update A
            call ZGEMM('N','N',cstrt-strtidx,cstp-cstrt,cstp-cstrt,cone,&
                      A(strtidx+1,cstrt+1),LDAB,ZM(1,1),LDZM,czero,CuM(1,1),LDC)
            A(strtidx+1:cstrt,cstrt+1:cstp) = CuM(1:cstrt-strtidx,1:cstp-cstrt)
            ! Column update B
            call ZGEMM('N','N',cstrt-strtidx,cstp-cstrt,cstp-cstrt,cone,&
                      B(strtidx+1,cstrt+1),LDAB,ZM(1,1),LDZM,czero,CuM(1,1),LDC)
            B(strtidx+1:cstrt,cstrt+1:cstp) = CuM(1:cstrt-strtidx,1:cstp-cstrt)
          end if
        end if

        if (compQ) then
          call ZGEMM('N','N',LDQ,cstp-cstrt,cstp-cstrt,cone,&
                     Q(1,qci+cstrt),LDQ,QM(1,1),LDQM,czero,CuM,LDC)
          Q(1:LDQ,qci+cstrt:qci+cstp-1) = CuM(1:LDQ,1:cstp-cstrt)
        end if

        if (compZ) then
          call ZGEMM('N','N',LDZ,cstp-cstrt,cstp-cstrt,cone,&
                     Z(1,zci+cstrt),LDZ,ZM(1,1),LDZM,czero,CuM,LDC)
          Z(1:LDZ,zci+cstrt:zci+cstp-1) = CuM(1:LDZ,1:cstp-cstrt)
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
        ! Use the Tightly-packed implementation on the current active part
        ! check for deflation at front
        call z_start_aed(defl, wsc, A, B, compAB, apCurrent, apTotal, &
                         Q, compQ, qci+cstrt, Z, compZ, zci+cstrt)
        if (cstp - cstrt .gt. mc) then
          ! Compute the shifts as the eigenvalues of the trailing mcxmc block
          ! Create copy (we do not want to change the pencil here)
          QM(1:mc,1:mc) = A(cstp-mc+1:cstp,cstp-mc+1:cstp)
          ZM(1:mc,1:mc) = B(cstp-mc+1:cstp,cstp-mc+1:cstp)
          call z_RQZ1(QM, ZM, .false., 0, mc, Q, .false., 1, Z, .false., 1, &
                      alpha, beta)
          ! Introduce the mc shifts
          call z_set_first_m_poles(alpha(1:mc),beta(1:mc),A,B,compAB,apCurrent,&
                                   apTotal,Q,compQ,qci+cstrt,Z,compZ,zci+cstrt)
          ! Check again for regular deflation at front
          call z_check_start_deflation(defl,A,B,compAB,apCurrent,apTotal,&
                                       Q,compQ,qci+cstrt)
          ! Deflations along subdiagonal:
          call z_check_interior_deflations(defl,A, B, apCurrent)
          ! apCurrent is updated; consequently cstrt, and cstp
          if (defl) then
            ! Make sure that we're handling the correct range
            call tApGet(apCurrent,cstrt,cstp)
          end if

          if (cstp - cstrt .gt. mc) then
            nit = nit + mc
            ! Chase
            ! Take as many steps of size kc as possible
            j = 0
            do i = cstrt+1, cstp-mc-kc, kc
              call z_swapmk(i, mc, kc,A, B, compAB, apCurrent, apTotal, &
                            Q, compQ, qci+i, Z, compZ, zci+i-1)
              nswaps = nswaps + (mc*kc)
              j = j + 1
            end do
            ! Take the remaining steps
            i = cstp-cstrt-mc-(j*kc)-2
            call z_swapmk(cstrt+(j*kc)+1, mc,i ,A, B, compAB, apCurrent,&
                          apTotal, Q, compQ, qci+cstrt+(j*kc)+1,&
                          Z, compZ, zci+cstrt+(j*kc))
            nswaps = nswaps + (mc*i)

            ! Check for deflations at the end
            call z_stop_aed(defl, wec, A, B, compAB, apCurrent, apTotal, &
                            Q, compQ, qci+cstp-wec, Z, compZ, zci+cstp-wec-1)
            ! Compute mc poles as the eigenvalues of the leading mcxmc block
            ! Create copy (we do not want to change the pencil here)
            QM(1:mc,1:mc) = A(cstrt+1:cstrt+mc,cstrt+1:cstrt+mc)
            ZM(1:mc,1:mc) = B(cstrt+1:cstrt+mc,cstrt+1:cstrt+mc)
            call z_RQZ1(QM, ZM, .false., 0, mc, Q, .false., 1, Z, .false., 1, &
                        alpha, beta)
            ! Introduce the mc poles
            call z_set_last_m_poles(alpha(1:mc),beta(1:mc),A,B,compAB,&
                            apCurrent,apTotal,Q,compQ,cstp-mc+1,Z,compZ,cstp-mc)
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

  ! return results
  do i = 1,stpidx-strtidx
    alpha(i) = A(i+strtidx,i+strtidx)
    beta(i) = B(i+strtidx,i+strtidx)
  end do

  ! free everything
  apCurrent => null()
  call tApFree(apHead)
  call tApFree(apTotal)

  call z_deallocateSwapMem
  call z_deallocateAedMem

  end subroutine

  subroutine z_getSettings(N,m,k,we,ws)
  ! Settings for tightly-packed chasing and aggressive
  ! deflation at the start and end
  !___________________________________________________________________________
    integer, intent(in)             ::  N
    integer, intent(out)            ::  m, k, we, ws

    if (N .lt. ssth) then
      k = 1
      m = 1
      we = 1
      ws = 1
    elseif (N .lt. 60) then
      k = 4
      m = 4
      we = 6
      ws = 2
    elseif (N .lt. 180) then
      k = 10
      m = 10
      we = 12
      ws = 4
    elseif (N .lt. 1001) then
      k = 32
      m = 32
      we = 34
      ws = 6
    elseif (N .lt. 3000) then
      k = 64
      m = 64
      we = 66
      ws = 12
    elseif (N .lt. 6000) then
      k = 128
      m = 128
      we = 130
      ws = 20
    else
      k = 256
      m = 256
      we = 266
      ws = 28
    end if
  end subroutine
end module
