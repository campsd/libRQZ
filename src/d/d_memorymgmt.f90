! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   auxilary memory for double multishift
! ___________________________________________________________________
module d_memorymgmt
  use u_parameters
  implicit none
  real(kind=dp), allocatable   :: QM(:,:), ZM(:,:), Q4(:,:), Z4(:,:), &
                                  QW(:,:), ZW(:,:), Ru(:,:), Cu(:,:)
  public d_allocateSwapMMem, d_allocateSwap4Mem, d_allocateAedMem, &
         d_allocateRCMem, d_deallocateSwapMMem, d_deallocateSwap4Mem, &
         d_deallocateAedMem, d_deallocateRCMem, d_getLeadingDim
contains

  subroutine d_allocateSwap4Mem
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for swapping  2x2 poles
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    integer             ::  ERR

    ERR = 0
    if (.not. allocated(Q4)) then
      allocate(Q4(4,4),STAT=ERR)
    end if

    if (.not. allocated(Z4)) then
      allocate(Z4(4,4),STAT=ERR)
    end if

    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE SWAP4 MEMORY'
    end if
  end subroutine

  subroutine d_allocateSwapMMem(m)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for swapping multiple poles
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! m       integer [IN]
  !           Batch size
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    integer, intent(in) ::  m
    integer             ::  ERR

    ERR = 0
    ! First check if everything is unallocated
    if (allocated(QM)) then
      deallocate(QM)
    end if

    if (allocated(ZM)) then
      deallocate(ZM)
    end if

    if (m .gt. ssth) then
      allocate(QM(m,m),ZM(m,m),STAT=ERR)
    else
      ! In this way they can be used for the single-shift call
      allocate(QM(ssth,ssth),ZM(ssth,ssth),STAT=ERR)
    end if
    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE SWAPM MEMORY'
    end if
  end subroutine

  subroutine d_allocateAedMem(w)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for aggressive deflation
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! w       integer [IN]
  !           Window size
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    integer, intent(in) ::  w
    integer             ::  ERR

    ERR = 0
    ! First check if everything is unallocated
    if (allocated(QW)) then
      deallocate(QW)
    end if

    if (allocated(ZW)) then
      deallocate(ZW)
    end if

    allocate(QW(w,w),ZW(w,w),STAT=ERR)

    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE AED MEMORY'
    end if
  end subroutine

  subroutine d_allocateRCMem(k,N)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for row and column updates
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! k       integer [IN]
  !           Update size
  ! N       integer [IN]
  !           Problem size
  !___________________________________________________________________________
  ! last edit: September 21, 2018
    integer, intent(in)   ::  k, N
    integer               ::  ERR

    ERR = 0
    if (allocated(Ru) .and. allocated(Cu)) then
      if (size(Ru,1) .lt. k) then
        deallocate(Cu, Ru,STAT=ERR)
        if (k .gt. ssth) then
          allocate(Cu(N,k),Ru(k,N),STAT=ERR)
        else
          allocate(Cu(N,ssth),Ru(ssth,N),STAT=ERR)
        end if
      else
        ERR = 0
      end if
    elseif (.not. allocated(Ru) .and. .not. allocated(Cu)) then
      if (k .gt. ssth) then
        allocate(Cu(N,k),Ru(k,N),STAT=ERR)
      else
        allocate(Cu(N,ssth),Ru(ssth,N),STAT=ERR)
      end if
    else
      stop 'UNEXPECTED CASE'
    end if

    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE R/C UPDATE MEMORY'
    end if
  end subroutine

  subroutine d_deallocateSwap4Mem
    integer ERR

    ERR = 0
    if (allocated(Q4)) then
      deallocate(Q4,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE Q4 SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(Z4)) then
      deallocate(Z4,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE Z4 SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine d_deallocateSwapMMem
    integer ERR

    ERR = 0
    if (allocated(QM)) then
      deallocate(QM,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE QM SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(ZM)) then
      deallocate(ZM,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE ZM SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine d_deallocateAedMem
    integer ERR

    ERR = 0
    if (allocated(QW)) then
      deallocate(QW, STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE QW AED MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(ZW)) then
      deallocate(ZW, STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE ZW AED MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine d_deallocateRCMem
    integer ERR

    ERR = 0
    if (allocated(Ru)) then
      deallocate(Ru, STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE Ru MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(Cu)) then
      deallocate(Cu,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE Cu MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine d_getLeadingDim(LD,M)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Gets the leading dimension of a 2D array (useful for LAPACK)
  ! ARGUMENTS
  !___________________________________________________________________________
  ! LD      integer [OUT]
  !           Leading dimension
  ! M       double array [IN]
  !           2D array
  !___________________________________________________________________________
  ! last edit: September 25, 2018
    integer, intent(out)          ::  LD
    real(kind=dp), intent(in)     ::  M(:,:)

    LD = size(M,1)
  end subroutine
end module
