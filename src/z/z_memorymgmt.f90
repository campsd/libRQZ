! libRQZ
! author: daan. camps@cs.kuleuven.be
! Description:
!   auxilary memory for complex double multishift
! ___________________________________________________________________
module z_memorymgmt
  use u_parameters
  implicit none
  complex(kind=dp), allocatable   :: QM(:,:), ZM(:,:), RuM(:,:), CuM(:,:), &
                                     QW(:,:), ZW(:,:), RuW(:,:), CuW(:,:)
  public z_allocateSwapMem, z_allocateAedMem, &
         z_deallocateSwapMem, z_deallocateAedMem, z_getLeadingDim
contains

  subroutine z_allocateSwapMem(m, N)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for swapping poles
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! m       integer [IN]
  !           Batch size
  ! N       integer [IN]
  !           Problem size
  !___________________________________________________________________________
  ! last edit: July 25, 2018
    integer, intent(in) ::  m, N
    integer             ::  ERR

    ! First check if everything is unallocated
    if (allocated(QM)) then
      deallocate(QM)
    end if

    if (allocated(ZM)) then
      deallocate(ZM)
    end if

    if (allocated(RuM)) then
      deallocate(RuM)
    end if

    if (allocated(CuM)) then
      deallocate(CuM)
    end if

    if (m+1 .gt. ssth) then
      allocate(QM(m+1,m+1),ZM(m+1,m+1),CuM(N,m+2),RuM(m+1,N),STAT=ERR)
    else
      ! In this way they can be used for the single-shift call
      allocate(QM(ssth,ssth),ZM(ssth,ssth),RuM(ssth,N),CuM(N,ssth),STAT=ERR)
    end if
    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE SWAP MEMORY'
    end if
  end subroutine

  subroutine z_allocateAedMem(w,N)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Allocates the auxilary arrays for aggressive deflation
  !
  ! ARGUMENTS
  !___________________________________________________________________________
  ! w       integer [IN]
  !           Window size
  ! N       integer [IN]
  !           Problem size
  !___________________________________________________________________________
  ! last edit: July 25, 2018
    integer, intent(in) ::  w, N
    integer             ::  ERR

    ! First check if everything is unallocated
    if (allocated(QW)) then
      deallocate(QW)
    end if

    if (allocated(ZW)) then
      deallocate(ZW)
    end if

    if (allocated(RuW)) then
      deallocate(RuW)
    end if

    if (allocated(CuW)) then
      deallocate(CuW)
    end if

    allocate(QW(w+1,w+1),ZW(w+1,w+1),RuW(w+1,N),CuW(N,w+1),STAT=ERR)

    if (ERR .gt. 0) then
      stop 'FAILED TO ALLOCATE AED MEMORY'
    end if
  end subroutine

  subroutine z_deallocateSwapMem
    integer ERR
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
    if (allocated(RuM)) then
      deallocate(RuM, STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE RuM SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(CuM)) then
      deallocate(CuM,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE CuM SWAP MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine z_deallocateAedMem
    integer ERR
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
    if (allocated(RuW)) then
      deallocate(RuW, STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE RuW AED MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
    if (allocated(CuW)) then
      deallocate(CuW,STAT=ERR)
      if (ERR .gt. 0) then
        write(*,*) 'FAILED TO DEALLOCATE CuW AED MEMORY, ERROR CODE', ERR
        stop
      end if
    end if
  end subroutine

  subroutine z_getLeadingDim(LD,M)
  ! DESCRIPTION
  !___________________________________________________________________________
  ! Gets the leading dimension of a 2D array (useful for LAPACK)
  ! ARGUMENTS
  !___________________________________________________________________________
  ! LD      integer [OUT]
  !           Leading dimension
  ! M       complex double array [IN]
  !           2D array
  !___________________________________________________________________________
  ! last edit: July 25, 2018
    integer, intent(out)          ::  LD
    complex(kind=dp), intent(in)  ::  M(:,:)

    LD = size(M,1)
  end subroutine
end module
