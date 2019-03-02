! libRQZ
! daan. camps@cs.kuleuven.be
! Description:
!   Handles the active parts during the iterations
! ___________________________________________________________________
module u_activeparts

  implicit none
  private

  public :: tAp, tApInit, tApFree, tApNext, tApPrevious, tApDelete, &
            tApInsertAfter, tApPut, tApGet, tApLength, tApI

  type :: tAp
    integer, pointer    ::  strt => null(), stp => null()
    type(tAp), pointer  ::  nextElem => null(), prevElem => null()
  end type tAp

contains

  ! Initialize the active parts list
  subroutine tApInit(self, strt, stp)
    type(tAp), pointer :: self
    integer, intent(in), optional :: strt, stp

    allocate(self)

    nullify(self%nextElem)
    nullify(self%prevElem)

    if (present(strt)) then
      allocate(self%strt)
      self%strt = strt
    else
      self%strt => null()
    end if

    if (present(stp)) then
      allocate(self%stp)
      self%stp = stp
    else
      self%stp => null()
    end if
  end subroutine tApInit

  ! Free the entire list, beginning at self
  subroutine tApFree(self)
    type(tAp), pointer :: self, current, next

    current => self
    do while (associated(current))
      next => current%nextElem
      if (associated(current%strt)) then
        deallocate(current%strt)
        nullify(current%strt)
      end if
      if (associated(current%stp)) then
        deallocate(current%stp)
        nullify(current%stp)
      end if
      deallocate(current)
      nullify(current)
      current => next
    end do
  end subroutine tApFree

  ! Return the next node after self
  function tApNext(self) result(next)
    type(tAp), pointer :: self, next
    if (.not. associated(self)) then
      next => null()
      return
    end if
    if (associated(self%nextElem)) then
      next => self%nextElem
    else
      next => null()
    end if
  end function tApNext

  ! Return the previous node before self
  function tApPrevious(self) result(previous)
    type(tAp), pointer :: self, previous

    if (associated(self%prevElem)) then
      previous => self%prevElem
    else
      previous => null()
    end if
  end function tApPrevious

  ! Insert node after self
  subroutine tApInsertAfter(self, strt, stp)
    type(tAp), intent(in), pointer :: self
    integer, intent(in), optional :: strt, stp
    type(tAp), pointer :: next

    allocate(next)
    nullify(next%nextElem)
    nullify(next%prevElem)

    if (present(strt)) then
      allocate(next%strt)
      next%strt = strt
    else
      next%strt => null()
    end if
    if (present(stp)) then
      allocate(next%stp)
      next%stp = stp
    else
      next%stp => null()
    end if

    ! make the connection
    if (associated(self%nextElem)) then
      next%nextElem => self%nextElem
      next%nextElem%prevElem => next
    end if

    self%nextElem => next
    next%prevElem => self

  end subroutine

  ! Delete node from chain (not a HEAD)
  subroutine tApDelete(self)
    type(tAp), pointer, intent(inout) :: self

    if (associated(self)) then
      ! change links
      if ((associated(self%nextElem)) .and. (associated(self%prevElem))) then
        ! middle element
        self%prevElem%nextElem => self%nextElem
        self%nextElem%prevElem => self%prevElem
        self%nextElem => null()
        self%prevElem => null()
      elseif (associated(self%nextElem)) then
        ! head element
        self%nextElem%prevElem => null()
        self%nextElem => null()
      elseif (associated(self%prevElem)) then
        ! tail element
        self%prevElem%nextElem => null()
        self%prevElem => null()
      end if

      ! delete self
      if (associated(self%strt)) then
        deallocate(self%strt)
        nullify(self%strt)
      end if
      if (associated(self%stp)) then
        deallocate(self%stp)
        nullify(self%stp)
      end if

      deallocate(self)
      nullify(self)
    end if
  end subroutine tApDelete

  ! Place the strt and stp index in self
  subroutine tApPut(self,strt, stp)
    type(tAp), pointer :: self
    integer, intent(in) :: strt, stp

    self%strt = strt
    self%stp = stp

  end subroutine tApPut

  ! Get the strt and stp index from self
  subroutine tApGet(self,strt,stp)
    type(tAp), pointer :: self
    integer, pointer, intent(out) :: strt, stp
    if (associated(self)) then
      strt => self%strt
      stp => self%stp
    else
      strt => null()
      stp => null()
    end if
  end subroutine tApGet

  ! Get the current number of active parts
  function tApLength(self) result(n)
    type(tAp), pointer :: self, current
    integer :: n

    if (associated(self)) then
      current => self
      ! go to the first element
      do while(associated(current%prevElem))
        current => current%prevElem
      end do
      n = 1
      do while(associated(current%nextElem))
        current => current%nextElem
        n = n + 1
      end do
    else
      n = 0
    end if

  end function tApLength

  ! Get the ith active part
  function tApI(self,i) result(ith)
    type(tAp), pointer :: self, ith
    integer, intent(in) :: i

    integer j

    if (i .gt. tApLength(self)) then
      ith => null()
    else
      ! go to the first
      ith => self
      do while(associated(ith%prevElem))
        ith => ith%prevElem
      end do
      ! take i-1 steps
      do j=1,i-1
        ith => ith%nextElem
      end do
    end if
  end function tApI

end module u_activeparts
