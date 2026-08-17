!
! Scratch space for the integral drivers.
!
! Replaces the C's MALLOC_INSTACK protocol:
!
!     #define MALLOC_INSTACK(var, n) \
!             var = (void *)(((uintptr_t)cache + 7) & (-(uintptr_t)8)); \
!             cache = (double *)(var + (n));
!
! which is a bump allocator over a buffer the caller supplies, plus the
! two-call convention where passing a NULL output means "tell me how big the
! buffer must be" instead of "compute the integral".
!
! That design is not an accident and is kept: it is what keeps malloc out of
! the shell-quartet loop, and it is what makes the library thread-safe without
! locks, since each thread brings its own buffer.
!
! ⚠️ WHAT IS NOT DONE HERE, AND WHY.  The obvious Fortran translation is an
! automatic array sized from the same formula.  That would be wrong: the C
! checks its own cache sizes against INT32_MAX, so these buffers reach tens of
! megabytes on large shell quartets, and an automatic array of that size is
! placed on the stack and will overflow it.  (PORT_TO_FORTRAN.md 3.5.)  Hence a
! persistent allocatable that grows to the high-water mark and is then reused.
!
! Two buffers, real and integer, rather than the C's single one.  The C stores
! both in the same `double *` and casts; keeping them apart is type-safe and
! removes the alignment arithmetic the macro above is doing.
!
module cint_workspace
   use cint_const,  only: dp
   use cint_screen, only: pair_data, cint_opt_t, OPT_NOVALUE
   implicit none
   private

   public :: cint_ws, ws_reset, ws_ensure, ws_ensure_pd, ws_alloc_d, ws_alloc_i, ws_release
   public :: ws_mark, ws_rewind
   public :: ws_opt_log_maxc, ws_opt_non0, ws_opt_idx

   type :: cint_ws
      real(dp), allocatable :: d(:)
      integer,  allocatable :: i(:)
      ! The primitive-pair table.  A third buffer rather than a slice of d
      ! because it holds a derived type; it lives here for the same reason
      ! the other two do -- one allocate per shell quartet was measurable.
      type(pair_data), allocatable :: p(:)
      integer :: dtop = 0     ! next free slot in d
      integer :: itop = 0     ! next free slot in i
      ! The optimizer, if the caller built one.
      !
      ! The C threads CINTOpt through the entry points, one argument per
      ! call.  Here it rides on the workspace instead, which is already the
      ! per-thread context and is already passed to everything -- so an
      ! optimizer costs no change to any of the 610 entry point signatures,
      ! and a caller that never builds one pays nothing.
      !
      ! POINTER, not allocatable: the caller owns it and typically keeps one
      ! per integral type, swapping which is current between calls.
      type(cint_opt_t), pointer :: opt => null()
   end type cint_ws

contains

   ! ---- reading the optimizer's caches into the workspace ---------------
   !
   ! The drivers all lay these out the same way: one contiguous block per
   ! shell in ws%d / ws%i, at a fixed offset.  The cache has them scattered
   ! by shell index instead, so they are copied rather than pointed at --
   ! a few dozen doubles against the logarithms, or a few dozen integers
   ! against the sort, that the copy replaces.

   subroutine ws_opt_log_maxc(ws, dst, sh, nprim)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in) :: dst, sh, nprim
      integer :: o
      o = ws%opt%log_maxc_off(sh)
      ws%d(dst:dst+nprim-1) = ws%opt%log_maxc(o:o+nprim-1)
   end subroutine ws_opt_log_maxc

   ! The cached Cartesian index map, if this quartet's key is in the table.
   ! Above the l_allow gen_idx was built with there is a hole, which the C
   ! indexes into regardless and dereferences NULL; here the miss is
   ! reported and the caller builds the map the usual way.
   logical function ws_opt_idx(ws, dst, key, n) result(hit)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in) :: dst, key, n
      integer :: o
      hit = .false.
      if (key < 0 .or. key > ubound(ws%opt%idx_off, 1)) return
      o = ws%opt%idx_off(key)
      if (o == OPT_NOVALUE) return
      ws%i(dst:dst+n-1) = ws%opt%idx(o:o+n-1)
      hit = .true.
   end function ws_opt_idx

   subroutine ws_opt_non0(ws, onon, osrt, sh, nprim, nctr, doff_n, doff_s)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in) :: onon, osrt, sh, nprim, nctr, doff_n, doff_s
      integer :: sn, ss
      sn = ws%opt%non0ctr_off(sh)
      ss = ws%opt%sortedidx_off(sh)
      ws%i(onon+doff_n : onon+doff_n+nprim-1) = ws%opt%non0ctr(sn:sn+nprim-1)
      ws%i(osrt+doff_s : osrt+doff_s+nprim*nctr-1) = &
         ws%opt%sortedidx(ss:ss+nprim*nctr-1)
   end subroutine ws_opt_non0

   ! Make sure there is room for nd doubles and ni integers, growing if there
   ! is not.
   !
   ! ⚠️ THIS MUST BE CALLED BEFORE THE BUFFERS ARE PASSED ANYWHERE.  ws_alloc_*
   ! used to grow on demand, which is fine in C -- realloc moves memory and the
   ! caller re-reads the pointer -- and is not fine here: a driver passes ws%d
   ! as an actual argument and then calls ws_alloc_d, whose move_alloc
   ! invalidates the dummy the callee is still writing through.  The symptom
   ! was components of a multi-component integral silently coming back zero,
   ! depending on whether that particular shell pair happened to make the
   ! buffer grow.
   !
   ! So: size once, up front, from the same formula the C uses for its cache,
   ! and let ws_alloc_* only ever hand out what is already there.
   subroutine ws_ensure(ws, nd, ni)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in) :: nd, ni
      if (.not. allocated(ws%d)) then
         allocate(ws%d(0:max(nd, 1024) - 1))
      else if (nd > size(ws%d)) then
         deallocate(ws%d)
         allocate(ws%d(0:nd - 1))
      end if
      if (.not. allocated(ws%i)) then
         allocate(ws%i(0:max(ni, 1024) - 1))
      else if (ni > size(ws%i)) then
         deallocate(ws%i)
         allocate(ws%i(0:ni - 1))
      end if
      ws%dtop = 0
      ws%itop = 0
   end subroutine ws_ensure

   ! Room for np primitive pairs.  Separate from ws_ensure because the pair
   ! count is known only after the shells are read, which is after the driver
   ! has already sized d and i.
   subroutine ws_ensure_pd(ws, np)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in) :: np
      if (.not. allocated(ws%p)) then
         allocate(ws%p(0:max(np, 256) - 1))
      else if (np > size(ws%p)) then
         deallocate(ws%p)
         allocate(ws%p(0:np - 1))
      end if
   end subroutine ws_ensure_pd

   ! Rewind to empty without giving the memory back, which is the point: the
   ! next shell quartet reuses it.
   subroutine ws_reset(ws)
      type(cint_ws), intent(inout) :: ws
      ws%dtop = 0
      ws%itop = 0
   end subroutine ws_reset

   ! Mark and rewind, for scratch that is used and finished with inside a
   ! loop.  The C gets this free: it passes `cache` by value to each c2s call,
   ! so every one of them bump-allocates from the same place and the space is
   ! implicitly reused.  A monotonically rising top does not, and a
   ! nine-component integral duly exhausted the buffer -- which the whole-
   ! catalogue sweep found and nothing narrower would have.
   pure subroutine ws_mark(ws, dmark, imark)
      type(cint_ws), intent(in)  :: ws
      integer,       intent(out) :: dmark, imark
      dmark = ws%dtop
      imark = ws%itop
   end subroutine ws_mark

   pure subroutine ws_rewind(ws, dmark, imark)
      type(cint_ws), intent(inout) :: ws
      integer,       intent(in)    :: dmark, imark
      ws%dtop = dmark
      ws%itop = imark
   end subroutine ws_rewind

   subroutine ws_release(ws)
      type(cint_ws), intent(inout) :: ws
      if (allocated(ws%d)) deallocate(ws%d)
      if (allocated(ws%i)) deallocate(ws%i)
      ws%dtop = 0
      ws%itop = 0
   end subroutine ws_release

   ! Hand out n doubles and return the offset of the first.  Callers index as
   ! ws%d(off + k) with k from 0, matching the C's pointer arithmetic.
   subroutine ws_alloc_d(ws, n, off)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in)  :: n
      integer, intent(out) :: off
      integer :: need

      off = ws%dtop
      need = ws%dtop + max(n, 0)
      if (need > size(ws%d)) then
         ! Never grow here -- see ws_ensure.  Reaching this means the size
         ! estimate is wrong, which is a bug to fix rather than paper over.
         error stop "cint workspace: double buffer too small; ws_ensure under-counted"
      end if
      ws%dtop = need
   end subroutine ws_alloc_d

   subroutine ws_alloc_i(ws, n, off)
      type(cint_ws), intent(inout) :: ws
      integer, intent(in)  :: n
      integer, intent(out) :: off
      integer :: need

      off = ws%itop
      need = ws%itop + max(n, 0)
      if (need > size(ws%i)) then
         error stop "cint workspace: integer buffer too small; ws_ensure under-counted"
      end if
      ws%itop = need
   end subroutine ws_alloc_i

end module cint_workspace
