!
! The optimizer: libcint's CINTOpt, which caches across shell quartets what
! the drivers would otherwise recompute for every one of them.
!
! Ported from src/optimizer.c.  cint_screen already holds the three routines
! that do the actual computing -- cint_log_max_pgto_coeff, cint_non0coeff_-
! byshell, cint_set_pairdata -- because the unoptimised drivers call them
! directly.  This module is the cache around them, plus the Cartesian index
! map, which is the fourth thing the C caches and the only one that needs a
! fake basis to build.
!
! Nothing here changes a number.  The C's own comment on the file says the
! optimizer cannot really speed things up if the driver is called only a few
! hundred times; the point is that a real calculation calls it millions of
! times, and then the four caches are most of the per-quartet setup.
!
! Ragged pointer arrays become one flat array plus offsets, as everywhere
! else in this port.  The C's NOVALUE sentinel -- a pointer to 0xffff... --
! becomes OPT_NOVALUE, a value no offset can take.
!
! WHAT IS DELIBERATELY NOT CACHED.  The C's five specialised 2e loops
! (CINT2e_n111_loop and friends), selected by which contraction counts are
! one, are a separate optimisation that happens to sit behind the same
! `opt != NULL` test.  They are not here: this module makes the generic loop
! cheaper, which is the part that applies to every arity.
!
module cint_opt
   use cint_const,     only: dp
   use cint_envs
   use cint_screen,    only: pair_data, cint_opt_t, cint_set_pairdata, &
                             cint_log_max_pgto_coeff, cint_non0coeff_byshell, &
                             OPT_NOVALUE, LMAX1, MAX_PGTO_FOR_PAIRDATA
   use cint_bas,       only: ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                             ATOM_OF, PTR_COORD
   implicit none
   private

   public :: cint_del_optimizer
   public :: opt_set_log_maxc, opt_set_non0coeff, opt_setij, opt_gen_idx, opt_finish

   ! The two procedures gen_idx needs.  The C passes them as bare `void (*)()`
   ! and casts at the call; Fortran wants the shapes written down, which is
   ! the only reason these exist.
   !
   ! They are arguments rather than a `use` because the init routines live in
   ! the modules that own each arity -- cint_3c1e, cint_1e_grids -- and those
   ! modules have optimizer entry points of their own.  Reaching for the init
   ! routines from here would make the dependency a cycle; taking them as
   ! arguments leaves this module underneath all of them.
   abstract interface
      subroutine opt_init_iface(envs, ng, shls, atm, natm, bas, nbas, env)
         import :: dp, cint_env_vars
         type(cint_env_vars), intent(inout) :: envs
         integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
         integer,  target :: atm(0:), bas(0:)
         real(dp), target :: env(0:)
      end subroutine opt_init_iface
      subroutine opt_index_iface(idx, envs)
         import :: cint_env_vars
         integer, intent(out) :: idx(0:)
         type(cint_env_vars), intent(in) :: envs
      end subroutine opt_index_iface
   end interface

contains

   subroutine cint_del_optimizer(opt)
      type(cint_opt_t), intent(inout) :: opt
      opt%active = .false.
      opt%nbas = 0; opt%order = 0; opt%ng = 0
      if (allocated(opt%log_maxc))      deallocate(opt%log_maxc)
      if (allocated(opt%log_maxc_off))  deallocate(opt%log_maxc_off)
      if (allocated(opt%non0ctr))       deallocate(opt%non0ctr)
      if (allocated(opt%non0ctr_off))   deallocate(opt%non0ctr_off)
      if (allocated(opt%sortedidx))     deallocate(opt%sortedidx)
      if (allocated(opt%sortedidx_off)) deallocate(opt%sortedidx_off)
      if (allocated(opt%pd))            deallocate(opt%pd)
      if (allocated(opt%pd_off))        deallocate(opt%pd_off)
      if (allocated(opt%idx))           deallocate(opt%idx)
      if (allocated(opt%idx_off))       deallocate(opt%idx_off)
   end subroutine cint_del_optimizer

   ! ---- the three caches cint_screen can fill --------------------------

   subroutine opt_set_log_maxc(opt, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: bas(0:), nbas
      real(dp), intent(in) :: env(0:)
      integer :: i, iprim, ictr, tot, off
      tot = 0
      do i = 0, nbas - 1
         tot = tot + bas(8*i + NPRIM_OF)
      end do
      if (tot == 0) return
      allocate(opt%log_maxc(0:tot-1), opt%log_maxc_off(0:nbas-1))
      off = 0
      do i = 0, nbas - 1
         iprim = bas(8*i + NPRIM_OF)
         ictr  = bas(8*i + NCTR_OF)
         opt%log_maxc_off(i) = off
         call cint_log_max_pgto_coeff(opt%log_maxc(off:), env, &
                                      bas(8*i + PTR_COEFF), iprim, ictr)
         off = off + iprim
      end do
   end subroutine opt_set_log_maxc

   subroutine opt_set_non0coeff(opt, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: bas(0:), nbas
      real(dp), intent(in) :: env(0:)
      integer :: i, iprim, ictr, tot, totc, offn, offs
      tot = 0; totc = 0
      do i = 0, nbas - 1
         tot  = tot  + bas(8*i + NPRIM_OF)
         totc = totc + bas(8*i + NPRIM_OF) * bas(8*i + NCTR_OF)
      end do
      if (tot == 0) return
      allocate(opt%non0ctr(0:tot-1), opt%non0ctr_off(0:nbas-1))
      allocate(opt%sortedidx(0:totc-1), opt%sortedidx_off(0:nbas-1))
      offn = 0; offs = 0
      do i = 0, nbas - 1
         iprim = bas(8*i + NPRIM_OF)
         ictr  = bas(8*i + NCTR_OF)
         opt%non0ctr_off(i)   = offn
         opt%sortedidx_off(i) = offs
         call cint_non0coeff_byshell(opt%sortedidx(offs:), opt%non0ctr(offn:), &
                                     env, bas(8*i + PTR_COEFF), iprim, ictr)
         offn = offn + iprim
         offs = offs + iprim * ictr
      end do
   end subroutine opt_set_non0coeff

   ! Primitive-pair data for every shell pair.  The C walks j <= i and fills
   ! the transpose itself rather than computing it twice; so does this, and
   ! the transpose is stored rather than indexed around because the driver
   ! reads pairs in the (ip, jp) order of whichever shell is the bra.
   subroutine opt_setij(opt, ng, atm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), atm(0:), bas(0:), nbas
      real(dp), intent(in) :: env(0:)
      integer  :: i, j, ip, jp, iprim, jprim, li, lj, tot, off, o_ij, ijkl_inc
      real(dp) :: expcutoff, rr, ri(0:2), rj(0:2)
      logical  :: empty

      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         expcutoff = EXPCUTOFF
      else
         expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if
      if (.not. allocated(opt%log_maxc)) call opt_set_log_maxc(opt, bas, nbas, env)

      tot = 0
      do i = 0, nbas - 1
         tot = tot + bas(8*i + NPRIM_OF)
      end do
      ! The C gives up on pairdata above MAX_PGTO_FOR_PAIRDATA because the
      ! table is quadratic in the primitive count; so does this, and the
      ! driver then falls back to computing pairs itself.
      if (tot == 0 .or. tot > MAX_PGTO_FOR_PAIRDATA) return

      allocate(opt%pd(0:tot*tot-1), opt%pd_off(0:nbas*nbas-1))
      opt%pd_off = OPT_NOVALUE

      if (ng(IINC) + ng(JINC) > ng(KINC) + ng(LINC)) then
         ijkl_inc = ng(IINC) + ng(JINC)
      else
         ijkl_inc = ng(KINC) + ng(LINC)
      end if

      off = 0
      do i = 0, nbas - 1
         ri = env(atm(6*bas(8*i + ATOM_OF) + PTR_COORD) : &
                  atm(6*bas(8*i + ATOM_OF) + PTR_COORD) + 2)
         iprim = bas(8*i + NPRIM_OF)
         li    = bas(8*i + ANG_OF)
         do j = 0, i
            rj = env(atm(6*bas(8*j + ATOM_OF) + PTR_COORD) : &
                     atm(6*bas(8*j + ATOM_OF) + PTR_COORD) + 2)
            jprim = bas(8*j + NPRIM_OF)
            lj    = bas(8*j + ANG_OF)
            rr = sum((ri - rj)**2)
            empty = cint_set_pairdata(opt%pd(off:), env, bas(8*i + PTR_EXP), &
                                      env, bas(8*j + PTR_EXP), ri, rj, &
                                      opt%log_maxc(opt%log_maxc_off(i):), &
                                      opt%log_maxc(opt%log_maxc_off(j):), &
                                      li + ijkl_inc, lj, iprim, jprim, rr, &
                                      expcutoff, env(PTR_RANGE_OMEGA))
            ! The C keeps shell pair (0,0) whether or not it is empty -- the
            ! buffer head has to point somewhere -- and that asymmetry is
            ! inherited rather than tidied, because the driver's fallback
            ! path depends on which pairs report NOVALUE.
            if (i == 0 .and. j == 0) then
               opt%pd_off(0) = off
               off = off + iprim * jprim
            else if (.not. empty) then
               o_ij = off
               opt%pd_off(i*nbas + j) = off
               off = off + iprim * jprim
               if (i /= j) then
                  opt%pd_off(j*nbas + i) = off
                  do ip = 0, iprim - 1
                  do jp = 0, jprim - 1
                     opt%pd(off) = opt%pd(o_ij + jp*iprim + ip)
                     off = off + 1
                  end do
                  end do
               end if
            end if
         end do
      end do
   end subroutine opt_setij

   ! ---- the Cartesian index map ----------------------------------------

   ! index_xyz depends on the angular momenta and nothing else, so the C
   ! builds it against a fake basis of one shell per l -- which is why this
   ! is the one cache that cannot be filled from the real shells.
   subroutine opt_gen_idx(opt, finit, findex, order, l_allow_in, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      procedure(opt_init_iface)  :: finit
      procedure(opt_index_iface) :: findex
      integer,  intent(in) :: order, l_allow_in, ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:)
      integer,  intent(in), target :: bas(0:)
      real(dp), intent(in) :: env(0:)

      integer, allocatable, target :: fakebas(:)
      integer, allocatable :: fatm(:)
      real(dp), allocatable :: fenv(:)
      type(cint_env_vars) :: envs
      integer :: i, j, k, l, max_l, fakenbas, l_allow, cumcart, nkey, nbuf, ptr, off
      integer :: shls(0:3)

      max_l = 0
      do i = 0, nbas - 1
         max_l = max(max_l, bas(8*i + ANG_OF))
      end do
      l_allow = min(max_l, l_allow_in)
      fakenbas = max_l + 1

      allocate(fakebas(0:8*fakenbas-1))
      fakebas = 0
      do i = 0, max_l
         fakebas(8*i + ANG_OF) = i
      end do
      ! index_xyz reads only ANG_OF out of bas, but the envs setup reads
      ! coordinates out of atm and env, so both have to be present and
      ! harmless.  One atom at the origin is enough.
      allocate(fatm(0:5)); fatm = 0; fatm(PTR_COORD) = PTR_ENV_START
      allocate(fenv(0:size(env)-1)); fenv = env
      fenv(PTR_ENV_START:PTR_ENV_START+2) = 0.0_dp

      cumcart = (l_allow+1) * (l_allow+2) * (l_allow+3) / 6
      nkey = LMAX1
      nbuf = cumcart
      do i = 2, order
         nkey = nkey * LMAX1
         nbuf = nbuf * cumcart
      end do
      allocate(opt%idx(0:3*nbuf-1), opt%idx_off(0:nkey-1))
      opt%idx_off = OPT_NOVALUE

      shls = 0
      off = 0
      select case (order)
      case (2)
         do i = 0, l_allow
         do j = 0, l_allow
            shls(0) = i; shls(1) = j
            call finit(envs, ng, shls, fatm, 1, fakebas, fakenbas, fenv)
            ptr = i*LMAX1 + j
            opt%idx_off(ptr) = off
            call findex(opt%idx(off:), envs)
            off = off + envs%nf * 3
         end do
         end do
      case (3)
         do i = 0, l_allow
         do j = 0, l_allow
         do k = 0, l_allow
            shls(0) = i; shls(1) = j; shls(2) = k
            call finit(envs, ng, shls, fatm, 1, fakebas, fakenbas, fenv)
            ptr = (i*LMAX1 + j)*LMAX1 + k
            opt%idx_off(ptr) = off
            call findex(opt%idx(off:), envs)
            off = off + envs%nf * 3
         end do
         end do
         end do
      case default
         do i = 0, l_allow
         do j = 0, l_allow
         do k = 0, l_allow
         do l = 0, l_allow
            shls(0) = i; shls(1) = j; shls(2) = k; shls(3) = l
            call finit(envs, ng, shls, fatm, 1, fakebas, fakenbas, fenv)
            ptr = ((i*LMAX1 + j)*LMAX1 + k)*LMAX1 + l
            opt%idx_off(ptr) = off
            call findex(opt%idx(off:), envs)
            off = off + envs%nf * 3
         end do
         end do
         end do
         end do
      end select
   end subroutine opt_gen_idx

   pure subroutine opt_finish(opt, order, ng, nbas)
      type(cint_opt_t), intent(inout) :: opt
      integer, intent(in) :: order, ng(0:), nbas
      opt%order = order
      opt%ng(0:min(7, size(ng)-1)) = ng(0:min(7, size(ng)-1))
      opt%nbas = nbas
      opt%active = .true.
   end subroutine opt_finish

end module cint_opt
