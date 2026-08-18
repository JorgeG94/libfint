!
! libcint's C entry points, as Fortran procedures with libcint's C names.
!
! WHY
! ---
! A Fortran program that uses libcint has two ways in.  The tidy one is the
! libcint_fortran module, which this library also provides.  The other is to
! write a `bind(C, name="int1e_drinv_sph")` interface and call the C entry
! point directly, which is what a caller does the moment it wants an integral
! the module never exposed.  metalquicha does both -- 25 procedures through
! the module and 14 raw bindings for multipoles and ESP -- and it is not
! unusual in that.
!
! Without this file, such a caller cannot link libfint at all, and the failure
! is a linker error listing symbols it never knew were C.  With it, the swap is
! a build option and nothing in the caller changes.
!
! WHAT IS AND IS NOT HERE
! -----------------------
! Not all 620 entry points: the ones a consumer actually drives, which now
! means everything a Fock build and its gradient need.  Overlap, kinetic and
! nuclear for the core Hamiltonian; int2e itself; int3c2e and int2c2e for
! density fitting; int2e_ip1 and the four one-electron first derivatives for a
! gradient; the optimizer lifecycle; and the shell-dimension helpers a driver
! needs to size its own buffers.  Both angular conventions throughout.
!
! Until this grew, every export here was one-electron.  There was no int2e in
! any form, so a C consumer could link the library and then fail to resolve the
! one integral it came for.
!
! Two shapes, both of which libcint exports for every integral:
!
!   the CINT3 form   name(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
!   the CINT2 form   cname(buf, shls, atm, natm, bas, nbas, env[, opt])
!
! dims and cache may be NULL, and the CINT2 form has no dims at all, so both
! compute the natural block shape from the shells -- exactly as libcint's own
! c2s layer does.  The CINT2 form carries opt for the two-electron family and
! not for the one-electron one, matching libcint: a one-electron integral has
! no optimizer to carry.
!
! THE OPTIMIZER IS HONOURED NOW, not ignored.  It used to be accepted and
! dropped, on the reasoning that libcint's CINTOpt is a foreign struct this
! library cannot read.  That is true and beside the point: the handle a caller
! passes in came from the *_optimizer entry points below, which are this
! library's own, so what it points at is a Fortran cint_opt_t and c_f_pointer
! reads it back exactly.  A program links libcint or libfint, never both, so a
! genuinely foreign CINTOpt cannot arrive here.  NULL still means "no
! optimizer", as it does in the C.
!
module cint_c_abi
   use iso_c_binding, only: c_int, c_double, c_ptr, c_f_pointer, c_associated, c_null_ptr
   use cint_const,     only: dp
   use cint_bas,       only: cint_cgto_cart, cint_cgto_spheric
   use cint_workspace, only: cint_ws
   use cint_gen_intor1, only: int1e_r_cart, int1e_r_sph, &
                              int1e_rr_cart, int1e_rr_sph, &
                              int1e_rrr_cart, int1e_rrr_sph, &
                              int1e_drinv_cart, int1e_drinv_sph
   use cint_gen_hess,   only: int1e_ipiprinv_cart, int1e_ipiprinv_sph, &
                              int1e_iprinvip_cart, int1e_iprinvip_sph
   use cint_1e_grids,   only: int1e_grids_cart, int1e_grids_sph
   ! The Fock-build set: overlap, kinetic and nuclear for the core Hamiltonian,
   ! the four-centre integral itself, the two- and three-centre ones density
   ! fitting needs, and the first derivatives a gradient needs.
   use cint_1e,         only: int1e_ovlp_cart, int1e_ovlp_sph, &
                              int1e_nuc_cart, int1e_nuc_sph
   use cint_gen_intor1, only: int1e_kin_cart, int1e_kin_sph
   use cint_2e,         only: int2e_cart, int2e_sph
   use cint_3c2e,       only: int3c2e_cart, int3c2e_sph, &
                              int2c2e_cart, int2c2e_sph
   use cint_gen_grad2,  only: int2e_ip1_cart, int2e_ip1_sph
   use cint_gen_grad1,  only: int1e_ipovlp_cart, int1e_ipovlp_sph, &
                              int1e_ipkin_cart, int1e_ipkin_sph, &
                              int1e_ipnuc_cart, int1e_ipnuc_sph, &
                              int1e_iprinv_cart, int1e_iprinv_sph
   ! Shell dimensions, and the optimizer lifecycle.  The pool, the claim and
   ! the c_loc handle already exist for the Fortran-facing module; exporting
   ! them under libcint's C names is all this adds, so both front doors hand
   ! out the same object and either can free it.
   use cint_envs,       only: cint_opt_t
   use libcint_fortran, only: libcint_2e_cart_optimizer, libcint_2e_sph_optimizer, &
                              libcint_3c2e_cart_optimizer, libcint_3c2e_sph_optimizer, &
                              libcint_2c2e_cart_optimizer, libcint_2c2e_sph_optimizer, &
                              libcint_del_optimizer, &
                              libcint_cgto_cart, libcint_cgto_sph, &
                              libcint_tot_cgto_cart, libcint_tot_cgto_sph
   implicit none
   private

   ! `integer`, not `integer(c_int)`, for the arrays that come back from
   ! c_f_pointer -- and the assertion below is what makes that legal.
   !
   ! These entry points used to declare them `integer(c_int)` and then write
   ! `int(pbas)` at every use.  The kinds are the same on every supported
   ! target, so the conversion could not change a value; what it did was make
   ! each use a whole-array expression, and gfortran materialises those.
   ! `-Warray-temporaries` counted 195 of them in this file, four per call on
   ! `bas` alone.  That put an O(nbas) copy in front of an integral whose own
   ! cost is O(1) in nbas -- so a C caller building a full matrix paid O(nbas^3)
   ! where libcint pays O(nbas^2), and the penalty grew with the molecule
   ! rather than staying a constant overhead.
   !
   ! Declaring the pointers as default integer removes the conversion instead
   ! of the copy, which is the same thing done properly: c_f_pointer requires
   ! an interoperable type, and default integer is interoperable exactly when
   ! it is c_int.  That is asserted rather than assumed.
   integer, parameter :: c_abi_int_is_c_int = 1 / merge(1, 0, kind(1) == c_int)
      !! Compile-time only: a division by zero here means default integer is
      !! not c_int on this target, and every c_f_pointer below is invalid.
      !! `include/libfint.f90` carries the same assertion for the same reason.

   ! ATM_SLOTS and BAS_SLOTS, and the two bas fields the extent needs.
   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: NPRIM_OF = 2, NCTR_OF = 3, PTR_EXP = 5, PTR_COEFF = 6
   integer, parameter :: PTR_ENV_START = 20

   ! One workspace per thread, reused.  The C reaches for its own stack on
   ! every call; this grows once and then never allocates again.
   !
   ! THREADPRIVATE is not optional.  Callers drive these from inside
   ! `!$omp parallel` -- metalquicha's AO evaluation does -- and a shared
   ! workspace means every thread writing over every other one's scratch.
   ! That does not fail cleanly: it segfaults somewhere unrelated, which is
   ! exactly how this was found.  The libcint_fortran module marks its own
   ! workspace the same way.
   type(cint_ws), save :: ws
   !$omp threadprivate(ws)

contains

   ! Point the workspace at the caller's optimizer for the duration of one call.
   !
   ! The C passes CINTOpt as an argument to every entry point; here it rides on
   ! the workspace, so the front door is where the two conventions meet.  A null
   ! pointer means no optimizer, which is what libcint's own API says a NULL
   ! CINTOpt means -- so a caller that never builds one still works, just
   ! without the shell-pair caching.
   !
   ! This is only sound because the pointer came from *this* library: the
   ! optimizer entry points below hand back a c_loc of a Fortran cint_opt_t.
   ! A CINTOpt built by the real libcint and passed in here would be a foreign
   ! struct and this would read it as the wrong type -- but that cannot happen,
   ! because a program links one of the two libraries, not both.
   subroutine bind_opt(opt)
      type(c_ptr), value :: opt
      type(cint_opt_t), pointer :: p
      ws%opt => null()
      if (.not. c_associated(opt)) return
      call c_f_pointer(opt, p)
      ws%opt => p
   end subroutine bind_opt

   ! How long env has to be for these shells: the largest offset any of them
   ! reaches, which is what lets a c_f_pointer give env a real extent.
   pure integer function env_len(shls, nsh, bas, nbas) result(n)
      integer,        intent(in) :: shls(:), bas(0:)
      integer,        intent(in) :: nsh, nbas
      integer :: i, sh, np, nc
      n = PTR_ENV_START
      do i = 1, nsh
         sh = shls(i)
         np = bas(BAS_SLOTS*sh + NPRIM_OF)
         nc = bas(BAS_SLOTS*sh + NCTR_OF)
         n = max(n, bas(BAS_SLOTS*sh + PTR_EXP)   + np)
         n = max(n, bas(BAS_SLOTS*sh + PTR_COEFF) + np*nc)
      end do
   end function env_len

   ! int1e_r_cart, 3 component(s) per shell pair
   function c_int1e_r_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_r_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_r_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_r_cart

   ! int1e_r_sph, 3 component(s) per shell pair
   function c_int1e_r_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_r_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_r_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_r_sph

   ! int1e_rr_cart, 9 component(s) per shell pair
   function c_int1e_rr_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rr_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rr_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rr_cart

   ! int1e_rr_sph, 9 component(s) per shell pair
   function c_int1e_rr_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rr_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rr_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rr_sph

   ! int1e_rrr_cart, 27 component(s) per shell pair
   function c_int1e_rrr_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rrr_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rrr_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rrr_cart

   ! int1e_rrr_sph, 27 component(s) per shell pair
   function c_int1e_rrr_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rrr_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rrr_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rrr_sph

   ! int1e_drinv_cart, 3 component(s) per shell pair
   function c_int1e_drinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_drinv_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_drinv_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_drinv_cart

   ! int1e_drinv_sph, 3 component(s) per shell pair
   function c_int1e_drinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_drinv_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_drinv_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_drinv_sph

   ! int1e_ipiprinv_cart, 9 component(s) per shell pair
   function c_int1e_ipiprinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipiprinv_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_ipiprinv_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipiprinv_cart

   ! int1e_ipiprinv_sph, 9 component(s) per shell pair
   function c_int1e_ipiprinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipiprinv_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_ipiprinv_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipiprinv_sph

   ! int1e_iprinvip_cart, 9 component(s) per shell pair
   function c_int1e_iprinvip_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinvip_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_iprinvip_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinvip_cart

   ! int1e_iprinvip_sph, 9 component(s) per shell pair
   function c_int1e_iprinvip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinvip_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_iprinvip_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinvip_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_r_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_r_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_r_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_r_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_r_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_r_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rr_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rr_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rr_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rr_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rrr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rrr_cart(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rrr_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rrr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      hv = int1e_rrr_sph(pout, d, fshls, patm, natm, pbas, &
                         nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rrr_sph

   ! int1e_grids_cart.  Not the shape of the others: shls carries FOUR
   ! entries -- two shells, then the half-open grid range -- and the block is
   ! di*dj*ngrids rather than a fixed component count.
   function c_int1e_grids_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                                 cache) result(ret) bind(C, name="int1e_grids_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), di, dj, ngrids
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      fshls(2) = pshls(3); fshls(3) = pshls(4)
      ngrids = fshls(3) - fshls(2)
      di = cint_cgto_cart(fshls(0), pbas); dj = cint_cgto_cart(fshls(1), pbas)
      d = [di, dj, ngrids, 1]
      call c_f_pointer(out, pout, [di*dj*ngrids])
      call c_f_pointer(env, penv, [env_len(pshls(1:2), 2, pbas, nbas)])
      hv = int1e_grids_cart(pout, d, fshls, patm, natm, pbas, &
                              nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_grids_cart

   ! int1e_grids_sph.  Not the shape of the others: shls carries FOUR
   ! entries -- two shells, then the half-open grid range -- and the block is
   ! di*dj*ngrids rather than a fixed component count.
   function c_int1e_grids_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                                 cache) result(ret) bind(C, name="int1e_grids_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), di, dj, ngrids
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = pshls(1); fshls(1) = pshls(2)
      fshls(2) = pshls(3); fshls(3) = pshls(4)
      ngrids = fshls(3) - fshls(2)
      di = cint_cgto_spheric(fshls(0), pbas); dj = cint_cgto_spheric(fshls(1), pbas)
      d = [di, dj, ngrids, 1]
      call c_f_pointer(out, pout, [di*dj*ngrids])
      call c_f_pointer(env, penv, [env_len(pshls(1:2), 2, pbas, nbas)])
      hv = int1e_grids_sph(pout, d, fshls, patm, natm, pbas, &
                              nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_grids_sph


   ! int1e_ovlp_cart, 2 shells, 1 component(s)
   function c_int1e_ovlp_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ovlp_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ovlp_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ovlp_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ovlp_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ovlp_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ovlp_cart

   ! int1e_ovlp_sph, 2 shells, 1 component(s)
   function c_int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ovlp_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ovlp_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ovlp_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ovlp_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ovlp_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ovlp_sph

   ! int1e_kin_cart, 2 shells, 1 component(s)
   function c_int1e_kin_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_kin_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_kin_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_kin_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_kin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_kin_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_kin_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_kin_cart

   ! int1e_kin_sph, 2 shells, 1 component(s)
   function c_int1e_kin_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_kin_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_kin_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_kin_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_kin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_kin_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_kin_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_kin_sph

   ! int1e_nuc_cart, 2 shells, 1 component(s)
   function c_int1e_nuc_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_nuc_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_nuc_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_nuc_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_nuc_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_nuc_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_nuc_cart

   ! int1e_nuc_sph, 2 shells, 1 component(s)
   function c_int1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_nuc_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_nuc_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_nuc_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_nuc_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_nuc_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_nuc_sph

   ! int2e_cart, 4 shells, 1 component(s)
   function c_int2e_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2e_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2e_cart

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2e_cart

   ! int2e_sph, 4 shells, 1 component(s)
   function c_int2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2e_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2e_sph

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2e_sph

   ! int3c2e_cart, 3 shells, 1 component(s)
   function c_int3c2e_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int3c2e_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [3])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 3
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 3, pbas, nbas)])
      call bind_opt(opt)
      hv = int3c2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int3c2e_cart

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int3c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint3c2e_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [3])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 3
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 3, pbas, nbas)])
      call bind_opt(opt)
      hv = int3c2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int3c2e_cart

   ! int3c2e_sph, 3 shells, 1 component(s)
   function c_int3c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int3c2e_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [3])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 3
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 3, pbas, nbas)])
      call bind_opt(opt)
      hv = int3c2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int3c2e_sph

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int3c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint3c2e_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [3])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 3
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 3, pbas, nbas)])
      call bind_opt(opt)
      hv = int3c2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int3c2e_sph

   ! int2c2e_cart, 2 shells, 1 component(s)
   function c_int2c2e_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2c2e_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int2c2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2c2e_cart

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2c2e_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int2c2e_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2c2e_cart

   ! int2c2e_sph, 2 shells, 1 component(s)
   function c_int2c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2c2e_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int2c2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2c2e_sph

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2c2e_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int2c2e_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2c2e_sph

   ! int2e_ip1_cart, 4 shells, 3 component(s)
   function c_int2e_ip1_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2e_ip1_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_ip1_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2e_ip1_cart

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_ip1_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2e_ip1_cart

   ! int2e_ip1_sph, 4 shells, 3 component(s)
   function c_int2e_ip1_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int2e_ip1_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_ip1_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int2e_ip1_sph

   ! the CINT2 spelling of the same thing: no dims, no cache
   function c2_int2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret) &
         bind(C, name="cint2e_ip1_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env, opt
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 4
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 4, pbas, nbas)])
      call bind_opt(opt)
      hv = int2e_ip1_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int2e_ip1_sph

   ! int1e_ipovlp_cart, 2 shells, 3 component(s)
   function c_int1e_ipovlp_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipovlp_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipovlp_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipovlp_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlp_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipovlp_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipovlp_cart

   ! int1e_ipovlp_sph, 2 shells, 3 component(s)
   function c_int1e_ipovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipovlp_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipovlp_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipovlp_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipovlp_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipovlp_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipovlp_sph

   ! int1e_ipkin_cart, 2 shells, 3 component(s)
   function c_int1e_ipkin_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipkin_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipkin_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipkin_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipkin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkin_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipkin_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipkin_cart

   ! int1e_ipkin_sph, 2 shells, 3 component(s)
   function c_int1e_ipkin_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipkin_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipkin_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipkin_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipkin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipkin_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipkin_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipkin_sph

   ! int1e_ipnuc_cart, 2 shells, 3 component(s)
   function c_int1e_ipnuc_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipnuc_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipnuc_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipnuc_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnuc_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipnuc_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipnuc_cart

   ! int1e_ipnuc_sph, 2 shells, 3 component(s)
   function c_int1e_ipnuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipnuc_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_ipnuc_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipnuc_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_ipnuc_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_ipnuc_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_ipnuc_sph

   ! int1e_iprinv_cart, 2 shells, 3 component(s)
   function c_int1e_iprinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinv_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_iprinv_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinv_cart

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_iprinv_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinv_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_cart(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_iprinv_cart(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_iprinv_cart

   ! int1e_iprinv_sph, 2 shells, 3 component(s)
   function c_int1e_iprinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinv_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(out, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(opt)
      hv = int1e_iprinv_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinv_sph

   ! the CINT2 spelling of the same thing: no dims, no cache, no opt
   function c2_int1e_iprinv_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_iprinv_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer,        pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), i, n
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls = 0
      d = 1
      do i = 1, 2
         fshls(i-1) = pshls(i)
         d(i-1) = cint_cgto_spheric(fshls(i-1), pbas)
      end do
      n = product(d)*3
      call c_f_pointer(buf, pout, [n])
      call c_f_pointer(env, penv, [env_len(pshls, 2, pbas, nbas)])
      call bind_opt(c_null_ptr)
      hv = int1e_iprinv_sph(pout, d, fshls, patm, natm, pbas, nbas, penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_iprinv_sph

   ! ---------------------------------------------------------------------
   ! The optimizer lifecycle, and the shell-dimension helpers.
   !
   ! libcint's contract: the caller declares `CINTOpt *opt = NULL`, passes its
   ! address to a *_optimizer routine, hands the resulting pointer to every
   ! integral call, and frees it with CINTdel_optimizer.  The handle is opaque,
   ! which is what lets this hand back a Fortran object instead of a C struct.
   ! ---------------------------------------------------------------------

   subroutine c_int2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint2e_sph_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_2e_sph_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int2e_sph_optimizer

   subroutine c_int2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint2e_cart_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_2e_cart_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int2e_cart_optimizer

   subroutine c_int3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint3c2e_sph_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_3c2e_sph_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int3c2e_sph_optimizer

   subroutine c_int3c2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint3c2e_cart_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_3c2e_cart_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int3c2e_cart_optimizer

   subroutine c_int2c2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint2c2e_sph_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_2c2e_sph_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int2c2e_sph_optimizer

   subroutine c_int2c2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
         bind(C, name="cint2c2e_cart_optimizer")
      type(c_ptr)           :: opt
      type(c_ptr),    value :: atm, bas, env
      integer(c_int), value :: natm, nbas
      integer, pointer :: patm(:), pbas(:)
      real(c_double), pointer :: penv(:)
      call c_f_pointer(atm, patm, [ATM_SLOTS*natm])
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      call c_f_pointer(env, penv, [env_len_all(patm, natm, pbas, nbas)])
      call libcint_2c2e_cart_optimizer(opt, patm, natm, pbas, nbas, penv)
   end subroutine c_int2c2e_cart_optimizer

   subroutine c_del_optimizer(opt) bind(C, name="CINTdel_optimizer")
      type(c_ptr) :: opt
      call libcint_del_optimizer(opt)
   end subroutine c_del_optimizer

   ! Shell dimensions.  A driver needs these to size its own buffers, and
   ! without them a C caller has to reimplement the l -> nf mapping.
   integer(c_int) function c_cgto_spheric(bas_id, bas) result(dim) &
         bind(C, name="CINTcgto_spheric")
      integer(c_int), value :: bas_id
      type(c_ptr),    value :: bas
      integer, pointer :: pbas(:)
      call c_f_pointer(bas, pbas, [BAS_SLOTS*(bas_id + 1)])
      dim = int(libcint_cgto_sph(bas_id, pbas), c_int)
   end function c_cgto_spheric

   integer(c_int) function c_cgto_cart(bas_id, bas) result(dim) &
         bind(C, name="CINTcgto_cart")
      integer(c_int), value :: bas_id
      type(c_ptr),    value :: bas
      integer, pointer :: pbas(:)
      call c_f_pointer(bas, pbas, [BAS_SLOTS*(bas_id + 1)])
      dim = int(libcint_cgto_cart(bas_id, pbas), c_int)
   end function c_cgto_cart

   integer(c_int) function c_tot_cgto_spheric(bas, nbas) result(dim) &
         bind(C, name="CINTtot_cgto_spheric")
      type(c_ptr),    value :: bas
      integer(c_int), value :: nbas
      integer, pointer :: pbas(:)
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      dim = int(libcint_tot_cgto_sph(pbas, nbas), c_int)
   end function c_tot_cgto_spheric

   integer(c_int) function c_tot_cgto_cart(bas, nbas) result(dim) &
         bind(C, name="CINTtot_cgto_cart")
      type(c_ptr),    value :: bas
      integer(c_int), value :: nbas
      integer, pointer :: pbas(:)
      call c_f_pointer(bas, pbas, [BAS_SLOTS*nbas])
      dim = int(libcint_tot_cgto_cart(pbas, nbas), c_int)
   end function c_tot_cgto_cart

   ! An optimizer names no shells, so its env bound is over every shell and
   ! every atom -- the same reasoning as env_extent_all in libcint_fortran.
   pure integer function env_len_all(atm, natm, bas, nbas) result(n)
      integer, intent(in) :: atm(0:), bas(0:), natm, nbas
      integer :: i, np, nc
      n = PTR_ENV_START
      do i = 0, natm - 1
         n = max(n, atm(ATM_SLOTS*i + 1) + 3)
      end do
      do i = 0, nbas - 1
         np = bas(BAS_SLOTS*i + NPRIM_OF)
         nc = bas(BAS_SLOTS*i + NCTR_OF)
         n = max(n, bas(BAS_SLOTS*i + PTR_EXP)   + np)
         n = max(n, bas(BAS_SLOTS*i + PTR_COEFF) + np*nc)
      end do
   end function env_len_all

end module cint_c_abi
