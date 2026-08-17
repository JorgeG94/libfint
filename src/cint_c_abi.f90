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
! The entry points a real consumer has been observed to bind, not all 620.
! Adding one is four lines and the list below says which; the generator could
! emit the whole set, and should if this list keeps growing.
!
! Two shapes, both of which libcint exports for every integral:
!
!   the CINT3 form   name(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
!   the CINT2 form   cname(buf, shls, atm, natm, bas, nbas, env)
!
! dims, opt and cache may all be NULL, and the CINT2 form has no dims at all,
! so both compute the natural block shape from the shells -- exactly as
! libcint's own c2s layer does.
!
! The optimizer argument is accepted and ignored.  libcint's CINTOpt is a C
! struct this library has no reading of; passing NULL to libcint means "no
! optimizer", and ignoring it here means the same thing, at the cost of the
! shell-pair caching a caller may have set up.  Correct, slower, and the
! alternative is decoding a foreign struct layout.
!
module cint_c_abi
   use iso_c_binding, only: c_int, c_double, c_ptr, c_f_pointer, c_associated
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
   implicit none
   private

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

   ! How long env has to be for these shells: the largest offset any of them
   ! reaches, which is what lets a c_f_pointer give env a real extent.
   pure integer function env_len(shls, nsh, bas, nbas) result(n)
      integer(c_int), intent(in) :: shls(:), bas(0:)
      integer,        intent(in) :: nsh, nbas
      integer :: i, sh, np, nc
      n = PTR_ENV_START
      do i = 1, nsh
         sh = int(shls(i))
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
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_r_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_r_cart

   ! int1e_r_sph, 3 component(s) per shell pair
   function c_int1e_r_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_r_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_r_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_r_sph

   ! int1e_rr_cart, 9 component(s) per shell pair
   function c_int1e_rr_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rr_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rr_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rr_cart

   ! int1e_rr_sph, 9 component(s) per shell pair
   function c_int1e_rr_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rr_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rr_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rr_sph

   ! int1e_rrr_cart, 27 component(s) per shell pair
   function c_int1e_rrr_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rrr_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rrr_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rrr_cart

   ! int1e_rrr_sph, 27 component(s) per shell pair
   function c_int1e_rrr_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_rrr_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rrr_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_rrr_sph

   ! int1e_drinv_cart, 3 component(s) per shell pair
   function c_int1e_drinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_drinv_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_drinv_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_drinv_cart

   ! int1e_drinv_sph, 3 component(s) per shell pair
   function c_int1e_drinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_drinv_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_drinv_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_drinv_sph

   ! int1e_ipiprinv_cart, 9 component(s) per shell pair
   function c_int1e_ipiprinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipiprinv_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_ipiprinv_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipiprinv_cart

   ! int1e_ipiprinv_sph, 9 component(s) per shell pair
   function c_int1e_ipiprinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_ipiprinv_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_ipiprinv_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_ipiprinv_sph

   ! int1e_iprinvip_cart, 9 component(s) per shell pair
   function c_int1e_iprinvip_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinvip_cart")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_iprinvip_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinvip_cart

   ! int1e_iprinvip_sph, 9 component(s) per shell pair
   function c_int1e_iprinvip_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, &
                            cache) result(ret) bind(C, name="int1e_iprinvip_sph")
      type(c_ptr),    value :: out, dims, shls, atm, bas, env, opt, cache
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(out, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_iprinvip_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_iprinvip_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_r_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_r_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_r_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_r_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_r_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*3])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_r_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_r_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rr_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rr_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rr_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*9])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rr_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rr_sph

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rrr_cart(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_cart")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rrr_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c2_int1e_rrr_cart

   ! the CINT2 spelling of the same thing: no dims, no opt, no cache
   function c2_int1e_rrr_sph(buf, shls, atm, natm, bas, nbas, env) result(ret) &
         bind(C, name="cint1e_rrr_sph")
      type(c_ptr),    value :: buf, shls, atm, bas, env
      integer(c_int), value :: natm, nbas
      integer(c_int)        :: ret
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:1), di, dj
      logical :: hv
      call c_f_pointer(shls, pshls, [2])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, 1, 1]
      call c_f_pointer(buf, pout, [di*dj*27])
      call c_f_pointer(env, penv, [env_len(pshls, 2, int(pbas), int(nbas))])
      hv = int1e_rrr_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                         int(nbas), penv, ws)
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
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), di, dj, ngrids
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      fshls(2) = int(pshls(3)); fshls(3) = int(pshls(4))
      ngrids = fshls(3) - fshls(2)
      di = cint_cgto_cart(fshls(0), int(pbas)); dj = cint_cgto_cart(fshls(1), int(pbas))
      d = [di, dj, ngrids, 1]
      call c_f_pointer(out, pout, [di*dj*ngrids])
      call c_f_pointer(env, penv, [env_len(pshls(1:2), 2, int(pbas), int(nbas))])
      hv = int1e_grids_cart(pout, d, fshls, int(patm), int(natm), int(pbas), &
                              int(nbas), penv, ws)
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
      integer(c_int), pointer :: pshls(:), patm(:), pbas(:)
      real(c_double), pointer :: pout(:), penv(:)
      integer :: d(0:3), fshls(0:3), di, dj, ngrids
      logical :: hv
      call c_f_pointer(shls, pshls, [4])
      call c_f_pointer(atm,  patm,  [ATM_SLOTS*natm])
      call c_f_pointer(bas,  pbas,  [BAS_SLOTS*nbas])
      fshls(0) = int(pshls(1)); fshls(1) = int(pshls(2))
      fshls(2) = int(pshls(3)); fshls(3) = int(pshls(4))
      ngrids = fshls(3) - fshls(2)
      di = cint_cgto_spheric(fshls(0), int(pbas)); dj = cint_cgto_spheric(fshls(1), int(pbas))
      d = [di, dj, ngrids, 1]
      call c_f_pointer(out, pout, [di*dj*ngrids])
      call c_f_pointer(env, penv, [env_len(pshls(1:2), 2, int(pbas), int(nbas))])
      hv = int1e_grids_sph(pout, d, fshls, int(patm), int(natm), int(pbas), &
                              int(nbas), penv, ws)
      ret = merge(1_c_int, 0_c_int, hv)
   end function c_int1e_grids_sph

end module cint_c_abi
