!
! libcint_fortran, backed by the Fortran port instead of the C library.
!
! Same module name and the same interfaces as include/libcint_fortran.f90, so
! a caller switches backends by choosing which of the two files to compile --
! WITH_FORTRAN_BACKEND in the top-level CMakeLists -- and changes no source.
! That is the point of the exercise: metalquicha says `use libcint_fortran`
! and never finds out which one it got.
!
! The constants come from libcint_interface, which is nothing but interface
! blocks and parameters: compiling and linking against it needs no C object.
! Taking the numbers from there rather than retyping them is what keeps the
! two backends' slot conventions identical by construction -- the 1-based slot
! names here are NOT the 0-based ones the port uses internally, and confusing
! the two has already cost this port one debugging session.
!
! NOTHING IS COPIED PER SHELL QUARTET.  The public interfaces take atm and bas
! as (SLOTS, *) tables, and the port wants flat 0-based arrays.  Fortran's
! column-major order means those are the same bytes, so the private helpers
! below declare rank-1 explicit-shape dummies and let sequence association do
! the reinterpretation.  Building a flat copy instead would reintroduce
! exactly the per-call cost the port was just measured free of.
!
! The optimizer is real here.  libcint's API hands back an opaque pointer and
! lets the caller free it when they choose, so the optimizers live in a fixed
! pool and the handle is a c_loc into it.  Underneath, the port carries the
! optimizer on the workspace rather than as an argument to every entry point,
! and this module is where the two conventions meet: ws_bind_opt points the
! workspace at the caller's optimizer for the duration of one call, and a
! null or absent handle means no optimizer -- which is exactly what a NULL
! CINTOpt means to libcint.
!
! An optimizer built for one integral must not be used for another.  The C
! does not check; the port does, comparing the ng and the arity it was built
! with, and falls back to computing rather than reading the wrong index map.
!
module libcint_fortran
    use iso_c_binding, only: c_int, c_double, c_double_complex, c_ptr, c_null_ptr, c_loc, c_f_pointer, c_associated
    use iso_fortran_env, only: real64, int32
    use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, CHARGE_OF, PTR_COORD, &
                                 NUC_MOD_OF, PTR_ZETA, ATOM_OF, ANG_OF, &
                                 NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, &
                                 PTR_EXPCUTOFF, PTR_COMMON_ORIG, PTR_RINV_ORIG, &
                                 PTR_RINV_ZETA, PTR_RANGE_OMEGA, PTR_ENV_START
    use cint_const,     only: port_dp => dp
    use cint_workspace, only: cint_ws
    use cint_envs,      only: cint_opt_t
    use cint_opt,       only: cint_del_optimizer
    use cint_g1e,       only: cint_all_1e_optimizer
    use cint_g2e,       only: cint_all_2e_optimizer, cint_all_3c2e_optimizer, &
                              cint_all_2c2e_optimizer
    ! The derivative and spinor families the port now has.  They were stubbed
    ! here while they were D9/D10 work; the stubs are gone.
    use cint_gen_grad1,   only: int1e_ipnuc_cart, int1e_ipnuc_sph, &
                                int1e_iprinv_cart, int1e_iprinv_sph
    use cint_gen_grad2,   only: int2e_ip1_cart, int2e_ip1_sph
    ! The scalar ECP integrals. libcint has no ECP code at all, so the C-backed
    ! twin of this module cannot offer these and a caller that wants them has
    ! to be built against this backend -- which is what `ECP_AVAILABLE` on the
    ! metalquicha side is for. Reached here rather than through a hand-written
    ! `bind(C)` interface so that the caller needs no C marshalling of its own.
    use cint_ecp_drv,     only: ecp_scalar_cart, ecp_scalar_sph
    use cint_ecp_num,     only: AS_ECPBAS_OFFSET, AS_NECPBAS
    use cint_gen_int3c2e, only: int3c2e_ip1_cart, int3c2e_ip1_sph, &
                                int3c2e_ip2_cart, int3c2e_ip2_sph, &
                                int2c2e_ip1_cart, int2c2e_ip1_sph
    use cint_gen_intor3,  only: int1e_spnucsp_spinor
    use cint_gen_intor4,  only: int2e_spsp1_spinor
    use cint_bas,       only: cint_cgto_cart, cint_cgto_spheric, cint_cgto_spinor, &
                              cint_tot_cgto_cart, cint_tot_cgto_spheric, &
                              cint_tot_pgto_spheric, cint_gto_norm
    use cint_1e,        only: int1e_nuc_cart, int1e_nuc_sph, &
                              int1e_ovlp_cart, int1e_ovlp_sph
    use cint_gen_intor1, only: int1e_kin_cart, int1e_kin_sph
    ! the ip-derivatives live with the other gradients, in grad1.c's module
    use cint_gen_grad1,  only: int1e_ipovlp_cart, int1e_ipovlp_sph, &
                               int1e_ipkin_cart,  int1e_ipkin_sph
    use cint_2e,        only: int2e_cart, int2e_sph
    use cint_3c2e,      only: int3c2e_cart, int3c2e_sph, int2c2e_cart, int2c2e_sph
    implicit none
    private

    public :: dp, ip, zp

    public :: LIBCINT_ATM_SLOTS, LIBCINT_BAS_SLOTS
    public :: LIBCINT_CHARGE_OF, LIBCINT_PTR_COORD, LIBCINT_NUC_MOD_OF, LIBCINT_PTR_ZETA
    public :: LIBCINT_ATOM_OF, LIBCINT_ANG_OF, LIBCINT_NPRIM_OF, LIBCINT_NCTR_OF
    public :: LIBCINT_KAPPA_OF, LIBCINT_PTR_EXP, LIBCINT_PTR_COEFF
    public :: LIBCINT_PTR_EXPCUTOFF, LIBCINT_PTR_COMMON_ORIG, LIBCINT_PTR_RINV_ORIG
    public :: LIBCINT_PTR_RINV_ZETA, LIBCINT_PTR_RANGE_OMEGA, LIBCINT_PTR_ENV_START

    public :: libcint_cgto_cart, libcint_cgto_sph, libcint_cgto_spinor
    public :: libcint_tot_cgto_sph, libcint_tot_cgto_cart, libcint_tot_pgto_sph
    public :: libcint_gto_norm

    public :: libcint_1e_ovlp_cart, libcint_1e_kin_cart, libcint_1e_nuc_cart
    public :: libcint_1e_ipovlp_cart
    public :: libcint_ecp_cart, libcint_ecp_sph
    public :: libcint_1e_ovlp_sph, libcint_1e_kin_sph, libcint_1e_nuc_sph
    public :: libcint_1e_ipovlp_sph
    public :: libcint_1e_ipkin_cart, libcint_1e_ipnuc_cart, libcint_1e_iprinv_cart
    public :: libcint_1e_ipkin_sph, libcint_1e_ipnuc_sph, libcint_1e_iprinv_sph
    public :: libcint_1e_spnucsp

    public :: libcint_2e_cart, libcint_2e_sph
    public :: libcint_3c2e_sph, libcint_2c2e_sph
    public :: libcint_3c2e_cart, libcint_2c2e_cart
    public :: libcint_2e_ip1_cart, libcint_2e_ip1_sph
    public :: libcint_3c2e_ip1_cart, libcint_3c2e_ip1_sph
    public :: libcint_3c2e_ip2_cart, libcint_3c2e_ip2_sph
    public :: libcint_2c2e_ip1_cart, libcint_2c2e_ip1_sph
    public :: libcint_2e_spsp1

    public :: libcint_2e_cart_optimizer, libcint_2e_sph_optimizer
    public :: libcint_3c2e_cart_optimizer, libcint_3c2e_sph_optimizer
    public :: libcint_2c2e_cart_optimizer, libcint_2c2e_sph_optimizer
    public :: libcint_2e_ip1_cart_optimizer, libcint_2e_ip1_sph_optimizer
    public :: libcint_del_optimizer
    public :: libcint_3c2e_ip1_cart_optimizer, libcint_3c2e_ip1_sph_optimizer
    public :: libcint_3c2e_ip2_cart_optimizer, libcint_3c2e_ip2_sph_optimizer
    public :: libcint_2c2e_ip1_cart_optimizer, libcint_2c2e_ip1_sph_optimizer
    public :: libcint_2e_spsp1_optimizer

    integer, parameter :: dp = c_double
    integer, parameter :: ip = c_int
    integer, parameter :: zp = c_double_complex

    integer, parameter :: libcint_real64_is_c_double = 1 / merge(1, 0, real64 == c_double)
    integer, parameter :: libcint_int32_is_c_int     = 1 / merge(1, 0, int32  == c_int)
    ! and the port's own kinds had better be the same ones, since the helpers
    ! below hand it these arrays without conversion
    integer, parameter :: libcint_port_dp_is_c_double = 1 / merge(1, 0, port_dp == c_double)
    integer, parameter :: libcint_port_int_is_c_int   = 1 / merge(1, 0, kind(1) == c_int)

    integer(c_int), parameter :: LIBCINT_ATM_SLOTS = ATM_SLOTS
    integer(c_int), parameter :: LIBCINT_BAS_SLOTS = BAS_SLOTS
    integer(c_int), parameter :: LIBCINT_CHARGE_OF = CHARGE_OF
    integer(c_int), parameter :: LIBCINT_PTR_COORD = PTR_COORD
    integer(c_int), parameter :: LIBCINT_NUC_MOD_OF = NUC_MOD_OF
    integer(c_int), parameter :: LIBCINT_PTR_ZETA = PTR_ZETA
    integer(c_int), parameter :: LIBCINT_ATOM_OF = ATOM_OF
    integer(c_int), parameter :: LIBCINT_ANG_OF = ANG_OF
    integer(c_int), parameter :: LIBCINT_NPRIM_OF = NPRIM_OF
    integer(c_int), parameter :: LIBCINT_NCTR_OF = NCTR_OF
    integer(c_int), parameter :: LIBCINT_KAPPA_OF = KAPPA_OF
    integer(c_int), parameter :: LIBCINT_PTR_EXP = PTR_EXP
    integer(c_int), parameter :: LIBCINT_PTR_COEFF = PTR_COEFF
    integer(c_int), parameter :: LIBCINT_PTR_EXPCUTOFF = PTR_EXPCUTOFF
    integer(c_int), parameter :: LIBCINT_PTR_COMMON_ORIG = PTR_COMMON_ORIG
    integer(c_int), parameter :: LIBCINT_PTR_RINV_ORIG = PTR_RINV_ORIG
    integer(c_int), parameter :: LIBCINT_PTR_RINV_ZETA = PTR_RINV_ZETA
    integer(c_int), parameter :: LIBCINT_PTR_RANGE_OMEGA = PTR_RANGE_OMEGA
    integer(c_int), parameter :: LIBCINT_PTR_ENV_START = PTR_ENV_START

    ! Which integral a dispatcher should evaluate.
    integer, parameter :: K_OVLP = 1, K_KIN = 2, K_NUC = 3, K_IPOVLP = 4, K_IPKIN = 5

    ! Scratch, reused across calls.  One per thread: metalquicha runs OpenMP
    ! over shell pairs, and a shared workspace would have two threads bump-
    ! allocating from one buffer.  With OpenMP off the directive is a comment
    ! and this is an ordinary saved variable, which is right for one thread.
    type(cint_ws), save :: ws
    !$omp threadprivate(ws)

    ! The optimizers the caller has built.  libcint's API hands back an opaque
    ! pointer and lets the caller decide when to release it, so these have to
    ! outlive the call that makes them; a fixed pool with c_loc handles is the
    ! smallest thing that does.
    !
    ! Shared, not threadprivate, unlike the workspace.  An optimizer is
    ! written once and only read afterwards, and the usual shape of a caller
    ! is to build it before a parallel region and use it inside -- which a
    ! threadprivate copy would break, because the thread that built it is not
    ! the one that reads it.
    integer, parameter :: LIBCINT_MAX_OPT = 64
    type(cint_opt_t), target, save :: opt_pool(LIBCINT_MAX_OPT)
    logical, save :: opt_taken(LIBCINT_MAX_OPT) = .false.

contains

    ! ========================================================================
    ! Bookkeeping.  Each public entry hands its (SLOTS, *) table to a helper
    ! whose dummy is the flat rank-1 equivalent, sized from the count it also
    ! receives; the association is the standard's sequence rule, not a copy.
    ! ========================================================================

    function libcint_cgto_cart(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = cgto_cart_(bas_id, bas)
    end function libcint_cgto_cart

    function cgto_cart_(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*(bas_id+1_ip) - 1)
        integer(ip) :: dim
        dim = cint_cgto_cart(bas_id, bas)
    end function cgto_cart_

    function libcint_cgto_sph(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = cgto_sph_(bas_id, bas)
    end function libcint_cgto_sph

    function cgto_sph_(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*(bas_id+1_ip) - 1)
        integer(ip) :: dim
        dim = cint_cgto_spheric(bas_id, bas)
    end function cgto_sph_

    function libcint_cgto_spinor(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip) :: dim
        dim = cgto_spinor_(bas_id, bas)
    end function libcint_cgto_spinor

    function cgto_spinor_(bas_id, bas) result(dim)
        integer(ip), intent(in) :: bas_id
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*(bas_id+1_ip) - 1)
        integer(ip) :: dim
        dim = cint_cgto_spinor(bas_id, bas)
    end function cgto_spinor_

    function libcint_tot_cgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = tot_cgto_sph_(bas, nbas)
    end function libcint_tot_cgto_sph

    function tot_cgto_sph_(bas, nbas) result(ntot)
        integer(ip), intent(in) :: nbas
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip) :: ntot
        ntot = cint_tot_cgto_spheric(bas, nbas)
    end function tot_cgto_sph_

    function libcint_tot_cgto_cart(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = tot_cgto_cart_(bas, nbas)
    end function libcint_tot_cgto_cart

    function tot_cgto_cart_(bas, nbas) result(ntot)
        integer(ip), intent(in) :: nbas
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip) :: ntot
        ntot = cint_tot_cgto_cart(bas, nbas)
    end function tot_cgto_cart_

    function libcint_tot_pgto_sph(bas, nbas) result(ntot)
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        integer(ip) :: ntot
        ntot = tot_pgto_sph_(bas, nbas)
    end function libcint_tot_pgto_sph

    function tot_pgto_sph_(bas, nbas) result(ntot)
        integer(ip), intent(in) :: nbas
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip) :: ntot
        ntot = cint_tot_pgto_spheric(bas, nbas)
    end function tot_pgto_sph_

    function libcint_gto_norm(n, alpha) result(norm)
        integer(ip), intent(in) :: n
        real(dp), intent(in) :: alpha
        real(dp) :: norm
        norm = cint_gto_norm(n, alpha)
    end function libcint_gto_norm

    ! ========================================================================
    ! How much of env the port can reach.
    !
    ! The port holds env by reference, so it needs an extent -- and an
    ! assumed-size dummy has none.  Rather than guess, bound it: the only
    ! offsets read are the fixed header slots, the exponent and coefficient
    ! blocks of the shells named in shls, and the coordinate and
    ! nuclear-model slots of every atom, which int1e_nuc walks in full.
    ! nbas does not enter, so this stays cheap inside the quartet loop.
    ! ========================================================================
    pure function env_extent(shls, nsh, atm, natm, bas, nbas) result(n)
        integer(ip), intent(in) :: nsh, natm, nbas
        integer(ip), intent(in) :: shls(nsh)
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip) :: n
        integer(ip) :: i, sh, np, nc
        n = LIBCINT_PTR_ENV_START
        do i = 1, nsh
            sh = shls(i)                          ! 0-based shell id
            np = bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_NPRIM_OF - 1)
            nc = bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_NCTR_OF  - 1)
            n = max(n, bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_PTR_EXP   - 1) + np)
            n = max(n, bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_PTR_COEFF - 1) + np*nc)
        end do
        do i = 0, natm - 1
            n = max(n, atm(LIBCINT_ATM_SLOTS*i + LIBCINT_PTR_COORD - 1) + 3_ip)
            ! PTR_ZETA and, one past it, PTR_FRAC_CHARGE, for which
            ! libcint_interface has no name
            n = max(n, atm(LIBCINT_ATM_SLOTS*i + LIBCINT_PTR_ZETA - 1) + 2_ip)
        end do
    end function env_extent

    ! `env_extent` for a call that also reads ECP shells.
    !
    ! An ECP row lives in the same `bas` array as the orbital shells but past
    ! them, at `ecpoff`, and its exponents and coefficients live in `env` like
    ! any other -- so the extent has to cover them or the callee reads past the
    ! end of what the caller passed. Slot 3 of an ECP row is RADI_POWER, the r
    ! exponent, sitting where NCTR_OF sits in an orbital row and meaning
    ! something else entirely: reading it as a contraction count gives an
    ! extent that is too short for an r^0 term and past the buffer for an r^2
    ! one. An ECP shell has one coefficient per primitive, so `np` alone.
    ! Not `pure`, and `bas` arrives already cut to length rather than assumed
    ! size. Both are for lfortran: it reports `LIBCINT_BAS_SLOTS` as a call to
    ! an impure getter when a PURE procedure uses it to slice an assumed-size
    ! dummy in an executable statement, which is how this was first written.
    ! `env_extent` beside it stays pure because it only uses the parameter in a
    ! specification expression, which lfortran accepts. Nothing here is hot --
    ! one call per shell pair against a kernel that integrates a radial grid --
    ! so the attribute bought nothing worth working around it for.
    function ecp_env_extent(shls, atm, natm, bas, nbas, ecpoff, necpbas) result(n)
        integer(ip), intent(in) :: natm, nbas, ecpoff, necpbas
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:)
        integer(ip) :: n
        integer(ip) :: i, sh, np

        n = env_extent(shls, 2_ip, atm, natm, bas(0:LIBCINT_BAS_SLOTS*nbas - 1), nbas)
        do i = 0, necpbas - 1
            sh = ecpoff + i
            np = bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_NPRIM_OF - 1)
            n = max(n, bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_PTR_EXP   - 1) + np)
            n = max(n, bas(LIBCINT_BAS_SLOTS*sh + LIBCINT_PTR_COEFF - 1) + np)
        end do
    end function ecp_env_extent

    ! ========================================================================
    ! One- and two-electron dispatch.  These take the flat views; the public
    ! entry points below are one line each.
    ! ========================================================================

    function run1e(which, cart, buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        type(c_ptr), intent(in), optional :: opt
        integer,     intent(in) :: which
        logical,     intent(in) :: cart
        integer(ip), intent(in) :: natm, nbas
        real(dp),    intent(out), target :: buf(*)
        integer(ip), intent(in),  target :: shls(2)
        integer(ip), intent(in),  target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in),  target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp),    intent(in),  target :: env(*)
        integer(ip) :: ret

        integer :: dims(0:3), fshls(0:1), di, dj, ncomp, nenv
        logical :: hv

        call ws_bind_opt(opt)
        fshls(0) = shls(1); fshls(1) = shls(2)
        if (cart) then
            di = cint_cgto_cart(fshls(0), bas);    dj = cint_cgto_cart(fshls(1), bas)
        else
            di = cint_cgto_spheric(fshls(0), bas); dj = cint_cgto_spheric(fshls(1), bas)
        end if
        ncomp = 1
        if (which == K_IPOVLP .or. which == K_IPKIN) ncomp = 3
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)

        select case (which)
        case (K_OVLP)
            if (cart) then
                hv = int1e_ovlp_cart(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                     bas, nbas, env(1:nenv), ws)
            else
                hv = int1e_ovlp_sph(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                    bas, nbas, env(1:nenv), ws)
            end if
        case (K_KIN)
            if (cart) then
                hv = int1e_kin_cart(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                    bas, nbas, env(1:nenv), ws)
            else
                hv = int1e_kin_sph(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                   bas, nbas, env(1:nenv), ws)
            end if
        case (K_NUC)
            if (cart) then
                hv = int1e_nuc_cart(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                    bas, nbas, env(1:nenv), ws)
            else
                hv = int1e_nuc_sph(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                   bas, nbas, env(1:nenv), ws)
            end if
        case (K_IPOVLP)
            if (cart) then
                hv = int1e_ipovlp_cart(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                       bas, nbas, env(1:nenv), ws)
            else
                hv = int1e_ipovlp_sph(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                      bas, nbas, env(1:nenv), ws)
            end if
        case (K_IPKIN)
            if (cart) then
                hv = int1e_ipkin_cart(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                      bas, nbas, env(1:nenv), ws)
            else
                hv = int1e_ipkin_sph(buf(1:di*dj*ncomp), dims, fshls, atm, natm, &
                                     bas, nbas, env(1:nenv), ws)
            end if
        case default
            error stop "libcint_fortran: unknown 1e integral"
        end select
        ret = merge(1_ip, 0_ip, hv)
    end function run1e

    ! ncen is 2, 3 or 4: the shorter forms name fewer shells and take a
    ! different envs setup, but one output convention.
    function run2e(ncen, cart, buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        integer(ip), intent(in) :: ncen, natm, nbas
        logical,     intent(in) :: cart
        real(dp),    intent(out), target :: buf(*)
        integer(ip), intent(in),  target :: shls(ncen)
        integer(ip), intent(in),  target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in),  target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp),    intent(in),  target :: env(*)
        integer(ip) :: ret

        type(c_ptr), intent(in), optional :: opt
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv

        call ws_bind_opt(opt)
        fshls = 0
        dims = 1
        do i = 1, ncen
            fshls(i-1) = shls(i)
            if (cart) then
                dims(i-1) = cint_cgto_cart(fshls(i-1), bas)
            else
                dims(i-1) = cint_cgto_spheric(fshls(i-1), bas)
            end if
        end do
        n = product(dims)
        nenv = env_extent(shls, ncen, atm, natm, bas, nbas)

        select case (ncen)
        case (4)
            if (cart) then
                hv = int2e_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            else
                hv = int2e_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            end if
        case (3)
            if (cart) then
                hv = int3c2e_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            else
                hv = int3c2e_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            end if
        case (2)
            if (cart) then
                hv = int2c2e_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            else
                hv = int2c2e_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
            end if
        case default
            error stop "libcint_fortran: unknown 2e centre count"
        end select
        ret = merge(1_ip, 0_ip, hv)
    end function run2e

    ! ========================================================================
    ! One-electron entry points
    ! ========================================================================

    ! ========================================================================
    ! Scalar ECP over a shell pair.
    !
    ! Signature deliberately the same shape as the 1e integrals above, so a
    ! caller that already has `libcint_1e_nuc_sph` working needs nothing new to
    ! call this. The ECP shells are not passed: they travel in `bas` past the
    ! orbital ones, and `env` slots 18 and 19 say where they start and how many
    ! there are, which is the convention PySCF writes and libcint reads.

    function run_ecp(cart, buf, shls, atm, natm, bas, nbas, env) result(ret)
        logical,     intent(in)  :: cart
        integer(ip), intent(in)  :: natm, nbas
        real(dp),    intent(out) :: buf(*)
        integer(ip), intent(in)  :: shls(2)
        integer(ip), intent(in)  :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        ! Plain assumed-size, one-based, which is what every other wrapper in
        ! this file takes. The ECP rows put its real length past `nbas`, so it
        ! cannot be sized by `nbas` the way `run1e`'s is.
        integer(ip), intent(in)  :: bas(*)
        real(dp),    intent(in)  :: env(*)
        integer(ip) :: ret

        integer(ip) :: fshls(0:1), di, dj, ecpoff, necpbas, nb_all, nenv

        ! Slots 18 and 19 sit below PTR_ENV_START, so they are readable before
        ! anything else about `env` is known.
        ecpoff  = int(env(AS_ECPBAS_OFFSET + 1), ip)
        necpbas = int(env(AS_NECPBAS + 1), ip)
        if (necpbas <= 0) then
            ! No ECP shells in this molecule. The block is zero rather than
            ! absent, so a caller can add it unconditionally.
            ret = 0_ip
            return
        end if
        ! Not `nbas + necpbas`: nothing requires the ECP rows to begin exactly
        ! where the orbital ones end. PySCF stacks them that way and the C never
        ! assumes it, so neither does this.
        nb_all = max(nbas, ecpoff + necpbas)

        fshls(0) = shls(1); fshls(1) = shls(2)
        ! Sliced rather than passed whole: `bas` is assumed-size here, because
        ! the ECP rows put its real length past `nbas`, and an assumed-size
        ! array cannot be handed to the assumed-shape dummies below.
        if (cart) then
            di = cint_cgto_cart(fshls(0), bas(1:LIBCINT_BAS_SLOTS*nb_all))
            dj = cint_cgto_cart(fshls(1), bas(1:LIBCINT_BAS_SLOTS*nb_all))
        else
            di = cint_cgto_spheric(fshls(0), bas(1:LIBCINT_BAS_SLOTS*nb_all))
            dj = cint_cgto_spheric(fshls(1), bas(1:LIBCINT_BAS_SLOTS*nb_all))
        end if
        nenv = ecp_env_extent(shls, atm, natm, bas(1:LIBCINT_BAS_SLOTS*nb_all), &
                              nbas, ecpoff, necpbas)

        if (cart) then
            ret = int(ecp_scalar_cart(buf(1:di*dj), fshls, atm, natm, &
                                      bas(1:LIBCINT_BAS_SLOTS*nb_all), nbas, &
                                      env(1:nenv), [di, dj]), ip)
        else
            ret = int(ecp_scalar_sph(buf(1:di*dj), fshls, atm, natm, &
                                     bas(1:LIBCINT_BAS_SLOTS*nb_all), nbas, &
                                     env(1:nenv), [di, dj]), ip)
        end if
    end function run_ecp

    function libcint_ecp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run_ecp(.true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_ecp_cart

    function libcint_ecp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run_ecp(.false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_ecp_sph

    function libcint_1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_OVLP, .true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ovlp_cart

    function libcint_1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_OVLP, .false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ovlp_sph

    function libcint_1e_kin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_KIN, .true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_kin_cart

    function libcint_1e_kin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_KIN, .false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_kin_sph

    function libcint_1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_NUC, .true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_nuc_cart

    function libcint_1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_NUC, .false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_nuc_sph

    function libcint_1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_IPOVLP, .true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipovlp_cart

    function libcint_1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_IPOVLP, .false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipovlp_sph

    function libcint_1e_ipkin_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_IPKIN, .true., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipkin_cart

    function libcint_1e_ipkin_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        integer(ip) :: ret
        ret = run1e(K_IPKIN, .false., buf, shls, atm, natm, bas, nbas, env)
    end function libcint_1e_ipkin_sph

    ! ========================================================================
    ! Two-electron entry points
    ! ========================================================================

    function libcint_2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(4_ip, .true., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_2e_cart

    function libcint_2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(4_ip, .false., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_2e_sph

    !> Three-centre ERI, the (mu nu | P) of density fitting.  shls has four
    !> entries though only three centres contribute; the fourth is ignored.
    function libcint_3c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(3_ip, .true., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_3c2e_cart

    function libcint_3c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(4)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(3_ip, .false., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_3c2e_sph

    function libcint_2c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(2_ip, .true., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_2c2e_cart

    function libcint_2c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out) :: buf(*)
        integer(ip), intent(in) :: shls(2)
        integer(ip), intent(in) :: atm(LIBCINT_ATM_SLOTS, *)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in) :: bas(LIBCINT_BAS_SLOTS, *)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in) :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        ret = run2e(2_ip, .false., buf, shls, atm, natm, bas, nbas, env, opt)
    end function libcint_2c2e_sph

    ! ========================================================================
    ! Optimizers.  The port has none; a null pointer is what every libcint
    ! entry point already reads as "no optimizer, use the general loop", so
    ! this is a supported state rather than a sentinel being abused.
    ! ========================================================================

    subroutine libcint_2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2e_cart_optimizer

    subroutine libcint_2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2e_sph_optimizer

    subroutine libcint_3c2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_cart_optimizer

    subroutine libcint_3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_sph_optimizer

    subroutine libcint_2c2e_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2c2e_cart_optimizer

    subroutine libcint_2c2e_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2c2e_sph_optimizer

    subroutine libcint_2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2e_ip1_cart_optimizer

    subroutine libcint_2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2e_ip1_sph_optimizer

    subroutine libcint_3c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_ip1_cart_optimizer

    subroutine libcint_3c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_ip1_sph_optimizer

    subroutine libcint_3c2e_ip2_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 1, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_ip2_cart_optimizer

    subroutine libcint_3c2e_ip2_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [0, 0, 1, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_3c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_3c2e_ip2_sph_optimizer

    subroutine libcint_2c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2c2e_ip1_cart_optimizer

    subroutine libcint_2c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2c2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2c2e_ip1_sph_optimizer

    subroutine libcint_2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env)
        type(c_ptr), intent(out) :: opt
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        real(dp), intent(in) :: env(*)
        integer, parameter :: ng(0:7) = [1, 1, 0, 0, 2, 4, 1, 1]
        integer :: slot, nenv
        slot = opt_claim()
        nenv = env_extent_all(atm, natm, bas, nbas)
        call cint_all_2e_optimizer(opt_pool(slot), ng, atm, natm, bas, nbas, env(1:nenv))
        opt = c_loc(opt_pool(slot))
    end subroutine libcint_2e_spsp1_optimizer

    ! Find the pool slot by comparing C addresses, not by ASSOCIATED.  A
    ! pointer obtained from C_F_POINTER is not required to compare associated
    ! with the module target it happens to address, and nvfortran indeed says
    ! it is not: the slot would never be found, every free would be a silent
    ! no-op, and the pool would fill after LIBCINT_MAX_OPT optimizers.
    subroutine libcint_del_optimizer(opt)
        type(c_ptr), intent(inout) :: opt
        integer :: i
        if (c_associated(opt)) then
            do i = 1, LIBCINT_MAX_OPT
                if (c_associated(opt, c_loc(opt_pool(i)))) then
                    call cint_del_optimizer(opt_pool(i))
                    opt_taken(i) = .false.
                    exit
                end if
            end do
        end if
        opt = c_null_ptr
    end subroutine libcint_del_optimizer

    ! Point the workspace at the caller's optimizer for the duration of one
    ! call.  The C passes CINTOpt as an argument to every entry point; here it
    ! rides on the workspace, so the front door is where the two conventions
    ! meet.  A null or absent pointer means no optimizer, which is exactly
    ! what libcint's own API says a NULL CINTOpt means.
    subroutine ws_bind_opt(opt)
        type(c_ptr), intent(in), optional :: opt
        type(cint_opt_t), pointer :: p
        ws%opt => null()
        if (.not. present(opt)) return
        if (.not. c_associated(opt)) return
        call c_f_pointer(opt, p)
        ws%opt => p
    end subroutine ws_bind_opt

    ! The next free slot in the pool.  Serial by construction: an optimizer is
    ! built before the work that uses it, not during.
    integer function opt_claim() result(slot)
        integer :: i
        do i = 1, LIBCINT_MAX_OPT
            if (.not. opt_taken(i)) then
                opt_taken(i) = .true.
                slot = i
                return
            end if
        end do
        error stop "libcint_fortran: more than LIBCINT_MAX_OPT live optimizers"
    end function opt_claim

    ! env_extent bounds env from the shells a call names; an optimizer names
    ! none, so it needs the bound over every shell as well as every atom.
    pure function env_extent_all(atm, natm, bas, nbas) result(n)
        integer(ip), intent(in) :: natm, nbas
        integer(ip), intent(in) :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip) :: n
        integer(ip) :: i, np, nc
        n = LIBCINT_PTR_ENV_START
        do i = 0, nbas - 1
            np = bas(LIBCINT_BAS_SLOTS*i + LIBCINT_NPRIM_OF - 1)
            nc = bas(LIBCINT_BAS_SLOTS*i + LIBCINT_NCTR_OF  - 1)
            n = max(n, bas(LIBCINT_BAS_SLOTS*i + LIBCINT_PTR_EXP   - 1) + np)
            n = max(n, bas(LIBCINT_BAS_SLOTS*i + LIBCINT_PTR_COEFF - 1) + np*nc)
        end do
        do i = 0, natm - 1
            n = max(n, atm(LIBCINT_ATM_SLOTS*i + LIBCINT_PTR_COORD - 1) + 3_ip)
            n = max(n, atm(LIBCINT_ATM_SLOTS*i + LIBCINT_PTR_ZETA  - 1) + 2_ip)
        end do
    end function env_extent_all

    ! ========================================================================
    ! The derivative and spinor families.  These were stubs while they were
    ! D9 and D10 work; they run on the port now, like everything else here.
    ! ========================================================================

    function libcint_1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), di, dj, nenv
        logical :: hv
        ws%opt => null()
        fshls = 0; fshls(0) = shls(1); fshls(1) = shls(2)
        di = cint_cgto_cart(fshls(0), bas);    dj = cint_cgto_cart(fshls(1), bas)
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int1e_ipnuc_cart(buf(1:di*dj*3), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_1e_ipnuc_cart

    function libcint_1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), di, dj, nenv
        logical :: hv
        ws%opt => null()
        fshls = 0; fshls(0) = shls(1); fshls(1) = shls(2)
        di = cint_cgto_spheric(fshls(0), bas); dj = cint_cgto_spheric(fshls(1), bas)
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int1e_ipnuc_sph(buf(1:di*dj*3), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_1e_ipnuc_sph

    function libcint_1e_iprinv_cart(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), di, dj, nenv
        logical :: hv
        ws%opt => null()
        fshls = 0; fshls(0) = shls(1); fshls(1) = shls(2)
        di = cint_cgto_cart(fshls(0), bas);    dj = cint_cgto_cart(fshls(1), bas)
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int1e_iprinv_cart(buf(1:di*dj*3), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_1e_iprinv_cart

    function libcint_1e_iprinv_sph(buf, shls, atm, natm, bas, nbas, env) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), di, dj, nenv
        logical :: hv
        ws%opt => null()
        fshls = 0; fshls(0) = shls(1); fshls(1) = shls(2)
        di = cint_cgto_spheric(fshls(0), bas); dj = cint_cgto_spheric(fshls(1), bas)
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int1e_iprinv_sph(buf(1:di*dj*3), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_1e_iprinv_sph

    function libcint_2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(4)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 4
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_cart(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 4_ip, atm, natm, bas, nbas)
        hv = int2e_ip1_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_2e_ip1_cart

    function libcint_2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(4)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 4
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_spheric(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 4_ip, atm, natm, bas, nbas)
        hv = int2e_ip1_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_2e_ip1_sph

    function libcint_3c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(3)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 3
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_cart(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 3_ip, atm, natm, bas, nbas)
        hv = int3c2e_ip1_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_3c2e_ip1_cart

    function libcint_3c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(3)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 3
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_spheric(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 3_ip, atm, natm, bas, nbas)
        hv = int3c2e_ip1_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_3c2e_ip1_sph

    function libcint_3c2e_ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(3)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 3
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_cart(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 3_ip, atm, natm, bas, nbas)
        hv = int3c2e_ip2_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_3c2e_ip2_cart

    function libcint_3c2e_ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(3)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 3
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_spheric(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 3_ip, atm, natm, bas, nbas)
        hv = int3c2e_ip2_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_3c2e_ip2_sph

    function libcint_2c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 2
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_cart(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int2c2e_ip1_cart(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_2c2e_ip1_cart

    function libcint_2c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) result(ret)
        real(dp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer(ip) :: ret
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 2
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_spheric(fshls(i-1), bas)
        end do
        n = product(dims) * 3
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int2c2e_ip1_sph(buf(1:n), dims, fshls, atm, natm, bas, nbas, env(1:nenv), ws)
        ret = merge(1_ip, 0_ip, hv)
    end function libcint_2c2e_ip1_sph

    ! A subroutine, not a function, because that is what the C-backed
    ! libcint_fortran declares.  A caller writing `call libcint_1e_spnucsp(...)`
    ! has to compile against either backend, and it did not.
    subroutine libcint_1e_spnucsp(buf, shls, atm, natm, bas, nbas, env)
        complex(zp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(2)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        integer :: dims(0:3), fshls(0:3), di, dj, nenv
        logical :: hv
        ws%opt => null()
        fshls = 0; fshls(0) = shls(1); fshls(1) = shls(2)
        di = cint_cgto_spinor(fshls(0), bas); dj = cint_cgto_spinor(fshls(1), bas)
        dims = [di, dj, 1, 1]
        nenv = env_extent(shls, 2_ip, atm, natm, bas, nbas)
        hv = int1e_spnucsp_spinor(buf(1:di*dj), dims, fshls, atm, natm, bas, nbas, &
                                  env(1:nenv), ws)
    end subroutine libcint_1e_spnucsp

    subroutine libcint_2e_spsp1(buf, shls, atm, natm, bas, nbas, env, opt)
        complex(zp), intent(out), target :: buf(*)
        integer(ip), intent(in), target :: shls(4)
        integer(ip), intent(in), target :: atm(0:LIBCINT_ATM_SLOTS*natm - 1)
        integer(ip), intent(in) :: natm
        integer(ip), intent(in), target :: bas(0:LIBCINT_BAS_SLOTS*nbas - 1)
        integer(ip), intent(in) :: nbas
        real(dp), intent(in), target :: env(*)
        type(c_ptr), intent(in), optional :: opt
        integer :: dims(0:3), fshls(0:3), i, n, nenv
        logical :: hv
        call ws_bind_opt(opt)
        fshls = 0; dims = 1
        do i = 1, 4
            fshls(i-1) = shls(i)
            dims(i-1) = cint_cgto_spinor(fshls(i-1), bas)
        end do
        n = product(dims)
        nenv = env_extent(shls, 4_ip, atm, natm, bas, nbas)
        hv = int2e_spsp1_spinor(buf(1:n), dims, fshls, atm, natm, bas, nbas, &
                                env(1:nenv), ws)
    end subroutine libcint_2e_spsp1

end module libcint_fortran
