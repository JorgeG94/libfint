!
! Low-level Fortran interface for libcint using iso_c_binding
!
! This module provides type-safe, direct bindings to libcint C functions.
! All types use C-compatible kinds (c_int, c_double, c_double_complex).
!
! For a higher-level API with native Fortran types, see libcint_fortran module.
!
module libcint_interface
    use iso_c_binding
    implicit none
    private

    ! ========================================================================
    ! Public API
    ! ========================================================================

    ! Constants
    public :: CHARGE_OF, PTR_COORD, NUC_MOD_OF, PTR_ZETA, ATM_SLOTS
    public :: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, BAS_SLOTS
    public :: PTR_EXPCUTOFF, PTR_COMMON_ORIG, PTR_RINV_ORIG, PTR_RINV_ZETA
    public :: PTR_RANGE_OMEGA, PTR_ENV_START

    ! Dimension and normalization functions
    public :: CINTcgto_cart, CINTcgto_spheric, CINTcgto_spinor
    public :: CINTtot_cgto_spheric, CINTtot_cgto_cart, CINTtot_pgto_spheric
    public :: CINTgto_norm

    ! One-electron integrals
    public :: cint1e_ovlp_cart, cint1e_nuc_cart, cint1e_kin_cart, cint1e_ipovlp_cart
    public :: cint1e_ovlp_sph, cint1e_nuc_sph, cint1e_kin_sph, cint1e_ipovlp_sph
    public :: cint1e_ipkin_cart, cint1e_ipnuc_cart, cint1e_iprinv_cart
    public :: cint1e_ipkin_sph, cint1e_ipnuc_sph, cint1e_iprinv_sph
    public :: cint1e_spnucsp

    ! Two-electron integrals
    public :: cint2e_cart, cint2e_sph
    public :: cint3c2e_sph, cint2c2e_sph
    public :: cint3c2e_cart, cint2c2e_cart
    public :: cint2e_ip1_cart, cint2e_ip1_sph
    public :: cint3c2e_ip1_cart, cint3c2e_ip1_sph
    public :: cint3c2e_ip2_cart, cint3c2e_ip2_sph
    public :: cint2c2e_ip1_cart, cint2c2e_ip1_sph
    public :: cint2e_spsp1

    ! Optimizers
    public :: cint2e_cart_optimizer, cint2e_sph_optimizer
    public :: cint3c2e_cart_optimizer, cint3c2e_sph_optimizer
    public :: cint2c2e_cart_optimizer, cint2c2e_sph_optimizer
    public :: cint2e_ip1_cart_optimizer, cint2e_ip1_sph_optimizer
    public :: cint3c2e_ip1_cart_optimizer, cint3c2e_ip1_sph_optimizer
    public :: cint3c2e_ip2_cart_optimizer, cint3c2e_ip2_sph_optimizer
    public :: cint2c2e_ip1_cart_optimizer, cint2c2e_ip1_sph_optimizer
    public :: cint2e_spsp1_optimizer
    public :: CINTdel_optimizer

    ! ========================================================================
    ! Constants from cint.h
    ! ========================================================================

    ! Slots of atm array (1-indexed in Fortran, but values are 0-indexed offsets for C)
    integer(c_int), parameter :: CHARGE_OF  = 1
    integer(c_int), parameter :: PTR_COORD  = 2
    integer(c_int), parameter :: NUC_MOD_OF = 3
    integer(c_int), parameter :: PTR_ZETA   = 4
    integer(c_int), parameter :: ATM_SLOTS  = 6

    ! Slots of bas array
    integer(c_int), parameter :: ATOM_OF    = 1
    integer(c_int), parameter :: ANG_OF     = 2
    integer(c_int), parameter :: NPRIM_OF   = 3
    integer(c_int), parameter :: NCTR_OF    = 4
    integer(c_int), parameter :: KAPPA_OF   = 5
    integer(c_int), parameter :: PTR_EXP    = 6
    integer(c_int), parameter :: PTR_COEFF  = 7
    integer(c_int), parameter :: BAS_SLOTS  = 8

    ! Global parameters in env array
    integer(c_int), parameter :: PTR_EXPCUTOFF    = 0
    integer(c_int), parameter :: PTR_COMMON_ORIG  = 1
    integer(c_int), parameter :: PTR_RINV_ORIG    = 4
    integer(c_int), parameter :: PTR_RINV_ZETA    = 7
    integer(c_int), parameter :: PTR_RANGE_OMEGA  = 8
    integer(c_int), parameter :: PTR_ENV_START    = 20

    ! ========================================================================
    ! Interface to C functions
    ! ========================================================================

    interface

        ! ====================================================================
        ! Basis set dimension functions
        ! ====================================================================

        ! Get number of Cartesian GTOs for a given shell
        function CINTcgto_cart(bas_id, bas) bind(C, name='CINTcgto_cart')
            import :: c_int
            integer(c_int), value, intent(in) :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_cart
        end function CINTcgto_cart

        ! Get number of spherical GTOs for a given shell
        function CINTcgto_spheric(bas_id, bas) bind(C, name='CINTcgto_spheric')
            import :: c_int
            integer(c_int), value, intent(in) :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_spheric
        end function CINTcgto_spheric

        ! Get number of spinor components for a given shell
        function CINTcgto_spinor(bas_id, bas) bind(C, name='CINTcgto_spinor')
            import :: c_int
            integer(c_int), value, intent(in) :: bas_id
            integer(c_int), intent(in) :: bas(*)
            integer(c_int) :: CINTcgto_spinor
        end function CINTcgto_spinor

        ! ====================================================================
        ! GTO normalization and total counts
        ! ====================================================================

        ! Compute GTO normalization factor
        function CINTgto_norm(n, a) bind(C, name='CINTgto_norm')
            import :: c_int, c_double
            integer(c_int), value, intent(in) :: n
            real(c_double), value, intent(in) :: a
            real(c_double) :: CINTgto_norm
        end function CINTgto_norm

        ! Total number of contracted spherical GTOs
        function CINTtot_cgto_spheric(bas, nbas) bind(C, name='CINTtot_cgto_spheric')
            import :: c_int
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            integer(c_int) :: CINTtot_cgto_spheric
        end function CINTtot_cgto_spheric

        ! Total contracted Cartesian functions across all shells.
        function CINTtot_cgto_cart(bas, nbas) bind(C, name='CINTtot_cgto_cart')
            import :: c_int
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            integer(c_int) :: CINTtot_cgto_cart
        end function CINTtot_cgto_cart

        ! Total number of primitive spherical GTOs
        function CINTtot_pgto_spheric(bas, nbas) bind(C, name='CINTtot_pgto_spheric')
            import :: c_int
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            integer(c_int) :: CINTtot_pgto_spheric
        end function CINTtot_pgto_spheric

        ! ====================================================================
        ! Optimizer functions
        ! ====================================================================

        ! Delete/free an optimizer object
        subroutine CINTdel_optimizer(opt) bind(C, name='CINTdel_optimizer')
            import :: c_ptr
            type(c_ptr), intent(inout) :: opt
        end subroutine CINTdel_optimizer

        ! ====================================================================
        ! One-electron integrals - Cartesian
        ! ====================================================================

        ! Overlap integral (Cartesian)
        function cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ovlp_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ovlp_cart
        end function cint1e_ovlp_cart

        ! Nuclear attraction integral (Cartesian)
        function cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_nuc_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_nuc_cart
        end function cint1e_nuc_cart

        ! Kinetic energy integral (Cartesian)
        function cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_kin_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_kin_cart
        end function cint1e_kin_cart

        ! Overlap gradient integral (Cartesian)
        function cint1e_ipovlp_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipovlp_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipovlp_cart
        end function cint1e_ipovlp_cart

        ! ====================================================================
        ! One-electron integrals - Spherical
        ! ====================================================================

        ! Overlap integral (Spherical)
        function cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ovlp_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ovlp_sph
        end function cint1e_ovlp_sph

        ! Nuclear attraction integral (Spherical)
        function cint1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_nuc_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_nuc_sph
        end function cint1e_nuc_sph

        ! Kinetic energy integral (Spherical)
        function cint1e_kin_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_kin_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_kin_sph
        end function cint1e_kin_sph

        ! Overlap gradient integral (Spherical)
        function cint1e_ipovlp_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipovlp_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipovlp_sph
        end function cint1e_ipovlp_sph

        ! ====================================================================
        ! One-electron integrals - Spinor
        ! ====================================================================

        ! Spinor nuclear attraction integral
        subroutine cint1e_spnucsp(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_spnucsp')
            import :: c_double_complex, c_int, c_double
            complex(c_double_complex), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint1e_spnucsp

        ! ====================================================================
        ! Two-electron integrals - Cartesian
        ! ====================================================================

        ! Electron repulsion integral (Cartesian)
        function cint2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2e_cart
        end function cint2e_cart

        ! Create optimizer for 2e integrals (Cartesian)
        subroutine cint2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_cart_optimizer

        ! ERI gradient (Cartesian)
        function cint2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_ip1_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2e_ip1_cart
        end function cint2e_ip1_cart

        ! Create optimizer for ERI gradients (Cartesian)
        subroutine cint2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_ip1_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_ip1_cart_optimizer

        ! ====================================================================
        ! Two-electron integrals - Spherical
        ! ====================================================================

        ! Electron repulsion integral (Spherical)
        function cint2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2e_sph
        end function cint2e_sph

        ! Three-centre two-electron integrals (Spherical), for density fitting.
        ! shls has four entries even though only three centres contribute: the
        ! fourth names a dummy s shell with a zero exponent, which is how
        ! libcint spells "no fourth centre".
        function cint3c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_sph
        end function cint3c2e_sph

        ! Two-centre two-electron integrals (Spherical): the fitting metric.
        function cint2c2e_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2c2e_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2c2e_sph
        end function cint2c2e_sph

        ! Three-centre two-electron integrals (Cartesian), for density
        ! fitting in a Cartesian basis. Same dummy-fourth-shell convention as
        ! the spherical form above.
        function cint3c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_cart
        end function cint3c2e_cart

        ! Two-centre two-electron integrals (Cartesian): the fitting metric.
        function cint2c2e_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2c2e_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2c2e_cart
        end function cint2c2e_cart

        ! Create optimizer for 2e integrals (Spherical)
        subroutine cint2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_sph_optimizer

        ! Create optimizer for three-centre 2e integrals (Spherical)
        subroutine cint3c2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_sph_optimizer

        ! Create optimizer for three-centre 2e integrals (Cartesian)
        subroutine cint3c2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_cart_optimizer

        ! Create optimizer for two-centre 2e integrals (Spherical)
        subroutine cint2c2e_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2c2e_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2c2e_sph_optimizer

        ! Create optimizer for two-centre 2e integrals (Cartesian)
        subroutine cint2c2e_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2c2e_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2c2e_cart_optimizer

        ! ERI gradient (Spherical)
        function cint2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_ip1_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2e_ip1_sph
        end function cint2e_ip1_sph

        ! Create optimizer for ERI gradients (Spherical)
        subroutine cint2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_ip1_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_ip1_sph_optimizer

        ! ====================================================================
        ! Two-electron integrals - Spinor
        ! ====================================================================

        ! Spinor electron repulsion integral
        subroutine cint2e_spsp1(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2e_spsp1')
            import :: c_double_complex, c_int, c_ptr, c_double
            complex(c_double_complex), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
        end subroutine cint2e_spsp1

        ! Create optimizer for spinor 2e integrals
        subroutine cint2e_spsp1_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2e_spsp1_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2e_spsp1_optimizer

        ! ====================================================================
        ! Nuclear-derivative integrals, for analytic gradients
        !
        ! Each returns three components per shell pair, the derivative with
        ! respect to the coordinates of the *first* shell's centre. Translational
        ! invariance gives the second centre's, which is why libcint generates
        ! only one of them.
        !
        ! iprinv is the exception and the reason gradients need it: it
        ! differentiates with respect to the origin of the 1/|r-R| operator
        ! rather than a basis centre, which is the Hellmann-Feynman term. The
        ! caller sets that origin in env(PTR_RINV_ORIG+1:+3) before calling.
        ! ====================================================================

        ! Gradient of the kinetic energy integral (Spherical)
        function cint1e_ipkin_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipkin_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipkin_sph
        end function cint1e_ipkin_sph

        ! Gradient of the kinetic energy integral (Cartesian)
        function cint1e_ipkin_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipkin_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipkin_cart
        end function cint1e_ipkin_cart

        ! Gradient of the nuclear attraction integral (Spherical)
        function cint1e_ipnuc_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipnuc_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipnuc_sph
        end function cint1e_ipnuc_sph

        ! Gradient of the nuclear attraction integral (Cartesian)
        function cint1e_ipnuc_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_ipnuc_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_ipnuc_cart
        end function cint1e_ipnuc_cart

        ! Gradient of the 1/|r-R| operator integral (Spherical)
        function cint1e_iprinv_sph(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_iprinv_sph')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_iprinv_sph
        end function cint1e_iprinv_sph

        ! Gradient of the 1/|r-R| operator integral (Cartesian)
        function cint1e_iprinv_cart(buf, shls, atm, natm, bas, nbas, env) &
                bind(C, name='cint1e_iprinv_cart')
            import :: c_double, c_int
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            integer(c_int) :: cint1e_iprinv_cart
        end function cint1e_iprinv_cart

        ! Gradient of the three-centre 2e, differentiated on the first index (Spherical)
        function cint3c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_ip1_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_ip1_sph
        end function cint3c2e_ip1_sph

        ! Gradient of the three-centre 2e, differentiated on the first index (Cartesian)
        function cint3c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_ip1_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_ip1_cart
        end function cint3c2e_ip1_cart

        ! Gradient of the three-centre 2e, differentiated on the auxiliary index (Spherical)
        function cint3c2e_ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_ip2_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_ip2_sph
        end function cint3c2e_ip2_sph

        ! Gradient of the three-centre 2e, differentiated on the auxiliary index (Cartesian)
        function cint3c2e_ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint3c2e_ip2_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint3c2e_ip2_cart
        end function cint3c2e_ip2_cart

        ! Gradient of the two-centre 2e (the fitting metric), differentiated on the first index (Spherical)
        function cint2c2e_ip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2c2e_ip1_sph')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2c2e_ip1_sph
        end function cint2c2e_ip1_sph

        ! Gradient of the two-centre 2e (the fitting metric), differentiated on the first index (Cartesian)
        function cint2c2e_ip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
                bind(C, name='cint2c2e_ip1_cart')
            import :: c_double, c_int, c_ptr
            real(c_double), intent(out) :: buf(*)
            integer(c_int), intent(in) :: shls(*)
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
            type(c_ptr), value, intent(in) :: opt
            integer(c_int) :: cint2c2e_ip1_cart
        end function cint2c2e_ip1_cart

        ! Optimizer for the three-centre 2e, differentiated on the first index (Spherical)
        subroutine cint3c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_ip1_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_ip1_sph_optimizer

        ! Optimizer for the three-centre 2e, differentiated on the first index (Cartesian)
        subroutine cint3c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_ip1_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_ip1_cart_optimizer

        ! Optimizer for the three-centre 2e, differentiated on the auxiliary index (Spherical)
        subroutine cint3c2e_ip2_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_ip2_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_ip2_sph_optimizer

        ! Optimizer for the three-centre 2e, differentiated on the auxiliary index (Cartesian)
        subroutine cint3c2e_ip2_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint3c2e_ip2_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint3c2e_ip2_cart_optimizer

        ! Optimizer for the two-centre 2e (the fitting metric), differentiated on the first index (Spherical)
        subroutine cint2c2e_ip1_sph_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2c2e_ip1_sph_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2c2e_ip1_sph_optimizer

        ! Optimizer for the two-centre 2e (the fitting metric), differentiated on the first index (Cartesian)
        subroutine cint2c2e_ip1_cart_optimizer(opt, atm, natm, bas, nbas, env) &
                bind(C, name='cint2c2e_ip1_cart_optimizer')
            import :: c_ptr, c_int, c_double
            type(c_ptr), intent(out) :: opt
            integer(c_int), intent(in) :: atm(*)
            integer(c_int), value, intent(in) :: natm
            integer(c_int), intent(in) :: bas(*)
            integer(c_int), value, intent(in) :: nbas
            real(c_double), intent(in) :: env(*)
        end subroutine cint2c2e_ip1_cart_optimizer

    end interface

end module libcint_interface
