!
! D8 acceptance: the libcint_fortran module, as a caller sees it.
!
! Only built with WITH_FORTRAN_BACKEND, where libcint_fortran is the port and
! libcint_interface is still the C -- so one program can put the two side by
! side through exactly the API metalquicha uses, 2D (SLOTS, *) tables and all.
!
! This is a different test from int1e_check/int2e_check/int3c2e_check, which
! call the port's own entry points directly.  What it adds is the shim: the
! sequence association that reinterprets the caller's tables, the env extent
! it works out for itself, the output dimensions it derives rather than being
! told.  Any of those could be wrong while every engine test still passes.
!
! THE BAR IS BIT-IDENTICAL UP TO FIVE RYS ROOTS AND 1e-11 ABOVE, and the split
! is the one D3 predicted rather than one chosen to make this pass.  Up to five
! roots the quadrature comes from the closed forms, which the port reproduces
! exactly; above that it comes from a tridiagonal eigensolve, and D3 replaced
! the C's vendored copy of DSTEMR with LAPACK's.  Asked to referee the two by
! a symmetry the integral itself obeys -- (ij|kl) = (ji|kl), which the exact
! answer satisfies and neither engine quite does -- the port comes out very
! slightly ahead:
!
!     roots   quartets       C          port
!       4        3720     1.13e-14    1.13e-14      (identical)
!       5        2528     3.17e-14    3.17e-14      (identical)
!       6         714     8.25e-14    7.43e-14
!       7          64     2.02e-13    1.92e-13
!
! so this is the accuracy floor of two different eigensolvers, not a defect on
! either side.  The tolerance sits two orders above that floor.
!
program backend_check
   use iso_c_binding
   use libcint_interface, only: cint1e_ovlp_cart, cint1e_ovlp_sph, &
                                cint1e_kin_cart,  cint1e_kin_sph,  &
                                cint1e_nuc_cart,  cint1e_nuc_sph,  &
                                cint2e_cart,      cint2e_sph,      &
                                cint3c2e_cart,    cint3c2e_sph,    &
                                cint2c2e_cart,    cint2c2e_sph,    &
                                cint1e_ipnuc_cart,  cint1e_ipnuc_sph,  &
                                cint1e_iprinv_cart, cint1e_iprinv_sph, &
                                cint2e_ip1_cart,    cint2e_ip1_sph,    &
                                cint3c2e_ip1_cart,  cint3c2e_ip1_sph,  &
                                cint3c2e_ip2_cart,  cint3c2e_ip2_sph,  &
                                cint2c2e_ip1_cart,  cint2c2e_ip1_sph,  &
                                CINTcgto_cart,    CINTcgto_spheric, &
                                ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, ANG_OF, &
                                PTR_COMMON_ORIG, PTR_RINV_ORIG
   use libcint_fortran
   use cint_test_systems
   implicit none

   ! libcint's own entry points for the integrals libcint_interface does not
   ! declare, here rather than there: other checks `use libcint_interface`
   ! whole and declare some of these names themselves.
   interface
      ! Dipole about PTR_COMMON_ORIG (Spherical)
      function cint1e_r_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_r_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_r_sph
      end function cint1e_r_sph

      ! Dipole about PTR_COMMON_ORIG (Cartesian)
      function cint1e_r_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_r_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_r_cart
      end function cint1e_r_cart

      ! Second moment about PTR_COMMON_ORIG (Spherical)
      function cint1e_rr_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_rr_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_rr_sph
      end function cint1e_rr_sph

      ! Second moment about PTR_COMMON_ORIG (Cartesian)
      function cint1e_rr_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_rr_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_rr_cart
      end function cint1e_rr_cart

      ! Third moment about PTR_COMMON_ORIG (Spherical)
      function cint1e_rrr_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_rrr_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_rrr_sph
      end function cint1e_rrr_sph

      ! Third moment about PTR_COMMON_ORIG (Cartesian)
      function cint1e_rrr_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_rrr_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_rrr_cart
      end function cint1e_rrr_cart

      ! <i| r nabla |j>, the dipole's nuclear derivative (Spherical)
      function cint1e_irp_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_irp_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_irp_sph
      end function cint1e_irp_sph

      ! <i| r nabla |j>, the dipole's nuclear derivative (Cartesian)
      function cint1e_irp_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_irp_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_irp_cart
      end function cint1e_irp_cart

      ! Overlap, second derivative on the bra (Spherical)
      function cint1e_ipipovlp_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipovlp_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipovlp_sph
      end function cint1e_ipipovlp_sph

      ! Overlap, second derivative on the bra (Cartesian)
      function cint1e_ipipovlp_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipovlp_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipovlp_cart
      end function cint1e_ipipovlp_cart

      ! Overlap, first derivative on bra and ket (Spherical)
      function cint1e_ipovlpip_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipovlpip_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipovlpip_sph
      end function cint1e_ipovlpip_sph

      ! Overlap, first derivative on bra and ket (Cartesian)
      function cint1e_ipovlpip_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipovlpip_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipovlpip_cart
      end function cint1e_ipovlpip_cart

      ! Kinetic, second derivative on the bra (Spherical)
      function cint1e_ipipkin_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipkin_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipkin_sph
      end function cint1e_ipipkin_sph

      ! Kinetic, second derivative on the bra (Cartesian)
      function cint1e_ipipkin_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipkin_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipkin_cart
      end function cint1e_ipipkin_cart

      ! Kinetic, first derivative on bra and ket (Spherical)
      function cint1e_ipkinip_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipkinip_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipkinip_sph
      end function cint1e_ipkinip_sph

      ! Kinetic, first derivative on bra and ket (Cartesian)
      function cint1e_ipkinip_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipkinip_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipkinip_cart
      end function cint1e_ipkinip_cart

      ! Nuclear attraction, second derivative on the bra (Spherical)
      function cint1e_ipipnuc_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipnuc_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipnuc_sph
      end function cint1e_ipipnuc_sph

      ! Nuclear attraction, second derivative on the bra (Cartesian)
      function cint1e_ipipnuc_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipipnuc_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipipnuc_cart
      end function cint1e_ipipnuc_cart

      ! Nuclear attraction, first derivative on bra and ket (Spherical)
      function cint1e_ipnucip_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipnucip_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipnucip_sph
      end function cint1e_ipnucip_sph

      ! Nuclear attraction, first derivative on bra and ket (Cartesian)
      function cint1e_ipnucip_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipnucip_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipnucip_cart
      end function cint1e_ipnucip_cart

      ! 1/|r-R| about PTR_RINV_ORIG, second derivative on the bra (Spherical)
      function cint1e_ipiprinv_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipiprinv_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipiprinv_sph
      end function cint1e_ipiprinv_sph

      ! 1/|r-R| about PTR_RINV_ORIG, second derivative on the bra (Cartesian)
      function cint1e_ipiprinv_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_ipiprinv_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_ipiprinv_cart
      end function cint1e_ipiprinv_cart

      ! 1/|r-R| about PTR_RINV_ORIG, first derivative on bra and ket (Spherical)
      function cint1e_iprinvip_sph(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_iprinvip_sph')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_iprinvip_sph
      end function cint1e_iprinvip_sph

      ! 1/|r-R| about PTR_RINV_ORIG, first derivative on bra and ket (Cartesian)
      function cint1e_iprinvip_cart(buf, shls, atm, natm, bas, nbas, env) &
              bind(C, name='cint1e_iprinvip_cart')
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         integer(c_int) :: cint1e_iprinvip_cart
      end function cint1e_iprinvip_cart

      ! ERI, second derivative on the first centre (Spherical)
      function cint2e_ipip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ipip1_sph')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ipip1_sph
      end function cint2e_ipip1_sph

      ! ERI, second derivative on the first centre (Cartesian)
      function cint2e_ipip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ipip1_cart')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ipip1_cart
      end function cint2e_ipip1_cart

      ! ERI, first derivative on the first two centres (Spherical)
      function cint2e_ipvip1_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ipvip1_sph')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ipvip1_sph
      end function cint2e_ipvip1_sph

      ! ERI, first derivative on the first two centres (Cartesian)
      function cint2e_ipvip1_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ipvip1_cart')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ipvip1_cart
      end function cint2e_ipvip1_cart

      ! ERI, first derivative on the first and on the second electron (Spherical)
      function cint2e_ip1ip2_sph(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ip1ip2_sph')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ip1ip2_sph
      end function cint2e_ip1ip2_sph

      ! ERI, first derivative on the first and on the second electron (Cartesian)
      function cint2e_ip1ip2_cart(buf, shls, atm, natm, bas, nbas, env, opt) &
              bind(C, name='cint2e_ip1ip2_cart')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         integer(c_int) :: cint2e_ip1ip2_cart
      end function cint2e_ip1ip2_cart

      ! Field of a point charge at PTR_RINV_ORIG, derivative on the charge (Spherical).
      ! libcint has no CINT2 spelling of this one, so it is the CINT3 form:
      ! dims, opt and cache may all be null.
      function int1e_drinv_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
              bind(C, name='int1e_drinv_sph')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         type(c_ptr), value, intent(in) :: dims
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         type(c_ptr), value, intent(in) :: cache
         integer(c_int) :: int1e_drinv_sph
      end function int1e_drinv_sph

      ! Field of a point charge at PTR_RINV_ORIG, derivative on the charge (Cartesian).
      ! libcint has no CINT2 spelling of this one, so it is the CINT3 form:
      ! dims, opt and cache may all be null.
      function int1e_drinv_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
              bind(C, name='int1e_drinv_cart')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         type(c_ptr), value, intent(in) :: dims
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         type(c_ptr), value, intent(in) :: cache
         integer(c_int) :: int1e_drinv_cart
      end function int1e_drinv_cart

      ! 1/|r-R_g| over the grid at env(PTR_GRIDS); shls(3:4) is the 0-based half-open range (Spherical).
      ! libcint has no CINT2 spelling of this one, so it is the CINT3 form:
      ! dims, opt and cache may all be null.
      function int1e_grids_sph(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
              bind(C, name='int1e_grids_sph')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         type(c_ptr), value, intent(in) :: dims
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         type(c_ptr), value, intent(in) :: cache
         integer(c_int) :: int1e_grids_sph
      end function int1e_grids_sph

      ! 1/|r-R_g| over the grid at env(PTR_GRIDS); shls(3:4) is the 0-based half-open range (Cartesian).
      ! libcint has no CINT2 spelling of this one, so it is the CINT3 form:
      ! dims, opt and cache may all be null.
      function int1e_grids_cart(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
              bind(C, name='int1e_grids_cart')
         import :: c_double, c_int, c_ptr
         real(c_double), intent(out) :: buf(*)
         type(c_ptr), value, intent(in) :: dims
         integer(c_int), intent(in) :: shls(*)
         integer(c_int), intent(in) :: atm(*)
         integer(c_int), value, intent(in) :: natm
         integer(c_int), intent(in) :: bas(*)
         integer(c_int), value, intent(in) :: nbas
         real(c_double), intent(in) :: env(*)
         type(c_ptr), value, intent(in) :: opt
         type(c_ptr), value, intent(in) :: cache
         integer(c_int) :: int1e_grids_cart
      end function int1e_grids_cart
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, di, dj, dk, dl, nq, stride, nc, nf
   integer :: ncmp, nbad, nexact, nloose
   real(dp) :: worst
   character(len=80) :: wat
   type(c_ptr) :: opt
   integer, parameter :: NGRID = 5
   integer :: g

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ! 15^4 * 9: a Cartesian g quartet of a second derivative
   allocate(cbuf(460000), fbuf(460000))
   ncmp = 0; nbad = 0; nexact = 0; nloose = 0; worst = 0; wat = "(none)"

   do ibasis = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)
      ! Off-centre origins, so a multipole or a point-charge field read from
      ! the wrong env slot -- or about zero -- cannot agree by accident.
      env(PTR_COMMON_ORIG+1:PTR_COMMON_ORIG+3) = [0.31_dp, -0.47_dp, 0.83_dp]
      env(PTR_RINV_ORIG+1:PTR_RINV_ORIG+3)     = [-0.62_dp, 0.25_dp, -0.39_dp]
      env(LIBCINT_PTR_GRIDS+1) = real(off, dp)
      do g = 0, NGRID-1
         env(off+3*g+1:off+3*g+3) = [0.4_dp*g - 0.7_dp, 0.9_dp - 0.3_dp*g, 0.2_dp*g]
      end do
      off = off + 3*NGRID

      ! The bookkeeping first: if these disagree the buffers are the wrong
      ! size and every comparison below is meaningless.
      do i = 0, nbas-1
         call check_int(int(CINTcgto_cart(i, bas)), int(libcint_cgto_cart(i, bas)), &
                        "cgto_cart", ibasis, i)
         call check_int(int(CINTcgto_spheric(i, bas)), int(libcint_cgto_sph(i, bas)), &
                        "cgto_sph", ibasis, i)
      end do

      ! An optimizer, requested and then ignored, because that is what a
      ! caller does and the null it gets back has to be harmless.
      call libcint_2e_sph_optimizer(opt, atm, 8_ip, bas, nbas, env)

      stride = max(1, int(nbas)/6)

      do sph = 0, 1
         ! ---- one electron, every pair ----
         do i = 0, nbas-1
         do j = 0, nbas-1
            call dims2(sph, i, j, bas, di, dj)
            nq = di*dj
            shls(1)=i; shls(2)=j
            if (sph == 1) then
               ret = cint1e_ovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ovlp_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ovlp_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "ovlp", ibasis, sph, i, j, 0, 0)

            if (sph == 1) then
               ret = cint1e_kin_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_kin_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_kin_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_kin_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "kin", ibasis, sph, i, j, 0, 0)

            if (sph == 1) then
               ret = cint1e_nuc_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_nuc_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_nuc_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_nuc_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "nuc", ibasis, sph, i, j, 0, 0)

            ! ---- two centre, every pair ----
            call dims2(sph, i, j, bas, di, dj)
            if (sph == 1) then
               ret = cint2c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint2c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj, "2c2e", ibasis, sph, i, j, 0, 0)

            ! ---- the derivative families, which were stubs here until the
            ! port grew them.  Three components each.
            call dims2(sph, i, j, bas, di, dj)
            if (sph == 1) then
               ret = cint1e_ipnuc_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnuc_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipnuc_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnuc_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "ipnuc", ibasis, sph, i, j, 0, 0, 1)

            if (sph == 1) then
               ret = cint1e_iprinv_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinv_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_iprinv_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinv_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "iprinv", ibasis, sph, i, j, 0, 0, 1)

            if (sph == 1) then
               ret = cint2c2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_ip1_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint2c2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_ip1_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "2c2e_ip1", ibasis, sph, i, j, 0, 0, 1)

            ! ---- multipoles, the point-charge field, second derivatives ----
            if (sph == 1) then
               ret = cint1e_r_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_r_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_r_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_r_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "r", ibasis, sph, i, j, 0, 0, 1)

            if (sph == 1) then
               ret = cint1e_rr_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_rr_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_rr_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_rr_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "rr", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_rrr_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_rrr_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_rrr_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_rrr_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*27, "rrr", ibasis, sph, i, j, 0, 0, 3)

            if (sph == 1) then
               ret = cint1e_irp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_irp_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_irp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_irp_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "irp", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipipovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipovlp_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipipovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipovlp_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipipovlp", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipovlpip_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipovlpip_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipovlpip_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipovlpip_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipovlpip", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipipkin_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipkin_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipipkin_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipkin_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipipkin", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipkinip_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipkinip_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipkinip_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipkinip_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipkinip", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipipnuc_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipnuc_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipipnuc_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipipnuc_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipipnuc", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipnucip_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnucip_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipnucip_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnucip_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipnucip", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_ipiprinv_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipiprinv_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipiprinv_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipiprinv_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "ipiprinv", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = cint1e_iprinvip_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinvip_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_iprinvip_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinvip_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*9, "iprinvip", ibasis, sph, i, j, 0, 0, 2)

            if (sph == 1) then
               ret = int1e_drinv_sph(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                                     c_null_ptr, c_null_ptr)
               ret = libcint_1e_drinv_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = int1e_drinv_cart(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                                      c_null_ptr, c_null_ptr)
               ret = libcint_1e_drinv_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "drinv", ibasis, sph, i, j, 0, 0, 1)

            ! a range that starts past the first point, so an offset dropped
            ! in either direction shows
            shls(3) = 1; shls(4) = NGRID
            if (sph == 1) then
               ret = int1e_grids_sph(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                                     c_null_ptr, c_null_ptr)
               ret = libcint_1e_grids_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = int1e_grids_cart(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                                      c_null_ptr, c_null_ptr)
               ret = libcint_1e_grids_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*(NGRID-1), "grids", ibasis, sph, i, j, 0, 0)
         end do
         end do

         ! ---- three and four centre, strided ----
         do i = 0, nbas-1, stride
         do j = 0, nbas-1, stride
         do k = 0, nbas-1, stride
            call dims2(sph, i, j, bas, di, dj)
            call dims1(sph, k, bas, dk)
            shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=0
            if (sph == 1) then
               ret = cint3c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk, "3c2e", ibasis, sph, i, j, k, 0)

            if (sph == 1) then
               ret = cint3c2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk*3, "3c2e_ip1", ibasis, sph, i, j, k, 0, 1)

            if (sph == 1) then
               ret = cint3c2e_ip2_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip2_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_ip2_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip2_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk*3, "3c2e_ip2", ibasis, sph, i, j, k, 0, 1)

            do l = 0, nbas-1, stride
               call dims1(sph, l, bas, dl)
               shls(4)=l
               if (sph == 1) then
                  ret = cint2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_sph(fbuf, shls, atm, 8_ip, bas, nbas, env, opt)
               else
                  ret = cint2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_cart(fbuf, shls, atm, 8_ip, bas, nbas, env, opt)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl, "2e", ibasis, sph, i, j, k, l)

               if (sph == 1) then
                  ret = cint2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
               else
                  ret = cint2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl*3, "2e_ip1", ibasis, sph, i, j, k, l, 1)

               if (sph == 1) then
                  ret = cint2e_ipip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ipip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
               else
                  ret = cint2e_ipip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ipip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl*9, "2e_ipip1", ibasis, sph, i, j, k, l, 2)

               if (sph == 1) then
                  ret = cint2e_ipvip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ipvip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
               else
                  ret = cint2e_ipvip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ipvip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl*9, "2e_ipvip1", ibasis, sph, i, j, k, l, 2)

               if (sph == 1) then
                  ret = cint2e_ip1ip2_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1ip2_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
               else
                  ret = cint2e_ip1ip2_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1ip2_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl*9, "2e_ip1ip2", ibasis, sph, i, j, k, l, 2)
            end do
         end do
         end do
         end do
      end do

      call libcint_del_optimizer(opt)
   end do

   print '(A,I0)',         "  values compared      : ", ncmp
   print '(A,ES10.2)',     "  worst diff / block   : ", worst
   print '(A,I0)',         "  blocks above 5 roots : ", nloose
   print '(A,I0)',         "  inexact at <=5 roots : ", nexact
   print '(A,I0)',         "  over tolerance       : ", nbad
   if (nbad > 0 .or. nexact > 0) then
      print '(A,A)',       "  worst at             : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS through libcint_fortran"

contains

   subroutine dims1(sp, sh, b, d)
      integer, intent(in) :: sp
      integer(c_int), intent(in) :: sh, b(BAS_SLOTS, *)
      integer, intent(out) :: d
      if (sp == 1) then
         d = int(CINTcgto_spheric(sh, b))
      else
         d = int(CINTcgto_cart(sh, b))
      end if
   end subroutine dims1

   subroutine dims2(sp, sa, sb, b, da, db)
      integer, intent(in) :: sp
      integer(c_int), intent(in) :: sa, sb, b(BAS_SLOTS, *)
      integer, intent(out) :: da, db
      call dims1(sp, sa, b, da)
      call dims1(sp, sb, b, db)
   end subroutine dims2

   subroutine check_int(a, b, what, ib, sh)
      integer, intent(in) :: a, b, ib, sh
      character(len=*), intent(in) :: what
      ncmp = ncmp + 1
      if (a /= b) then
         nbad = nbad + 1
         if (worst < 1.0_dp) then
            worst = 1.0_dp
            write(wat,'(A,A,I0,A,I0,A,I0,A,I0)') trim(what), " basis=", ib, &
               " shell=", sh, " C=", a, " F=", b
         end if
      end if
   end subroutine check_int

   ! Rys order for a set of shells: the expression the envs setup uses.
   !
   ! `inc` is the integral's own ng increment total.  It is not optional
   ! bookkeeping: a derivative raises the angular ceiling on the shell it acts
   ! on, so int2e_ip1 over four p shells needs six roots where int2e over the
   ! same four needs five.  Counting from the bare l undercounts, and the
   ! block then gets judged at tol=0 when it went through the eigensolver --
   ! which is where the port and the C are allowed to differ.  That mistake
   ! is what 169,353 "failures" at a worst of 2.4e-14 looked like.
   integer function rys_order(b, s1, s2, s3, s4, inc)
      integer(c_int), intent(in) :: b(BAS_SLOTS, *)
      integer, intent(in) :: s1, s2, s3, s4
      integer, intent(in), optional :: inc
      integer :: extra
      extra = 0
      if (present(inc)) extra = inc
      rys_order = (int(b(ANG_OF, s1+1)) + int(b(ANG_OF, s2+1)) &
                 + int(b(ANG_OF, s3+1)) + int(b(ANG_OF, s4+1)) + extra) / 2 + 1
   end function rys_order

   subroutine cmp(c, f, n, what, ib, sp, ii, jj, kk, ll, inc)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, sp, ii, jj, kk, ll
      character(len=*), intent(in) :: what
      integer, intent(in), optional :: inc
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(bas, ii, jj, kk, ll, inc)
      if (nr <= 5) then
         tol = 0.0_dp
      else
         tol = 1.0e-11_dp
         nloose = nloose + 1
      end if
      ! Scale by the largest integral in the block rather than element by
      ! element: an entry near zero divided by itself measures nothing.
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         r = abs(f(m)-c(m)) / scal
         ncmp = ncmp + 1
         if (r > worst) then
            worst = r
            write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ib, " sph=", sp, " shls (",ii,",",jj,",",kk,",",ll, &
               ") m=", m, " roots=", nr
         end if
         if (r > tol) nbad = nbad + 1
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmp

end program backend_check
