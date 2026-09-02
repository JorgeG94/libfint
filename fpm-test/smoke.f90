! A smoke test for the fpm build: one helium 1s function, whose self-overlap
! is 1 when the contraction coefficient is normalised.  The real verification
! suite is under test/ and is driven by CMake, because it needs libcint to
! compare against.
program smoke
   use libcint_fortran
   implicit none
   integer(ip) :: atm(LIBCINT_ATM_SLOTS, 1), bas(LIBCINT_BAS_SLOTS, 1), shls(2)
   real(dp)    :: env(64), buf(16)
   integer(ip) :: off, ret

   atm = 0; bas = 0; env = 0.0_dp; buf = 0.0_dp
   off = LIBCINT_PTR_ENV_START

   atm(LIBCINT_CHARGE_OF, 1) = 2
   atm(LIBCINT_PTR_COORD, 1) = off
   env(off+1) = 0.0_dp; env(off+2) = 0.0_dp; env(off+3) = 0.0_dp
   off = off + 3

   bas(LIBCINT_ATOM_OF,   1) = 0
   bas(LIBCINT_ANG_OF,    1) = 0
   bas(LIBCINT_NPRIM_OF,  1) = 1
   bas(LIBCINT_NCTR_OF,   1) = 1
   bas(LIBCINT_PTR_EXP,   1) = off
   env(off+1) = 1.2_dp
   off = off + 1
   bas(LIBCINT_PTR_COEFF, 1) = off
   env(off+1) = libcint_gto_norm(bas(LIBCINT_ANG_OF,1), env(bas(LIBCINT_PTR_EXP,1)+1))
   off = off + 1

   shls = [0, 0]
   ret  = libcint_1e_ovlp_cart(buf, shls, atm, 1, bas, 1, env)

   print '(A,F18.14)', " <1s|1s> = ", buf(1)
   if (abs(buf(1) - 1.0_dp) > 1.0e-12_dp) then
      print *, "FAIL: expected 1.0"
      stop 1
   end if
   print *, "OK -- no libcint linked"
end program smoke
