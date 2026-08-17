!
! The one-electron spinor driver, and the integrals that sit on it.
!
! Ported from CINT1e_spinor_drv (src/cint1e.c:239) and the four c2s_s[fi]_1e[i]
! transforms (src/cart2sph.c:4869-5040).
!
! The primitive loop underneath is the same cint_1e_loop the real drivers use:
! a spinor integral is a Cartesian integral with a different transform on the
! way out, and nothing below the transform knows the difference.  What changes
! is the output type, which is complex, and the dimensions, which count spinor
! components rather than Cartesian or spherical ones.
!
! Four transforms, two axes:
!
!   sf / si   spin free or spin included.  The spin-included ones read four
!             Cartesian blocks -- the three sigma components and the scalar --
!             where the spin-free ones read only the scalar.
!   plain / i whether the ket transform carries a factor of i, which the
!             operators with an odd number of momentum factors need.
!
! A separate module from cint_1e rather than more arms on cint_1e_drv: the
! driver's return type differs, and Fortran has no way to make that a runtime
! choice without either a union or two entry points.  Two entry points is the
! honest one.
!
module cint_1e_spinor
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d, ws_mark, ws_rewind
   use cint_g1e,       only: cint_init_int1e_envvars
   use cint_1e,        only: cint_1e_loop, cint_gout1e, int1e_cache_size
   use cint_bas,       only: cint_cgto_spinor, KAPPA_OF
   use cint_cart2spinor
   implicit none
   private

   public :: cint_1e_spinor_drv
   public :: int1e_ovlp_spinor, int1e_nuc_spinor
   ! int2c2e_spinor wants exactly this transform.  Its envs already substitutes
   ! the k shell into the j slots -- j_l, nfj, x_ctr(1), shls(1) -- so c2s_sf_1e
   ! applies to a two-centre integral unchanged, which is why the C reaches for
   ! the same function rather than writing a two-centre one.
   public :: apply_c2s_spinor_1e, zcopy_ij

contains

   ! Merge a split real/imaginary block into the caller's complex array.
   pure subroutine zcopy_ij(out, ooff, gR, gI, ni, mi, mj)
      complex(dp), intent(inout) :: out(0:)
      real(dp),    intent(in)    :: gR(0:), gI(0:)
      integer,     intent(in)    :: ooff, ni, mi, mj
      integer :: i, j
      do j = 0, mj - 1
      do i = 0, mi - 1
         out(ooff + j*ni + i) = cmplx(gR(j*mi+i), gI(j*mi+i), dp)
      end do
      end do
   end subroutine zcopy_ij

   ! The four transforms.  They differ only in which bra kernel runs and
   ! whether the ket carries the i, so they share one body.
   subroutine apply_c2s_spinor_1e(opij, gctr, goff, dims, envs, ws, spin_included, with_i)
      complex(dp), intent(inout) :: opij(0:)
      real(dp),    intent(in)    :: gctr(0:)
      integer,     intent(in)    :: goff, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      type(cint_ws), intent(inout) :: ws
      logical, intent(in) :: spin_included, with_i

      integer :: i_kp, j_kp, i_ctr, j_ctr, di, dj, ni, ofj, nfj, nf2j, nf
      integer :: ic, jc, o1r, o1i, o2r, o2i, gb, blk

      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      ni = dims(0)
      ofj = ni * dj
      nfj = envs%nfj
      nf2j = nfj + nfj
      nf = envs%nf

      call ws_alloc_d(ws, di*nf2j, o1r)
      call ws_alloc_d(ws, di*nf2j, o1i)
      call ws_alloc_d(ws, di*dj,   o2r)
      call ws_alloc_d(ws, di*dj,   o2i)

      ! The spin-included transform reads four blocks laid end to end, each
      ! nf*i_ctr*j_ctr long, in the order x, y, z, scalar.
      blk = nf * i_ctr * j_ctr
      gb = goff
      do jc = 0, j_ctr - 1
      do ic = 0, i_ctr - 1
         if (spin_included) then
            call a_bra_cart2spinor_si(ws%d(o1r:), ws%d(o1i:), &
                                      gctr(gb:), gctr(gb+blk:), gctr(gb+2*blk:), &
                                      gctr(gb+3*blk:), nfj, i_kp, envs%i_l)
         else
            call a_bra_cart2spinor_sf(ws%d(o1r:), ws%d(o1i:), gctr(gb:), &
                                      nfj, i_kp, envs%i_l)
         end if
         if (with_i) then
            call a_iket_cart2spinor(ws%d(o2r:), ws%d(o2i:), ws%d(o1r:), ws%d(o1i:), &
                                    di, j_kp, envs%j_l)
         else
            call a_ket_cart2spinor(ws%d(o2r:), ws%d(o2i:), ws%d(o1r:), ws%d(o1i:), &
                                   di, j_kp, envs%j_l)
         end if
         call zcopy_ij(opij, ofj*jc + di*ic, ws%d(o2r:), ws%d(o2i:), ni, di, dj)
         gb = gb + nf
      end do
      end do
   end subroutine apply_c2s_spinor_1e

   function cint_1e_spinor_drv(out, dims, envs, ws, spin_included, with_i, int1e_type) &
         result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      logical,     intent(in)    :: spin_included, with_i
      integer,     intent(in)    :: int1e_type
      logical :: has_value
      integer :: dmark, imark

      integer :: nc, ogctr, n, nout, counts(0:3), nd, ni

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%ncomp_e1

      call int1e_cache_size(envs, nd, ni)
      ! the four split scratch blocks the spinor transform wants, which the
      ! real driver's estimate has no reason to carry
      nd = nd + 4 * envs%nfi * envs%nfj * 4 * (envs%i_l*2 + 2) + 64
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*envs%ncomp_tensor, ogctr)
      ws%d(ogctr:ogctr+nc*envs%ncomp_tensor-1) = 0.0_dp

      has_value = cint_1e_loop(ws%d, ogctr, envs, ws, int1e_type)

      counts(0) = cint_cgto_spinor(envs%shls(0), envs%bas)
      counts(1) = cint_cgto_spinor(envs%shls(1), envs%bas)
      counts(2) = 1
      counts(3) = 1
      nout = dims(0) * dims(1)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, envs%ncomp_tensor - 1
            call ws_rewind(ws, dmark, imark)
            call apply_c2s_spinor_1e(out(nout*n:), ws%d, ogctr + nc*n, dims, envs, ws, &
                                     spin_included, with_i)
         end do
      else
         do n = 0, envs%ncomp_tensor - 1
            call c2s_zset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_1e_spinor_drv

   ! ---- entry points -------------------------------------------------

   function int1e_ovlp_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e
      has_value = cint_1e_spinor_drv(out, dims, envs, ws, .false., .false., INT1E_OVLP)
   end function int1e_ovlp_spinor

   function int1e_nuc_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_nuc_spinor
      has_value = cint_1e_spinor_drv(out, dims, envs, ws, .false., .false., INT1E_NUC)
   end function int1e_nuc_spinor

   ! int1e_nuc's gout sums over the Rys roots; the spinor path needs the same
   ! kernel the real one uses, re-exported here because cint_1e keeps it
   ! private to its own driver.
   subroutine cint_gout1e_nuc_spinor(gout, g, idx, envs, empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: empty
      integer  :: n, i, ix, iy, iz
      real(dp) :: s
      if (empty /= 0) then
         do n = 0, envs%nf - 1
            ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
            s = 0.0_dp
            do i = 0, envs%nrys_roots - 1
               s = s + g(ix+i) * g(iy+i) * g(iz+i)
            end do
            gout(n) = s
         end do
      else
         do n = 0, envs%nf - 1
            ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
            s = 0.0_dp
            do i = 0, envs%nrys_roots - 1
               s = s + g(ix+i) * g(iy+i) * g(iz+i)
            end do
            gout(n) = gout(n) + s
         end do
      end if
   end subroutine cint_gout1e_nuc_spinor

end module cint_1e_spinor
