!
! The one-electron G tensor: setting up the environment, the index map from
! Cartesian components into it, and the overlap recursion itself.
!
! Ported from src/g1e.c.  The g array is a flat (nrys_roots, i, j) block per
! Cartesian direction, three of them end to end, and every stride is kept as
! the C computes it -- see PORT_TO_FORTRAN.md 3.6 for why the index arithmetic
! is copied rather than tidied.
!
module cint_g1e
   use cint_const,     only: dp
   use cint_envs
   use cint_bas,       only: cint_cart_comp, cint_square_dist, &
                             ANG_OF, NCTR_OF, ATOM_OF, PTR_COORD, &
                             CHARGE_OF, NUC_MOD_OF, PTR_ZETA, PTR_FRAC_CHARGE, &
                             ATM_SLOTS, GAUSSIAN_NUC, FRAC_CHARGE_NUC
   use cint_rys_roots, only: cint_rys_roots_lr
   use cint_opt,       only: cint_del_optimizer, opt_set_log_maxc, &
                             opt_set_non0coeff, opt_setij, opt_gen_idx, opt_finish
   implicit none
   private

   public :: cint_init_int1e_envvars, cint_g1e_index_xyz, cint_g1e_ovlp
   public :: cint_g1e_nuc, cint_nuc_mod
   public :: cint_common_fac_sp, cint_prim_to_ctr_0, cint_prim_to_ctr_1
   ! The G1E_* operations the generated gout kernels reach through.  In the C
   ! these hide behind macros in g1e.h; the macros split into two kinds --
   ! calls like G1E_D_J, and offset assignments like G1E_R_I -- which is the
   ! distinction the D1 spike found and the generator now dispatches on.
   public :: cint_nabla1i_1e, cint_nabla1j_1e, cint_nabla1k_1e
   public :: cint_x1i_1e, cint_x1j_1e, cint_x1k_1e
   public :: cint_all_1e_optimizer

   real(dp), parameter :: SQRTPI = 1.7724538509055160272981674833411451_dp
   real(dp), parameter :: PI     = 3.1415926535897932384626433832795029_dp

contains

   subroutine cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      ! no INTENT(IN) on the three tables: envs points at them, and F2018
      ! 8.5.10 forbids an INTENT(IN) dummy as a pointer target.  They are
      ! never written -- see the note on cint_env_vars.
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)

      integer :: i_sh, j_sh, dli, dlj, ip
      logical :: ibase

      envs%natm = natm
      envs%nbas = nbas
      envs%atm(0:) => atm
      envs%bas(0:) => bas
      envs%env(0:) => env
      envs%shls(0:size(shls)-1) = shls
      envs%ng(0:min(7,size(ng)-1)) = ng(0:min(7,size(ng)-1))
      envs%omega = env(PTR_RANGE_OMEGA)

      i_sh = shls(0)
      j_sh = shls(1)
      envs%i_l = bas_of(envs, ANG_OF, i_sh)
      envs%j_l = bas_of(envs, ANG_OF, j_sh)
      envs%x_ctr(0) = bas_of(envs, NCTR_OF, i_sh)
      envs%x_ctr(1) = bas_of(envs, NCTR_OF, j_sh)
      envs%nfi = (envs%i_l+1)*(envs%i_l+2)/2
      envs%nfj = (envs%j_l+1)*(envs%j_l+2)/2
      envs%nf = envs%nfi * envs%nfj
      envs%common_factor = 1.0_dp
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = envs%j_l + ng(JINC)

      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, i_sh))
      envs%ri(0) = env(ip); envs%ri(1) = env(ip+1); envs%ri(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, j_sh))
      envs%rj(0) = env(ip); envs%rj(1) = env(ip+1); envs%rj(2) = env(ip+2)

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_tensor = ng(TENSOR)
      if (ng(SLOT_RYS_ROOTS) > 0) then
         envs%nrys_roots = ng(SLOT_RYS_ROOTS)
      else
         envs%nrys_roots = (envs%li_ceil + envs%lj_ceil)/2 + 1
      end if

      ! The recursion is built along whichever centre carries more angular
      ! momentum, so the strides and the sign of rirj both follow from that.
      ibase = envs%li_ceil > envs%lj_ceil
      if (ibase) then
         dli = envs%li_ceil + envs%lj_ceil + 1
         dlj = envs%lj_ceil + 1
         envs%rirj = envs%ri - envs%rj
      else
         dli = envs%li_ceil + 1
         dlj = envs%li_ceil + envs%lj_ceil + 1
         envs%rirj = envs%rj - envs%ri
      end if
      envs%g_stride_i = envs%nrys_roots
      envs%g_stride_j = envs%nrys_roots * dli
      envs%g_size     = envs%nrys_roots * dli * dlj
      envs%g_stride_k = envs%g_size
      envs%g_stride_l = envs%g_size
   end subroutine cint_init_int1e_envvars

   ! For each (i,j) Cartesian component pair, the three offsets into g at
   ! which its x, y and z factors live.
   subroutine cint_g1e_index_xyz(idx, envs)
      integer, intent(out) :: idx(0:)
      type(cint_env_vars), intent(in) :: envs

      integer :: i, j, n, di, dj, ofx, ofy, ofz, ofjx, ofjy, ofjz
      integer :: i_nx(0:CART_MAX-1), i_ny(0:CART_MAX-1), i_nz(0:CART_MAX-1)
      integer :: j_nx(0:CART_MAX-1), j_ny(0:CART_MAX-1), j_nz(0:CART_MAX-1)

      di = envs%g_stride_i
      dj = envs%g_stride_j
      call cint_cart_comp(i_nx, i_ny, i_nz, envs%i_l)
      call cint_cart_comp(j_nx, j_ny, j_nz, envs%j_l)

      ofx = 0
      ofy = envs%g_size
      ofz = envs%g_size * 2
      n = 0
      do j = 0, envs%nfj - 1
         ofjx = ofx + dj * j_nx(j)
         ofjy = ofy + dj * j_ny(j)
         ofjz = ofz + dj * j_nz(j)
         do i = 0, envs%nfi - 1
            idx(n+0) = ofjx + di * i_nx(i)
            idx(n+1) = ofjy + di * i_ny(i)
            idx(n+2) = ofjz + di * i_nz(i)
            n = n + 3
         end do
      end do
   end subroutine cint_g1e_index_xyz

   ! Overlap G tensor by the Obara-Saika recursion: build up in i along the
   ! dominant centre, then transfer to j.
   function cint_g1e_ovlp(g, goff, envs) result(has_value)
      ! TARGET for the three block pointers -- see cint_g1e_nuc.
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: goff
      type(cint_env_vars), intent(in) :: envs
      integer :: has_value

      integer  :: gx, gy, gz, nmax, lj, di, dj, i, j, n, ptr
      real(dp) :: aij, aij2, rijrx(0:2), rx(0:2)
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      real(dp) :: rirj0, rirj1, rirj2

      gx = goff
      gy = goff + envs%g_size
      gz = goff + envs%g_size * 2
      aij = envs%ai + envs%aj
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      gxp(0) = 1.0_dp
      gyp(0) = 1.0_dp
      gzp(0) = envs%fac * SQRTPI * PI / (aij * sqrt(aij))

      nmax = envs%li_ceil + envs%lj_ceil
      has_value = 1
      if (nmax == 0) return

      if (envs%li_ceil > envs%lj_ceil) then
         lj = envs%lj_ceil
         di = envs%g_stride_i
         dj = envs%g_stride_j
         rx = envs%ri
      else
         lj = envs%li_ceil
         di = envs%g_stride_j
         dj = envs%g_stride_i
         rx = envs%rj
      end if
      rijrx = envs%rij - rx
      rirj0 = envs%rirj(0); rirj1 = envs%rirj(1); rirj2 = envs%rirj(2)

      gxp(di) = rijrx(0) * gxp(0)
      gyp(di) = rijrx(1) * gyp(0)
      gzp(di) = rijrx(2) * gzp(0)

      aij2 = 0.5_dp / aij
      do i = 1, nmax - 1
         gxp((i+1)*di) = i * aij2 * gxp((i-1)*di) + rijrx(0) * gxp(i*di)
         gyp((i+1)*di) = i * aij2 * gyp((i-1)*di) + rijrx(1) * gyp(i*di)
         gzp((i+1)*di) = i * aij2 * gzp((i-1)*di) + rijrx(2) * gzp(i*di)
      end do

      do j = 1, lj
         ptr = dj * j
         n = ptr
         do i = 0, nmax - j
            gxp(n) = gxp(n+di-dj) + rirj0 * gxp(n-dj)
            gyp(n) = gyp(n+di-dj) + rirj1 * gyp(n-dj)
            gzp(n) = gzp(n+di-dj) + rirj2 * gzp(n-dj)
            n = n + di
         end do
      end do
   end function cint_g1e_ovlp

   ! The s and p normalisation the cart2sph tables leave out, so that those
   ! transforms can be pure copies for l <= 1.
   pure function cint_common_fac_sp(l) result(f)
      integer, intent(in) :: l
      real(dp) :: f
      select case (l)
      case (0); f = 0.282094791773878143_dp
      case (1); f = 0.488602511902919921_dp
      case default; f = 1.0_dp
      end select
   end function cint_common_fac_sp

   ! Primitive -> contracted, overwriting.  sortedidx and non0ctr are unused
   ! here, exactly as in the C: the zero-coefficient contractions still have
   ! to be written, because this is the call that initialises the block.
   ! buf holds both the contracted target and the primitive source: they are
   ! always the same workspace, so one array and two offsets rather than two
   ! dummies aliased onto it, which Fortran does not allow.
   pure subroutine cint_prim_to_ctr_0(buf, gcoff, gpoff, coeff, coff, nf, nprim, nctr)
      real(dp), intent(inout) :: buf(0:)
      real(dp), intent(in)    :: coeff(0:)
      integer,  intent(in)    :: gcoff, gpoff, coff, nf, nprim, nctr
      integer  :: i, n
      real(dp) :: c0
      do i = 0, nctr - 1
         c0 = coeff(coff + nprim*i)
         do n = 0, nf - 1
            buf(gcoff + nf*i + n) = c0 * buf(gpoff + n)
         end do
      end do
   end subroutine cint_prim_to_ctr_0

   ! Primitive -> contracted, accumulating.  Only the contractions with a
   ! non-zero coefficient are touched, which is what sortedidx is for.
   pure subroutine cint_prim_to_ctr_1(buf, gcoff, gpoff, coeff, coff, nf, &
                                      nprim, nctr, non0ctr, sortedidx, sioff)
      real(dp), intent(inout) :: buf(0:)
      real(dp), intent(in)    :: coeff(0:)
      integer,  intent(in)    :: sortedidx(0:)
      integer,  intent(in)    :: gcoff, gpoff, coff, nf, nprim, nctr, non0ctr, sioff
      integer  :: i, j, n, gcb
      real(dp) :: c0
      do i = 0, non0ctr - 1
         j = sortedidx(sioff + i)
         c0 = coeff(coff + nprim*j)
         ! The row base out of the inner loop.  The C carries gc and gp as
         ! two pointers, so its body is base+n on both sides; leaving
         ! gcoff + nf*j inside cost 73.4 million instructions against the
         ! C's 50.8 million for the same axpy.
         gcb = gcoff + nf*j
         do n = 0, nf - 1
            buf(gcb + n) = buf(gcb + n) + c0 * buf(gpoff + n)
         end do
      end do
   end subroutine cint_prim_to_ctr_1

   ! d/dx on the bra: f(i) = i*g(i-1) - 2*ai*g(i+1)
   subroutine cint_nabla1i_1e(g, foff, goff, li, lj, lk, envs)
      ! ONE array, two offsets.  The C writes these as f and g, distinct
      ! pointers into the same buffer -- which is fine there, and is illegal
      ! in Fortran: associating the same actual argument with both an
      ! intent(in) and an intent(inout) dummy lets the compiler assume they
      ! cannot alias, and gfortran duly miscompiled int1e_kin.  Saying it is
      ! one array is both true and safe.
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      real(dp) :: ai2
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      ai2 = -2.0_dp * envs%ai
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do k = 0, lk
         do j = 0, lj
            ptr = dj*j + dk*k
            g(fx+ptr) = ai2 * g(gx+ptr+1)
            g(fy+ptr) = ai2 * g(gy+ptr+1)
            g(fz+ptr) = ai2 * g(gz+ptr+1)
            do i = 1, li
               g(fx+ptr+i) = i * g(gx+ptr+i-1) + ai2 * g(gx+ptr+i+1)
               g(fy+ptr+i) = i * g(gy+ptr+i-1) + ai2 * g(gy+ptr+i+1)
               g(fz+ptr+i) = i * g(gz+ptr+i-1) + ai2 * g(gz+ptr+i+1)
            end do
         end do
      end do
   end subroutine cint_nabla1i_1e

   ! d/dx on the ket
   subroutine cint_nabla1j_1e(g, foff, goff, li, lj, lk, envs)
      ! ONE array, two offsets.  The C writes these as f and g, distinct
      ! pointers into the same buffer -- which is fine there, and is illegal
      ! in Fortran: associating the same actual argument with both an
      ! intent(in) and an intent(inout) dummy lets the compiler assume they
      ! cannot alias, and gfortran duly miscompiled int1e_kin.  Saying it is
      ! one array is both true and safe.
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      real(dp) :: aj2
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      aj2 = -2.0_dp * envs%aj
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do k = 0, lk
         ptr = dk*k
         do i = ptr, ptr + li
            g(fx+i) = aj2 * g(gx+i+dj)
            g(fy+i) = aj2 * g(gy+i+dj)
            g(fz+i) = aj2 * g(gz+i+dj)
         end do
         do j = 1, lj
            ptr = dj*j + dk*k
            do i = ptr, ptr + li
               g(fx+i) = j * g(gx+i-dj) + aj2 * g(gx+i+dj)
               g(fy+i) = j * g(gy+i-dj) + aj2 * g(gy+i+dj)
               g(fz+i) = j * g(gz+i-dj) + aj2 * g(gz+i+dj)
            end do
         end do
      end do
   end subroutine cint_nabla1j_1e

   ! d/dx on the third centre (3c1e)
   subroutine cint_nabla1k_1e(g, foff, goff, li, lj, lk, envs)
      ! ONE array, two offsets.  The C writes these as f and g, distinct
      ! pointers into the same buffer -- which is fine there, and is illegal
      ! in Fortran: associating the same actual argument with both an
      ! intent(in) and an intent(inout) dummy lets the compiler assume they
      ! cannot alias, and gfortran duly miscompiled int1e_kin.  Saying it is
      ! one array is both true and safe.
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      real(dp) :: ak2
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      ak2 = -2.0_dp * envs%ak
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do j = 0, lj
         ptr = dj*j
         do i = ptr, ptr + li
            g(fx+i) = ak2 * g(gx+i+dk)
            g(fy+i) = ak2 * g(gy+i+dk)
            g(fz+i) = ak2 * g(gz+i+dk)
         end do
      end do
      do k = 1, lk
         do j = 0, lj
            ptr = dj*j + dk*k
            do i = ptr, ptr + li
               g(fx+i) = k * g(gx+i-dk) + ak2 * g(gx+i+dk)
               g(fy+i) = k * g(gy+i-dk) + ak2 * g(gy+i+dk)
               g(fz+i) = k * g(gz+i-dk) + ak2 * g(gz+i+dk)
            end do
         end do
      end do
   end subroutine cint_nabla1k_1e

   ! Shift the origin of the bra index: f(i) = g(i+1) + r*g(i)
   subroutine cint_x1i_1e(g, foff, goff, r, li, lj, lk, envs)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: r(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do k = 0, lk
         do j = 0, lj
            ptr = dj*j + dk*k
            do i = ptr, ptr + li
               g(fx+i) = g(gx+i+1) + r(0) * g(gx+i)
               g(fy+i) = g(gy+i+1) + r(1) * g(gy+i)
               g(fz+i) = g(gz+i+1) + r(2) * g(gz+i)
            end do
         end do
      end do
   end subroutine cint_x1i_1e

   subroutine cint_x1j_1e(g, foff, goff, r, li, lj, lk, envs)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: r(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do k = 0, lk
         do j = 0, lj
            ptr = dj*j + dk*k
            do i = ptr, ptr + li
               g(fx+i) = g(gx+i+dj) + r(0) * g(gx+i)
               g(fy+i) = g(gy+i+dj) + r(1) * g(gy+i)
               g(fz+i) = g(gz+i+dj) + r(2) * g(gz+i)
            end do
         end do
      end do
   end subroutine cint_x1j_1e

   subroutine cint_x1k_1e(g, foff, goff, r, li, lj, lk, envs)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: r(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk
      type(cint_env_vars), intent(in) :: envs
      integer :: i, j, k, ptr, dj, dk, gs, fx, fy, fz, gx, gy, gz
      dj = envs%g_stride_j; dk = envs%g_stride_k; gs = envs%g_size
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2
      do k = 0, lk
         do j = 0, lj
            ptr = dj*j + dk*k
            do i = ptr, ptr + li
               g(fx+i) = g(gx+i+dk) + r(0) * g(gx+i)
               g(fy+i) = g(gy+i+dk) + r(1) * g(gy+i)
               g(fz+i) = g(gz+i+dk) + r(2) * g(gz+i)
            end do
         end do
      end do
   end subroutine cint_x1k_1e

   ! The width of the nuclear charge distribution, as a factor on the Rys
   ! argument.  Zero zeta means a point nucleus and the factor is 1.
   pure function cint_nuc_mod(aij, nuc_id, atm, env) result(tau)
      real(dp), intent(in) :: aij
      integer,  intent(in) :: nuc_id, atm(0:)
      real(dp), intent(in) :: env(0:)
      real(dp) :: tau, zeta
      if (nuc_id < 0) then
         zeta = env(PTR_RINV_ZETA)
      else if (atm(ATM_SLOTS*nuc_id + NUC_MOD_OF) == GAUSSIAN_NUC) then
         zeta = env(atm(ATM_SLOTS*nuc_id + PTR_ZETA))
      else
         zeta = 0.0_dp
      end if
      if (zeta > 0.0_dp) then
         tau = sqrt(zeta / (aij + zeta))
      else
         tau = 1.0_dp
      end if
   end function cint_nuc_mod

   ! Nuclear attraction (or 1/|r-R| with the origin from PTR_RINV_ORIG when
   ! nuc_id < 0).  Unlike the overlap this needs Rys quadrature, because the
   ! Coulomb operator is not separable in x, y, z.
   function cint_g1e_nuc(g, goff, envs, nuc_id) result(has_value)
      ! TARGET for the three block pointers below.
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: goff, nuc_id
      type(cint_env_vars), intent(in) :: envs
      integer :: has_value

      integer  :: nrys_roots, i, j, n, nmax, lj, di, dj, gx, gy, gz, cr, err
      real(dp) :: u(0:MXRYSROOTS-1)
      real(dp) :: crij(0:2), rx(0:2), aij, tau, x, fac1
      real(dp) :: rijrx, rijry, rijrz, aij2, ru, rt, r0, r1, r2
      ! One pointer per coordinate block, as the C keeps `double *gz = g +
      ! envs->g_size*2` -- the recursion below is three writes and six reads
      ! per iteration, and every one of them was paying gx+ or gy+ or gz+.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      real(dp) :: rirj0, rirj1, rirj2

      nrys_roots = envs%nrys_roots
      gx = goff; gy = goff + envs%g_size; gz = goff + envs%g_size*2
      aij = envs%ai + envs%aj
      tau = cint_nuc_mod(aij, nuc_id, envs%atm, envs%env)

      if (nuc_id < 0) then
         fac1 = 2.0_dp*PI * envs%fac * tau / aij
         crij = envs%env(PTR_RINV_ORIG:PTR_RINV_ORIG+2) - envs%rij
      else if (envs%atm(ATM_SLOTS*nuc_id + NUC_MOD_OF) == FRAC_CHARGE_NUC) then
         fac1 = 2.0_dp*PI * (-envs%env(envs%atm(ATM_SLOTS*nuc_id + PTR_FRAC_CHARGE))) &
              * envs%fac * tau / aij
         cr = envs%atm(ATM_SLOTS*nuc_id + PTR_COORD)
         crij = envs%env(cr:cr+2) - envs%rij
      else
         fac1 = 2.0_dp*PI * (-real(abs(envs%atm(ATM_SLOTS*nuc_id + CHARGE_OF)), dp)) &
              * envs%fac * tau / aij
         cr = envs%atm(ATM_SLOTS*nuc_id + PTR_COORD)
         crij = envs%env(cr:cr+2) - envs%rij
      end if

      x = aij * tau * tau * sum(crij*crij)
      err = cint_rys_roots_lr(nrys_roots, x, u, g(gz:))
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)
      do i = 0, nrys_roots - 1
         gxp(i) = 1.0_dp
         gyp(i) = 1.0_dp
         gzp(i) = gzp(i) * fac1
      end do

      has_value = 1
      nmax = envs%li_ceil + envs%lj_ceil
      if (nmax == 0) return

      if (envs%li_ceil > envs%lj_ceil) then
         lj = envs%lj_ceil; di = envs%g_stride_i; dj = envs%g_stride_j; rx = envs%ri
      else
         lj = envs%li_ceil; di = envs%g_stride_j; dj = envs%g_stride_i; rx = envs%rj
      end if
      rijrx = envs%rij(0) - rx(0)
      rijry = envs%rij(1) - rx(1)
      rijrz = envs%rij(2) - rx(2)
      aij2 = 0.5_dp / aij
      rirj0 = envs%rirj(0); rirj1 = envs%rirj(1); rirj2 = envs%rirj(2)

      do n = 0, nrys_roots - 1
         ru = tau * tau * u(n) / (1.0_dp + u(n))
         rt = aij2 - aij2 * ru
         r0 = rijrx + ru * crij(0)
         r1 = rijry + ru * crij(1)
         r2 = rijrz + ru * crij(2)
         gxp(di+n) = r0 * gxp(n)
         gyp(di+n) = r1 * gyp(n)
         gzp(di+n) = r2 * gzp(n)
         do i = 1, nmax - 1
            gxp(di+n+i*di) = i * rt * gxp(-di+n+i*di) + r0 * gxp(n+i*di)
            gyp(di+n+i*di) = i * rt * gyp(-di+n+i*di) + r1 * gyp(n+i*di)
            gzp(di+n+i*di) = i * rt * gzp(-di+n+i*di) + r2 * gzp(n+i*di)
         end do
      end do

      do j = 1, lj
         do i = 0, nmax - j
            do n = 0, nrys_roots - 1
               gxp(j*dj+n+i*di) = gxp((j-1)*dj+di+n+i*di) &
                                + rirj0 * gxp((j-1)*dj+n+i*di)
               gyp(j*dj+n+i*di) = gyp((j-1)*dj+di+n+i*di) &
                                + rirj1 * gyp((j-1)*dj+n+i*di)
               gzp(j*dj+n+i*di) = gzp((j-1)*dj+di+n+i*di) &
                                + rirj2 * gzp((j-1)*dj+n+i*di)
            end do
         end do
      end do
   end function cint_g1e_nuc

   ! ---- the optimizer builder for this arity ----------------------------
   !
   ! It lives here rather than in cint_opt because it names this module's
   ! envs init, and cint_opt has to stay underneath every arity.

   subroutine cint_all_1e_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_set_log_maxc(opt, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int1e_envvars, cint_g1e_index_xyz, 2, ANG_MAX, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 2, ng, nbas)
   end subroutine cint_all_1e_optimizer


end module cint_g1e
