!
! The Breit gauge term, and the two p2 integrals it is built from.
!
! Ported from src/breit.c.  Three things live here, and the first is the only
! reason the file is hand-written rather than generated:
!
!   * CINTgout2e_int2e_breit_r{1,2}p2.  Upstream marks these "based on" a
!     description rather than generated from one, and the distinction is real:
!     feeding those descriptions to the generator gives a nine-component
!     tensor (ng tensor slot 9), while these contract it to one scalar in a
!     single accumulator.  That contraction is not something the generator
!     expresses, so inventing a feature to reproduce one hand-written pair
!     would be the tail wagging the dog.  Transcribed instead, term for term,
!     with the accumulation order kept because nine additions into one
!     running total is not the same number as nine into nine.
!
!   * the gauge assembly, _int2e_breit_drv: the Breit term is
!     [gaunt - gauge_r1 + gauge_r2] / 2, evaluated as the C evaluates it.
!
!   * the four int2e_breit_* entry points the assembly serves.
!
module cint_breit_gauge
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws
   use cint_g2e,       only: cint_init_int2e_envvars
   use cint_g2e_ops
   use cint_2e,        only: cint_2e_drv
   use cint_2e_spinor, only: cint_2e_spinor_drv
   use cint_bas,       only: cint_cgto_spinor
   use cint_gen_breit1
   use cint_gen_gaunt1
      use cint_g2e,       only: cint_all_2e_optimizer
   use cint_opt,       only: cint_del_optimizer
   implicit none
   private

   public :: int2e_breit_r1p2_cart, int2e_breit_r1p2_sph, int2e_breit_r1p2_spinor
   public :: int2e_breit_r2p2_cart, int2e_breit_r2p2_sph, int2e_breit_r2p2_spinor
   public :: int2e_breit_ssp1ssp2_spinor, int2e_breit_ssp1sps2_spinor
   public :: int2e_breit_sps1ssp2_spinor, int2e_breit_sps1sps2_spinor
   public :: int2e_breit_r1p2_optimizer
   public :: int2e_breit_r2p2_optimizer
   public :: int2e_breit_sps1sps2_optimizer
   public :: int2e_breit_sps1ssp2_optimizer
   public :: int2e_breit_ssp1sps2_optimizer
   public :: int2e_breit_ssp1ssp2_optimizer

contains

   ! (NABLA i R0 j|DOT NABLA-R12 |k NABLA l), contracted.
   subroutine CINTgout2e_int2e_breit_r1p2(gout, g, idx, envs, gout_empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: gout_empty
      integer  :: nf, nroots, ix, iy, iz, i, n
      integer  :: g0, g1, g3, g4, g5, g6, g7, g12, g15, gs3
      real(dp) :: s

      nf = envs%nf
      nroots = envs%nrys_roots
      gs3 = envs%g_size * 3
      ! sixteen intermediates, though only these are read
      g0 = 0
      g1 = g0 + gs3;   g3 = g1 + 2*gs3;  g4 = g3 + gs3
      g5 = g4 + gs3;   g6 = g5 + gs3;    g7 = g6 + gs3
      g12 = g7 + 5*gs3
      g15 = g12 + 3*gs3

      call cint_nabla1l_2e(g, g1, g0, envs%i_l+2, envs%j_l+2, envs%k_l+0, envs%l_l+0, envs)
      call cint_x1j_2e(g, g3, g1, envs%rj, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      call cint_nabla1j_2e(g, g4, g0, envs%i_l+1, envs%j_l+1, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g5, g0, envs%i_l+1, envs%j_l+1, envs%k_l, envs%l_l, envs)
      do ix = 0, gs3 - 1
         g(g4+ix) = g(g4+ix) + g(g5+ix)
      end do
      call cint_nabla1j_2e(g, g5, g1, envs%i_l+1, envs%j_l+1, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g6, g1, envs%i_l+1, envs%j_l+1, envs%k_l, envs%l_l, envs)
      do ix = 0, gs3 - 1
         g(g5+ix) = g(g5+ix) + g(g6+ix)
      end do
      call cint_x1j_2e(g, g7, g5, envs%rj, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g12, g4, envs%i_l+0, envs%j_l, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g15, g7, envs%i_l+0, envs%j_l, envs%k_l, envs%l_l, envs)

      do n = 0, nf - 1
         ix = idx(0+n*3); iy = idx(1+n*3); iz = idx(2+n*3)
         s = 0.0_dp
         do i = 0, nroots - 1
            s = s + g(g15+ix+i) * g(g0+iy+i)  * g(g0+iz+i)
            s = s + g(g12+ix+i) * g(g3+iy+i)  * g(g0+iz+i)
            s = s + g(g12+ix+i) * g(g0+iy+i)  * g(g3+iz+i)
            s = s + g(g3+ix+i)  * g(g12+iy+i) * g(g0+iz+i)
            s = s + g(g0+ix+i)  * g(g15+iy+i) * g(g0+iz+i)
            s = s + g(g0+ix+i)  * g(g12+iy+i) * g(g3+iz+i)
            s = s + g(g3+ix+i)  * g(g0+iy+i)  * g(g12+iz+i)
            s = s + g(g0+ix+i)  * g(g3+iy+i)  * g(g12+iz+i)
            s = s + g(g0+ix+i)  * g(g0+iy+i)  * g(g15+iz+i)
         end do
         if (gout_empty /= 0) then
            gout(n) = s
         else
            gout(n) = gout(n) + s
         end if
      end do
   end subroutine CINTgout2e_int2e_breit_r1p2

   ! The same with the origin shift on the ket rather than the bra.
   subroutine CINTgout2e_int2e_breit_r2p2(gout, g, idx, envs, gout_empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: gout_empty
      integer  :: nf, nroots, ix, iy, iz, i, n
      integer  :: g0, g2, g3, g4, g5, g7, g8, g12, g15, gs3
      real(dp) :: s

      nf = envs%nf
      nroots = envs%nrys_roots
      gs3 = envs%g_size * 3
      g0 = 0
      g2 = g0 + 2*gs3; g3 = g2 + gs3;  g4 = g3 + gs3
      g5 = g4 + gs3;   g7 = g5 + 2*gs3; g8 = g7 + gs3
      g12 = g8 + 4*gs3
      g15 = g12 + 3*gs3

      call cint_x1l_2e(g, g2, g0, envs%rl, envs%i_l+2, envs%j_l+1, envs%k_l+0, envs%l_l+1, envs)
      call cint_nabla1l_2e(g, g3, g2, envs%i_l+2, envs%j_l+1, envs%k_l+0, envs%l_l+0, envs)
      call cint_nabla1j_2e(g, g4, g0, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g5, g0, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      do ix = 0, gs3 - 1
         g(g4+ix) = g(g4+ix) + g(g5+ix)
      end do
      call cint_nabla1j_2e(g, g7, g3, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g8, g3, envs%i_l+1, envs%j_l+0, envs%k_l, envs%l_l, envs)
      do ix = 0, gs3 - 1
         g(g7+ix) = g(g7+ix) + g(g8+ix)
      end do
      call cint_nabla1i_2e(g, g12, g4, envs%i_l+0, envs%j_l, envs%k_l, envs%l_l, envs)
      call cint_nabla1i_2e(g, g15, g7, envs%i_l+0, envs%j_l, envs%k_l, envs%l_l, envs)

      do n = 0, nf - 1
         ix = idx(0+n*3); iy = idx(1+n*3); iz = idx(2+n*3)
         s = 0.0_dp
         do i = 0, nroots - 1
            s = s + g(g15+ix+i) * g(g0+iy+i)  * g(g0+iz+i)
            s = s + g(g12+ix+i) * g(g3+iy+i)  * g(g0+iz+i)
            s = s + g(g12+ix+i) * g(g0+iy+i)  * g(g3+iz+i)
            s = s + g(g3+ix+i)  * g(g12+iy+i) * g(g0+iz+i)
            s = s + g(g0+ix+i)  * g(g15+iy+i) * g(g0+iz+i)
            s = s + g(g0+ix+i)  * g(g12+iy+i) * g(g3+iz+i)
            s = s + g(g3+ix+i)  * g(g0+iy+i)  * g(g12+iz+i)
            s = s + g(g0+ix+i)  * g(g3+iy+i)  * g(g12+iz+i)
            s = s + g(g0+ix+i)  * g(g0+iy+i)  * g(g15+iz+i)
         end do
         if (gout_empty /= 0) then
            gout(n) = s
         else
            gout(n) = gout(n) + s
         end if
      end do
   end subroutine CINTgout2e_int2e_breit_r2p2

   ! ---- the r1p2 / r2p2 entry points ----------------------------------

   function int2e_breit_r1p2_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 2, 0, 1, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r1p2
      hv = cint_2e_drv(out, dims, envs, ws, C2S_CART_2E1)
   end function int2e_breit_r1p2_cart

   function int2e_breit_r1p2_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 2, 0, 1, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r1p2
      hv = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_breit_r1p2_sph

   function int2e_breit_r1p2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 2, 0, 1, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r1p2
      hv = cint_2e_spinor_drv(out, dims, envs, ws, .false., .true., .false., .true.)
   end function int2e_breit_r1p2_spinor

   function int2e_breit_r2p2_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 1, 0, 2, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r2p2
      hv = cint_2e_drv(out, dims, envs, ws, C2S_CART_2E1)
   end function int2e_breit_r2p2_cart

   function int2e_breit_r2p2_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 1, 0, 2, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r2p2
      hv = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_breit_r2p2_sph

   function int2e_breit_r2p2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 1, 0, 2, 4, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_breit_r2p2
      hv = cint_2e_spinor_drv(out, dims, envs, ws, .false., .true., .false., .true.)
   end function int2e_breit_r2p2_spinor

   ! ---- the gauge assembly --------------------------------------------
   !
   ! [1/2 gaunt] - [1/2 xxx sigma1.r1] - [-1/2 xxx sigma1.(-r2)], which the C
   ! evaluates as (-(gaunt) - r1 + r2)/2 in two passes over the buffer.  The
   ! sign and the halving are kept in exactly that shape: the alternative,
   ! collecting all three and combining once, rounds differently.
   !
   ! `which` selects which of the four (ssp/sps) x (ssp/sps) combinations to
   ! assemble, so the four entry points below are one line each.
   function breit_assemble(out, dims, shls, atm, natm, bas, nbas, env, ws, which) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas, which
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv

      integer :: counts(0:3), nop, i
      complex(dp), allocatable :: buf(:), buf1(:)
      logical :: h1, h2

      counts(0) = cint_cgto_spinor(shls(0), bas)
      counts(1) = cint_cgto_spinor(shls(1), bas)
      counts(2) = cint_cgto_spinor(shls(2), bas)
      counts(3) = cint_cgto_spinor(shls(3), bas)
      nop = counts(0) * counts(1) * counts(2) * counts(3)
      allocate(buf(0:nop-1), buf1(0:nop-1))
      buf = (0.0_dp, 0.0_dp); buf1 = (0.0_dp, 0.0_dp)

      hv = gaunt_of(which, buf1, counts, shls, atm, natm, bas, nbas, env, ws)
      h1 = gauge_r1_of(which, buf, counts, shls, atm, natm, bas, nbas, env, ws)
      hv = h1 .or. hv
      if (hv) then
         do i = 0, nop - 1
            buf1(i) = -buf1(i) - buf(i)
         end do
      end if
      h2 = gauge_r2_of(which, buf, counts, shls, atm, natm, bas, nbas, env, ws)
      hv = h2 .or. hv
      if (hv) then
         do i = 0, nop - 1
            buf1(i) = (buf1(i) + buf(i)) * 0.5_dp
         end do
      end if

      call copy_to_out(out, buf1, dims, counts)
   end function breit_assemble

   ! The three pieces, dispatched.  Written as three small selects rather
   ! than a table of procedure pointers: the four combinations are fixed and
   ! a pointer table would need four abstract interfaces to say the same.
   function gaunt_of(which, b, counts, shls, atm, natm, bas, nbas, env, ws) result(hv)
      integer,  intent(in) :: which, counts(0:3), shls(0:), natm, nbas
      complex(dp), intent(inout) :: b(0:)
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      select case (which)
      case (1); hv = int2e_ssp1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (2); hv = int2e_ssp1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (3); hv = int2e_sps1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case default
                hv = int2e_sps1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      end select
   end function gaunt_of

   function gauge_r1_of(which, b, counts, shls, atm, natm, bas, nbas, env, ws) result(hv)
      integer,  intent(in) :: which, counts(0:3), shls(0:), natm, nbas
      complex(dp), intent(inout) :: b(0:)
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      select case (which)
      case (1); hv = int2e_gauge_r1_ssp1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (2); hv = int2e_gauge_r1_ssp1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (3); hv = int2e_gauge_r1_sps1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case default
                hv = int2e_gauge_r1_sps1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      end select
   end function gauge_r1_of

   function gauge_r2_of(which, b, counts, shls, atm, natm, bas, nbas, env, ws) result(hv)
      integer,  intent(in) :: which, counts(0:3), shls(0:), natm, nbas
      complex(dp), intent(inout) :: b(0:)
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      select case (which)
      case (1); hv = int2e_gauge_r2_ssp1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (2); hv = int2e_gauge_r2_ssp1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case (3); hv = int2e_gauge_r2_sps1ssp2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      case default
                hv = int2e_gauge_r2_sps1sps2_spinor(b, counts, shls, atm, natm, bas, nbas, env, ws)
      end select
   end function gauge_r2_of

   ! The assembled block is packed at `counts`; the caller may want it at
   ! `dims`, which can be larger.
   pure subroutine copy_to_out(out, inb, dims, counts)
      complex(dp), intent(inout) :: out(0:)
      complex(dp), intent(in)    :: inb(0:)
      integer,     intent(in)    :: dims(0:), counts(0:3)
      integer :: i, j, k, l, ni, nj, nk, nij, nijk, di, dj, dk, dl, dij, dijk
      ni = dims(0); nj = dims(1); nk = dims(2)
      nij = ni*nj; nijk = nij*nk
      di = counts(0); dj = counts(1); dk = counts(2); dl = counts(3)
      dij = di*dj; dijk = dij*dk
      do l = 0, dl - 1
         do k = 0, dk - 1
            do j = 0, dj - 1
            do i = 0, di - 1
               out(l*nijk + k*nij + j*ni + i) = inb(l*dijk + k*dij + j*di + i)
            end do
            end do
         end do
      end do
   end subroutine copy_to_out

   ! ---- entry points ---------------------------------------------------

   function int2e_breit_ssp1ssp2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      hv = breit_assemble(out, dims, shls, atm, natm, bas, nbas, env, ws, 1)
   end function int2e_breit_ssp1ssp2_spinor

   function int2e_breit_ssp1sps2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      hv = breit_assemble(out, dims, shls, atm, natm, bas, nbas, env, ws, 2)
   end function int2e_breit_ssp1sps2_spinor

   function int2e_breit_sps1ssp2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      hv = breit_assemble(out, dims, shls, atm, natm, bas, nbas, env, ws, 3)
   end function int2e_breit_sps1ssp2_spinor

   function int2e_breit_sps1sps2_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      hv = breit_assemble(out, dims, shls, atm, natm, bas, nbas, env, ws, 4)
   end function int2e_breit_sps1sps2_spinor

   ! ---- optimizers ------------------------------------------------------

   subroutine int2e_breit_r1p2_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [2, 2, 0, 1, 4, 1, 1, 1]
      call cint_all_2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_breit_r1p2_optimizer

   subroutine int2e_breit_r2p2_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [2, 1, 0, 2, 4, 1, 1, 1]
      call cint_all_2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_breit_r2p2_optimizer

   subroutine int2e_breit_sps1sps2_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int2e_breit_sps1sps2_optimizer

   subroutine int2e_breit_sps1ssp2_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int2e_breit_sps1ssp2_optimizer

   subroutine int2e_breit_ssp1sps2_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int2e_breit_ssp1sps2_optimizer

   subroutine int2e_breit_ssp1ssp2_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int2e_breit_ssp1ssp2_optimizer


end module cint_breit_gauge
