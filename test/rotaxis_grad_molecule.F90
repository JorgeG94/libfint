!
! The rotated-axis gradient on a real molecule, validated by finite
! difference of the energy it is the derivative of.
!
! WHY THIS EXISTS ALONGSIDE rotaxis_grad_check.  That check compares the
! path against libfint's own int2e_ip1, which catches an error in the
! expansion.  It cannot catch an error that both paths share, and it cannot
! see the three things a caller actually has to get right: the sign, the
! layout of the three derivative blocks, and which atom a block belongs to.
! A finite difference of the energy sees all of them, because it never
! touches a derivative integral at all.
!
! The quantity is the two-electron energy of a fixed density,
!
!     E2 = 1/2 sum_{ijkl} D_ij D_kl (ij|kl)
!
! with D held CONSTANT as the nuclei move -- so dE2/dX is exactly the
! integral derivative and nothing else, which is what is under test.  Each
! quartet contributes to all four centres, reached by permuting the shells
! rather than by calling ip2: the integral is symmetric under those
! exchanges, and permuting also exercises the driver on orderings the
! canonical energy path would never produce.
!
! A water dimer in 6-31G by default, through the PUBLIC wrappers rather
! than the module entries, because the wrappers are what a caller sees.
! The finite difference costs 6*natm full energy evaluations over a naive
! quartet loop, so the default is small; `rotaxis_grad_molecule 24` runs
! the same thing on caffeine when there is time for it.
!
program rotaxis_grad_molecule
   use iso_c_binding, only: c_int, c_ptr
   use cint_const, only: dp
   use cint_bas, only: ATM_SLOTS, BAS_SLOTS, ANG_OF, NCTR_OF, NPRIM_OF, KAPPA_OF, &
                       PTR_EXP, PTR_COEFF, ATOM_OF, PTR_COORD, CHARGE_OF, &
                       cint_gto_norm
   use cint_envs, only: PTR_ENV_START
   use libcint_fortran, only: libcint_2e_sph, libcint_2e_ip1_rotaxis_sph, &
                              libcint_rotaxis_grad_supported
   implicit none

   integer, parameter :: MAXA = 40, MAXB = 300, MAXE = 40000
   integer(c_int) :: atm(ATM_SLOTS, MAXA), bas(BAS_SLOTS, MAXB), natm, nbas
   real(dp) :: env(MAXE)
   integer  :: envoff, ao_loc(0:MAXB), nao, coord_ptr(MAXA)
   real(dp), allocatable :: dm(:,:)
   real(dp) :: ganal(3, MAXA), gfd(3, MAXA)
   integer  :: ia, ax, nskip, nuse
   character(len=8) :: arg
   real(dp) :: h, ep, em, worst, scal
   real(dp), parameter :: TOL = 2.0e-7_dp   ! a central difference on ~1e2 Ha

   nuse = 6
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg); read(arg,*) nuse
   end if
   call build()
   print '(A,I0,A,I0,A,I0)', "  6-31G: atoms ", natm, "  shells ", nbas, &
        "  basis functions ", nao
   allocate(dm(nao, nao))
   call fake_density()

   call analytic(ganal, nskip)
   print '(A,I0,A)', "  gradient assembled through the public wrappers (", nskip, &
        " quartets fell back to Rys)"

   ! central difference of the same energy expression, same density
   h = 2.0e-4_dp
   gfd = 0.0_dp
   do ia = 1, natm
      do ax = 1, 3
         call shift(ia, ax, +h); ep = energy()
         call shift(ia, ax, -2.0_dp*h); em = energy()
         call shift(ia, ax, +h)
         gfd(ax, ia) = (ep - em)/(2.0_dp*h)
      end do
   end do

   scal = max(maxval(abs(gfd(:,1:natm))), 1.0e-30_dp)
   worst = maxval(abs(ganal(:,1:natm) - gfd(:,1:natm)))/scal
   print '(A,ES12.4)', "  largest gradient component : ", scal
   print '(A,ES12.4)', "  worst scaled difference    : ", worst
   print '(A)', "  atom   analytic dE2/dz      finite difference"
   do ia = 1, min(natm, 5)
      print '(I6,2ES22.12)', ia, ganal(3, ia), gfd(3, ia)
   end do
   if (worst > TOL) then
      print '(A)', "  RESULT: FAIL"
      stop 1
   end if
   print '(A)', "  RESULT: PASS (rotated-axis gradient matches finite difference)"

contains

   ! Every shell here is one contraction and not an L shell, so the
   ! spherical dimension is just 2l+1.
   integer function nsph(ish) result(n)
      integer, intent(in) :: ish
      n = 2*bas(ANG_OF+1, ish+1) + 1
   end function nsph

   ! dE2/dX with D fixed: every quartet contributes to all four centres.
   subroutine analytic(g, nfall)
      real(dp), intent(out) :: g(3, MAXA)
      integer,  intent(out) :: nfall
      integer  :: i, j, k, l, p, shq(4), sh4(4), dims(4), n, ret
      integer  :: a, b, c, d, di, dj, dk, dl, at, m
      real(dp) :: buf(3*20**4), v, s8
      integer  :: ord(4, 4)
      g = 0.0_dp; nfall = 0
      ! the four permutations that put each slot first, as exchanges the
      ! integral is invariant under
      ord(:,1) = [1,2,3,4]; ord(:,2) = [2,1,3,4]
      ord(:,3) = [3,4,1,2]; ord(:,4) = [4,3,1,2]
      do i = 0, nbas-1
      do j = 0, i
      do k = 0, i
      do l = 0, k
         if (i == k .and. j < l) cycle
         shq = [i, j, k, l]
         s8 = 1.0_dp
         if (i /= j) s8 = s8*2.0_dp
         if (k /= l) s8 = s8*2.0_dp
         if (i /= k .or. j /= l) s8 = s8*2.0_dp
         do p = 1, 4
            sh4 = [shq(ord(1,p)), shq(ord(2,p)), shq(ord(3,p)), shq(ord(4,p))]
            if (.not. libcint_rotaxis_grad_supported(sh4, bas, nbas)) then
               nfall = nfall + 1
               cycle
            end if
            dims = [nsph(sh4(1)), nsph(sh4(2)), nsph(sh4(3)), nsph(sh4(4))]
            n = dims(1)*dims(2)*dims(3)*dims(4)
            ret = libcint_2e_ip1_rotaxis_sph(buf, sh4, atm, natm, bas, nbas, env)
            if (ret == 0) cycle
            at = bas(ATOM_OF+1, sh4(1)+1) + 1
            di = dims(1); dj = dims(2); dk = dims(3); dl = dims(4)
            do m = 0, 2
            do d = 0, dl-1
            do c = 0, dk-1
            do b = 0, dj-1
            do a = 0, di-1
               v = buf(a + di*(b + dj*(c + dk*d)) + m*n + 1)
               g(m+1, at) = g(m+1, at) + s8*v &
                  * dm(ao_loc(sh4(1))+a+1, ao_loc(sh4(2))+b+1) &
                  * dm(ao_loc(sh4(3))+c+1, ao_loc(sh4(4))+d+1)
            end do
            end do
            end do
            end do
            end do
         end do
      end do
      end do
      end do
      end do
      ! ip1 is MINUS the derivative with respect to the centre, and E2
      ! carries the 1/2
      g = -0.5_dp*g
   end subroutine analytic

   real(dp) function energy() result(e)
      integer  :: i, j, k, l, shq(4), dims(4), ret, a, b, c, d, di, dj, dk, dl
      real(dp) :: buf(20**4), s8
      e = 0.0_dp
      do i = 0, nbas-1
      do j = 0, i
      do k = 0, i
      do l = 0, k
         if (i == k .and. j < l) cycle
         shq = [i, j, k, l]
         s8 = 1.0_dp
         if (i /= j) s8 = s8*2.0_dp
         if (k /= l) s8 = s8*2.0_dp
         if (i /= k .or. j /= l) s8 = s8*2.0_dp
         dims = [nsph(i), nsph(j), nsph(k), nsph(l)]
         ret = libcint_2e_sph(buf, shq, atm, natm, bas, nbas, env)
         if (ret == 0) cycle
         di = dims(1); dj = dims(2); dk = dims(3); dl = dims(4)
         do d = 0, dl-1
         do c = 0, dk-1
         do b = 0, dj-1
         do a = 0, di-1
            e = e + s8*buf(a + di*(b + dj*(c + dk*d)) + 1) &
               * dm(ao_loc(i)+a+1, ao_loc(j)+b+1) * dm(ao_loc(k)+c+1, ao_loc(l)+d+1)
         end do
         end do
         end do
         end do
      end do
      end do
      end do
      end do
      e = 0.5_dp*e
   end function energy

   subroutine shift(iat, axis, d)
      integer,  intent(in) :: iat, axis
      real(dp), intent(in) :: d
      env(coord_ptr(iat) + axis) = env(coord_ptr(iat) + axis) + d
   end subroutine shift

   subroutine fake_density()
      integer :: a, b
      integer(8) :: s
      s = 918273_8
      do a = 1, nao
         do b = 1, a
            s = mod(s*6364136223846793005_8 + 1442695040888963407_8, 2147483647_8)
            dm(a,b) = (real(s,dp)/2147483647.0_dp - 0.5_dp)*0.4_dp
            dm(b,a) = dm(a,b)
         end do
      end do
   end subroutine fake_density

   subroutine build()
      real(dp) :: g(3, 24)
      integer  :: z(24), i
      ! a water dimer in the first six slots -- a real hydrogen-bonded
      ! geometry, bohr -- then caffeine's atoms after it, so the argument
      ! selects how much of the list to use.
      real(dp), parameter :: wd(3,6) = reshape([ &
          0.0000_dp,  0.0000_dp,  0.0000_dp, &
          1.4200_dp,  0.0000_dp,  1.0900_dp, &
         -1.4200_dp,  0.0000_dp,  1.0900_dp, &
          0.0000_dp,  0.0000_dp,  5.5900_dp, &
          0.0000_dp,  1.4300_dp,  6.6900_dp, &
          0.0000_dp, -1.4300_dp,  6.6900_dp], [3,6])
      integer, parameter :: wz(6) = [8,1,1,8,1,1]
      z = [6,6,6,6,7,7,7,7,8,8,6,6,6,1,1,1,1,1,1,1,1,1,1,1]
      g = reshape([ &
         0.00_dp, 0.00_dp, 0.00_dp,   2.62_dp, 0.10_dp, 0.05_dp, &
         4.08_dp, 2.30_dp, 0.00_dp,   2.30_dp, 4.55_dp, 0.10_dp, &
        -0.30_dp, 2.50_dp, 0.05_dp,   5.10_dp,-1.80_dp, 0.00_dp, &
        -1.90_dp, 4.60_dp, 0.10_dp,   6.60_dp, 3.00_dp,-0.10_dp, &
         3.50_dp,-2.00_dp, 0.05_dp,   3.20_dp, 6.70_dp, 0.15_dp, &
        -4.20_dp, 4.10_dp, 0.20_dp,   7.70_dp,-2.60_dp, 0.10_dp, &
         8.60_dp, 4.00_dp,-0.20_dp,  -1.10_dp,-1.60_dp, 0.00_dp, &
         7.80_dp, 2.10_dp, 1.50_dp,  -4.60_dp, 6.00_dp, 0.30_dp, &
        -5.30_dp, 3.20_dp, 1.60_dp,  -4.70_dp, 3.10_dp,-1.40_dp, &
         7.40_dp,-4.50_dp, 0.20_dp,   9.00_dp,-2.10_dp, 1.80_dp, &
         9.00_dp,-2.20_dp,-1.50_dp,   9.90_dp, 2.90_dp,-0.30_dp, &
         8.40_dp, 5.20_dp, 1.40_dp,   8.30_dp, 5.00_dp,-1.90_dp], [3,24])
      if (nuse <= 6) then
         g(:,1:6) = wd
         z(1:6) = wz
      end if
      atm = 0; bas = 0; env = 0.0_dp; natm = 0; nbas = 0; envoff = PTR_ENV_START
      do i = 1, nuse
         natm = natm + 1
         atm(CHARGE_OF+1, natm) = z(i)
         atm(PTR_COORD+1, natm) = envoff
         coord_ptr(natm) = envoff
         env(envoff+1:envoff+3) = g(:, i)
         envoff = envoff + 3
      end do
      do i = 1, nuse
         call add_6_31g(i-1, z(i))
      end do
      ao_loc(0) = 0
      do i = 0, nbas-1
         ao_loc(i+1) = ao_loc(i) + nsph(i)
      end do
      nao = ao_loc(nbas)
   end subroutine build

   ! 6-31G, split-valence: the s/p shells are separate here rather than
   ! fused, which keeps this test independent of the L-shell handling that
   ! l_shell_check already covers.
   subroutine add_6_31g(iat, z)
      integer, intent(in) :: iat, z
      real(dp) :: e6(6), c6(6), e3(3), c3(3), e1(1)
      if (z == 1) then
         e3 = [18.7311370_dp, 2.8253937_dp, 0.6401217_dp]
         c3 = [0.03349460_dp, 0.23472695_dp, 0.81375733_dp]
         call put(iat, 0, 3, e3, c3)
         e1 = [0.1612778_dp]; call put(iat, 0, 1, e1, [1.0_dp])
      else
         ! one tight core s, then a valence s and p pair (inner + outer)
         e6 = [3047.5249_dp, 457.36951_dp, 103.94869_dp, 29.210155_dp, 9.2866630_dp, 3.1639270_dp]
         c6 = [0.0018347_dp, 0.0140373_dp, 0.0688426_dp, 0.2321844_dp, 0.4679413_dp, 0.3623120_dp]
         if (z == 7) e6 = e6*1.17_dp
         if (z == 8) e6 = e6*1.36_dp
         call put(iat, 0, 6, e6, c6)
         e3 = [7.8682724_dp, 1.8812885_dp, 0.5442493_dp]
         if (z == 7) e3 = e3*1.17_dp
         if (z == 8) e3 = e3*1.36_dp
         c3 = [-0.1193324_dp, -0.1608542_dp, 1.1434564_dp]
         call put(iat, 0, 3, e3, c3)
         c3 = [0.0689991_dp, 0.3164240_dp, 0.7443083_dp]
         call put(iat, 1, 3, e3, c3)
         e1 = [0.1687144_dp]
         if (z == 7) e1 = e1*1.17_dp
         if (z == 8) e1 = e1*1.36_dp
         call put(iat, 0, 1, e1, [1.0_dp])
         call put(iat, 1, 1, e1, [1.0_dp])
      end if
   end subroutine add_6_31g

   subroutine put(iat, l, np, e, c)
      integer,  intent(in) :: iat, l, np
      real(dp), intent(in) :: e(:), c(:)
      integer :: ip_, eptr, cptr
      eptr = envoff; env(eptr+1:eptr+np) = e(1:np); envoff = envoff + np
      cptr = envoff
      do ip_ = 1, np
         env(cptr+ip_) = c(ip_)*cint_gto_norm(l, e(ip_))
      end do
      envoff = envoff + np
      nbas = nbas + 1
      bas(ATOM_OF+1, nbas) = iat; bas(ANG_OF+1, nbas) = l
      bas(NPRIM_OF+1, nbas) = np; bas(NCTR_OF+1, nbas) = 1
      bas(KAPPA_OF+1, nbas) = 0
      bas(PTR_EXP+1, nbas) = eptr; bas(PTR_COEFF+1, nbas) = cptr
   end subroutine put

end program rotaxis_grad_molecule
