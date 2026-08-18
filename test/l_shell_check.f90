!
! L shells: the packed answer against the split one.
!
! An L shell carries s and p on one set of exponents, with a contraction
! column for each.  libcint's `bas` has one ANG_OF slot and cannot say that,
! so every packer splits such a shell into an s shell and a p shell -- which
! is correct, and recomputes the exponential prefactors, the pair data, the
! Rys roots and the b/c setup twice for exponents that are the same.
!
! THIS TEST CANNOT BE BIT-IDENTICAL AGAINST THE C, and not because anything
! here is approximate: libcint cannot be given the input.  The bar it is held
! to instead is that the packed calculation reproduces the split one, which
! is the same integral computed a different way -- so the comparison is
! against libfint's own split-shell path, which int2e_check already holds to
! the last bit against libcint.  Transitively that is a comparison with the C.
!
! The two runs differ in the last bits for two structural reasons, both
! deliberate:
!
!   * nrys_roots comes from the highest l in the quartet, so the ss
!     sub-block of an L quartet is computed with more roots than a genuine
!     ss shell pair would use.  Rys quadrature is exact to degree 2n-1, so
!     the extra roots contribute exactly zero -- in exact arithmetic.  In
!     floating point they contribute rounding.
!
!   * the shell normalisation factor moves.  For a plain shell it rides in
!     common_factor; an L shell's two sub-blocks want different ones, so it
!     rides on the contraction coefficient instead.  Same product, different
!     association order.
!
! Both are a few ulp, and both are ABSOLUTE: the difference is set by the
! size of the largest term that went into the sum, not by the size of the
! answer.  So the scaled measure is the one the tolerance is on, and the
! per-element one is reported with a loose gate because it necessarily grows
! as the element shrinks.  Measured on the basis below: scaled to the largest
! integral in the set the two agree to 4e-15 spherical and 1e-14 Cartesian;
! element by element, over the entries above 1e-3, to 5e-13.
!
! The two bases share ONE env.  A split shell points at the same exponents
! and at one column of the same coefficient block as the L shell it came
! from, so the two runs cannot disagree because a number was typed twice.
!
program l_shell_check
   use cint_const, only: dp
   use cint_bas, only: ATM_SLOTS, BAS_SLOTS, ATOM_OF, ANG_OF, NPRIM_OF, &
                       NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, &
                       CHARGE_OF, PTR_COORD, KAPPA_SP_SHELL, &
                       cint_gto_norm, cint_cgto_spheric, cint_cgto_cart
   use cint_envs, only: PTR_ENV_START
   use cint_workspace, only: cint_ws
   use cint_2e, only: int2e_sph, int2e_cart
   implicit none

   integer, parameter :: MAXBAS = 64, MAXENV = 4096, NATM = 3
   integer  :: atm(0:ATM_SLOTS*NATM-1)
   integer  :: pbas(0:BAS_SLOTS*MAXBAS-1), sbas(0:BAS_SLOTS*MAXBAS-1)
   real(dp) :: env(0:MAXENV-1)
   real(dp) :: coords(3, 0:NATM-1)
   integer  :: npbas, nsbas, envoff

   real(dp), allocatable :: eri_packed(:), eri_split(:)
   integer  :: nao_p, nao_s, sph
   real(dp) :: worst_rel, worst_abs
   integer  :: nbad
   real(dp), parameter :: TOL = 1.0e-13_dp
   real(dp), parameter :: TOL_ELEM = 1.0e-11_dp

   call build_basis()

   nbad = 0
   do sph = 1, 0, -1
      nao_p = nao_of(pbas, npbas, sph)
      nao_s = nao_of(sbas, nsbas, sph)
      if (nao_p /= nao_s) then
         print "(A,I0,A,I0)", "  FAIL: packed nao ", nao_p, " /= split nao ", nao_s
         stop 1
      end if
      allocate(eri_packed(0:nao_p**4-1), eri_split(0:nao_s**4-1))
      eri_packed = 0.0_dp; eri_split = 0.0_dp
      call all_eri(pbas, npbas, sph, nao_p, eri_packed)
      call all_eri(sbas, nsbas, sph, nao_s, eri_split)
      call compare(eri_packed, eri_split, nao_p**4, sph)
      deallocate(eri_packed, eri_split)
   end do

   if (nbad > 0) stop 1
   print "(A)", "  RESULT: PASS (L-packed reproduces split-shell)"

contains

   ! Three atoms, and a basis with an L shell of every shape the packing can
   ! produce: segmented, generally contracted, and one primitive.  The d
   ! shell is there so the cart-to-spherical stages run on a quartet that
   ! also has an L shell in it -- an L shell needs no transform, and the code
   ! that decides that has to leave the d alone.
   subroutine build_basis()
      real(dp) :: e3(3), c3(3), e2(2), c22(2,4), e1(1), c1(1)
      real(dp) :: cL3(3,2), cL1(1,2)

      atm = 0
      call set_atom(0, 6, [0.0_dp, 0.0_dp, 0.0_dp])
      call set_atom(1, 1, [0.0_dp, 0.0_dp, 1.9_dp])
      call set_atom(2, 8, [1.6_dp, 0.0_dp, -1.1_dp])

      npbas = 0; nsbas = 0
      envoff = PTR_ENV_START + 3*NATM
      env = 0.0_dp
      call set_atom_coords()

      ! carbon: a core s, a segmented L shell, a one-primitive L shell
      e3 = [71.6168370_dp, 13.0450960_dp, 3.5305122_dp]
      c3 = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      call add_shell(0, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)

      e3 = [2.9412494_dp, 0.6834831_dp, 0.2222899_dp]
      cL3(:,1) = [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp]
      cL3(:,2) = [ 0.15591627_dp, 0.60768372_dp, 0.39195739_dp]
      call add_shell(0, 1, 3, 1, e3, cL3, .true.)

      e1 = [0.1687144_dp]
      cL1(:,1) = [1.0_dp]
      cL1(:,2) = [1.0_dp]
      call add_shell(0, 1, 1, 1, e1, cL1, .true.)

      ! hydrogen: two plain s shells
      e3 = [18.7311370_dp, 2.8253937_dp, 0.6401217_dp]
      c3 = [0.03349460_dp, 0.23472695_dp, 0.81375733_dp]
      call add_shell(1, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)
      e1 = [0.1612778_dp]
      c1 = [1.0_dp]
      call add_shell(1, 0, 1, 1, e1, reshape(c1, [1,1]), .false.)

      ! oxygen: a core s, a GENERALLY CONTRACTED L shell -- two contractions
      ! sharing the same two primitives -- and a d polarisation shell
      e3 = [130.7093200_dp, 23.8088610_dp, 6.4436083_dp]
      c3 = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      call add_shell(2, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)

      e2 = [5.0331513_dp, 1.1695961_dp]
      c22(:,1) = [-0.09996723_dp, 0.39951283_dp]   ! s, contraction 0
      c22(:,2) = [ 0.21000000_dp, 0.83000000_dp]   ! s, contraction 1
      c22(:,3) = [ 0.15591627_dp, 0.60768372_dp]   ! p, contraction 0
      c22(:,4) = [ 0.44000000_dp,-0.27000000_dp]   ! p, contraction 1
      call add_shell(2, 1, 2, 2, e2, c22, .true.)

      e1 = [0.8000000_dp]
      c1 = [1.0_dp]
      call add_shell(2, 2, 1, 1, e1, reshape(c1, [1,1]), .false.)
   end subroutine build_basis

   subroutine set_atom(ia, z, r)
      integer,  intent(in) :: ia, z
      real(dp), intent(in) :: r(3)
      atm(ATM_SLOTS*ia + CHARGE_OF) = z
      atm(ATM_SLOTS*ia + PTR_COORD) = PTR_ENV_START + 3*ia
      coords(:, ia) = r
   end subroutine set_atom

   subroutine set_atom_coords()
      integer :: ia
      do ia = 0, NATM - 1
         env(PTR_ENV_START + 3*ia : PTR_ENV_START + 3*ia + 2) = coords(:, ia)
      end do
   end subroutine set_atom_coords

   ! One shell into the packed basis, and its equivalent into the split one.
   !
   ! `c` is (nprim, ncol) with ncol = nctr for a plain shell and 2*nctr for
   ! an L shell -- the s columns then the p columns, which is the layout a
   ! plain shell with 2*nctr contractions already has, and the reason no
   ! packer has to move a number to mark one.
   !
   ! The split emission is s(ctr 0), p(ctr 0), s(ctr 1), p(ctr 1)... rather
   ! than every s and then every p.  That is what makes the two bases share
   ! an AO ordering: the driver lays a shell's output out contraction-major,
   ! so the packed L shell gives [s p p p][s p p p] and the split shells give
   ! exactly the same sequence.  With the s shells grouped first they would
   ! not, and the test would need a permutation to say something it can
   ! instead say directly.
   subroutine add_shell(ia, l, nprim, nctr, e, c, is_sp)
      integer,  intent(in) :: ia, l, nprim, nctr
      real(dp), intent(in) :: e(:), c(:,:)
      logical,  intent(in) :: is_sp
      integer :: ncol, m, ip, eptr, cptr, ic, lm

      ncol = nctr
      if (is_sp) ncol = 2*nctr
      eptr = envoff
      env(eptr:eptr+nprim-1) = e(1:nprim)
      envoff = envoff + nprim
      cptr = envoff
      do m = 1, ncol
         lm = l
         if (is_sp) then
            lm = 1
            if (m <= nctr) lm = 0
         end if
         do ip = 1, nprim
            env(cptr + (m-1)*nprim + ip-1) = c(ip, m) * cint_gto_norm(lm, e(ip))
         end do
      end do
      envoff = envoff + ncol*nprim

      call put(pbas, npbas, ia, l, nprim, nctr, eptr, cptr, is_sp)
      if (is_sp) then
         do ic = 0, nctr - 1
            call put(sbas, nsbas, ia, 0, nprim, 1, eptr, cptr + ic*nprim, .false.)
            call put(sbas, nsbas, ia, 1, nprim, 1, eptr, cptr + (nctr+ic)*nprim, .false.)
         end do
      else
         call put(sbas, nsbas, ia, l, nprim, nctr, eptr, cptr, .false.)
      end if
   end subroutine add_shell

   subroutine put(bas, nbas, ia, l, nprim, nctr, eptr, cptr, is_sp)
      integer, intent(inout) :: bas(0:), nbas
      integer, intent(in)    :: ia, l, nprim, nctr, eptr, cptr
      logical, intent(in)    :: is_sp
      integer :: b
      b = BAS_SLOTS*nbas
      bas(b + ATOM_OF)   = ia
      bas(b + ANG_OF)    = l
      bas(b + NPRIM_OF)  = nprim
      bas(b + NCTR_OF)   = nctr
      bas(b + KAPPA_OF)  = 0
      if (is_sp) bas(b + KAPPA_OF) = KAPPA_SP_SHELL
      bas(b + PTR_EXP)   = eptr
      bas(b + PTR_COEFF) = cptr
      nbas = nbas + 1
   end subroutine put

   integer function shell_dim(bas, ish, sph) result(d)
      integer, intent(in) :: bas(0:), ish, sph
      if (sph == 1) then
         d = cint_cgto_spheric(ish, bas)
      else
         d = cint_cgto_cart(ish, bas)
      end if
   end function shell_dim

   integer function nao_of(bas, nbas, sph) result(n)
      integer, intent(in) :: bas(0:), nbas, sph
      integer :: i
      n = 0
      do i = 0, nbas - 1
         n = n + shell_dim(bas, i, sph)
      end do
   end function nao_of

   ! Every quartet, scattered into one nao^4 block so the two bases can be
   ! compared without knowing anything about each other's shell structure.
   subroutine all_eri(bas, nbas, sph, nao, out)
      integer,  intent(in)    :: nbas, sph, nao
      integer,  intent(inout), target :: bas(0:)
      real(dp), intent(inout) :: out(0:)
      integer  :: ao_loc(0:MAXBAS), i, j, k, l, di, dj, dk, dl
      integer  :: shls(0:3), dims(0:3), a, b, c, d
      real(dp), allocatable :: buf(:)
      logical  :: hv
      type(cint_ws) :: ws

      ao_loc(0) = 0
      do i = 0, nbas - 1
         ao_loc(i+1) = ao_loc(i) + shell_dim(bas, i, sph)
      end do

      do i = 0, nbas - 1
      do j = 0, nbas - 1
      do k = 0, nbas - 1
      do l = 0, nbas - 1
         di = shell_dim(bas, i, sph); dj = shell_dim(bas, j, sph)
         dk = shell_dim(bas, k, sph); dl = shell_dim(bas, l, sph)
         allocate(buf(0:di*dj*dk*dl-1))
         buf = 0.0_dp
         shls = [i, j, k, l]
         dims = [di, dj, dk, dl]
         if (sph == 1) then
            hv = int2e_sph(buf, dims, shls, atm, NATM, bas, nbas, env, ws)
         else
            hv = int2e_cart(buf, dims, shls, atm, NATM, bas, nbas, env, ws)
         end if
         if (hv) then
            do d = 0, dl - 1
            do c = 0, dk - 1
            do b = 0, dj - 1
            do a = 0, di - 1
               out(((ao_loc(l)+d)*nao + ao_loc(k)+c)*nao*nao &
                   + (ao_loc(j)+b)*nao + ao_loc(i)+a) &
                  = buf(((d*dk + c)*dj + b)*di + a)
            end do
            end do
            end do
            end do
         end if
         deallocate(buf)
      end do
      end do
      end do
      end do
   end subroutine all_eri

   ! TWO MEASURES, BOTH REPORTED.
   !
   ! Scaled to the largest integral in the set is the one the tolerance is on.
   ! An ERI tensor is mostly entries that are zero or nearly zero by symmetry
   ! and distance, and a few ulp of rounding on one of those is a relative
   ! error of one against itself and of nothing against the integral it will
   ! be summed with.
   !
   ! Per element is the stronger statement and is reported alongside, over
   ! the entries big enough for it to mean anything.  It is what would catch
   ! a component permuted, a coefficient applied to the wrong sub-block, or a
   ! normalisation factor left in common_factor -- none of which show up as a
   ! small error.
   subroutine compare(a, b, n, sph)
      real(dp), intent(in) :: a(0:), b(0:)
      integer,  intent(in) :: n, sph
      integer  :: m, bad, nelem
      real(dp) :: d, r, scal, worst_elem
      character(len=4) :: what
      worst_rel = 0.0_dp; worst_abs = 0.0_dp; worst_elem = 0.0_dp
      bad = 0; nelem = 0
      scal = maxval(abs(b(0:n-1)))
      do m = 0, n - 1
         d = abs(a(m) - b(m))
         if (d > worst_abs) worst_abs = d
         r = d / scal
         if (r > worst_rel) worst_rel = r
         if (r > TOL) bad = bad + 1
         if (abs(b(m)) > 1.0e-3_dp) then
            nelem = nelem + 1
            r = d / abs(b(m))
            if (r > worst_elem) worst_elem = r
            if (r > TOL_ELEM) bad = bad + 1
         end if
      end do
      what = "cart"
      if (sph == 1) what = "sph "
      print "(A,A,1X,I0,A,ES10.2,A,ES10.2,A,I0,A,ES10.2,A,I0)", "  ", what, n, &
         "  worst abs ", worst_abs, "  scaled rel ", worst_rel, &
         "   |x|>1e-3: ", nelem, " worst rel ", worst_elem, "  over tol: ", bad
      nbad = nbad + bad
   end subroutine compare

end program l_shell_check
