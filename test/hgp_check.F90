!
! Obara-Saika/HGP against the Rys path.
!
! int2e_hgp_cart/sph compute the same integrals int2e_cart/sph do, by a
! different algorithm: the Obara-Saika vertical recurrence inside the
! primitive loop, contraction, then Head-Gordon-Pople's horizontal
! transfers once per contracted quartet.  THIS TEST CANNOT BE BIT-IDENTICAL, and not
! because either side is approximate: two exact algorithms round
! differently.  The bar is instead the one l_shell_check sets -- agreement
! to a scaled tolerance, where the scale is the largest integral in the
! quartet, because both algorithms' rounding is set by the size of the
! largest term that went into a sum and not by the size of the answer.
!
! The reference is libfint's own Rys path, which int2e_check holds to the
! last bit against libcint; transitively this compares with the C, and it
! needs no C to run.
!
! Coverage: every quartet of the reference bases up to the generated l,
! Cartesian and spherical, strided as int2e_check strides; then the
! adversarial systems, which push the Boys argument to both ends; then the
! L-shell basis of l_shell_check, packed, because an L shell is where the
! driver's column and component bookkeeping is most easily wrong.
!
! The recurrences themselves are checked before this, in Python, against a
! McMurchie-Davidson reference: scripts/hgp/check_numeric.py.
!
program hgp_check
   use cint_const, only: dp
   use cint_bas, only: ATM_SLOTS, BAS_SLOTS, ANG_OF, NCTR_OF, NPRIM_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, &
                       ATOM_OF, PTR_COORD, CHARGE_OF, KAPPA_SP_SHELL, &
                       cint_cgto_cart, cint_cgto_spheric, cint_gto_norm
   use cint_envs, only: PTR_ENV_START
   use cint_workspace, only: cint_ws
   use cint_test_systems
   use cint_2e, only: int2e_cart, int2e_sph
   use cint_hgp_2e, only: int2e_hgp_cart, int2e_hgp_sph, hgp_supported
   use cint_hgp_kernels, only: HGP_LMAX
   implicit none

   integer, parameter :: MAXBAS = 200, NATM = 8
   integer,  allocatable :: atm(:), bas(:)
   real(dp), allocatable :: env(:)
   integer  :: nbas
   type(cint_ws) :: ws

   ! Scaled to the largest value in the quartet, which is the size of the
   ! largest term either algorithm summed.  Below ABS_FLOOR the largest
   ! value is itself rounding -- a quartet that vanishes by symmetry, where
   ! Rys leaves 1e-18 and the rotated frame an exact zero -- and the
   ! comparison is absolute.
   real(dp), parameter :: TOL = 1.0e-12_dp
   real(dp), parameter :: ABS_FLOOR = 1.0e-14_dp
   integer  :: ncmp, nbad, nq, nskip, nshown
   real(dp) :: worst
   character(len=96) :: wat

   ncmp = 0; nbad = 0; nq = 0; nskip = 0; worst = 0.0_dp; wat = "(none)"

   call reference_bases()
   call adversarial()
   call l_shells()

   print '(A,I0,A,I0,A)', "  quartets        : ", nq, "  (", nskip, " beyond the generated l, skipped)"
   print '(A,I0)',        "  values compared : ", ncmp
   print '(A,ES10.2,A,A)',"  worst scaled    : ", worst, "   at ", trim(wat)
   print '(A,I0)',        "  over tolerance  : ", nbad
   if (nbad > 0) stop 1
   print '(A)', "  RESULT: PASS (Obara-Saika/HGP reproduces Rys)"

contains

   subroutine reference_bases()
      integer(4), allocatable :: catm(:,:), cbas(:,:)
      real(dp),   allocatable :: cenv(:)
      integer(4) :: coff, cnbas
      integer :: ibasis, i, stride

      allocate(catm(ATM_SLOTS, NATM), cbas(BAS_SLOTS, MAXBAS), cenv(20000))
      do ibasis = 1, n_reference_basis
         catm = 0; cbas = 0; cenv = 0.0_dp
         coff = PTR_ENV_START
         call setup_c2h6_geometry(catm, cenv, coff)
         call setup_reference_basis(ibasis, cbas, cenv, cnbas, coff)
         nbas = cnbas
         allocate(atm(0:ATM_SLOTS*NATM-1), bas(0:BAS_SLOTS*nbas-1), env(0:size(cenv)-1))
         do i = 1, NATM
            atm(ATM_SLOTS*(i-1):ATM_SLOTS*i-1) = int(catm(:, i))
         end do
         do i = 1, nbas
            bas(BAS_SLOTS*(i-1):BAS_SLOTS*i-1) = int(cbas(:, i))
         end do
         env = cenv
         stride = max(1, nbas/5)
         call sweep(stride, "ref:" // trim(reference_basis_name(ibasis)))
         deallocate(atm, bas, env)
      end do
   end subroutine reference_bases

   subroutine adversarial()
      integer(4), allocatable :: catm(:,:), cbas(:,:)
      real(dp),   allocatable :: cenv(:)
      integer(4) :: cnatm, cnbas
      integer :: which, i

      allocate(catm(ATM_SLOTS, NATM), cbas(BAS_SLOTS, MAXBAS), cenv(20000))
      do which = 1, n_adversarial
         catm = 0; cbas = 0; cenv = 0.0_dp
         call setup_adversarial(which, catm, cnatm, cbas, cnbas, cenv)
         nbas = cnbas
         allocate(atm(0:ATM_SLOTS*NATM-1), bas(0:BAS_SLOTS*nbas-1), env(0:size(cenv)-1))
         do i = 1, NATM
            atm(ATM_SLOTS*(i-1):ATM_SLOTS*i-1) = int(catm(:, i))
         end do
         do i = 1, nbas
            bas(BAS_SLOTS*(i-1):BAS_SLOTS*i-1) = int(cbas(:, i))
         end do
         env = cenv
         call sweep(1, "adv:" // trim(adversarial_name(which)))
         deallocate(atm, bas, env)
      end do
   end subroutine adversarial

   ! The packed L-shell basis of l_shell_check: three atoms, L shells that
   ! are segmented, generally contracted and single-primitive, and a d.
   subroutine l_shells()
      integer :: envoff
      real(dp) :: e3(3), c3(3), e2(2), c22(2,4), e1(1), cL3(3,2), cL1(1,2)

      allocate(atm(0:ATM_SLOTS*3-1), bas(0:BAS_SLOTS*MAXBAS-1), env(0:4000))
      atm = 0; bas = 0; env = 0.0_dp
      nbas = 0
      envoff = PTR_ENV_START
      call set_atom(0, 6, [0.0_dp, 0.0_dp, 0.0_dp], envoff)
      call set_atom(1, 1, [0.0_dp, 0.0_dp, 1.9_dp], envoff)
      call set_atom(2, 8, [1.6_dp, 0.0_dp, -1.1_dp], envoff)

      e3 = [71.6168370_dp, 13.0450960_dp, 3.5305122_dp]
      c3 = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      call add_shell(0, 0, 3, 1, e3, reshape(c3, [3,1]), .false., envoff)
      e3 = [2.9412494_dp, 0.6834831_dp, 0.2222899_dp]
      cL3(:,1) = [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp]
      cL3(:,2) = [ 0.15591627_dp, 0.60768372_dp, 0.39195739_dp]
      call add_shell(0, 1, 3, 1, e3, cL3, .true., envoff)
      e1 = [0.1687144_dp]
      cL1(:,1) = [1.0_dp]; cL1(:,2) = [1.0_dp]
      call add_shell(0, 1, 1, 1, e1, cL1, .true., envoff)
      e2 = [0.8_dp, 0.3_dp]
      c22(:,1) = [0.7_dp, 0.4_dp]; c22(:,2) = [-0.2_dp, 0.9_dp]
      c22(:,3) = [0.5_dp, 0.6_dp]; c22(:,4) = [ 0.3_dp, -0.8_dp]
      call add_shell(1, 1, 2, 2, e2, c22, .true., envoff)
      e1 = [1.1_dp]
      call add_shell(1, 0, 1, 1, e1, reshape([1.0_dp], [1,1]), .false., envoff)
      e2 = [1.3_dp, 0.45_dp]
      call add_shell(2, 2, 2, 1, e2, reshape([0.6_dp, 0.5_dp], [2,1]), .false., envoff)
      e3 = [5.0_dp, 1.2_dp, 0.35_dp]
      call add_shell(2, 1, 3, 1, e3, reshape([0.2_dp, 0.5_dp, 0.4_dp], [3,1]), .false., envoff)

      call sweep(1, "Lshells")
      deallocate(atm, bas, env)
   end subroutine l_shells

   subroutine set_atom(ia, z, r, envoff)
      integer,  intent(in) :: ia, z
      real(dp), intent(in) :: r(3)
      integer,  intent(inout) :: envoff
      atm(ATM_SLOTS*ia + CHARGE_OF) = z
      atm(ATM_SLOTS*ia + PTR_COORD) = envoff
      env(envoff:envoff+2) = r
      envoff = envoff + 3
   end subroutine set_atom

   subroutine add_shell(ia, l, nprim, nctr, e, c, is_sp, envoff)
      integer,  intent(in) :: ia, l, nprim, nctr
      real(dp), intent(in) :: e(:), c(:,:)
      logical,  intent(in) :: is_sp
      integer,  intent(inout) :: envoff
      integer :: ncol, m, ip, eptr, cptr, lm, b
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
      b = BAS_SLOTS*nbas
      bas(b + ATOM_OF) = ia
      bas(b + ANG_OF) = l
      bas(b + NPRIM_OF) = nprim
      bas(b + NCTR_OF) = nctr
      bas(b + KAPPA_OF) = 0
      if (is_sp) bas(b + KAPPA_OF) = KAPPA_SP_SHELL
      bas(b + PTR_EXP) = eptr
      bas(b + PTR_COEFF) = cptr
      nbas = nbas + 1
   end subroutine add_shell

   ! Every quartet of the current basis at the given stride, both layouts.
   subroutine sweep(stride, label)
      integer, intent(in) :: stride
      character(len=*), intent(in) :: label
      integer :: i, j, k, l, sph, shls(0:3), dims(0:3), n, natm
      real(dp), allocatable :: rys(:), rot(:)
      logical :: hv1, hv2

      natm = size(atm) / ATM_SLOTS
      nshown = 0
      do sph = 0, 1
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
      do l = 0, nbas-1, stride
         shls = [i, j, k, l]
         if (.not. hgp_supported(shls, bas)) then
            nskip = nskip + 1
            cycle
         end if
         if (sph == 1) then
            dims = [cint_cgto_spheric(i, bas), cint_cgto_spheric(j, bas), &
                    cint_cgto_spheric(k, bas), cint_cgto_spheric(l, bas)]
         else
            dims = [cint_cgto_cart(i, bas), cint_cgto_cart(j, bas), &
                    cint_cgto_cart(k, bas), cint_cgto_cart(l, bas)]
         end if
         n = product(dims)
         allocate(rys(0:n-1), rot(0:n-1))
         rys = 0.0_dp; rot = 0.0_dp
         if (sph == 1) then
            hv1 = int2e_sph(rys, dims, shls, atm, natm, bas, nbas, env, ws)
            hv2 = int2e_hgp_sph(rot, dims, shls, atm, natm, bas, nbas, env, ws)
         else
            hv1 = int2e_cart(rys, dims, shls, atm, natm, bas, nbas, env, ws)
            hv2 = int2e_hgp_cart(rot, dims, shls, atm, natm, bas, nbas, env, ws)
         end if
         nq = nq + 1
         call compare(rys, rot, n, label, sph, shls)
         deallocate(rys, rot)
      end do
      end do
      end do
      end do
      end do
   end subroutine sweep

   subroutine compare(a, b, n, label, sph, shls)
      real(dp), intent(in) :: a(0:), b(0:)
      integer,  intent(in) :: n, sph, shls(0:3)
      character(len=*), intent(in) :: label
      real(dp) :: scale, r
      integer :: m
      scale = maxval(abs(a(0:n-1)))
      do m = 0, n - 1
         r = abs(a(m) - b(m)) / max(scale, ABS_FLOOR)
         ncmp = ncmp + 1
         if (abs(a(m) - b(m)) <= ABS_FLOOR) r = 0.0_dp
         if (r > worst) then
            worst = r
            write(wat, '(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(label), " sph=", sph, &
               " quartet (", shls(0), ",", shls(1), ",", shls(2), ",", shls(3), ") m=", m
         end if
         if (r > TOL) then
            nbad = nbad + 1
            nshown = nshown + 1
            if (nshown <= 8) then
               print '(A,A,I0,A,4(I0,1X),A,I0,A,2ES23.15)', "  MISMATCH ", trim(label), sph, &
                  " quartet ", shls, " m=", m, " rys/rot ", a(m), b(m)
            end if
         end if
      end do
   end subroutine compare

end program hgp_check
