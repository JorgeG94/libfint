!
! The scalar ECP integrals themselves.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp.c.  See cint_ecp_num for why
! the upstream is PySCF rather than libcint, and cint_ecp_rad for why the
! agreement below this point is a tolerance rather than an equality.
!
! An ECP splits into two terms and they are computed by separate routines
! because they have different structure, not merely different parameters:
!
!   type 1, the local term, is the l = -1 channel.  It multiplies the
!           density by a radial function and nothing else, so the angular
!           work is a single contraction.
!   type 2, the semi-local term, projects onto each l channel in turn.  It
!           carries a projector on both sides, which is where the lambda
!           ladder and most of the cost come from.
!
! **The spin-orbit term is not ported.**  `ECPtype_so_cart` and the angular
! momentum matrices it needs are about a third of the C's algorithm and serve
! `ECPso_*`, which is a different integral with a different output shape.
! Nothing here silently drops it: it was never claimed.
!
module cint_ecp
   use cint_const, only: dp
   use cint_bas, only: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       BAS_SLOTS, PTR_COORD, ATM_SLOTS, cint_len_cart
   use cint_g1e, only: cint_common_fac_sp
   use cint_ecp_num, only: ecp_gauss_chebyshev, ECP_EXPCUTOFF, ECP_LEVEL0, &
                           ECP_LEVEL_MAX, SO_TYPE_OF, ECP_LMAX
   use cint_ecp_ang, only: ecp_type1_rad_ang, ecp_type2_facs_ang, CART_POW_Y, CART_POW_Z
   use cint_ecp_rad, only: ecp_rad_part, ecp_type1_rad_part, ecp_type1_static_facs, &
                           ecp_type2_facs_rad
   implicit none
   private

   public :: ecp_loc_ecpbas
   public :: ecp_check_3c_overlap
   public :: ecp_type1_cart
   public :: ecp_type2_cart
   public :: ecp_type_scalar_cart

contains

   !> Group consecutive ECP shells that share an atom, an l and a spin-orbit flag
   !>
   !> Shells differing only in the r exponent are evaluated together by
   !> `ecp_rad_part`, which sums them, so the grouping is what lets one pass
   !> over the radial grid serve all of them.  `ecploc` comes back holding
   !> nslots+1 boundaries, the last being `necpbas`.
   pure subroutine ecp_loc_ecpbas(ecploc, nslots, ecpbas, necpbas)
      integer, intent(out) :: ecploc(0:)
      integer, intent(out) :: nslots
      integer, intent(in) :: ecpbas(0:)
      integer, intent(in) :: necpbas

      integer :: i, l, so, atm_id, atm_last, l_last, so_last

      ecploc(0) = 0
      nslots = 0
      if (necpbas == 0) return

      atm_last = ecpbas(ATOM_OF)
      l_last = ecpbas(ANG_OF)
      so_last = ecpbas(SO_TYPE_OF)
      do i = 1, necpbas - 1
         atm_id = ecpbas(ATOM_OF + i*BAS_SLOTS)
         l = ecpbas(ANG_OF + i*BAS_SLOTS)
         so = ecpbas(SO_TYPE_OF + i*BAS_SLOTS)
         if (atm_id /= atm_last .or. l /= l_last .or. so /= so_last) then
            nslots = nslots + 1
            ecploc(nslots) = i
            atm_last = atm_id
            l_last = l
            so_last = so
         end if
      end do
      nslots = nslots + 1
      ecploc(nslots) = necpbas
   end subroutine ecp_loc_ecpbas

   pure function distance_square(r1, r2) result(d)
      real(dp), intent(in) :: r1(3), r2(3)
      real(dp) :: d
      d = (r1(1) - r2(1))**2 + (r1(2) - r2(2))**2 + (r1(3) - r2(3))**2
   end function distance_square

   !> Whether a shell pair and an ECP centre overlap enough to bother with
   !>
   !> The screen that makes an ECP affordable on a real molecule: a potential
   !> is short ranged, and most shell pairs are nowhere near it.
   !>
   !> Only the last primitive of each shell is tested, which is exact rather
   !> than approximate given libcint's convention that primitives are sorted
   !> with the smallest exponent last -- the most diffuse one, and so the one
   !> that reaches furthest. Testing any other would be the conservative test
   !> done wrong.
   pure function ecp_check_3c_overlap(shls, atm, bas, env, rc, ecpshls, ecpbas) result(keep)
      integer, intent(in) :: shls(0:1)
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      real(dp), intent(in) :: rc(3)
      integer, intent(in) :: ecpshls(0:1)
      integer, intent(in) :: ecpbas(0:)
      logical :: keep

      integer :: ish, jsh, iprim, jprim, csh, kprim, pk
      real(dp) :: ai, aj, ri(3), rj(3)
      real(dp) :: rrab, rrca, rrcb, aiaj, aij, ajak, aiak, aijk, eijk

      ish = shls(0)
      jsh = shls(1)
      iprim = bas(NPRIM_OF + ish*BAS_SLOTS)
      jprim = bas(NPRIM_OF + jsh*BAS_SLOTS)
      ai = env(bas(PTR_EXP + ish*BAS_SLOTS) + iprim - 1)
      aj = env(bas(PTR_EXP + jsh*BAS_SLOTS) + jprim - 1)
      ri = env(atm(PTR_COORD + bas(ATOM_OF + ish*BAS_SLOTS)*ATM_SLOTS): &
               atm(PTR_COORD + bas(ATOM_OF + ish*BAS_SLOTS)*ATM_SLOTS) + 2)
      rj = env(atm(PTR_COORD + bas(ATOM_OF + jsh*BAS_SLOTS)*ATM_SLOTS): &
               atm(PTR_COORD + bas(ATOM_OF + jsh*BAS_SLOTS)*ATM_SLOTS) + 2)

      rrab = distance_square(ri, rj)
      rrca = distance_square(rc, ri)
      rrcb = distance_square(rc, rj)
      aiaj = ai*aj*rrab
      aij = ai + aj

      keep = .false.
      do csh = ecpshls(0), ecpshls(1) - 1
         kprim = ecpbas(csh*BAS_SLOTS + NPRIM_OF)
         pk = ecpbas(csh*BAS_SLOTS + PTR_EXP)
         aijk = aij + env(pk + kprim - 1)
         aiak = ai*env(pk + kprim - 1)*rrca
         ajak = aj*env(pk + kprim - 1)*rrcb
         eijk = (aiaj + aiak + ajak)/aijk
         if (eijk < ECP_EXPCUTOFF) then
            keep = .true.
            return
         end if
      end do
   end function ecp_check_3c_overlap

   !> Contraction coefficients scaled by the Gaussian prefactor at the ECP
   pure subroutine scale_coeff(cei, ci, ai, r2ca, npi, nci, li)
      real(dp), intent(out) :: cei(0:)
      real(dp), intent(in) :: ci(0:), ai(0:)
      real(dp), intent(in) :: r2ca
      integer, intent(in) :: npi, nci, li

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp) :: tmp, common_fac
      integer :: ip, ic

      common_fac = cint_common_fac_sp(li)*4.0_dp*PI
      do ip = 0, npi - 1
         tmp = exp(-ai(ip)*r2ca)*common_fac
         do ic = 0, nci - 1
            cei(ic*npi + ip) = ci(ic*npi + ip)*tmp
         end do
      end do
   end subroutine scale_coeff

   pure function close_enough(x, y) result(c)
      real(dp), intent(in) :: x, y
      logical :: c
      c = abs(x - y) <= 1.0e-12_dp*max(abs(x), abs(y))
   end function close_enough

   !> The local (type 1) ECP term, over one Cartesian shell pair
   !>
   !> Accumulates into `gctr`, which the caller has zeroed: both ECP types add
   !> into the same block.  Returns whether anything was added at all, which
   !> the driver uses to skip the transform on a pair the screen rejected.
   !>
   !> **The radial grid is adaptive**, and that is the part worth reading
   !> slowly.  It starts at 31 points and doubles, each level reusing the
   !> previous points by striding a 2047-point grid rather than rebuilding
   !> one, and stops when every primitive pair has repeated itself twice.
   !> Twice and not once: one agreement can be a coincidence of a quadrature
   !> that has not yet resolved the integrand at all, and the C's own comment
   !> is that this converges from below.
   function ecp_type1_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env) &
      result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer, intent(in) :: shls(0:1)
      integer, intent(in) :: ecpbas(0:)
      integer, intent(in) :: necpbas
      integer, intent(in) :: atm(0:), bas(0:)
      integer, intent(in) :: natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      integer :: ish, jsh, li, lj, npi, npj, nci, ncj, nfi, nfj
      integer :: pai, paj, pci, pcj, pri, prj
      integer :: lilj1, d1, d2, d3, di1, di2, di3, dj1, dj2, dj3
      integer :: nrs_full, nslots, iloc, atm_id, prc
      integer :: i, n, ip, jp, ic, jc, mi, mj
      integer :: i1, i2, i3, j1, j2, j3, ix, iy, iz, jx, jy, jz
      integer :: level, nrs0, start, step, nrs, pr, po
      real(dp) :: rca(3), rcb(3), rij(3), ri(3), rj(3), rc(3)
      real(dp) :: fac, wtscale
      logical :: all_conv, pair_close
      integer, allocatable :: ecploc(:), converged(:)
      real(dp), allocatable :: rs(:), ws(:), ur(:), rad_all(:), rad_ang(:)
      real(dp), allocatable :: rad_ang_all(:), ifac(:), jfac(:), cei(:), cej(:), plast(:)

      has_value = 0
      if (necpbas == 0) return

      ish = shls(0)
      jsh = shls(1)
      li = bas(ANG_OF + ish*BAS_SLOTS)
      lj = bas(ANG_OF + jsh*BAS_SLOTS)
      npi = bas(NPRIM_OF + ish*BAS_SLOTS)
      npj = bas(NPRIM_OF + jsh*BAS_SLOTS)
      nci = bas(NCTR_OF + ish*BAS_SLOTS)
      ncj = bas(NCTR_OF + jsh*BAS_SLOTS)
      nfi = cint_len_cart(li)
      nfj = cint_len_cart(lj)
      pai = bas(PTR_EXP + ish*BAS_SLOTS)
      paj = bas(PTR_EXP + jsh*BAS_SLOTS)
      pci = bas(PTR_COEFF + ish*BAS_SLOTS)
      pcj = bas(PTR_COEFF + jsh*BAS_SLOTS)
      pri = atm(PTR_COORD + bas(ATOM_OF + ish*BAS_SLOTS)*ATM_SLOTS)
      prj = atm(PTR_COORD + bas(ATOM_OF + jsh*BAS_SLOTS)*ATM_SLOTS)
      ri = env(pri:pri + 2)
      rj = env(prj:prj + 2)

      lilj1 = li + lj + 1
      d1 = lilj1
      d2 = d1*d1
      d3 = d2*d1
      di1 = li + 1; di2 = di1*di1; di3 = di2*di1
      dj1 = lj + 1; dj2 = dj1*dj1; dj3 = dj2*dj1

      nrs_full = 2**ECP_LEVEL_MAX
      allocate (rs(0:nrs_full - 1), ws(0:nrs_full - 1))
      call ecp_gauss_chebyshev(rs, ws, nrs_full)

      allocate (ecploc(0:necpbas))
      call ecp_loc_ecpbas(ecploc, nslots, ecpbas, necpbas)

      allocate (ur(0:nrs_full - 1))
      allocate (rad_all(0:npi*npj*d2 - 1))
      allocate (rad_ang(0:d3 - 1))
      allocate (rad_ang_all(0:nci*ncj*d3 - 1))
      allocate (ifac(0:nfi*di3 - 1), jfac(0:nfj*dj3 - 1))
      allocate (cei(0:npi*nci - 1), cej(0:npj*ncj - 1))
      allocate (plast(0:d2 - 1), converged(0:npi*npj - 1))

      do iloc = 0, nslots - 1
         ! The local term is the l = -1 channel, and only that channel.
         if (ecpbas(ANG_OF + ecploc(iloc)*BAS_SLOTS) /= -1 .or. &
             ecpbas(SO_TYPE_OF + ecploc(iloc)*BAS_SLOTS) == 1) cycle

         atm_id = ecpbas(ATOM_OF + ecploc(iloc)*BAS_SLOTS)
         prc = atm(PTR_COORD + atm_id*ATM_SLOTS)
         rc = env(prc:prc + 2)

         if (.not. ecp_check_3c_overlap(shls, atm, bas, env, rc, &
                                        ecploc(iloc:iloc + 1), ecpbas)) cycle

         has_value = 1
         rca = rc - ri
         rcb = rc - rj
         call scale_coeff(cei, env(pci:), env(pai:), &
                          rca(1)**2 + rca(2)**2 + rca(3)**2, npi, nci, li)
         call scale_coeff(cej, env(pcj:), env(paj:), &
                          rcb(1)**2 + rcb(2)**2 + rcb(3)**2, npj, ncj, lj)

         converged = 0
         rad_all = 0.0_dp
         nrs0 = 2**ECP_LEVEL0 - 1
         step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
         start = step - 1
         wtscale = real(step, dp)

         do level = ECP_LEVEL0, ECP_LEVEL_MAX
            call ecp_rad_part(ur, nrs, rs, start, nrs0, step, &
                              ecploc(iloc), ecploc(iloc + 1), ecpbas, env)
            if (nrs == 0) exit
            do n = 0, nrs - 1
               ur(n) = ur(n)*ws(start + n*step)*wtscale
            end do

            all_conv = .true.
            do ip = 0, npi - 1
               do jp = 0, npj - 1
                  if (converged(ip*npj + jp) >= 2) cycle
                  pr = (ip*npj + jp)*d2
                  do i = 0, d2 - 1
                     plast(i) = rad_all(pr + i)
                     ! Halved because the previous level's points reappear in
                     ! this one at twice the weight: the refinement reuses the
                     ! grid rather than starting over.
                     rad_all(pr + i) = rad_all(pr + i)*0.5_dp
                  end do
                  rij = env(pai + ip)*rca + env(paj + jp)*rcb
                  call ecp_type1_rad_part(rad_all(pr:), li + lj, &
                                          sqrt(rij(1)**2 + rij(2)**2 + rij(3)**2)*2.0_dp, &
                                          env(pai + ip) + env(paj + jp), &
                                          ur, rs(start:), nrs, step)
                  pair_close = .true.
                  do i = 0, d2 - 1
                     if (.not. close_enough(plast(i), rad_all(pr + i))) then
                        pair_close = .false.
                        exit
                     end if
                  end do
                  if (pair_close) then
                     converged(ip*npj + jp) = converged(ip*npj + jp) + 1
                     if (converged(ip*npj + jp) < 2) all_conv = .false.
                  else
                     converged(ip*npj + jp) = 0
                     all_conv = .false.
                  end if
               end do
            end do

            if (all_conv) exit
            nrs0 = 2**level - 1
            step = 2**(ECP_LEVEL_MAX - level)
            start = (start - 1)/2
            wtscale = wtscale*0.5_dp
         end do

         rad_ang_all = 0.0_dp
         do ip = 0, npi - 1
            do jp = 0, npj - 1
               rij = env(pai + ip)*rca + env(paj + jp)*rcb
               call ecp_type1_rad_ang(rad_ang, li + lj, rij, rad_all((ip*npj + jp)*d2:))
               do ic = 0, nci - 1
                  do jc = 0, ncj - 1
                     fac = cei(ic*npi + ip)*cej(jc*npj + jp)
                     pr = (ic*ncj + jc)*d3
                     do n = 0, d3 - 1
                        rad_ang_all(pr + n) = rad_ang_all(pr + n) + fac*rad_ang(n)
                     end do
                  end do
               end do
            end do
         end do

         call ecp_type1_static_facs(ifac, li, rca)
         call ecp_type1_static_facs(jfac, lj, rcb)

         do ic = 0, nci - 1
            do jc = 0, ncj - 1
               pr = (ic*ncj + jc)*d3
               do mi = 0, nfi - 1
                  iy = CART_POW_Y(mi); iz = CART_POW_Z(mi); ix = li - iy - iz
                  do mj = 0, nfj - 1
                     jy = CART_POW_Y(mj); jz = CART_POW_Z(mj); jx = lj - jy - jz
                     po = (jc*nfj + mj)*nci*nfi + ic*nfi + mi
                     do i1 = 0, ix
                        do i2 = 0, iy
                           do i3 = 0, iz
                              do j1 = 0, jx
                                 do j2 = 0, jy
                                    do j3 = 0, jz
                                       gctr(po) = gctr(po) &
                                          + ifac(mi*di3 + i1*di2 + i2*di1 + i3) &
                                          *jfac(mj*dj3 + j1*dj2 + j2*dj1 + j3) &
                                          *rad_ang_all(pr + (i1 + j1)*d2 + (i2 + j2)*d1 + i3 + j3)
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
      end do

      deallocate (rs, ws, ecploc, ur, rad_all, rad_ang, rad_ang_all)
      deallocate (ifac, jfac, cei, cej, plast, converged)
   end function ecp_type1_cart

   !> The semi-local (type 2) ECP term, over one Cartesian shell pair
   !>
   !> Accumulates into `gctr` beside the type-1 term.
   !>
   !> This is the expensive half.  Where type 1 has one radial integral per
   !> primitive pair, this has a projector on each side, so the radial factors
   !> carry a lambda index for the bra and another for the ket and the
   !> convergence test runs over the pair of them.  The adaptive refinement is
   !> the same scheme as type 1's and converges per (contraction, contraction,
   !> total degree) rather than per primitive pair, because the radial factors
   !> here are already contracted when they are formed.
   function ecp_type2_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env) &
      result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer, intent(in) :: shls(0:1)
      integer, intent(in) :: ecpbas(0:)
      integer, intent(in) :: necpbas
      integer, intent(in) :: atm(0:), bas(0:)
      integer, intent(in) :: natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      real(dp), parameter :: PI = acos(-1.0_dp)
      integer :: ish, jsh, li, lj, nci, ncj, nfi, nfj, di, pri, prj, prc
      integer :: nrs_full, nslots, iloc, atm_id, lc, dlc, lilj1, lilc1, ljlc1
      integer :: d2, d3, im, mq, i, j, n, ic, jc, lab, ijl, nrs
      integer :: level, nrs0, start, step, pradi, pradj, prur, prad, mi, mj, po
      real(dp) :: rca(3), rcb(3), ri(3), rj(3), rc(3)
      real(dp) :: dca, dcb, s, wtscale, common_fac, acc
      logical :: all_conv, pair_close
      integer, allocatable :: ecploc(:), converged(:)
      real(dp), allocatable :: rs(:), ws(:), rur(:), radi(:), radj(:)
      real(dp), allocatable :: angi(:), angj(:), buf(:), rad_all(:), plast(:)

      has_value = 0
      if (necpbas == 0) return

      ish = shls(0)
      jsh = shls(1)
      li = bas(ANG_OF + ish*BAS_SLOTS)
      lj = bas(ANG_OF + jsh*BAS_SLOTS)
      nci = bas(NCTR_OF + ish*BAS_SLOTS)
      ncj = bas(NCTR_OF + jsh*BAS_SLOTS)
      nfi = cint_len_cart(li)
      nfj = cint_len_cart(lj)
      di = nfi*nci
      pri = atm(PTR_COORD + bas(ATOM_OF + ish*BAS_SLOTS)*ATM_SLOTS)
      prj = atm(PTR_COORD + bas(ATOM_OF + jsh*BAS_SLOTS)*ATM_SLOTS)
      ri = env(pri:pri + 2)
      rj = env(prj:prj + 2)
      common_fac = cint_common_fac_sp(li)*cint_common_fac_sp(lj)*16.0_dp*PI*PI

      lilj1 = li + lj + 1
      nrs_full = 2**ECP_LEVEL_MAX
      allocate (rs(0:nrs_full - 1), ws(0:nrs_full - 1))
      call ecp_gauss_chebyshev(rs, ws, nrs_full)

      allocate (ecploc(0:necpbas))
      call ecp_loc_ecpbas(ecploc, nslots, ecpbas, necpbas)

      ! Sized for the largest lc this can meet, then indexed with the actual
      ! one: the buffers outlive the loop over ECP channels and each channel
      ! has its own lc.
      allocate (rad_all(0:nci*ncj*lilj1*(li + ECP_LMAX + 1)*(lj + ECP_LMAX + 1) - 1))
      allocate (angi(0:(li + 1)*nfi*(ECP_LMAX*2 + 1)*(li + ECP_LMAX + 1) - 1))
      allocate (angj(0:(lj + 1)*nfj*(ECP_LMAX*2 + 1)*(lj + ECP_LMAX + 1) - 1))
      allocate (buf(0:nfi*(ECP_LMAX*2 + 1)*(lj + ECP_LMAX + 1) - 1))
      allocate (rur(0:nrs_full*lilj1 - 1))
      allocate (radi(0:nci*(li + ECP_LMAX + 1)*nrs_full - 1))
      allocate (radj(0:ncj*(lj + ECP_LMAX + 1)*nrs_full - 1))
      allocate (plast(0:(li + ECP_LMAX + 1)*(lj + ECP_LMAX + 1) - 1))
      allocate (converged(0:nci*ncj*lilj1 - 1))

      do iloc = 0, nslots - 1
         lc = ecpbas(ANG_OF + ecploc(iloc)*BAS_SLOTS)
         ! l = -1 is the local channel, which type 1 owns.
         if (lc == -1 .or. ecpbas(SO_TYPE_OF + ecploc(iloc)*BAS_SLOTS) == 1) cycle

         atm_id = ecpbas(ATOM_OF + ecploc(iloc)*BAS_SLOTS)
         prc = atm(PTR_COORD + atm_id*ATM_SLOTS)
         rc = env(prc:prc + 2)

         if (.not. ecp_check_3c_overlap(shls, atm, bas, env, rc, &
                                        ecploc(iloc:iloc + 1), ecpbas)) cycle

         has_value = 1
         rca = rc - ri
         rcb = rc - rj
         dca = sqrt(rca(1)**2 + rca(2)**2 + rca(3)**2)
         dcb = sqrt(rcb(1)**2 + rcb(2)**2 + rcb(3)**2)

         dlc = lc*2 + 1
         lilc1 = li + lc + 1
         ljlc1 = lj + lc + 1
         d2 = lilc1*ljlc1
         d3 = lilj1*d2
         im = nfi*dlc
         mq = dlc*ljlc1

         converged = 0
         rad_all(0:nci*ncj*d3 - 1) = 0.0_dp
         nrs0 = 2**ECP_LEVEL0 - 1
         step = 2**(ECP_LEVEL_MAX - ECP_LEVEL0)
         start = step - 1
         wtscale = real(step, dp)

         do level = ECP_LEVEL0, ECP_LEVEL_MAX
            call ecp_rad_part(rur, nrs, rs, start, nrs0, step, &
                              ecploc(iloc), ecploc(iloc + 1), ecpbas, env)
            do i = 0, nrs - 1
               rur(i) = rur(i)*ws(start + i*step)*wtscale
               ! The r^lab ladder, built once and indexed by total degree.
               do lab = 1, li + lj
                  rur(nrs*lab + i) = rur(nrs*(lab - 1) + i)*rs(start + i*step)
               end do
            end do

            call ecp_type2_facs_rad(radi, ish, lc, dca, rs(start:), nrs, step, bas, env)
            call ecp_type2_facs_rad(radj, jsh, lc, dcb, rs(start:), nrs, step, bas, env)

            all_conv = .true.
            ijl = 0
            do ic = 0, nci - 1
               do jc = 0, ncj - 1
                  pradi = ic*nrs*lilc1
                  pradj = jc*nrs*ljlc1
                  do lab = 0, li + lj
                     if (converged(ijl) < 2) then
                        prur = lab*nrs
                        prad = ijl*d2
                        do i = 0, d2 - 1
                           plast(i) = rad_all(prad + i)
                           rad_all(prad + i) = rad_all(prad + i)*0.5_dp
                        end do
                        do i = 0, lilc1 - 1
                           do j = 0, ljlc1 - 1
                              s = rad_all(prad + i*ljlc1 + j)
                              do n = 0, nrs - 1
                                 s = s + rur(prur + n)*radi(pradi + n*lilc1 + i) &
                                     *radj(pradj + n*ljlc1 + j)
                              end do
                              rad_all(prad + i*ljlc1 + j) = s
                           end do
                        end do
                        pair_close = .true.
                        do i = 0, d2 - 1
                           if (.not. close_enough(plast(i), rad_all(prad + i))) then
                              pair_close = .false.
                              exit
                           end if
                        end do
                        if (pair_close) then
                           converged(ijl) = converged(ijl) + 1
                           if (converged(ijl) < 2) all_conv = .false.
                        else
                           converged(ijl) = 0
                           all_conv = .false.
                        end if
                     end if
                     ijl = ijl + 1
                  end do
               end do
            end do

            if (all_conv) exit
            nrs0 = 2**level - 1
            step = 2**(ECP_LEVEL_MAX - level)
            start = (start - 1)/2
            wtscale = wtscale*0.5_dp
         end do

         call ecp_type2_facs_ang(angi, li, lc, rca)
         call ecp_type2_facs_ang(angj, lj, lc, rcb)

         do ic = 0, nci - 1
            do jc = 0, ncj - 1
               prad = (ic*ncj + jc)*d3
               do i = 0, li
                  do j = 0, lj
                     ! buf(ljlc1, im) = rad_block(ljlc1, lilc1) x angi(lilc1, im)
                     do n = 0, im - 1
                        do mj = 0, ljlc1 - 1
                           acc = 0.0_dp
                           do mi = 0, lilc1 - 1
                              acc = acc + rad_all(prad + (i + j)*d2 + mi*ljlc1 + mj) &
                                    *angi(i*nfi*dlc*lilc1 + n*lilc1 + mi)
                           end do
                           buf(n*ljlc1 + mj) = acc
                        end do
                     end do
                     ! gctr(nfi, nfj) += fac * buf^T(nfi, mq) x angj(mq, nfj)
                     do mj = 0, nfj - 1
                        do mi = 0, nfi - 1
                           acc = 0.0_dp
                           do n = 0, mq - 1
                              acc = acc + buf(mi*mq + n) &
                                    *angj(j*nfj*dlc*ljlc1 + mj*mq + n)
                           end do
                           po = jc*nfj*di + ic*nfi + mj*di + mi
                           gctr(po) = gctr(po) + common_fac*acc
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (rs, ws, ecploc, rad_all, angi, angj, buf, rur, radi, radj)
      deallocate (plast, converged)
   end function ecp_type2_cart

   !> Both ECP terms, accumulated into the same Cartesian block
   function ecp_type_scalar_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env) &
      result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer, intent(in) :: shls(0:1), ecpbas(0:), necpbas, atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value
      integer :: h1, h2

      h1 = ecp_type1_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
      h2 = ecp_type2_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
      has_value = ior(h1, h2)
   end function ecp_type_scalar_cart

end module cint_ecp
