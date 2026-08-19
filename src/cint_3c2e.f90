!
! Three- and two-centre two-electron integrals: (ij|k) and (i|k).
!
! Ported from src/cint3c2e.c and src/cint2c2e.c.  These are the density-fitting
! integrals; between them and int2e they cover every two-electron path a
! Coulomb-metric code takes.
!
! Both reuse the four-centre engine wholesale.  The trick is entirely in the
! envs setup (cint_init_int3c2e_envvars / cint_init_int2c2e_envvars, in
! cint_g2e.f90): the auxiliary shell goes into the *l* slot with k left at
! zero, and the missing centres get zero exponents and coincident coordinates
! so the Rys quadrature degenerates to the right thing.  Downstream, g0_2e and
! the 2D-to-4D transfers are the same code the four-centre case runs.
!
! As in cint_2e.f90, only the unoptimised loops are here; the C's five
! contraction-count specialisations produce the same numbers faster.
!
module cint_3c2e
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_ensure_pd, ws_alloc_d, ws_alloc_i, ws_mark, ws_rewind, &
                             ws_opt_log_maxc, ws_opt_non0, ws_opt_idx
   use cint_screen,    only: pair_data, cint_set_pairdata, &
                             cint_log_max_pgto_coeff, cint_non0coeff_byshell
   use cint_g1e,       only: cint_prim_to_ctr_0, cint_prim_to_ctr_1, &
                             cint_g1e_index_xyz
   use cint_g2e
   ! The optimizer builder comes in with the bare `use cint_g2e` above; naming
   ! it again here adds no visibility and ifx says so, remark #6536.
   use cint_2e,        only: cint_gout2e
   use cint_1e,        only: int1e_cache_size, apply_c2s_cart_1e, &
                             apply_c2s_sph_1e, apply_c2s_dset0
   use cint_bas,       only: NPRIM_OF, PTR_EXP, PTR_COEFF
   use cint_blas,      only: cint_dmat_transpose, cint_dplus_transpose
   use cint_cart2sph,  only: cint_c2s_bra_sph, cint_c2s_ket_sph, &
                                RESULT_IN_GCART, RESULT_IN_GSPH
   implicit none
   private

   public :: cint_3c2e_loop_nopt, cint_3c2e_drv
   public :: cint_2c2e_loop_nopt, cint_2c2e_drv
   public :: int3c2e_cart, int3c2e_sph, int2c2e_cart, int2c2e_sph
   ! cint_3c1e reuses the three-index output layer verbatim
   public :: apply_c2s_cart_3c2e1, apply_c2s_sph_3c2e1, apply_c2s_dset0_3c
   public :: int3c2e_optimizer
   public :: int3c2e_ssc_optimizer
   public :: int2c2e_optimizer

   real(dp), parameter :: R_GUESS = 8.0_dp

contains

   ! ---- three centre --------------------------------------------------

   function cint_3c2e_loop_nopt(gctr, gctroff, envs, ws, empty) result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET for the block pointers bound below
      type(cint_ws), intent(inout), target :: ws
      logical,  intent(inout) :: empty
      logical :: has_value

      integer :: i_sh, j_sh, k_sh, i_ctr, j_ctr, k_ctr
      logical :: use_opt
      logical :: idx_cached
      integer :: i_prim, j_prim, k_prim
      integer :: ai_o, aj_o, ak_o, ci_o, cj_o, ck_o
      integer :: n_comp, nc, nf, lij
      integer :: leng, lenk, lenj, leni, len0
      integer :: og, ogout, ogctri, ogctrj, ogctrk
      integer :: oidx, olm, onon, osrt
      integer :: ip, jp, kp, kij
      real(dp) :: expcutoff, rr_ij, cutoff, fac1i, fac1j, fac1k
      real(dp) :: omega, omega2, theta, aij, dist_ij, rkl(0:2)
      ! Same alias chain as the four-centre loop, one level shorter: the C's
      ! `#define gctrm gctr / mempty empty / m_ctr n_comp` makes the top link
      ! terminate on the caller's flag when n_comp == 1.
      integer, parameter :: E_I = 0, E_J = 1, E_K = 2, E_G = 3, E_M = 4
      logical :: flag(0:4)
      integer :: ii, ij, ik, ig, im
      ! The same hoists the four-centre loop makes, and for the same reason:
      ! this is a three-deep primitive nest whose innermost body is one
      ! cint_g0_2e and one gout, so every descriptor rebuilt per iteration
      ! shows up.  See the note in cint_2e_loop_nopt.
      real(dp), pointer, contiguous :: aip(:), ajp(:), akp(:)
      real(dp), pointer, contiguous :: cip(:), cjp(:), ckp(:)
      type(pair_data), pointer, contiguous :: pdp(:)
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:)
      integer,  pointer, contiguous :: nonp(:)

      has_value = .false.
      i_sh = envs%shls(0); j_sh = envs%shls(1); k_sh = envs%shls(2)
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1); k_ctr = envs%x_ctr(2)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      j_prim = bas_of(envs, NPRIM_OF, j_sh)
      k_prim = bas_of(envs, NPRIM_OF, k_sh)
      ai_o = bas_of(envs, PTR_EXP, i_sh); aj_o = bas_of(envs, PTR_EXP, j_sh)
      ak_o = bas_of(envs, PTR_EXP, k_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh); cj_o = bas_of(envs, PTR_COEFF, j_sh)
      ck_o = bas_of(envs, PTR_COEFF, k_sh)

      expcutoff = envs%expcutoff
      rr_ij = sum(envs%rirj * envs%rirj)
      nf = envs%nf
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 3, envs)

      call ws_alloc_d(ws, i_prim + j_prim, olm)
      if (use_opt) then
         call ws_opt_log_maxc(ws, olm,          i_sh, i_prim)
         call ws_opt_log_maxc(ws, olm + i_prim, j_sh, j_prim)
      else
      call cint_log_max_pgto_coeff(ws%d(olm:),        envs%env, ci_o, i_prim, i_ctr)
      call cint_log_max_pgto_coeff(ws%d(olm+i_prim:), envs%env, cj_o, j_prim, j_ctr)
      end if

      call ws_ensure_pd(ws, i_prim*j_prim)
      if (cint_set_pairdata(ws%p, envs%env, ai_o, envs%env, aj_o, envs%ri, envs%rj, &
                            ws%d(olm:), ws%d(olm+i_prim:), envs%li_ceil, envs%lj_ceil, &
                            i_prim, j_prim, rr_ij, expcutoff, &
                            envs%env(PTR_RANGE_OMEGA))) then
         return
      end if

      omega = envs%env(PTR_RANGE_OMEGA)
      if (omega < 0.0_dp .and. envs%rys_order > 1) then
         omega2 = omega * omega
         lij = envs%li_ceil + envs%lj_ceil
         if (lij > 0) then
            dist_ij = sqrt(rr_ij)
            aij = envs%env(ai_o + i_prim - 1) + envs%env(aj_o + j_prim - 1)
            theta = omega2 / (omega2 + aij)
            expcutoff = expcutoff + lij * &
               log((dist_ij + theta*R_GUESS + 1.0_dp) / (dist_ij + 1.0_dp))
         end if
         if (envs%lk_ceil > 0) then
            theta = omega2 / (omega2 + envs%env(ak_o + k_prim - 1))
            expcutoff = expcutoff + envs%lk_ceil * log(theta*R_GUESS + 1.0_dp)
         end if
      end if

      call ws_alloc_i(ws, nf*3, oidx)
      ! NOT `use_opt .and. ws_opt_idx(...)`: Fortran does not
      ! guarantee short-circuit evaluation, so the compiler is free to
      ! call ws_opt_idx even when use_opt is false -- and that
      ! dereferences ws%opt, which is null precisely when use_opt is
      ! false.  gfortran short-circuits at -O2 and does not at -O0, so
      ! this segfaulted only in a debug build.
      idx_cached = .false.
      if (use_opt) idx_cached = ws_opt_idx(ws, oidx, cint_opt_idx_key(envs, 3), nf*3)
      if (.not. idx_cached) then
         call cint_g2e_index_xyz(ws%i(oidx:), envs)
      end if

      call ws_alloc_i(ws, i_prim + j_prim + k_prim, onon)
      call ws_alloc_i(ws, i_prim*i_ctr + j_prim*j_ctr + k_prim*k_ctr, osrt)
      if (use_opt) then
         call ws_opt_non0(ws, onon, osrt, i_sh, i_prim, i_ctr, 0, 0)
         call ws_opt_non0(ws, onon, osrt, j_sh, j_prim, j_ctr, i_prim, i_prim*i_ctr)
         call ws_opt_non0(ws, onon, osrt, k_sh, k_prim, k_ctr, i_prim+j_prim, &
                          i_prim*i_ctr+j_prim*j_ctr)
      else
      call cint_non0coeff_byshell(ws%i(osrt:), ws%i(onon:), envs%env, ci_o, i_prim, i_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr:), ws%i(onon+i_prim:), &
                                  envs%env, cj_o, j_prim, j_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr+j_prim*j_ctr:), &
                                  ws%i(onon+i_prim+j_prim:), envs%env, ck_o, k_prim, k_ctr)
      end if

      nc = i_ctr * j_ctr * k_ctr
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenk = nf * nc * n_comp
      lenj = nf * i_ctr * j_ctr * n_comp
      leni = nf * i_ctr * n_comp
      len0 = nf * n_comp

      call ws_alloc_d(ws, leng, og)
      flag = .true.
      flag(E_M) = empty
      ii = E_I; ij = E_J; ik = E_K; ig = E_G; im = E_M
      if (n_comp == 1) then
         ogctrk = gctroff; ik = im
      else
         call ws_alloc_d(ws, lenk, ogctrk)
      end if
      if (k_ctr == 1) then
         ogctrj = ogctrk; ij = ik
      else
         call ws_alloc_d(ws, lenj, ogctrj)
      end if
      if (j_ctr == 1) then
         ogctri = ogctrj; ii = ij
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (i_ctr == 1) then
         ogout = ogctri; ig = ii
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      aip(0:) => envs%env(ai_o:); ajp(0:) => envs%env(aj_o:)
      akp(0:) => envs%env(ak_o:)
      cip(0:) => envs%env(ci_o:); cjp(0:) => envs%env(cj_o:)
      ckp(0:) => envs%env(ck_o:)
      pdp(0:)   => ws%p
      gp(0:)    => ws%d(og:og+leng-1)
      goutp(0:) => ws%d(ogout:ogout+len0-1)
      idxp(0:)  => ws%i(oidx:oidx+nf*3-1)
      nonp(0:)  => ws%i(onon:onon+i_prim+j_prim+k_prim-1)

      ! the auxiliary centre is its own charge centre; rkl never moves
      rkl = envs%rk

      ! prim2ctr is expanded at all three sites rather than called.  gfortran
      ! will not inline it (it takes ws by reference and writes flag), and at
      ! the innermost site the call was costing more than the two-line body.
      do kp = 0, k_prim - 1
         envs%ak = akp(kp)
         if (k_ctr == 1) then
            fac1k = envs%common_factor * ckp(kp)
         else
            fac1k = envs%common_factor
            flag(ij) = .true.
         end if

         do jp = 0, j_prim - 1
            envs%aj = ajp(jp)
            if (j_ctr == 1) then
               fac1j = fac1k * cjp(jp)
            else
               fac1j = fac1k
               flag(ii) = .true.
            end if

            do ip = 0, i_prim - 1
               kij = jp * i_prim + ip
               if (pdp(kij)%cceij > expcutoff) cycle
               envs%ai = aip(ip)
               cutoff = expcutoff - pdp(kij)%cceij
               if (i_ctr == 1) then
                  fac1i = fac1j * cip(ip) * pdp(kij)%eij
               else
                  fac1i = fac1j * pdp(kij)%eij
               end if
               envs%fac = fac1i
               if (cint_g0_2e(gp, pdp(kij)%rij, rkl, cutoff, envs) /= 0) then
                  call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig)))
                  if (i_ctr > 1) then
                     if (flag(ii)) then
                        call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, &
                                                ci_o+ip, len0, i_prim, i_ctr)
                     else
                        call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, &
                                                ci_o+ip, len0, i_prim, i_ctr, &
                                                nonp(ip), ws%i, osrt + ip*i_ctr)
                     end if
                  end if
                  flag(ii) = .false.
               end if
            end do

            if (.not. flag(ii)) then
               if (j_ctr > 1) then
                  if (flag(ij)) then
                     call cint_prim_to_ctr_0(ws%d, ogctrj, ogctri, envs%env, &
                                             cj_o+jp, leni, j_prim, j_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctrj, ogctri, envs%env, &
                                             cj_o+jp, leni, j_prim, j_ctr, &
                                             nonp(i_prim+jp), ws%i, &
                                             osrt + i_prim*i_ctr + jp*j_ctr)
                  end if
               end if
               flag(ij) = .false.
            end if
         end do

         if (.not. flag(ij)) then
            if (k_ctr > 1) then
               if (flag(ik)) then
                  call cint_prim_to_ctr_0(ws%d, ogctrk, ogctrj, envs%env, &
                                          ck_o+kp, lenj, k_prim, k_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctrk, ogctrj, envs%env, &
                                          ck_o+kp, lenj, k_prim, k_ctr, &
                                          nonp(i_prim+j_prim+kp), ws%i, &
                                          osrt + i_prim*i_ctr + j_prim*j_ctr + kp*k_ctr)
               end if
            end if
            flag(ik) = .false.
         end if
      end do

      if (n_comp > 1 .and. .not. flag(ik)) then
         if (flag(im)) then
            call cint_dmat_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
            flag(im) = .false.
         else
            call cint_dplus_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
         end if
      end if
      empty = flag(im)
      has_value = .not. empty
   end function cint_3c2e_loop_nopt

   ! ---- two centre ----------------------------------------------------

   function cint_2c2e_loop_nopt(gctr, gctroff, envs, ws, empty) result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET for the block pointers bound below
      type(cint_ws), intent(inout), target :: ws
      logical,  intent(inout) :: empty
      logical :: has_value

      integer :: i_sh, k_sh, i_ctr, k_ctr, i_prim, k_prim
      logical :: use_opt
      logical :: idx_cached
      integer :: ai_o, ak_o, ci_o, ck_o
      integer :: n_comp, nc, nf, leng, lenk, leni, len0
      integer :: og, ogout, ogctri, ogctrk, oidx, onon, osrt
      integer :: ip, kp
      real(dp) :: expcutoff, fac1i, fac1k
      integer, parameter :: E_I = 0, E_K = 1, E_G = 2, E_M = 3
      logical :: flag(0:3)
      integer :: ii, ik, ig, im
      ! Same hoists as the three-centre loop above.
      real(dp), pointer, contiguous :: aip(:), akp(:), cip(:), ckp(:)
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:), nonp(:)

      has_value = .false.
      i_sh = envs%shls(0); k_sh = envs%shls(1)
      i_ctr = envs%x_ctr(0); k_ctr = envs%x_ctr(1)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      k_prim = bas_of(envs, NPRIM_OF, k_sh)
      ai_o = bas_of(envs, PTR_EXP, i_sh); ak_o = bas_of(envs, PTR_EXP, k_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh); ck_o = bas_of(envs, PTR_COEFF, k_sh)

      expcutoff = envs%expcutoff
      nf = envs%nf
      ! note: ncomp_tensor alone, not the product -- the C's CINT2c2e_loop_nopt
      ! reads only that field, though its driver uses the full product.
      n_comp = envs%ncomp_tensor

      nc = i_ctr * k_ctr
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenk = nf * nc * n_comp
      leni = nf * i_ctr * n_comp
      len0 = nf * n_comp

      call ws_alloc_d(ws, leng, og)
      flag = .true.
      flag(E_M) = empty
      ii = E_I; ik = E_K; ig = E_G; im = E_M
      if (n_comp == 1) then
         ogctrk = gctroff; ik = im
      else
         call ws_alloc_d(ws, lenk, ogctrk)
      end if
      if (k_ctr == 1) then
         ogctri = ogctrk; ii = ik
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (i_ctr == 1) then
         ogout = ogctri; ig = ii
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 2, envs)

      ! int2c2e goes through the 1e index map, not the 2e one: with j and l
      ! absent the g array has the 1e shape.
      call ws_alloc_i(ws, nf*3, oidx)
      ! NOT `use_opt .and. ws_opt_idx(...)`: Fortran does not
      ! guarantee short-circuit evaluation, so the compiler is free to
      ! call ws_opt_idx even when use_opt is false -- and that
      ! dereferences ws%opt, which is null precisely when use_opt is
      ! false.  gfortran short-circuits at -O2 and does not at -O0, so
      ! this segfaulted only in a debug build.
      idx_cached = .false.
      if (use_opt) idx_cached = ws_opt_idx(ws, oidx, cint_opt_idx_key(envs, 2), nf*3)
      if (.not. idx_cached) then
         call cint_g1e_index_xyz(ws%i(oidx:), envs)
      end if

      call ws_alloc_i(ws, i_prim + k_prim, onon)
      call ws_alloc_i(ws, i_prim*i_ctr + k_prim*k_ctr, osrt)
      if (i_ctr > 1) then
         if (use_opt) then
            call ws_opt_non0(ws, onon, osrt, i_sh, i_prim, i_ctr, 0, 0)
         else
         call cint_non0coeff_byshell(ws%i(osrt:), ws%i(onon:), envs%env, ci_o, i_prim, i_ctr)
         end if
      end if
      if (k_ctr > 1) then
         if (use_opt) then
            call ws_opt_non0(ws, onon, osrt, k_sh, k_prim, k_ctr, i_prim, i_prim*i_ctr)
         else
         call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr:), ws%i(onon+i_prim:), &
                                     envs%env, ck_o, k_prim, k_ctr)
         end if
      end if

      aip(0:) => envs%env(ai_o:); akp(0:) => envs%env(ak_o:)
      cip(0:) => envs%env(ci_o:); ckp(0:) => envs%env(ck_o:)
      gp(0:)    => ws%d(og:og+leng-1)
      goutp(0:) => ws%d(ogout:ogout+len0-1)
      idxp(0:)  => ws%i(oidx:oidx+nf*3-1)
      nonp(0:)  => ws%i(onon:onon+i_prim+k_prim-1)

      do kp = 0, k_prim - 1
         envs%ak = akp(kp)
         envs%al = 0.0_dp
         if (k_ctr == 1) then
            fac1k = envs%common_factor * ckp(kp)
         else
            fac1k = envs%common_factor
            flag(ii) = .true.
         end if
         do ip = 0, i_prim - 1
            envs%ai = aip(ip)
            envs%aj = 0.0_dp
            if (i_ctr == 1) then
               fac1i = fac1k * cip(ip)
            else
               fac1i = fac1k
            end if
            envs%fac = fac1i
            if (cint_g0_2e(gp, envs%ri, envs%rk, expcutoff, envs) /= 0) then
               call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig)))
               if (i_ctr > 1) then
                  if (flag(ii)) then
                     call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, &
                                             ci_o+ip, len0, i_prim, i_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, &
                                             ci_o+ip, len0, i_prim, i_ctr, &
                                             nonp(ip), ws%i, osrt + ip*i_ctr)
                  end if
               end if
               flag(ii) = .false.
            end if
         end do
         if (.not. flag(ii)) then
            if (k_ctr > 1) then
               if (flag(ik)) then
                  call cint_prim_to_ctr_0(ws%d, ogctrk, ogctri, envs%env, &
                                          ck_o+kp, leni, k_prim, k_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctrk, ogctri, envs%env, &
                                          ck_o+kp, leni, k_prim, k_ctr, &
                                          nonp(i_prim+kp), ws%i, &
                                          osrt + i_prim*i_ctr + kp*k_ctr)
               end if
            end if
            flag(ik) = .false.
         end if
      end do

      if (n_comp > 1 .and. .not. flag(ik)) then
         if (flag(im)) then
            call cint_dmat_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
            flag(im) = .false.
         else
            call cint_dplus_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
         end if
      end if
      empty = flag(im)
      has_value = .not. empty
   end function cint_2c2e_loop_nopt

   ! Shared by both loops; the same body as cint_2e's, kept local rather than
   ! exported because it is four lines and depends on the workspace type.

   ! ---- drivers -------------------------------------------------------

   pure subroutine int3c2e_cache_size(envs, nd, ni)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(out) :: nd, ni
      integer :: ip, jp, kp, ic, jc, kc, nc, n_comp, nf
      integer :: leng, lenk, lenj, leni, len0, buflen

      ip = bas_of(envs, NPRIM_OF, envs%shls(0))
      jp = bas_of(envs, NPRIM_OF, envs%shls(1))
      kp = bas_of(envs, NPRIM_OF, envs%shls(2))
      ic = envs%x_ctr(0); jc = envs%x_ctr(1); kc = envs%x_ctr(2)
      nf = envs%nf
      nc = ic*jc*kc
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenk = nf * nc * n_comp
      lenj = nf * ic*jc * n_comp
      leni = nf * ic * n_comp
      len0 = nf * n_comp
      ! the three c2s_sph_3c2e1 scratch blocks
      buflen = envs%nfi * envs%nfk * (2*envs%j_l + 1)

      nd = nf*nc*n_comp + ip+jp &
         + leng + lenk + lenj + leni + len0 + 3*buflen + 64
      ni = nf*3 + ip+jp+kp + ip*ic + jp*jc + kp*kc + 64
   end subroutine int3c2e_cache_size

   function cint_3c2e_drv(out, dims, envs, ws, c2s_kind) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind
      logical :: has_value
      integer :: dmark, imark

      integer :: nc, n_comp, ogctr, n, nout, counts(0:3), nd, ni
      logical :: empty, dummy

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      call int3c2e_cache_size(envs, nd, ni)
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! NOT zeroed.  No C driver zeroes gctr -- CINT1e_drv, CINT3c2e_drv,
      ! CINT2c2e_drv, CINT3c1e_drv and CINT1e_grids_drv all just take the
      ! buffer and thread an empty flag through the loop, so the first
      ! contraction sets rather than accumulates.  The port threads the same
      ! flag, so the memset was pure extra writes.

      empty = .true.
      dummy = cint_3c2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      select case (c2s_kind)
      case (C2S_SPH_3C2E1)
         counts(0) = (envs%i_l*2+1) * envs%x_ctr(0)
         counts(1) = (envs%j_l*2+1) * envs%x_ctr(1)
         counts(2) = (envs%k_l*2+1) * envs%x_ctr(2)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfj * envs%x_ctr(1)
         counts(2) = envs%nfk * envs%x_ctr(2)
      end select
      counts(3) = 1
      nout = dims(0) * dims(1) * dims(2)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_CART_3C2E1)
               call apply_c2s_cart_3c2e1(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs)
            case (C2S_SPH_3C2E1)
               call apply_c2s_sph_3c2e1(out(nout*n:), ogctr+nc*n, dims, envs, ws)
            case default
               error stop "cint_3c2e_drv: spinor c2s is D9"
            end select
         end do
      else
         do n = 0, n_comp - 1
            call apply_c2s_dset0_3c(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_3c2e_drv

   function cint_2c2e_drv(out, dims, envs, ws, c2s_kind) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind
      logical :: has_value
      integer :: dmark, imark

      integer :: nc, n_comp, ogctr, n, nout, counts(0:3), nd, ni
      logical :: empty, dummy

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1)
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor

      ! int1e_cache_size is exactly right here: it sizes leng from envs%g_size
      ! (already the 2e one) and reads shls(1)/x_ctr(1), which for int2c2e are
      ! the k shell -- the same substitution the envs setup makes.
      call int1e_cache_size(envs, nd, ni)
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! NOT zeroed.  No C driver zeroes gctr -- CINT1e_drv, CINT3c2e_drv,
      ! CINT2c2e_drv, CINT3c1e_drv and CINT1e_grids_drv all just take the
      ! buffer and thread an empty flag through the loop, so the first
      ! contraction sets rather than accumulates.  The port threads the same
      ! flag, so the memset was pure extra writes.

      empty = .true.
      dummy = cint_2c2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      select case (c2s_kind)
      case (C2S_SPH_2C2E1)
         counts(0) = (envs%i_l*2+1) * envs%x_ctr(0)
         counts(1) = (envs%k_l*2+1) * envs%x_ctr(1)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfk * envs%x_ctr(1)
      end select
      counts(2) = 1
      counts(3) = 1
      nout = dims(0) * dims(1)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_CART_2C2E1)
               call apply_c2s_cart_1e(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs)
            case (C2S_SPH_2C2E1)
               call apply_c2s_sph_1e(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs, ws)
            case default
               error stop "cint_2c2e_drv: spinor c2s is D9"
            end select
         end do
      else
         do n = 0, n_comp - 1
            call apply_c2s_dset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_2c2e_drv

   ! ---- output layout -------------------------------------------------

   ! The three-index analogue of dcopy_iklj: the g array is laid out (i,k,j)
   ! and the caller wants (i,j,k).
   pure subroutine dcopy_ikj(fijk, foff, gctr, goff, ni, nj, mi, mj, mk)
      real(dp), intent(inout) :: fijk(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: foff, goff, ni, nj, mi, mj, mk
      integer :: i, j, k, nij, mik
      nij = ni*nj; mik = mi*mk
      do k = 0, mk - 1
         do j = 0, mj - 1
            do i = 0, mi - 1
               fijk(foff + k*nij + ni*j + i) = gctr(goff + k*mi + mik*j + i)
            end do
         end do
      end do
   end subroutine dcopy_ikj

   subroutine apply_c2s_cart_3c2e1(fijk, gctr, dims, envs)
      real(dp), intent(inout) :: fijk(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: ic, jc, kc, ni, nj, ofj, ofk, base, gb
      ni = dims(0); nj = dims(1)
      ofj = ni * envs%nfj
      ofk = ni * nj * envs%nfk
      gb = 0
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         base = ofk*kc + ofj*jc + envs%nfi*ic
         call dcopy_ikj(fijk, base, gctr, gb, ni, nj, envs%nfi, envs%nfj, envs%nfk)
         gb = gb + envs%nf
      end do
      end do
      end do
   end subroutine apply_c2s_cart_3c2e1

   ! ONE ARRAY, TWO OFFSETS.  Both the source and the destination are the
   ! workspace, and passing the whole of it as an intent(in) dummy and again
   ! as an intent(inout) one is the aliasing violation this port keeps
   ! rediscovering.  The slices do not overlap; saying so with offsets is
   ! what makes that a fact rather than a hope.
   subroutine sph2e_inner(buf, soff, coff, l, nbra, ncall, sizsph, sizcart, loc)
      real(dp), intent(inout) :: buf(0:)
      integer,  intent(in)    :: soff, coff, l, nbra, ncall, sizsph, sizcart
      integer,  intent(out)   :: loc
      integer :: n, r
      if (l <= 1) then
         loc = RESULT_IN_GCART
         return
      end if
      do n = 0, ncall - 1
         r = cint_c2s_ket_sph(buf(soff + n*sizsph:), buf(coff + n*sizcart:), &
                              nbra, nbra, l)
      end do
      loc = RESULT_IN_GSPH
   end subroutine sph2e_inner

   ! As in cint_2e.f90: an offset into the workspace, not a section of it,
   ! because sph2e_inner indexes the same buffer.
   subroutine apply_c2s_sph_3c2e1(out, goff, dims, envs, ws)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: goff, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      type(cint_ws), intent(inout) :: ws

      integer :: i_l, j_l, k_l, di, dj, dk, ni, nj, nfi, nfk, nfik

      integer :: icp
      integer :: ofj, ofk, ic, jc, kc, buflen
      integer :: o1, o2, o3, gb, base, cur, loc

      i_l = envs%i_l; j_l = envs%j_l; k_l = envs%k_l
      di = i_l*2+1; dj = j_l*2+1; dk = k_l*2+1
      ni = dims(0); nj = dims(1)
      nfi = envs%nfi; nfk = envs%nfk; nfik = nfi*nfk
      ofj = ni*dj; ofk = ni*nj*dk
      buflen = nfik*dj
      call ws_alloc_d(ws, buflen, o1)
      call ws_alloc_d(ws, buflen, o2)
      call ws_alloc_d(ws, buflen, o3)

      gb = goff
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         ! j, then k, then i -- the order the (i,k,j) layout makes cheapest.
         ! Each transform is the identity for l <= 1 and then reports its
         ! result as still being at the source; normalise into o1 so the rest
         ! of the chain has one place to look, as cint_2e.f90 does.
         loc = cint_c2s_ket_sph(ws%d(o1:), ws%d(gb:), nfik, nfik, j_l)
         if (loc == RESULT_IN_GCART) then
            do icp = 0, buflen - 1
               ws%d(o1+icp) = ws%d(gb+icp)
            end do
         end if

         call sph2e_inner(ws%d, o2, o1, k_l, nfi, dj, nfi*dk, nfik, loc)
         if (loc == RESULT_IN_GCART) then
            cur = o1
         else
            cur = o2
         end if

         base = ofk*kc + ofj*jc + di*ic
         loc = cint_c2s_bra_sph(ws%d(o3:), dk*dj, ws%d(cur:), i_l)
         if (loc == RESULT_IN_GCART) then
            call dcopy_ikj(out, base, ws%d, cur, ni, nj, di, dj, dk)
         else
            call dcopy_ikj(out, base, ws%d, o3, ni, nj, di, dj, dk)
         end if
         gb = gb + envs%nf
      end do
      end do
      end do
   end subroutine apply_c2s_sph_3c2e1

   pure subroutine apply_c2s_dset0_3c(out, dims, counts)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), counts(0:)
      integer :: i, j, k, ni, nij
      ni = dims(0); nij = ni*dims(1)
      do k = 0, counts(2) - 1
      do j = 0, counts(1) - 1
      do i = 0, counts(0) - 1
         out(k*nij + j*ni + i) = 0.0_dp
      end do
      end do
      end do
   end subroutine apply_c2s_dset0_3c

   ! ---- entry points --------------------------------------------------

   function int3c2e_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int3c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_3c2e_drv(out, dims, envs, ws, C2S_CART_3C2E1)
   end function int3c2e_cart

   function int3c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int3c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_3c2e_drv(out, dims, envs, ws, C2S_SPH_3C2E1)
   end function int3c2e_sph

   function int2c2e_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int2c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2c2e_drv(out, dims, envs, ws, C2S_CART_2C2E1)
   end function int2c2e_cart

   function int2c2e_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int2c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2c2e_drv(out, dims, envs, ws, C2S_SPH_2C2E1)
   end function int2c2e_sph

   ! ---- optimizers ------------------------------------------------------

   subroutine int3c2e_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int3c2e_optimizer

   subroutine int3c2e_ssc_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call int3c2e_optimizer(opt, atm, natm, bas, nbas, env)
   end subroutine int3c2e_ssc_optimizer

   subroutine int2c2e_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2c2e_optimizer


end module cint_3c2e
