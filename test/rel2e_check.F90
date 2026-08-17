!
! D9 acceptance: the relativistic two-electron catalogue against the C.
!
! Twelve spinor integrals, all from auto_intor.cl's own descriptions and all
! emitted by the same backend as the non-relativistic ones:
!
!   int2e_srsr1, int2e_srsr1srsr2      sigma.r rather than sigma.p -- the
!                                      x-operations rather than the nablas
!   the four Gaunt integrals           a different operator, not just a
!                                      different bra and ket
!   the eight Breit gauge integrals    both halves, r1 and r2, which is what
!                                      the gauge term is assembled from
!   int2e_breit_r1p2, r2p2             the two hand-written kernels, the only
!                                      part of the relativistic 2e catalogue
!                                      the generator does not produce
!   the four assembled Breit terms     [gaunt - r1 + r2] / 2, which is what a
!                                      caller actually asks for
!
! The two Breit halves are the harder case in every dimension: ng raises three
! angular ceilings at once, tot_bits reaches four, and the g array carries
! sixteen intermediates.
!
! Same bar as gen2e_check -- bit-identical up to five Rys roots, 1e-11 above,
! with the root count computed from each integral's own ng increments.
!
program rel2e_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                KAPPA_OF, ANG_OF, CINTcgto_spinor
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_gen_intor4, only: int2e_srsr1_spinor, int2e_srsr1srsr2_spinor
   use cint_gen_gaunt1
   use cint_gen_breit1
   use cint_breit_gauge
   implicit none

   ! One interface shape for all of them; the name is the only difference.
   abstract interface
      function cfun(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C) result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: r
      end function
   end interface

   interface
      function c_srsr1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_srsr1_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_srsr12(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_srsr1srsr2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_g_sspssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_ssp1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_g_sspsps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_ssp1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_g_spsssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_sps1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_g_spssps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_sps1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b1_sspssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r1_ssp1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b1_sspsps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r1_ssp1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b1_spsssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r1_sps1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b1_spssps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r1_sps1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b2_sspssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r2_ssp1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b2_sspsps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r2_ssp1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b2_spsssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r2_sps1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_b2_spssps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_gauge_r2_sps1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_r1p2(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_r1p2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_r2p2(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_r2p2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_br_sspssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_ssp1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_br_sspsps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_ssp1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_br_spsssp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_sps1ssp2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_br_spssps(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_breit_sps1sps2_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   complex(c_double_complex), allocatable :: cb(:), fb(:)
   integer :: ibasis, kmode, di, dj, dk, dl, nq, dims(0:3), fshls(0:3)
   integer :: ncmp, nbad, nexact, lang, stride
   real(dp) :: worst
   character(len=80) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   ! Each integral's ng increments -- the four numbers the generator emits,
   ! and what decides where the eigensolver takes over.
   integer, parameter :: NG_SRSR1(4)   = [1, 1, 0, 0]
   integer, parameter :: NG_SRSR12(4)  = [1, 1, 1, 1]
   integer, parameter :: NG_SSPSSP(4)  = [0, 1, 0, 1]
   integer, parameter :: NG_SSPSPS(4)  = [0, 1, 1, 0]
   integer, parameter :: NG_SPSSSP(4)  = [1, 0, 0, 1]
   integer, parameter :: NG_SPSSPS(4)  = [1, 0, 1, 0]
   integer, parameter :: NG_B1_SSPSSP(4) = [1, 3, 0, 1]
   integer, parameter :: NG_B1_SSPSPS(4) = [1, 3, 1, 0]
   integer, parameter :: NG_B1_SPSSSP(4) = [2, 2, 0, 1]
   integer, parameter :: NG_B1_SPSSPS(4) = [2, 2, 1, 0]
   integer, parameter :: NG_R1P2(4)      = [2, 2, 0, 1]
   integer, parameter :: NG_R2P2(4)      = [2, 1, 0, 2]

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; nexact = 0; worst = 0.0_dp; wat = "(none)"

   do ibasis = 1, n_reference_basis
   do kmode = 0, 2
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)
      do i = 1, nbas
         lang = int(bas(ANG_OF, i))
         select case (kmode)
         case (1); bas(KAPPA_OF, i) = int(-lang - 1, c_int)
         case (2); if (lang > 0) bas(KAPPA_OF, i) = int(lang, c_int)
         end select
      end do

      allocate(fatm(0:6*8-1), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      ! These carry up to four angular increments and sixteen g
      ! intermediates, so the quartets are expensive; stride widely.
      stride = max(1, int(nbas)/2)
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
      do l = 0, nbas-1, stride
         di = int(CINTcgto_spinor(i, bas)); dj = int(CINTcgto_spinor(j, bas))
         dk = int(CINTcgto_spinor(k, bas)); dl = int(CINTcgto_spinor(l, bas))
         if (di <= 0 .or. dj <= 0 .or. dk <= 0 .or. dl <= 0) cycle
         nq = di*dj*dk*dl
         shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=l
         fshls(0)=i; fshls(1)=j; fshls(2)=k; fshls(3)=l
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=dl
         allocate(cb(nq), fb(nq))

         cb=0; fb=0
         ret = c_srsr1(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_srsr1_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "srsr1", ibasis, kmode, i, j, k, l, NG_SRSR1)

         cb=0; fb=0
         ret = c_srsr12(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_srsr1srsr2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "srsr1srsr2", ibasis, kmode, i, j, k, l, NG_SRSR12)

         cb=0; fb=0
         ret = c_g_sspssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_ssp1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "gaunt ssp1ssp2", ibasis, kmode, i, j, k, l, NG_SSPSSP)

         cb=0; fb=0
         ret = c_g_sspsps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_ssp1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "gaunt ssp1sps2", ibasis, kmode, i, j, k, l, NG_SSPSPS)

         cb=0; fb=0
         ret = c_g_spsssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_sps1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "gaunt sps1ssp2", ibasis, kmode, i, j, k, l, NG_SPSSSP)

         cb=0; fb=0
         ret = c_g_spssps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_sps1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "gaunt sps1sps2", ibasis, kmode, i, j, k, l, NG_SPSSPS)

         cb=0; fb=0
         ret = c_b1_sspssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r1_ssp1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r1 ssp1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SSPSSP)

         cb=0; fb=0
         ret = c_b1_sspsps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r1_ssp1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r1 ssp1sps2", ibasis, kmode, i, j, k, l, NG_B1_SSPSPS)

         cb=0; fb=0
         ret = c_b1_spsssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r1_sps1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r1 sps1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SPSSSP)

         cb=0; fb=0
         ret = c_b1_spssps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r1_sps1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r1 sps1sps2", ibasis, kmode, i, j, k, l, NG_B1_SPSSPS)

         cb=0; fb=0
         ret = c_b2_sspssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r2_ssp1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r2 ssp1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SSPSSP)

         cb=0; fb=0
         ret = c_b2_sspsps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r2_ssp1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r2 ssp1sps2", ibasis, kmode, i, j, k, l, NG_B1_SSPSPS)

         cb=0; fb=0
         ret = c_b2_spsssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r2_sps1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r2 sps1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SPSSSP)

         cb=0; fb=0
         ret = c_b2_spssps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_gauge_r2_sps1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit r2 sps1sps2", ibasis, kmode, i, j, k, l, NG_B1_SPSSPS)

         cb=0; fb=0
         ret = c_r1p2(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_r1p2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit_r1p2", ibasis, kmode, i, j, k, l, NG_R1P2)

         cb=0; fb=0
         ret = c_r2p2(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_r2p2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit_r2p2", ibasis, kmode, i, j, k, l, NG_R2P2)

         ! the assembled gauge term: [gaunt - r1 + r2] / 2
         cb=0; fb=0
         ret = c_br_sspssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_ssp1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit ssp1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SSPSSP)

         cb=0; fb=0
         ret = c_br_sspsps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_ssp1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit ssp1sps2", ibasis, kmode, i, j, k, l, NG_B1_SSPSPS)

         cb=0; fb=0
         ret = c_br_spsssp(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_sps1ssp2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit sps1ssp2", ibasis, kmode, i, j, k, l, NG_B1_SPSSSP)

         cb=0; fb=0
         ret = c_br_spssps(cb, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr, c_null_ptr)
         hv = int2e_breit_sps1sps2_spinor(fb, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cb, fb, nq, "breit sps1sps2", ibasis, kmode, i, j, k, l, NG_B1_SPSSPS)

         deallocate(cb, fb)
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do
   end do

   print '(A,I0)',     "  values compared      : ", ncmp
   print '(A,ES10.2)', "  worst diff / block   : ", worst
   print '(A,I0)',     "  inexact at <=5 roots : ", nexact
   print '(A,I0)',     "  over tolerance       : ", nbad
   if (nbad > 0 .or. nexact > 0) then
      print '(A,A)',   "  worst at             : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS"

contains

   integer function rys_order(s1, s2, s3, s4, inc)
      integer, intent(in) :: s1, s2, s3, s4, inc(4)
      rys_order = (int(bas(ANG_OF, s1+1)) + inc(1) + int(bas(ANG_OF, s2+1)) + inc(2) &
                 + int(bas(ANG_OF, s3+1)) + inc(3) + int(bas(ANG_OF, s4+1)) + inc(4)) / 2 + 1
   end function rys_order

   subroutine cmp(c, f, n, what, ib, km, ii, jj, kk, ll, inc)
      complex(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, km, ii, jj, kk, ll, inc(4)
      character(len=*), intent(in) :: what
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(ii, jj, kk, ll, inc)
      if (nr <= 5) then
         tol = 0.0_dp
      else
         tol = 1.0e-11_dp
      end if
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         ncmp = ncmp + 1
         r = abs(c(m) - f(m)) / scal
         if (r > worst) worst = r
         if (r > tol) then
            nbad = nbad + 1
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ib, " kmode=", km, " shls (",ii,",",jj,",",kk,",",ll,") m=", m, &
               " roots=", nr
         end if
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmp

end program rel2e_check
