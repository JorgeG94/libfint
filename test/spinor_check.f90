!
! D9, first light: the Cartesian-to-spinor transforms against the C.
!
! The ten CINTc2s_*_spinor* entry points are exported from libcint.so, so they
! can be driven directly with pseudorandom input rather than through an
! integral.  That is the right shape for this layer: it is pure linear algebra
! over a fixed coefficient table, the table was already verified bitwise in D2,
! and an integral would only obscure which of the ten disagreed.
!
! Every l from 0 to 12 -- the table's own ceiling -- and both signs of kappa,
! which select the two coefficient blocks.  The bar is bitwise: the port
! accumulates real and imaginary parts separately in the C's term order for
! exactly this reason.
!
program spinor_check
   use iso_c_binding
   use cint_const, only: dp
   use cint_cart2spinor
   implicit none

   interface
      subroutine c_bra_e1sf(gsp, nket, gcart, kappa, l) bind(C, name='CINTc2s_bra_spinor_e1sf')
         import :: c_int, c_double_complex, c_double
         complex(c_double_complex) :: gsp(*)
         integer(c_int), value :: nket, kappa, l
         real(c_double) :: gcart(*)
      end subroutine
      subroutine c_bra_sf(gsp, nket, gcart, kappa, l) bind(C, name='CINTc2s_bra_spinor_sf')
         import :: c_int, c_double_complex
         complex(c_double_complex) :: gsp(*), gcart(*)
         integer(c_int), value :: nket, kappa, l
      end subroutine
      subroutine c_bra(gsp, nket, gcart, kappa, l) bind(C, name='CINTc2s_bra_spinor')
         import :: c_int, c_double_complex
         complex(c_double_complex) :: gsp(*), gcart(*)
         integer(c_int), value :: nket, kappa, l
      end subroutine
      subroutine c_bra_si(gsp, nket, gcart, kappa, l) bind(C, name='CINTc2s_bra_spinor_si')
         import :: c_int, c_double_complex
         complex(c_double_complex) :: gsp(*), gcart(*)
         integer(c_int), value :: nket, kappa, l
      end subroutine
      subroutine c_ket(gsp, nbra, gcart, kappa, l) bind(C, name='CINTc2s_ket_spinor')
         import :: c_int, c_double_complex
         complex(c_double_complex) :: gsp(*), gcart(*)
         integer(c_int), value :: nbra, kappa, l
      end subroutine
      subroutine c_iket(gsp, nbra, gcart, kappa, l) bind(C, name='CINTc2s_iket_spinor')
         import :: c_int, c_double_complex
         complex(c_double_complex) :: gsp(*), gcart(*)
         integer(c_int), value :: nbra, kappa, l
      end subroutine
      subroutine c_ket_sf1(a, b, gc, lds, ldc, nctr, kappa, l) bind(C, name='CINTc2s_ket_spinor_sf1')
         import :: c_int, c_double_complex, c_double
         complex(c_double_complex) :: a(*), b(*)
         real(c_double) :: gc(*)
         integer(c_int), value :: lds, ldc, nctr, kappa, l
      end subroutine
      subroutine c_iket_sf1(a, b, gc, lds, ldc, nctr, kappa, l) bind(C, name='CINTc2s_iket_spinor_sf1')
         import :: c_int, c_double_complex, c_double
         complex(c_double_complex) :: a(*), b(*)
         real(c_double) :: gc(*)
         integer(c_int), value :: lds, ldc, nctr, kappa, l
      end subroutine
      subroutine c_ket_si1(a, b, gc, lds, ldc, nctr, kappa, l) bind(C, name='CINTc2s_ket_spinor_si1')
         import :: c_int, c_double_complex, c_double
         complex(c_double_complex) :: a(*), b(*)
         real(c_double) :: gc(*)
         integer(c_int), value :: lds, ldc, nctr, kappa, l
      end subroutine
      subroutine c_iket_si1(a, b, gc, lds, ldc, nctr, kappa, l) bind(C, name='CINTc2s_iket_spinor_si1')
         import :: c_int, c_double_complex, c_double
         complex(c_double_complex) :: a(*), b(*)
         real(c_double) :: gc(*)
         integer(c_int), value :: lds, ldc, nctr, kappa, l
      end subroutine
   end interface

   integer, parameter :: LMAX = 12
   integer :: l, ks, kappa, nf, nd, n, ncmp, nbad
   integer, parameter :: NKET = 5, NBRA = 5, NCTR = 3
   real(dp),    allocatable :: rc(:), rc4(:)
   complex(dp), allocatable :: zc(:), zc2(:)
   complex(dp), allocatable :: oc(:), of(:), oc2(:), of2(:)
   character(len=40) :: wat

   ncmp = 0; nbad = 0; wat = "(none)"

   do l = 0, LMAX
      nf = (l+1)*(l+2)/2
      do ks = 1, 2
         ! kappa < 0 selects the j = l + 1/2 block, kappa > 0 the j = l - 1/2
         ! one.  l = 0 has no j = l - 1/2, so nd would be zero; skip it.
         if (ks == 1) then
            kappa = -l - 1
         else
            kappa = l
         end if
         nd = cint_len_spinor_kl(kappa, l)
         if (nd <= 0) cycle

         ! --- bra_spinor_e1sf: real cartesian in, two spinor halves out
         call alloc_r(rc, nf*NKET)
         call alloc_z(oc, 2*NKET*nd); call alloc_z(of, 2*NKET*nd)
         call c_bra_e1sf(oc, NKET, rc, kappa, l)
         call c2s_bra_spinor_e1sf(of, NKET, rc, kappa, l)
         call cmpz(oc, of, 2*NKET*nd, "bra_e1sf", l, kappa)

         ! --- bra_spinor_sf: complex cartesian in
         call alloc_z(zc, nf*NKET)
         call c_bra_sf(oc, NKET, zc, kappa, l)
         call c2s_bra_spinor_sf(of, NKET, zc, kappa, l)
         call cmpz(oc, of, 2*NKET*nd, "bra_sf", l, kappa)

         ! --- bra_spinor: the full 2*nf cartesian block, one spinor out
         call alloc_z(zc, nf*2*NKET)
         call c_bra(oc, NKET, zc, kappa, l)
         call c2s_bra_spinor(of, NKET, zc, kappa, l)
         call cmpz(oc, of, NKET*nd, "bra", l, kappa)

         ! --- bra_spinor_si: two cartesian halves in
         call alloc_z(zc, 2*nf*NKET)
         call c_bra_si(oc, NKET, zc, kappa, l)
         call c2s_bra_spinor_si(of, NKET, zc, kappa, l)
         call cmpz(oc, of, NKET*nd, "bra_si", l, kappa)

         ! --- ket_spinor / iket_spinor
         call alloc_z(zc, nf*2*NBRA)
         call alloc_z(oc, NBRA*nd); call alloc_z(of, NBRA*nd)
         call c_ket(oc, NBRA, zc, kappa, l)
         call c2s_ket_spinor(of, NBRA, zc, kappa, l)
         call cmpz(oc, of, NBRA*nd, "ket", l, kappa)
         call c_iket(oc, NBRA, zc, kappa, l)
         call c2s_iket_spinor(of, NBRA, zc, kappa, l)
         call cmpz(oc, of, NBRA*nd, "iket", l, kappa)

         ! --- ket_spinor_sf1 / iket_spinor_sf1, two output arrays
         call alloc_r(rc, nf*NBRA*NCTR)
         call alloc_z(oc,  nd*NBRA*NCTR); call alloc_z(of,  nd*NBRA*NCTR)
         call alloc_z(oc2, nd*NBRA*NCTR); call alloc_z(of2, nd*NBRA*NCTR)
         call c_ket_sf1(oc, oc2, rc, NBRA, NBRA, NCTR, kappa, l)
         call c2s_ket_spinor_sf1(of, of2, rc, NBRA, NBRA, NCTR, kappa, l)
         call cmpz(oc, of, nd*NBRA*NCTR, "ket_sf1(a)", l, kappa)
         call cmpz(oc2, of2, nd*NBRA*NCTR, "ket_sf1(b)", l, kappa)
         call c_iket_sf1(oc, oc2, rc, NBRA, NBRA, NCTR, kappa, l)
         call c2s_iket_spinor_sf1(of, of2, rc, NBRA, NBRA, NCTR, kappa, l)
         call cmpz(oc, of, nd*NBRA*NCTR, "iket_sf1(a)", l, kappa)
         call cmpz(oc2, of2, nd*NBRA*NCTR, "iket_sf1(b)", l, kappa)

         ! --- ket_spinor_si1 / iket_spinor_si1, four cartesian blocks in
         call alloc_r(rc4, 4*nf*NBRA*NCTR)
         call c_ket_si1(oc, oc2, rc4, NBRA, NBRA, NCTR, kappa, l)
         call c2s_ket_spinor_si1(of, of2, rc4, NBRA, NBRA, NCTR, kappa, l)
         call cmpz(oc, of, nd*NBRA*NCTR, "ket_si1(a)", l, kappa)
         call cmpz(oc2, of2, nd*NBRA*NCTR, "ket_si1(b)", l, kappa)
         call c_iket_si1(oc, oc2, rc4, NBRA, NBRA, NCTR, kappa, l)
         call c2s_iket_spinor_si1(of, of2, rc4, NBRA, NBRA, NCTR, kappa, l)
         call cmpz(oc, of, nd*NBRA*NCTR, "iket_si1(a)", l, kappa)
         call cmpz(oc2, of2, nd*NBRA*NCTR, "iket_si1(b)", l, kappa)
      end do
   end do

   print '(A,I0)',  "  values compared : ", ncmp
   print '(A,I0)',  "  differing       : ", nbad
   if (nbad > 0) then
      print '(A,A)', "  first at        : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS (bit-identical)"

contains

   ! A cheap reproducible generator: no intrinsic RNG, so the two sides get
   ! the same input on any compiler and the test is rerunnable.
   pure real(dp) function sample(i)
      integer, intent(in) :: i
      sample = sin(real(i, dp) * 0.7139_dp) + 0.25_dp * cos(real(i, dp) * 2.311_dp)
   end function sample

   subroutine alloc_r(a, n)
      real(dp), allocatable, intent(inout) :: a(:)
      integer, intent(in) :: n
      integer :: i
      if (allocated(a)) deallocate(a)
      allocate(a(0:n-1))
      do i = 0, n-1
         a(i) = sample(i)
      end do
   end subroutine alloc_r

   subroutine alloc_z(a, n)
      complex(dp), allocatable, intent(inout) :: a(:)
      integer, intent(in) :: n
      integer :: i
      if (allocated(a)) deallocate(a)
      allocate(a(0:n-1))
      do i = 0, n-1
         a(i) = cmplx(sample(i), sample(i+7919), dp)
      end do
   end subroutine alloc_z

   subroutine cmpz(a, b, n, what, ll, kk)
      complex(dp), intent(in) :: a(0:), b(0:)
      integer, intent(in) :: n, ll, kk
      character(len=*), intent(in) :: what
      integer :: m
      do m = 0, n-1
         ncmp = ncmp + 1
         if (a(m) /= b(m)) then
            nbad = nbad + 1
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0)') trim(what), &
               " l=", ll, " kappa=", kk, " i=", m
         end if
      end do
   end subroutine cmpz

end program spinor_check
