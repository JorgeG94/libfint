!
! The two-electron recursion operations: nabla and x on each of the four
! centres.
!
! Ported from src/g2e.c:4574-4925.  These are what the G2E_* macros expand to,
! and so what every generated derivative or relativistic 2e integral is built
! out of -- the four-centre counterparts of cint_nabla1?_1e and cint_x1?_1e in
! cint_g1e.f90.
!
! ONE ARRAY, TWO OFFSETS, for the same reason as the 1e ones: the C passes f
! and g as distinct pointers into the same buffer, which is legal there and is
! illegal aliasing in Fortran.  Associating one actual with an intent(in) and
! an intent(inout) dummy lets the compiler assume they cannot overlap, and
! gfortran duly miscompiled the 1e case when it was written that way.  Saying
! it is one array is both true and safe.
!
! The x-operations take their origin as an explicit vector rather than reading
! it from envs, because the G2E_R0? and G2E_RC? macros differ only in which
! vector they pass -- envs->r? for the first, a locally computed dr? for the
! second.
!
module cint_g2e_ops
   use cint_const, only: dp
   use cint_envs
   implicit none
   private

   public :: cint_nabla1i_2e, cint_nabla1j_2e, cint_nabla1k_2e, cint_nabla1l_2e
   public :: cint_x1i_2e, cint_x1j_2e, cint_x1k_2e, cint_x1l_2e

contains

   ! ( \nabla i j | kl )
   subroutine cint_nabla1i_2e(g, foff, goff, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, l, n, ptr, di, dk, dl, dj, nroots, gs
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: ai2
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      nroots = envs%nrys_roots; gs = envs%g_size
      ai2 = -2.0_dp * envs%ai
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
      do l = 0, ll
      do k = 0, lk
         ptr = dj*j + dl*l + dk*k
         ! f(...,0,...) = -2*ai*g(...,1,...)
         do n = ptr, ptr + nroots - 1
            g(fx+n) = ai2 * g(gx+di+n)
            g(fy+n) = ai2 * g(gy+di+n)
            g(fz+n) = ai2 * g(gz+di+n)
         end do
         ptr = ptr + di
         ! f(...,i,...) = i*g(...,i-1,...) - 2*ai*g(...,i+1,...)
         do i = 1, li
            do n = ptr, ptr + nroots - 1
               g(fx+n) = i*g(gx-di+n) + ai2*g(gx+di+n)
               g(fy+n) = i*g(gy-di+n) + ai2*g(gy+di+n)
               g(fz+n) = i*g(gz-di+n) + ai2*g(gz+di+n)
            end do
            ptr = ptr + di
         end do
      end do
      end do
      end do
   end subroutine cint_nabla1i_2e

   ! ( i \nabla j | kl )
   subroutine cint_nabla1j_2e(g, foff, goff, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, l, n, ptr, di, dk, dl, dj, nroots, gs
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: aj2
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      nroots = envs%nrys_roots; gs = envs%g_size
      aj2 = -2.0_dp * envs%aj
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do l = 0, ll
      do k = 0, lk
         ptr = dl*l + dk*k
         do i = 0, li
            do n = ptr, ptr + nroots - 1
               g(fx+n) = aj2 * g(gx+dj+n)
               g(fy+n) = aj2 * g(gy+dj+n)
               g(fz+n) = aj2 * g(gz+dj+n)
            end do
            ptr = ptr + di
         end do
      end do
      end do
      do j = 1, lj
         do l = 0, ll
         do k = 0, lk
            ptr = dj*j + dl*l + dk*k
            do i = 0, li
               do n = ptr, ptr + nroots - 1
                  g(fx+n) = j*g(gx-dj+n) + aj2*g(gx+dj+n)
                  g(fy+n) = j*g(gy-dj+n) + aj2*g(gy+dj+n)
                  g(fz+n) = j*g(gz-dj+n) + aj2*g(gz+dj+n)
               end do
               ptr = ptr + di
            end do
         end do
         end do
      end do
   end subroutine cint_nabla1j_2e

   ! ( ij | \nabla k l )
   subroutine cint_nabla1k_2e(g, foff, goff, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, l, n, ptr, di, dk, dl, dj, nroots, gs
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: ak2
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      nroots = envs%nrys_roots; gs = envs%g_size
      ak2 = -2.0_dp * envs%ak
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
      do l = 0, ll
         ptr = dj*j + dl*l
         do i = 0, li
            do n = ptr, ptr + nroots - 1
               g(fx+n) = ak2 * g(gx+dk+n)
               g(fy+n) = ak2 * g(gy+dk+n)
               g(fz+n) = ak2 * g(gz+dk+n)
            end do
            ptr = ptr + di
         end do
         do k = 1, lk
            ptr = dj*j + dl*l + dk*k
            do i = 0, li
               do n = ptr, ptr + nroots - 1
                  g(fx+n) = k*g(gx-dk+n) + ak2*g(gx+dk+n)
                  g(fy+n) = k*g(gy-dk+n) + ak2*g(gy+dk+n)
                  g(fz+n) = k*g(gz-dk+n) + ak2*g(gz+dk+n)
               end do
               ptr = ptr + di
            end do
         end do
      end do
      end do
   end subroutine cint_nabla1k_2e

   ! ( ij | k \nabla l )
   subroutine cint_nabla1l_2e(g, foff, goff, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      type(cint_env_vars), intent(in) :: envs
      integer  :: i, j, k, l, n, ptr, di, dk, dl, dj, nroots, gs
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: al2
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      nroots = envs%nrys_roots; gs = envs%g_size
      al2 = -2.0_dp * envs%al
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
         do k = 0, lk
            ptr = dj*j + dk*k
            do i = 0, li
               do n = ptr, ptr + nroots - 1
                  g(fx+n) = al2 * g(gx+dl+n)
                  g(fy+n) = al2 * g(gy+dl+n)
                  g(fz+n) = al2 * g(gz+dl+n)
               end do
               ptr = ptr + di
            end do
         end do
         do l = 1, ll
            do k = 0, lk
               ptr = dj*j + dl*l + dk*k
               do i = 0, li
                  do n = ptr, ptr + nroots - 1
                     g(fx+n) = l*g(gx-dl+n) + al2*g(gx+dl+n)
                     g(fy+n) = l*g(gy-dl+n) + al2*g(gy+dl+n)
                     g(fz+n) = l*g(gz-dl+n) + al2*g(gz+dl+n)
                  end do
                  ptr = ptr + di
               end do
            end do
         end do
      end do
   end subroutine cint_nabla1l_2e

   ! The four x-operations.  All the same body with a different shift stride;
   ! written out rather than parameterised because the stride is what makes
   ! them different, and a routine whose only content is a stride is harder to
   ! check against the C than four that read like it.
   !
   ! ( x^1 i j | kl ).  r is the shift from the origin to the centre of |i>.
   subroutine cint_x1i_2e(g, foff, goff, r, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_2e(g, foff, goff, r, li, lj, lk, ll, envs, envs%g_stride_i)
   end subroutine cint_x1i_2e

   ! ( i x^1 j | kl )
   subroutine cint_x1j_2e(g, foff, goff, r, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_2e(g, foff, goff, r, li, lj, lk, ll, envs, envs%g_stride_j)
   end subroutine cint_x1j_2e

   ! ( ij | x^1 k l )
   subroutine cint_x1k_2e(g, foff, goff, r, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_2e(g, foff, goff, r, li, lj, lk, ll, envs, envs%g_stride_k)
   end subroutine cint_x1k_2e

   ! ( ij | k x^1 l )
   subroutine cint_x1l_2e(g, foff, goff, r, li, lj, lk, ll, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_2e(g, foff, goff, r, li, lj, lk, ll, envs, envs%g_stride_l)
   end subroutine cint_x1l_2e

   ! f(...) = g(shifted by one along `stride`) + r * g(...).  The loop nest is
   ! identical in all four C routines; only `stride` and the vector differ.
   subroutine x1_2e(g, foff, goff, r, li, lj, lk, ll, envs, stride)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, lk, ll, stride
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      integer :: i, j, k, l, n, ptr, di, dk, dl, dj, nroots, gs
      integer :: fx, fy, fz, gx, gy, gz
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      nroots = envs%nrys_roots; gs = envs%g_size
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
      do l = 0, ll
      do k = 0, lk
         ptr = dj*j + dl*l + dk*k
         do i = 0, li
            do n = ptr, ptr + nroots - 1
               g(fx+n) = g(gx+stride+n) + r(0) * g(gx+n)
               g(fy+n) = g(gy+stride+n) + r(1) * g(gy+n)
               g(fz+n) = g(gz+stride+n) + r(2) * g(gz+n)
            end do
            ptr = ptr + di
         end do
      end do
      end do
      end do
   end subroutine x1_2e

end module cint_g2e_ops
