!
! Public entry points for the scalar ECP integrals.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp.c.
!
! These are the routines an SCF calls, and they have the shape every other
! one-electron integral in this library has: a shell pair in, a block out,
! `dims` saying where in a larger matrix that block belongs.  What is
! different is where the potential comes from -- not an argument, but two
! slots of `env` pointing into `bas`:
!
!     ecpbas  = bas + int(env(AS_ECPBAS_OFFSET)) * BAS_SLOTS
!     necpbas = int(env(AS_NECPBAS))
!
! That convention is PySCF's rather than libcint's, and following it is what
! lets an `env` built by `gto.Mole` be handed to this library unchanged.  A
! caller building its own arrays has to set both slots; leaving them zero is
! not an error, it means "no ECP", and the integral comes back zero.
!
! There is no cache argument.  The C threads a `cache` pointer through every
! routine and sizes it with `ECPscalar_cache_size`; here the working arrays
! are allocatables owned by the routine that needs them, so the sizing
! routine has no counterpart and none is provided.
!
module cint_ecp_drv
   use cint_const, only: dp
   use cint_bas, only: ANG_OF, NCTR_OF, BAS_SLOTS, cint_len_cart
   use cint_cart2sph, only: cint_c2s_bra_sph, cint_c2s_ket_sph_copy
   use cint_ecp_num, only: AS_ECPBAS_OFFSET, AS_NECPBAS
   use cint_ecp, only: ecp_type1_cart, ecp_type2_cart, ecp_type_scalar_cart
   implicit none
   private

   public :: ECP_TERM_TYPE1, ECP_TERM_TYPE2, ECP_TERM_SCALAR
   public :: ecp_scalar_cart_block, ecp_scalar_sph_block
   public :: ecp_scalar_cart, ecp_scalar_sph

   !> Which ECP term a driver should evaluate
   integer, parameter :: ECP_TERM_TYPE1 = 1
   integer, parameter :: ECP_TERM_TYPE2 = 2
   integer, parameter :: ECP_TERM_SCALAR = 3

contains

   !> Dispatch to one of the three Cartesian kernels
   function cart_kernel(term, gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env) &
      result(has_value)
      integer, intent(in) :: term
      real(dp), intent(inout) :: gctr(0:)
      integer, intent(in) :: shls(0:1), ecpbas(0:), necpbas, atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      select case (term)
      case (ECP_TERM_TYPE1)
         has_value = ecp_type1_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
      case (ECP_TERM_TYPE2)
         has_value = ecp_type2_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
      case default
         has_value = ecp_type_scalar_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
      end select
   end function cart_kernel

   !> One Cartesian shell-pair block, zeroed and filled
   function ecp_scalar_cart_block(term, gctr, shls, ecpbas, necpbas, atm, natm, &
                                  bas, nbas, env) result(has_value)
      integer, intent(in) :: term
      real(dp), intent(out) :: gctr(0:)
      integer, intent(in) :: shls(0:1), ecpbas(0:), necpbas, atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      integer :: ish, jsh, nfi, nfj, nci, ncj

      ish = shls(0); jsh = shls(1)
      nfi = cint_len_cart(bas(ANG_OF + ish*BAS_SLOTS))
      nfj = cint_len_cart(bas(ANG_OF + jsh*BAS_SLOTS))
      nci = bas(NCTR_OF + ish*BAS_SLOTS)
      ncj = bas(NCTR_OF + jsh*BAS_SLOTS)

      gctr(0:nfi*nfj*nci*ncj - 1) = 0.0_dp
      has_value = cart_kernel(term, gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
   end function ecp_scalar_cart_block

   !> One spherical shell-pair block
   !>
   !> s and p shells are already spherical in libcint's ordering, so a pair of
   !> them is filled in place and no transform runs at all -- which is the
   !> common case and the reason the C special-cases it rather than
   !> transforming an identity.
   function ecp_scalar_sph_block(term, gctr, shls, ecpbas, necpbas, atm, natm, &
                                 bas, nbas, env) result(has_value)
      integer, intent(in) :: term
      real(dp), intent(out) :: gctr(0:)
      integer, intent(in) :: shls(0:1), ecpbas(0:), necpbas, atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      integer :: ish, jsh, li, lj, nfi, nfj, nci, ncj, ngcart, di, dji, nij, j, loc
      real(dp), allocatable :: gcart(:), gtmp(:)

      ish = shls(0); jsh = shls(1)
      li = bas(ANG_OF + ish*BAS_SLOTS)
      lj = bas(ANG_OF + jsh*BAS_SLOTS)
      nfi = cint_len_cart(li); nfj = cint_len_cart(lj)
      nci = bas(NCTR_OF + ish*BAS_SLOTS)
      ncj = bas(NCTR_OF + jsh*BAS_SLOTS)
      ngcart = nfi*nfj*nci*ncj

      if (li < 2 .and. lj < 2) then
         gctr(0:ngcart - 1) = 0.0_dp
         has_value = cart_kernel(term, gctr, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)
         return
      end if

      di = nfi*nci
      dji = di*(lj*2 + 1)
      nij = (li*2 + 1)*(lj*2 + 1)*nci*ncj

      allocate (gcart(0:ngcart - 1))
      gcart = 0.0_dp
      has_value = cart_kernel(term, gcart, shls, ecpbas, necpbas, atm, natm, bas, nbas, env)

      if (has_value /= 0) then
         if (li < 2) then
            do j = 0, ncj - 1
               loc = cint_c2s_ket_sph_copy(gctr(j*dji:), gcart(j*nfj*di:), di, di, lj)
            end do
         else if (lj < 2) then
            loc = cint_c2s_bra_sph(gctr, (lj*2 + 1)*nci*ncj, gcart, li)
         else
            allocate (gtmp(0:di*(lj*2 + 1)*ncj - 1))
            do j = 0, ncj - 1
               loc = cint_c2s_ket_sph_copy(gtmp(j*dji:), gcart(j*nfj*di:), di, di, lj)
            end do
            loc = cint_c2s_bra_sph(gctr, (lj*2 + 1)*nci*ncj, gtmp, li)
            deallocate (gtmp)
         end if
      else
         gctr(0:nij - 1) = 0.0_dp
      end if

      deallocate (gcart)
   end function ecp_scalar_sph_block

   !> Copy a block into a larger matrix, or straight out when `dims` is absent
   subroutine distribute(out, gctr, dims, di, dj, zero_only)
      real(dp), intent(inout) :: out(0:)
      real(dp), intent(in) :: gctr(0:)
      integer, intent(in), optional :: dims(0:1)
      integer, intent(in) :: di, dj
      logical, intent(in) :: zero_only

      integer :: i, j, ni

      if (.not. present(dims)) then
         if (zero_only) then
            out(0:di*dj - 1) = 0.0_dp
         else
            out(0:di*dj - 1) = gctr(0:di*dj - 1)
         end if
      else
         ni = dims(0)
         do j = 0, dj - 1
            do i = 0, di - 1
               if (zero_only) then
                  out(i + j*ni) = 0.0_dp
               else
                  out(i + j*ni) = gctr(i + j*di)
               end if
            end do
         end do
      end if
   end subroutine distribute

   pure subroutine ecpbas_from_env(env, bas, offset, necpbas)
      real(dp), intent(in) :: env(0:)
      integer, intent(in) :: bas(0:)
      integer, intent(out) :: offset, necpbas
      offset = int(env(AS_ECPBAS_OFFSET))*BAS_SLOTS
      necpbas = int(env(AS_NECPBAS))
   end subroutine ecpbas_from_env

   !> Scalar ECP over a shell pair, spherical basis
   function ecp_scalar_sph(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      integer :: li, lj, di, dj, offset, necpbas
      real(dp), allocatable :: buf(:)

      li = bas(ANG_OF + shls(0)*BAS_SLOTS)
      lj = bas(ANG_OF + shls(1)*BAS_SLOTS)
      di = (li*2 + 1)*bas(NCTR_OF + shls(0)*BAS_SLOTS)
      dj = (lj*2 + 1)*bas(NCTR_OF + shls(1)*BAS_SLOTS)

      call ecpbas_from_env(env, bas, offset, necpbas)
      allocate (buf(0:max(di*dj, cint_len_cart(li)*cint_len_cart(lj) &
                          *bas(NCTR_OF + shls(0)*BAS_SLOTS) &
                          *bas(NCTR_OF + shls(1)*BAS_SLOTS)) - 1))

      has_value = ecp_scalar_sph_block(ECP_TERM_SCALAR, buf, shls, bas(offset:), &
                                       necpbas, atm, natm, bas, nbas, env)
      call distribute(out, buf, dims, di, dj, has_value == 0)
      deallocate (buf)
   end function ecp_scalar_sph

   !> Scalar ECP over a shell pair, Cartesian basis
   function ecp_scalar_cart(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      integer :: di, dj, offset, necpbas
      real(dp), allocatable :: buf(:)

      di = cint_len_cart(bas(ANG_OF + shls(0)*BAS_SLOTS))*bas(NCTR_OF + shls(0)*BAS_SLOTS)
      dj = cint_len_cart(bas(ANG_OF + shls(1)*BAS_SLOTS))*bas(NCTR_OF + shls(1)*BAS_SLOTS)

      call ecpbas_from_env(env, bas, offset, necpbas)
      allocate (buf(0:di*dj - 1))

      has_value = ecp_scalar_cart_block(ECP_TERM_SCALAR, buf, shls, bas(offset:), &
                                        necpbas, atm, natm, bas, nbas, env)
      call distribute(out, buf, dims, di, dj, has_value == 0)
      deallocate (buf)
   end function ecp_scalar_cart

end module cint_ecp_drv
