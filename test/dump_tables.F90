!
! Writes every generated table as an unformatted stream, in the order
! scripts/check_tables_bitwise.py expects, so the bytes can be compared
! against the ones gcc emitted into the C object files.
!
program dump_tables
   use cint_const, only: dp
#ifdef HAVE_SR_TABLES
   use cint_tab_polyfits
#endif
   use cint_tab_roots_x0
   use cint_tab_jacobi
   use cint_tab_cart2sph
#ifdef HAVE_SR_TABLES
   use cint_tab_sr_roots
#endif
   use cint_tab_jacobi_ext
   implicit none

   character(len=256) :: path
   integer :: u

   if (command_argument_count() < 1) then
      write(*,'(A)') "usage: dump_tables <outfile>"
      stop 2
   end if
   call get_command_argument(1, path)

   open(newunit=u, file=trim(path), access='stream', form='unformatted', &
        status='replace', action='write')

#ifdef HAVE_SR_TABLES
   write(u) COS_14_14
#endif
   write(u) POLY_SMALLX_R0, POLY_SMALLX_R1, POLY_SMALLX_W0, POLY_SMALLX_W1
   write(u) POLY_LARGEX_RT, POLY_LARGEX_WW
   write(u) JACOBI_ALPHA, JACOBI_BETA
   write(u) JACOBI_RN_PART2, JACOBI_SN, JACOBI_COEF, JACOBI_COEF_ORDER
   write(u) g_trans_cart2sph, g_trans_cart2jR, g_trans_cart2jI, len_cart
#ifdef HAVE_SR_TABLES
   ! Only present when WITH_POLYNOMIAL_FIT is on, exactly as in the C.
   write(u) SR_DATA0_X, SR_DATA0_W, SR_DATA1_X, SR_DATA1_W
   write(u) SR_DATA2_X, SR_DATA2_W, SR_DATAL_X, SR_DATAL_W
#endif
   write(u) lJACOBI_ALPHA, lJACOBI_BETA, lJACOBI_RN_PART2, lJACOBI_SN, lJACOBI_COEF
   write(u) qJACOBI_ALPHA, qJACOBI_BETA, qJACOBI_RN_PART2, qJACOBI_SN, qJACOBI_COEF

   close(u)
end program dump_tables
