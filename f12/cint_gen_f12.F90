! Guarded by COMPILE_F12.  Without it this file compiles to nothing at all,
! which is what lets every source in fortran/src be listed unconditionally --
! and the reason it has to be a guard rather than "just always build it" is
! that the chain below ends at cint_tab_stg_roots, a 100 MB generated table
! that is deliberately not in the repository.  See
! scripts/stg_roots_to_fortran.py.
!
! No #else stub: fpm scans for module names lexically, before the
! preprocessor runs, so a second `module cint_gen_f12` in this file reads to it
! as a duplicate definition.  An absent module is also the better error for
! anyone who USEs it without asking for F12.
#ifdef COMPILE_F12
!
! The F12 entry points: the Yukawa potential and the Slater-type geminal.
!
! Ported from src/cint2e_f12.c.  Only built with WITH_F12, as in the C.
!
! Ten integrals, each spherical only -- the C emits no Cartesian and no
! spinor form for this family -- over the ordinary two-electron driver.
! Every one of them reuses a gout the catalogue already has; the whole of
! F12 is a different envs setup and a different quadrature underneath the
! same machinery, which is why this file is as short as it is.
!
module cint_gen_f12
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws
   use cint_2e,        only: cint_2e_drv, cint_gout2e
   use cint_g2e_f12,   only: cint_init_int2e_yp_envvars, &
                             cint_init_int2e_stg_envvars, &
                             cint_all_2e_stg_optimizer
   use cint_gen_grad2, only: CINTgout2e_int2e_ip1
   use cint_gen_hess,  only: CINTgout2e_int2e_ipip1, CINTgout2e_int2e_ipvip1, &
                             CINTgout2e_int2e_ip1ip2
   implicit none
   private

   public :: int2e_yp_sph
   public :: int2e_yp_optimizer
   public :: int2e_yp_ip1_sph
   public :: int2e_yp_ip1_optimizer
   public :: int2e_yp_ipip1_sph
   public :: int2e_yp_ipip1_optimizer
   public :: int2e_yp_ipvip1_sph
   public :: int2e_yp_ipvip1_optimizer
   public :: int2e_yp_ip1ip2_sph
   public :: int2e_yp_ip1ip2_optimizer
   public :: int2e_stg_sph
   public :: int2e_stg_optimizer
   public :: int2e_stg_ip1_sph
   public :: int2e_stg_ip1_optimizer
   public :: int2e_stg_ipip1_sph
   public :: int2e_stg_ipip1_optimizer
   public :: int2e_stg_ipvip1_sph
   public :: int2e_stg_ipvip1_optimizer
   public :: int2e_stg_ip1ip2_sph
   public :: int2e_stg_ip1ip2_optimizer

contains

   function int2e_yp_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]

      call cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_yp_sph

   subroutine int2e_yp_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_yp_optimizer

   function int2e_yp_ip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]

      call cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_yp_ip1_sph

   subroutine int2e_yp_ip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_yp_ip1_optimizer

   function int2e_yp_ipip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 0, 0, 0, 2, 1, 1, 9]

      call cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ipip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_yp_ipip1_sph

   subroutine int2e_yp_ipip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [2, 0, 0, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_yp_ipip1_optimizer

   function int2e_yp_ipvip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 1, 0, 0, 2, 1, 1, 9]

      call cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ipvip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_yp_ipvip1_sph

   subroutine int2e_yp_ipvip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 1, 0, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_yp_ipvip1_optimizer

   function int2e_yp_ip1ip2_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 0, 1, 0, 2, 1, 1, 9]

      call cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ip1ip2
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_yp_ip1ip2_sph

   subroutine int2e_yp_ip1ip2_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 0, 1, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_yp_ip1ip2_optimizer

   function int2e_stg_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]

      call cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_stg_sph

   subroutine int2e_stg_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_stg_optimizer

   function int2e_stg_ip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]

      call cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_stg_ip1_sph

   subroutine int2e_stg_ip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 0, 0, 0, 1, 1, 1, 3]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_stg_ip1_optimizer

   function int2e_stg_ipip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [2, 0, 0, 0, 2, 1, 1, 9]

      call cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ipip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_stg_ipip1_sph

   subroutine int2e_stg_ipip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [2, 0, 0, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_stg_ipip1_optimizer

   function int2e_stg_ipvip1_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 1, 0, 0, 2, 1, 1, 9]

      call cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ipvip1
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_stg_ipvip1_sph

   subroutine int2e_stg_ipvip1_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 1, 0, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_stg_ipvip1_optimizer

   function int2e_stg_ip1ip2_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) &
         result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [1, 0, 1, 0, 2, 1, 1, 9]

      call cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => CINTgout2e_int2e_ip1ip2
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_stg_ip1ip2_sph

   subroutine int2e_stg_ip1ip2_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [1, 0, 1, 0, 2, 1, 1, 9]
      call cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_stg_ip1ip2_optimizer

end module cint_gen_f12
#endif
