!
! The integral environment, and the dispatch interfaces that take it.
!
! Both live here together, deliberately.  The three procedure pointers in
! CINTEnvVars have it as an argument, and the routines they point at live in
! modules that need the type -- so putting the type in one module and the
! abstract interfaces in another makes a cycle Fortran cannot resolve.  One
! base module breaks it, and this is the module every other one depends on.
! (PORT_TO_FORTRAN.md 3.3.)
!
! Differences from the C's struct, all forced and none behavioural:
!
!   * The unions are gone.  `union {FINT nfk; FINT grids_offset;}` and its two
!     siblings were aliases for readability -- one field used under two names
!     depending on the integral -- and Fortran has no unions.  Each slot keeps
!     one name and a comment saying what else it means.  (3.2.)
!
!   * `double ai[1]` and friends become scalars.  They are arrays of one in the
!     C only so that a pointer can be taken to them.  (3.7.)
!
!   * ri, rj, rk, rl are pointers into `env` in the C; here they are copies of
!     the three coordinates.  Three doubles is cheaper than carrying an offset
!     and indexing through it, and it removes the aliasing.
!
!   * The K&R-declared function pointers become typed procedure pointers.  The
!     C writes `FINT (*f_g0_2e)()` with no prototype, so nothing checks the
!     calls; Fortran demands one exact interface per slot, which is a gain.
!
module cint_envs
   use cint_const, only: dp
   implicit none
   public

   integer, parameter :: MXRYSROOTS = 32
   integer, parameter :: ANG_MAX    = 15
   integer, parameter :: CART_MAX   = 136
   integer, parameter :: SHLS_MAX   = 1048576
   integer, parameter :: LMAX1 = 16                  ! > ANG_MAX, the C's
   integer, parameter :: MAX_PGTO_FOR_PAIRDATA = 2048

   ! ng[] slot names, from src/cint_config.h.in.  Note SLOT_RYS_ROOTS and
   ! POS_E2 are the same slot -- that is the C's doing, not a transcription
   ! error: a 1e integral has no second electron to position.
   integer, parameter :: IINC = 0, JINC = 1, KINC = 2, LINC = 3
   integer, parameter :: GSHIFT = 4, POS_E1 = 5, POS_E2 = 6
   integer, parameter :: SLOT_RYS_ROOTS = 6, TENSOR = 7

   real(dp), parameter :: EXPCUTOFF     = 60.0_dp
   real(dp), parameter :: MIN_EXPCUTOFF = 40.0_dp

   ! env[] slots, 0-based as the C has them.
   integer, parameter :: PTR_EXPCUTOFF    = 0
   integer, parameter :: PTR_COMMON_ORIG  = 1
   integer, parameter :: PTR_RINV_ORIG    = 4
   integer, parameter :: PTR_RINV_ZETA    = 7
   integer, parameter :: PTR_RANGE_OMEGA  = 8
   integer, parameter :: PTR_F12_ZETA     = 9
   integer, parameter :: PTR_GRIDS         = 12
   integer, parameter :: PTR_ENV_START    = 20

   ! Which cart-to-spherical routine a driver should apply.  The C selects
   ! this by passing a function pointer; an enum is clearer and, unlike a
   ! procedure pointer, can be compared for equality the way CINT1e_drv does
   ! when it works out the output dimensions.
   integer, parameter :: C2S_CART_1E = 1
   integer, parameter :: C2S_SPH_1E  = 2
   integer, parameter :: C2S_SF_1E   = 3   ! spinor, spin-free   (D9)
   integer, parameter :: C2S_SI_1E   = 4   ! spinor, spin-included (D9)

   ! Which 2D-to-4D transfer a 2e integral uses.  The C stores a function
   ! pointer (f_g0_2d4d); an enum does the same job here without a second
   ! abstract interface, because unlike f_gout this slot is chosen once in
   ! CINTinit_int2e_EnvVars from the angular momenta and never by a caller.
   ! Which 2D-recursion kernel a two-electron integral uses.
   integer, parameter :: G0_COULOMB = 0, G0_F12 = 1

   integer, parameter :: G2D4D_UNROLLED    = 1
   integer, parameter :: G2D4D_SR_UNROLLED = 2
   integer, parameter :: G2D4D_IK          = 3
   integer, parameter :: G2D4D_KJ          = 4
   integer, parameter :: G2D4D_IL          = 5
   integer, parameter :: G2D4D_LJ          = 6
   ! int2c2e has no j or l index, so the 2D array *is* the answer and there is
   ! nothing to transfer -- the C spells that by pointing f_g0_2d4d straight
   ! at CINTg0_2e_2d.
   integer, parameter :: G2D4D_NONE        = 7

   integer, parameter :: C2S_CART_2E1 = 11
   integer, parameter :: C2S_SPH_2E1  = 12
   integer, parameter :: C2S_CART_3C2E1 = 13
   integer, parameter :: C2S_SPH_3C2E1  = 14
   integer, parameter :: C2S_CART_2C2E1 = 15
   integer, parameter :: C2S_SPH_2C2E1  = 16

   ! int1e_type, as make_g1e_gout switches on it.
   integer, parameter :: INT1E_OVLP = 0
   integer, parameter :: INT1E_RINV = 1
   integer, parameter :: INT1E_NUC  = 2

   type :: cint_env_vars
      ! The caller's tables, held by reference the way the C does.  These
      ! were copies until they were measured: env alone is the whole
      ! environment array, and copying it per shell quartet cost more than
      ! the integral.  shls stays a copy -- it is four integers, and keeping
      ! it one means the entry points need no `target` on their smallest
      ! argument.
      integer,  pointer     :: atm(:), bas(:)
      real(dp), pointer     :: env(:)
      ! four at most, so a fixed array rather than an allocatable: this is
      ! per shell quartet, and the allocation showed up in the profile
      integer :: shls(0:3)
      integer :: natm, nbas

      integer :: i_l, j_l, k_l, l_l
      integer :: nfi, nfj
      ! nfk doubles as grids_offset in int1e_grids; nfl as ngrids
      integer :: nfk, nfl
      integer :: nf
      integer :: rys_order
      integer :: x_ctr(0:3)

      integer :: gbits
      integer :: ncomp_e1, ncomp_e2, ncomp_tensor

      integer :: li_ceil, lj_ceil, lk_ceil, ll_ceil
      integer :: g_stride_i, g_stride_k, g_stride_l, g_stride_j
      integer :: nrys_roots
      integer :: g_size

      integer :: g2d_ijmax, g2d_klmax
      real(dp) :: common_factor
      real(dp) :: expcutoff
      real(dp) :: rirj(0:2)
      real(dp) :: rkrl(0:2)

      real(dp) :: ri(0:2)
      real(dp) :: rj(0:2)
      real(dp) :: rk(0:2)
      ! rl doubles as the grids pointer in int1e_grids
      real(dp) :: rl(0:2)

      ! assigned during the primitive loop; scalars, not length-1 arrays
      real(dp) :: ai, aj, ak, al
      real(dp) :: fac
      real(dp) :: rij(0:2)
      real(dp) :: rkl(0:2)

      ! Where this integral's grid coordinates start in env.  The C keeps a
      ! `double *grids` pointer in the union with rl; an offset says the same
      ! thing and keeps rl a plain vector for every other integral.
      integer :: grids_off

      ! rx_in_rijrx / rx_in_rklrx are pointers in the C, aliasing ri or rj
      ! (resp. rk or rl) depending on which centre the 2D recursion is built
      ! along.  Copies here, like ri..rl.
      real(dp) :: rx_in_rijrx(0:2)
      real(dp) :: rx_in_rklrx(0:2)

      integer :: f_g0_2d4d
      ! The ng[] the integral was built with, kept verbatim.  The optimizer
      ! caches an index map per (arity, ng), and the only honest way to ask
      ! whether a given optimizer belongs to a given integral is to compare
      ! the ng they were each built from.  Reconstructing it from the
      ! ceilings almost works and then does not: int2c2e overwrites j_l with
      ! k_l after init, so lj_ceil - j_l stops meaning ng(JINC).
      integer :: ng(0:7)
      ! env(PTR_RANGE_OMEGA), copied in at setup.  Reaching through the
      ! env pointer component for it cost 16.4 million instructions against
      ! the C's 4.7 million, and it is read once per primitive quartet.
      real(dp) :: omega
      ! Which 2D-recursion kernel this integral wants.  Zero is the
      ! ordinary Coulomb one and the driver calls it directly; the F12
      ! kernels set it and the driver goes through the hook in cint_2e.
      !
      ! An integer and not a procedure pointer, which is what this was
      ! first.  A pointer component costs 1.9% on the two-electron path --
      ! measured with the dispatch removed, so it is the component and not
      ! the branch -- and every one of the 620 entry points declares one of
      ! these on the stack.
      ! THE ONE COMPONENT THAT KEEPS A DEFAULT, and the reason the rest do
      ! not.  Every other component mirrors a CINTEnvVars field that the C
      ! leaves uninitialised and the init routine sets before anything reads
      ! it; carrying a default here only paid gfortran to zero 52 components
      ! at every declaration -- 170 instructions per shell pair, 3.8% of a
      ! cheap integral like int2c2e, before any work happened.  g0_kind is
      ! different: it has no C counterpart, it selects the F12 kernel, and
      ! only the F12 setup assigns it, so every other path reads the default.
      integer :: g0_kind = G0_COULOMB
      ! WHICH OF THE FOUR SHELLS IS AN L SHELL -- s and p on shared
      ! exponents, see KAPPA_SP_SHELL in cint_bas -- as a bitmask, bit x for
      ! shls(x).
      !
      ! THE SECOND COMPONENT THAT KEEPS A DEFAULT, and for g0_kind's reason:
      ! only cint_init_int2e_envvars ever assigns it, so every other init
      ! routine -- 1e, 3c1e, 3c2e, 2c2e, the spinor ones, the 620 generated
      ! entry points that go through them -- reads the default and behaves
      ! exactly as it did before L shells existed.
      !
      ! One integer rather than four logicals because the question the hot
      ! paths ask is "is any of them an L shell?", and that is one compare
      ! against zero instead of four loads and three ors.
      integer :: sp_mask = 0

      procedure(gout_iface), nopass, pointer :: f_gout => null()
   end type cint_env_vars

   ! The Rys intermediate coefficients, one set per root.  The C calls this
   ! Rys2eT and passes a pointer to it through f_g0_2d4d.
   ! NO DEFAULT INITIALISERS HERE, deliberately.  This type is a local
   ! variable of cint_g0_2e, which runs once per primitive quartet, and
   ! initialisers would make every one of those calls zero 2.3 kB before
   ! writing over it.  When they were there the function was 55% of the
   ! profile; the C's Rys2eT is an uninitialised stack struct and so is this.
   ! Only the first nrys_roots entries of each array are ever read, and they
   ! are all written before that happens.
   type :: rys_2e_t
      real(dp) :: c00x(0:MXRYSROOTS-1)
      real(dp) :: c00y(0:MXRYSROOTS-1)
      real(dp) :: c00z(0:MXRYSROOTS-1)
      real(dp) :: c0px(0:MXRYSROOTS-1)
      real(dp) :: c0py(0:MXRYSROOTS-1)
      real(dp) :: c0pz(0:MXRYSROOTS-1)
      real(dp) :: b01(0:MXRYSROOTS-1)
      real(dp) :: b00(0:MXRYSROOTS-1)
      real(dp) :: b10(0:MXRYSROOTS-1)
   end type rys_2e_t

   abstract interface
      ! The gout kernels: contract the g intermediates into the output for one
      ! primitive pair.  Uniform across 1e and 2e, despite the C declaring the
      ! slot without a prototype.
      ! The 2D-recursion entry.  The C stores this as a function pointer too
      ! (cint.h.in:210) and calls through it in every 2e loop.  Here it is
      ! null for the ordinary Coulomb kernel and the driver calls cint_g0_2e
      ! directly, so the common path keeps a call LTO can see through; only
      ! the F12 kernels, which live behind WITH_F12 and cannot be named from
      ! cint_2e without a build-time dependency on an optional module, set it.
      function g0_2e_iface(g, rij, rkl, cutoff, envs) result(ok)
         import :: dp, cint_env_vars
         real(dp), intent(inout) :: g(0:)
         real(dp), intent(in)    :: rij(0:2), rkl(0:2)
         real(dp), intent(in)    :: cutoff
         type(cint_env_vars), intent(in) :: envs
         integer :: ok
      end function g0_2e_iface

      subroutine gout_iface(gout, g, idx, envs, empty)
         import :: dp, cint_env_vars
         ! ASSUMED-SIZE, not assumed-shape, for gout and idx.  The C passes
         ! three bare pointers here; an assumed-shape dummy makes the caller
         ! build an array descriptor and the callee unpack one, on a call
         ! that happens once per primitive quartet.  The extents are not
         ! lost, only unenforced: gout is nf*ncomp long, idx is nf*3, and
         ! both are derivable from the envs that travels with them.
         real(dp), intent(inout) :: gout(0:*)
         ! g is INOUT, not IN: a gout kernel for anything with a derivative in
         ! it builds its intermediates inside the same array, through the
         ! G1E_* operations.  The C passes a bare `double *` and the mutation
         ! is invisible; declaring it here is the honest version.
         !
         ! TARGET so a kernel can bind one pointer per g block and read
         ! g1p(ix) rather than g(g1+ix), which is what the C's own generated
         ! kernels do with `double *g1 = g0 + envs->g_size*3`.
         real(dp), intent(inout), target :: g(0:)
         integer,  intent(in)    :: idx(0:*)
         type(cint_env_vars), intent(in) :: envs
         integer,  intent(in)    :: empty
      end subroutine gout_iface
   end interface

   ! The C's PairData, one per primitive pair.
   type :: pair_data
      real(dp) :: rij(0:2) = 0.0_dp
      real(dp) :: eij = 0.0_dp
      real(dp) :: cceij = 0.0_dp
   end type pair_data


   ! The C spells "this shell pair is entirely negligible" as a pointer to
   ! 0xffff...; an offset table says it with a value no offset can take.
   integer, parameter :: OPT_NOVALUE = -1

   ! The C's CINTOpt (include/cint.h.in:145), which is five ragged pointer
   ! arrays.  Each becomes one flat array plus a table of offsets into it --
   ! the idiom the rest of the port already uses for the workspace, and the
   ! reason none of this needs a Fortran pointer.
   !
   ! It caches exactly what cint_2e_loop_nopt recomputes for every shell
   ! quartet: the logged coefficient bounds, the non-zero contraction map,
   ! the primitive-pair data, and the Cartesian index map.  Nothing here
   ! changes a result; the C's own comment is that it is not worth building
   ! unless the driver is called more than a few hundred times.
   !
   ! `ng` and `order` record what it was built for.  The C trusts the caller
   ! to pass the optimizer that belongs to the integral; a mismatch there
   ! reads the wrong index map and returns wrong numbers quietly, so the
   ! drivers here check before using one.
   type :: cint_opt_t
      logical :: active = .false.
      integer :: nbas = 0
      integer :: order = 0
      integer :: ng(0:7)
      real(dp), allocatable :: log_maxc(:)
      integer,  allocatable :: log_maxc_off(:)
      integer,  allocatable :: non0ctr(:), non0ctr_off(:)
      integer,  allocatable :: sortedidx(:), sortedidx_off(:)
      type(pair_data), allocatable :: pd(:)
      integer,  allocatable :: pd_off(:)
      integer,  allocatable :: idx(:), idx_off(:)
   end type cint_opt_t


contains

   ! Is this optimizer the one for this integral?  The C does not ask -- it
   ! trusts the caller to pass the CINTOpt built for the entry point being
   ! called, and passing the wrong one reads the wrong index map and returns
   ! wrong numbers with no complaint.  The check is a handful of integer
   ! comparisons against what the envs already knows, so it is made.
   !
   ! ng is not stored on the envs, but every slot of it is recoverable: the
   ! increments are the difference between each ceiling and its l, and the
   ! rest are copied across verbatim by the init routines.
   pure logical function cint_opt_usable(opt, order, envs) result(ok)
      type(cint_opt_t),    intent(in) :: opt
      integer,             intent(in) :: order
      type(cint_env_vars), intent(in) :: envs
      ok = opt%active .and. opt%order == order .and. opt%nbas == envs%nbas &
           .and. all(opt%ng == envs%ng)
   end function cint_opt_usable

   ! Where this quartet's Cartesian index map sits in the cache.  The key is
   ! the C's, base LMAX1 in the angular momenta, truncated to the arity.
   pure integer function cint_opt_idx_key(envs, order) result(key)
      type(cint_env_vars), intent(in) :: envs
      integer,             intent(in) :: order
      select case (order)
      case (2);      key = envs%i_l*LMAX1 + envs%j_l
      case (3);      key = (envs%i_l*LMAX1 + envs%j_l)*LMAX1 + envs%k_l
      case default;  key = ((envs%i_l*LMAX1 + envs%j_l)*LMAX1 + envs%k_l)*LMAX1 &
                         + envs%l_l
      end select
   end function cint_opt_idx_key

   ! Is shls(x) an L shell?  x is 0..3 in the i, j, k, l order of shls.
   pure logical function envs_is_sp(envs, x) result(yes)
      type(cint_env_vars), intent(in) :: envs
      integer,             intent(in) :: x
      yes = btest(envs%sp_mask, x)
   end function envs_is_sp

   ! bas(SLOT, I) and atm(SLOT, I), in the C's spelling.
   pure function bas_of(envs, slot, i) result(v)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(in) :: slot, i
      integer :: v
      v = envs%bas(8 * i + slot)
   end function bas_of

   pure function atm_of(envs, slot, i) result(v)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(in) :: slot, i
      integer :: v
      v = envs%atm(6 * i + slot)
   end function atm_of

end module cint_envs
