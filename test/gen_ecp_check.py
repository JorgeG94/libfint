#!/usr/bin/env python3
"""Generate test/ecp_check.F90 from PySCF.

    python3 test/gen_ecp_check.py

Every other check in this suite compares against libcint at run time, by
linking the fetched C.  This one cannot: libcint has no ECP code, so there is
nothing to link.  The reference is PySCF's `mol.intor('ECPscalar_sph')`, and it
has to be baked in.

That makes the choice of system load-bearing, because a committed reference is
only ever as good as what it covers -- and this check has already been caught
claiming coverage it did not have.  It began as Cu/H alone, with a docstring
asserting three projected channels and an overlap screen that "has to reject as
well as accept".  Cu/LANL2DZ carries two projected channels, and every one of
its 135 screens accepts, the largest exponent argument reaching 17.7 against a
cutoff of 39.  So `lc >= 2`, the `has_value == 0` path and the zero-fill were
all dead while the file said otherwise.

Hence two systems:

  Cu/H at 2.0 bohr in LANL2DZ
      local (l = -1) and two projected channels; r exponents 0, 1 and 2, so
      three of the four RADI_POWER branches; s, p and d orbital shells, so all
      four branches of the spherical transform; contracted shells; a second
      atom with no ECP.  Every screen accepts.

  I/H at 25 bohr in def2-SVP
      lc = 2, which is the only way to reach the `case default` in
      type2_facs_ang that transforms the projector.  With H that far out, 144
      of 364 screens reject, which reaches `has_value == 0`, the zero-fill and
      `distribute`'s zeroing branch.

Still unreached, and unreachable rather than untested: `RADI_POWER` case 3 and
the general branch in cint_ecp_rad.  No ECP set in PySCF's library -- crenbl,
lanl2dz, def2, sbkjc, stuttgart -- carries an r exponent above 2.

The tolerance is not bit-identity.  cint_ecp_rad explains the main reason: the
C reaches its modified spherical Bessel function through a 400-entry
interpolation table and this port evaluates the series.

There is a second reason, smaller and irreducible.  PySCF integrates on a
*static table*, `rs_gauss_chebyshev2047`, whose literals carry 16 significant
digits -- one short of a double round trip.  893 of its 2047 abscissae and 1274
of its weights differ by one ULP from what its own generator produces.  This
port calls the generator, so it matches PySCF's `ECPgauss_chebyshev` exactly and
PySCF's production grid to about a ULP per point.

The tolerance is stated against the largest element of the matrix rather than
element by element, because that is how the error behaves.  Across seven ECP
systems -- Cu, Au and W in LANL2DZ, I and Ba in def2-SVP, Pt and Hg in
def2-TZVP -- it runs between 7.8e-17 and 1.8e-14 of the largest element, which
is itself between 3 and 308.  TOL is 1e-11, about three orders of slack,
because the radial refinement stops on a convergence test and a pair sitting
near that threshold can take a different number of levels here than in the C.

A purely element-wise relative bound would be the wrong instrument: on Au2 an
element of magnitude 5e-6 carries 5e-14 of absolute error, which reads as
1.1e-08 and means nothing.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# tag, geometry, basis (and ECP), label, whether the screen should reject
SYSTEMS = [
    ("cu", "Cu 0 0 0; H 0 0 2.0", "lanl2dz", "Cu/H LANL2DZ", 0),
    ("i", "I 0 0 0; H 0 0 25.0", "def2-svp", "I/H(far) def2-SVP", 1),
]

# The C ABI check gets its own system, and the choice is not arbitrary.
# W in SBKJC is the one arrangement to hand where a fault in `ecp_env_len`
# is *detectable*: its highest-ptr_coeff ECP shell is an r^0 term, so
# mistaking slot 3 (RADI_POWER) for a contraction count makes the declared
# env extent one element too SHORT, which -fcheck=bounds traps. On Cu the same
# fault makes it two too long, which reads nothing and traps nothing -- the
# original bug was found by review precisely because no test could see it.
# W also carries lc up to 3, higher than either system above.
C_SYSTEM = ("W 0 0 0; H 0 0 3.2", "sbkjc", 1)


def wrap(values, indent="      ", per_line=3):
    out, line = [], []
    for i, v in enumerate(values):
        line.append(v)
        if len(line) == per_line or i == len(values) - 1:
            out.append(indent + ", ".join(line))
            line = []
    return ", &\n".join(out)


def cwrap(values, per_line=6):
    """Comma-separated, wrapped, for a C initialiser."""
    out = []
    for i in range(0, len(values), per_line):
        out.append("    " + ", ".join(values[i:i + per_line]))
    return ",\n".join(out)


def emit_system(tag, geom, basis, label, expect_reject, spin=0):
    import numpy as np
    from pyscf import gto

    mol = gto.M(atom=geom, basis=basis, ecp=basis, spin=spin, verbose=0)
    ref = mol.intor("ECPscalar_sph")          # also fills env[18], env[19]
    bas = np.vstack((mol._bas, mol._ecpbas))
    env = mol._env
    atm = np.asarray(mol._atm)
    ao_loc = mol.ao_loc_nr()

    d = {
        "u": tag.upper(), "label": label, "geom": geom, "basis": basis,
        "expect_reject": expect_reject,
        "natm": mol.natm, "nbas_all": bas.shape[0], "nbas": mol.nbas,
        "nenv": len(env), "nao": mol.nao_nr(),
        "atm": wrap([str(int(x)) for x in atm.ravel()], per_line=12),
        "bas": wrap([str(int(x)) for x in bas.ravel()], per_line=12),
        "env": wrap(["%.17e_dp" % v for v in env]),
        "ao_loc": wrap([str(int(x)) for x in ao_loc], per_line=12),
        "ref": wrap(["%.17e_dp" % v for v in ref.ravel(order="F")]),
        "lmax": int(max(r[1] for r in np.asarray(mol._ecpbas))),
        "_atm": [int(x) for x in atm.ravel()],
        "_bas": [int(x) for x in bas.ravel()],
        "_ao_loc": [int(x) for x in ao_loc],
        "_env": [float(x) for x in env],
        "_ref": [float(x) for x in ref.ravel(order="F")],
    }
    return DATA_BLOCK % d, d


def main():
    blocks, metas = [], []
    for tag, geom, basis, label, rej in SYSTEMS:
        blk, meta = emit_system(tag, geom, basis, label, rej)
        blocks.append(blk)
        metas.append(meta)

    calls = "\n".join(
        '   call check_system("%(label)s", %(u)s_NATM, %(u)s_NBAS, %(u)s_NAO, &\n'
        '                     %(u)s_ATM, %(u)s_BAS, %(u)s_AO_LOC, %(u)s_ENV, %(u)s_REF, &\n'
        '                     %(expect_reject)d, bad)' % m for m in metas)

    src = PROGRAM % {
        "pyscf_version": __import__("pyscf").__version__,
        "data": "\n".join(blocks),
        "calls": calls,
        "systems": "\n".join(
            "!   %-22s lc up to %d, %d AOs%s"
            % (m["label"], m["lmax"], m["nao"],
               ", screen rejects" if m["expect_reject"] else "")
            for m in metas),
    }
    with open(os.path.join(HERE, "ecp_check.F90"), "w") as fh:
        fh.write(src)

    # The C ABI gets its own check, from the same data, because the Fortran
    # one cannot reach it: ecp_check calls the module directly, and the
    # bind(C) wrappers are a separate layer. Both bugs fixed in 3042cff lived
    # in that layer and neither was visible to any test.
    geom, basis, spin = C_SYSTEM
    _, m = emit_system("c", geom, basis, "C ABI", 0, spin=spin)
    cm = dict(m)
    cm["pyscf_version"] = __import__("pyscf").__version__
    cm["c_atm"] = cwrap([str(x) for x in m["_atm"]], per_line=12)
    cm["c_bas"] = cwrap([str(x) for x in m["_bas"]], per_line=12)
    cm["c_ao_loc"] = cwrap([str(x) for x in m["_ao_loc"]], per_line=12)
    cm["c_env"] = cwrap(["%.17e" % v for v in m["_env"]], per_line=4)
    cm["c_ref"] = cwrap(["%.17e" % v for v in m["_ref"]], per_line=4)
    with open(os.path.join(HERE, "ecp_abi_check.c"), "w") as fh:
        fh.write(C_PROGRAM % cm)
    print("%d systems -> ecp_check.F90, ecp_abi_check.c" % len(metas))


C_PROGRAM = r"""/*
 * The ECP C ABI, driven the way a C consumer drives it.
 *
 * GENERATED by test/gen_ecp_check.py -- do not edit.
 * Reference: PySCF %(pyscf_version)s, mol.intor('ECPscalar_sph')
 * System:    %(geom)s in %(basis)s with the matching ECP
 *
 * ecp_check exercises the Fortran module directly.  That leaves
 * ECPscalar_sph and ECPscalar_cart in src/cint_c_abi.f90 -- the layer a C
 * program actually links -- called by nothing, and it is the layer most able
 * to fail quietly.  Both faults fixed in "Fix the radial grid size, and two
 * latent faults in the C ABI" lived there: an env extent computed by reading
 * an ecpbas row's RADI_POWER slot as though it were a contraction count, and
 * a bas extent assuming the ECP rows begin exactly at nbas.  Neither was
 * visible to any test, because no test called the layer.
 *
 * What is asserted here:
 *   - the symbols exist under libcint's names, so a rename fails at link;
 *   - values against PySCF, through the dims-present path;
 *   - the dims-absent path agrees with it bitwise -- one computation, two
 *     ways of placing the answer;
 *   - ECPscalar_cart agrees about has_value.
 *
 * bas below is the orbital shells followed by the ECP shells, and nbas is the
 * orbital count alone; env[18] and env[19] say where the ECP rows start and
 * how many there are.  That is PySCF's convention and the reason the two
 * extents above are subtle.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ECPscalar_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, void *opt, double *cache);
int ECPscalar_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, void *opt, double *cache);

#define BAS_SLOTS 8
#define ANG_OF    1
#define NCTR_OF   3

#define NATM %(natm)d
#define NBAS_ALL %(nbas_all)d
#define NBAS %(nbas)d
#define NENV %(nenv)d
#define NAO  %(nao)d

/* relative to the largest element; see gen_ecp_check.py */
#define TOL 1e-11

static int atm[NATM*6] = {
%(c_atm)s};

static int bas[NBAS_ALL*BAS_SLOTS] = {
%(c_bas)s};

static int ao_loc[NBAS+1] = {
%(c_ao_loc)s};

static double env[NENV] = {
%(c_env)s};

static double ref[NAO*NAO] = {
%(c_ref)s};

int main(void)
{
    static double strided[NAO*NAO];
    static double contig[NAO*NAO];
    double blk[4096], cblk[4096];
    int dims[2] = {NAO, NAO};
    int ish, jsh, i, j, bad = 0;
    double scale = 0.0, wabs = 0.0, d;

    for (i = 0; i < NAO*NAO; i++) { strided[i] = 0.0; contig[i] = 0.0; }

    for (ish = 0; ish < NBAS; ish++) {
        for (jsh = 0; jsh < NBAS; jsh++) {
            int shls[2]; int hv, hvc, di, dj, dic, djc;
            int li = bas[ish*BAS_SLOTS+ANG_OF], lj = bas[jsh*BAS_SLOTS+ANG_OF];
            shls[0] = ish; shls[1] = jsh;
            di = (2*li+1)*bas[ish*BAS_SLOTS+NCTR_OF];
            dj = (2*lj+1)*bas[jsh*BAS_SLOTS+NCTR_OF];

            /* dims present: the driver writes into the caller's matrix */
            hv = ECPscalar_sph(strided + ao_loc[jsh]*NAO + ao_loc[ish],
                               dims, shls, atm, NATM, bas, NBAS, env, NULL, NULL);

            /* dims absent: a natural block the caller places itself */
            for (i = 0; i < di*dj; i++) blk[i] = 1.0;   /* poisoned */
            hvc = ECPscalar_sph(blk, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL);
            if (hvc != hv) {
                printf("FAIL has_value differs between the dims paths, shells %%d %%d\n", ish, jsh);
                bad++;
            }
            for (j = 0; j < dj; j++)
                for (i = 0; i < di; i++)
                    contig[(ao_loc[jsh]+j)*NAO + ao_loc[ish]+i] = blk[j*di+i];

            /* the Cartesian entry point: same contract, different shape */
            dic = ((li+1)*(li+2)/2)*bas[ish*BAS_SLOTS+NCTR_OF];
            djc = ((lj+1)*(lj+2)/2)*bas[jsh*BAS_SLOTS+NCTR_OF];
            for (i = 0; i < dic*djc; i++) cblk[i] = 1.0;
            hvc = ECPscalar_cart(cblk, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL);
            if (hvc != hv) {
                printf("FAIL ECPscalar_cart disagrees about has_value, shells %%d %%d\n", ish, jsh);
                bad++;
            }
        }
    }

    for (i = 0; i < NAO*NAO; i++) if (fabs(ref[i]) > scale) scale = fabs(ref[i]);
    if (scale <= 0.0) { printf("FAIL the reference matrix is identically zero\n"); return 1; }

    for (i = 0; i < NAO*NAO; i++) {
        d = fabs(strided[i] - ref[i]);
        if (d > wabs) wabs = d;
        /* the two placements are one computation: they must agree exactly */
        if (strided[i] != contig[i]) {
            printf("FAIL dims-present and dims-absent differ at %%d: %%.17e vs %%.17e\n",
                   i, strided[i], contig[i]);
            bad++;
            break;
        }
    }

    printf("ecp_abi_check: %%dx%%d  largest |ref| %%.3e  worst/largest %%.3e\n",
           NAO, NAO, scale, wabs/scale);
    if (wabs/scale > TOL) { printf("FAIL outside tolerance\n"); bad++; }
    if (bad) { printf("ecp_abi_check: FAILED\n"); return 1; }
    printf("ecp_abi_check: OK\n");
    return 0;
}
"""


DATA_BLOCK = '''
   ! ---- %(label)s : %(geom)s in %(basis)s ----
   integer, parameter :: %(u)s_NATM = %(natm)d
   integer, parameter :: %(u)s_NBAS_ALL = %(nbas_all)d
   integer, parameter :: %(u)s_NBAS = %(nbas)d
   integer, parameter :: %(u)s_NENV = %(nenv)d
   integer, parameter :: %(u)s_NAO = %(nao)d

   integer, parameter :: %(u)s_ATM(0:%(u)s_NATM*6 - 1) = [ &
%(atm)s]

   integer, parameter :: %(u)s_BAS(0:%(u)s_NBAS_ALL*BAS_SLOTS - 1) = [ &
%(bas)s]

   integer, parameter :: %(u)s_AO_LOC(0:%(u)s_NBAS) = [ &
%(ao_loc)s]

   real(dp), parameter :: %(u)s_ENV(0:%(u)s_NENV - 1) = [ &
%(env)s]

   real(dp), parameter :: %(u)s_REF(%(u)s_NAO*%(u)s_NAO) = [ &
%(ref)s]
'''


PROGRAM = '''!
! The scalar ECP integrals against PySCF.
!
! GENERATED by test/gen_ecp_check.py -- do not edit.
! Reference: PySCF %(pyscf_version)s, mol.intor('ECPscalar_sph')
!
! Systems:
%(systems)s
!
! The reference is PySCF and not libcint because libcint has no ECP code --
! cint_ecp and its supporting modules are translated from
! pyscf/lib/gto/nr_ecp.c.  There is nothing to link at run time, so unlike the
! rest of this suite the reference is committed.
!
! Two systems, and the reason is a lesson rather than thoroughness.  This began
! as Cu/H alone, claiming in its own header to cover three projected channels
! and to force the overlap screen to reject.  It covered two channels and
! rejected nothing.  I/H at 25 bohr is here to reach lc = 2 and the rejection
! path; see gen_ecp_check.py for the numbers.
!
! The tolerance is stated against the largest element rather than element by
! element, and is not bit-identity.  gen_ecp_check.py explains both.
!
program ecp_check
   use cint_const, only: dp
   use cint_bas, only: ANG_OF, NCTR_OF, BAS_SLOTS
   use cint_ecp_num, only: AS_ECPBAS_OFFSET, AS_NECPBAS
   use cint_ecp_drv, only: ecp_scalar_sph, ecp_scalar_cart
   implicit none

   !> Relative to the largest element of the matrix, not to each element
   real(dp), parameter :: TOL = 1.0e-11_dp
%(data)s
   integer :: bad

   bad = 0
%(calls)s
   call check_no_ecp(bad)

   if (bad > 0) then
      write (*, '(a,i0,a)') "ecp_check: FAILED (", bad, " problem(s))"
      stop 1
   end if
   write (*, '(a)') "ecp_check: OK"

contains

   subroutine check_system(label, natm, nbas, nao, atm, bas, ao_loc, env, ref, &
                           expect_reject, bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: natm, nbas, nao
      integer, intent(in) :: atm(0:), bas(0:), ao_loc(0:)
      real(dp), intent(in) :: env(0:), ref(:)
      integer, intent(in) :: expect_reject
      integer, intent(inout) :: bad

      real(dp), allocatable :: mat(:, :), blk(:), cblk(:)
      real(dp) :: worst, wabs, d, denom, scale
      integer :: ish, jsh, di, dj, dic, djc, i0, j0, i, li, lj
      integer :: r, cc, hv, hvc, nzero, nscreened, nstale

      allocate (mat(nao, nao))
      mat = 0.0_dp
      nscreened = 0
      nstale = 0
      worst = 0.0_dp

      do ish = 0, nbas - 1
         do jsh = 0, nbas - 1
            li = bas(ANG_OF + ish*BAS_SLOTS)
            lj = bas(ANG_OF + jsh*BAS_SLOTS)
            di = (2*li + 1)*bas(NCTR_OF + ish*BAS_SLOTS)
            dj = (2*lj + 1)*bas(NCTR_OF + jsh*BAS_SLOTS)

            allocate (blk(0:di*dj - 1))
            ! Poisoned rather than zeroed. A driver that returns without
            ! writing would otherwise read as a correct zero block, which is
            ! exactly what the screen-rejection path does.
            blk = 1.0_dp
            hv = ecp_scalar_sph(blk, [ish, jsh], atm, natm, bas, nbas, env)

            if (hv == 0) then
               nscreened = nscreened + 1
               if (any(blk /= 0.0_dp)) nstale = nstale + 1
            end if

            ! The Cartesian entry point, which nothing else in the suite
            ! calls. Its values go through the same kernels as the spherical
            ! ones and are checked there; what is checked here is that it
            ! agrees about has_value and honours the same zeroing contract.
            dic = ((li + 1)*(li + 2))/2*bas(NCTR_OF + ish*BAS_SLOTS)
            djc = ((lj + 1)*(lj + 2))/2*bas(NCTR_OF + jsh*BAS_SLOTS)
            allocate (cblk(0:dic*djc - 1))
            cblk = 1.0_dp
            hvc = ecp_scalar_cart(cblk, [ish, jsh], atm, natm, bas, nbas, env)
            if (hvc /= hv) then
               write (*, '(a,a,a,i0,1x,i0)') "FAIL ", label, &
                  ": cart and sph disagree on has_value, shells ", ish, jsh
               bad = bad + 1
            end if
            if (hvc == 0 .and. any(cblk /= 0.0_dp)) then
               write (*, '(a,a,a,i0,1x,i0)') "FAIL ", label, &
                  ": cart left a rejected block unwritten, shells ", ish, jsh
               bad = bad + 1
            end if
            deallocate (cblk)

            i0 = ao_loc(ish)
            j0 = ao_loc(jsh)
            do i = 0, dj - 1
               mat(i0 + 1:i0 + di, j0 + 1 + i) = blk(i*di:i*di + di - 1)
            end do
            deallocate (blk)
         end do
      end do

      if (nstale > 0) then
         write (*, '(a,a,a,i0,a)') "FAIL ", label, ": ", nstale, &
            " block(s) reported has_value = 0 but were left unwritten"
         bad = bad + 1
      end if
      if (expect_reject /= 0 .and. nscreened == 0) then
         write (*, '(a,a,a)') "FAIL ", label, &
            ": the overlap screen rejected nothing, so that path went untested"
         bad = bad + 1
      end if

      scale = 0.0_dp
      do cc = 1, nao
         do r = 1, nao
            scale = max(scale, abs(ref((cc - 1)*nao + r)))
         end do
      end do
      ! An all-zero reference would make this vacuous rather than passing,
      ! which is worth failing on: asking for an ECP the basis does not define
      ! gives exactly that, and it reads as perfect agreement.
      if (scale <= 0.0_dp) then
         write (*, '(a,a,a)') "FAIL ", label, ": the reference matrix is identically zero"
         bad = bad + 1
         deallocate (mat)
         return
      end if

      wabs = 0.0_dp
      nzero = 0
      do cc = 1, nao
         do r = 1, nao
            denom = ref((cc - 1)*nao + r)
            d = abs(mat(r, cc) - denom)
            wabs = max(wabs, d)
            if (abs(denom) > 1.0e-12_dp) then
               if (d/abs(denom) > worst) worst = d/abs(denom)
            else if (d > TOL*scale) then
               nzero = nzero + 1
            end if
         end do
      end do

      write (*, '(a,a,a,i0,a,i0)') "  ", label, ": ", nao, "x", nao
      write (*, '(a,es11.3,a,es11.3)') "     largest |ref| ", scale, &
         "   worst/largest ", wabs/scale
      write (*, '(a,i0,a,i0)') "     screened pairs ", nscreened, &
         "   non-zero where PySCF is zero ", nzero

      if (wabs/scale > TOL) then
         write (*, '(a,a,a)') "FAIL ", label, ": outside tolerance"
         bad = bad + 1
      end if
      if (nzero > 0) then
         write (*, '(a,a,a)') "FAIL ", label, ": non-zero where PySCF is zero"
         bad = bad + 1
      end if
      deallocate (mat)
   end subroutine check_system

   !> env slots 18 and 19 left at zero: the documented "no ECP" contract
   !>
   !> cint_ecp_drv says an unset pair of slots is not an error but an empty
   !> potential, and that the integral comes back zero.  Nothing asserted it,
   !> and the early returns it exercises -- in ecp_loc_ecpbas and both
   !> Cartesian kernels -- are the ones a caller meets before it has an ECP.
   subroutine check_no_ecp(bad)
      integer, intent(inout) :: bad
      integer, parameter :: NATM = 1, NBAS = 1
      integer :: atm(0:5), bas(0:BAS_SLOTS - 1)
      real(dp) :: env(0:63), blk(0:0)
      integer :: hv

      atm = 0
      bas = 0
      env = 0.0_dp
      atm(1) = 20                       !! PTR_COORD
      bas(2) = 1                        !! NPRIM_OF
      bas(3) = 1                        !! NCTR_OF
      bas(5) = 23                       !! PTR_EXP
      bas(6) = 24                       !! PTR_COEFF
      env(23) = 1.0_dp
      env(24) = 1.0_dp
      env(AS_ECPBAS_OFFSET) = 0.0_dp
      env(AS_NECPBAS) = 0.0_dp

      blk = 1.0_dp                      !! poisoned
      hv = ecp_scalar_sph(blk, [0, 0], atm, NATM, bas, NBAS, env)
      if (hv /= 0) then
         write (*, '(a)') "FAIL no-ECP: has_value should be 0 when necpbas = 0"
         bad = bad + 1
      end if
      if (blk(0) /= 0.0_dp) then
         write (*, '(a)') "FAIL no-ECP: the block was not zeroed"
         bad = bad + 1
      end if
      write (*, '(a)') "  no ECP (necpbas = 0): zeroed, has_value = 0"
   end subroutine check_no_ecp

end program ecp_check
'''


if __name__ == "__main__":
    sys.exit(main())
