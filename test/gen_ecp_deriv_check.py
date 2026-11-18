#!/usr/bin/env python3
"""Generate test/ecp_deriv_check.F90 from PySCF.

    python3 test/gen_ecp_deriv_check.py

The first derivatives of the scalar ECP integrals, `ECPscalar_ipnuc` and
`ECPscalar_iprinv`, against `mol.intor`.  As with gen_ecp_check.py the
reference is committed, because libcint has no ECP code to link.

Every geometry is off the axes on purpose.  With the atoms on z, most x and y
components of a derivative vanish by symmetry, and a transposed or
mis-addressed component then agrees with the reference by being zero.

What each system is for:

  Cu/H in LANL2DZ, spherical and Cartesian
      s, p and d orbital shells, so the raised shell reaches f and all three
      branches of the s/p normalisation in `l_down` are taken;
      contracted shells; the one-atom integral on the atom that has an ECP and
      on the atom that has none.  The second has to come back zero and say
      so, and PySCF is not the reference for it: its entry point returns
      early there without writing its output.

  Na/Cl in LANL2DZ
      two atoms with an ECP, so the one-atom integral has to find a block of
      shells that does not start at row zero, and the two one-atom integrals
      have to add up to the whole.  An f function is added on chlorine, which no
      LANL2DZ set has: it is the only bra here above d, and so the only thing
      that reaches the general branch of `l_up` and a raised shell of g.

  I/H at 25 bohr in def2-SVP
      projected channels up to lc = 2, and an overlap screen that rejects.

The tolerance and how it is stated follow gen_ecp_check.py: relative to the
largest element of the matrix, 1e-11.

A reference array is written in pieces of at most CHUNK_LINES continuation
lines and joined in a constant expression.  One statement per array ran to 961
continuation lines, and a compiler is entitled to stop at 255; Intel's stops
at 511.
"""
import os
import sys

if os.environ.get("PYTHONHASHSEED") != "0":
    os.execve(
        sys.executable,
        [sys.executable, *sys.argv],
        {**os.environ, "PYTHONHASHSEED": "0"},
    )

HERE = os.path.dirname(os.path.abspath(__file__))

# tag, geometry (bohr), basis, label, spin, and the integrals to reference:
# (name, cart, atom) with atom None for ipnuc.
SYSTEMS = [
    ("cu", "Cu 0 0 0; H 0.4 -0.3 1.9", "lanl2dz", "Cu/H LANL2DZ", 0,
     [("ipnuc", False, None), ("ipnuc", True, None),
      ("iprinv", False, 0), ("iprinv", True, 0), ("iprinv", False, 1)]),
    ("nacl", "Na 0 0 0; Cl 0.5 0.3 4.4", "lanl2dz+f", "Na/Cl LANL2DZ+f", 0,
     [("ipnuc", False, None), ("iprinv", False, 0), ("iprinv", False, 1)]),
    ("i", "I 0 0 0; H 3.0 -4.0 24.5", "def2-svp", "I/H(far) def2-SVP", 0,
     [("ipnuc", False, None)]),
]


#: Continuation lines per statement, kept under the 255 the standard promises.
CHUNK_LINES = 200

#: The f primitive added to chlorine where a basis name ends in "+f".
EXTRA_F = [[3, [0.8, 1.0]]]


def reference_lines(name, values, per_line=3):
    """Declarations for one reference array, split so no statement is long."""
    size = CHUNK_LINES * per_line
    chunks = [values[i:i + size] for i in range(0, len(values), size)]
    if len(chunks) == 1:
        return ["   real(dp), parameter :: %s(%d) = [ &" % (name, len(values)),
                wrap(values, per_line=per_line) + "]"]
    lines = []
    for k, chunk in enumerate(chunks):
        lines += ["   real(dp), parameter :: %s_P%d(%d) = [ &" % (name, k, len(chunk)),
                  wrap(chunk, per_line=per_line) + "]"]
    lines.append("   real(dp), parameter :: %s(%d) = [ &" % (name, len(values)))
    lines.append(wrap(["%s_P%d" % (name, k) for k in range(len(chunks))], per_line=4) + "]")
    return lines


def wrap(values, indent="      ", per_line=3):
    out = []
    for i in range(0, len(values), per_line):
        out.append(indent + ", ".join(values[i:i + per_line]))
    return ", &\n".join(out)


def emit_system(tag, geom, basis, label, spin, integrals):
    import numpy as np
    from pyscf import gto

    ecp = basis
    orbital = basis
    if basis.endswith("+f"):
        ecp = basis[:-2]
        orbital = {"Na": ecp, "Cl": gto.basis.load(ecp, "Cl") + EXTRA_F}
    mol = gto.M(atom=geom, basis=orbital, ecp=ecp, spin=spin, unit="bohr", verbose=0)
    mol.intor("ECPscalar_sph")                # fills env[18], env[19]
    bas = np.vstack((mol._bas, mol._ecpbas))
    u = tag.upper()
    lines = [
        "",
        "   ! ---- %s : %s (bohr) in %s ----" % (label, geom, basis),
        "   integer, parameter :: %s_NATM = %d" % (u, mol.natm),
        "   integer, parameter :: %s_NBAS_ALL = %d" % (u, bas.shape[0]),
        "   integer, parameter :: %s_NBAS = %d" % (u, mol.nbas),
        "   integer, parameter :: %s_NENV = %d" % (u, len(mol._env)),
        "",
        "   integer, parameter :: %s_ATM(0:%s_NATM*6 - 1) = [ &" % (u, u),
        wrap([str(int(x)) for x in np.asarray(mol._atm).ravel()], per_line=12) + "]",
        "",
        "   integer, parameter :: %s_BAS(0:%s_NBAS_ALL*BAS_SLOTS - 1) = [ &" % (u, u),
        wrap([str(int(x)) for x in bas.ravel()], per_line=12) + "]",
        "",
        "   real(dp), parameter :: %s_ENV(0:%s_NENV - 1) = [ &" % (u, u),
        wrap(["%.17e_dp" % v for v in mol._env]) + "]",
    ]
    calls = []
    for k, (name, cart, atom) in enumerate(integrals):
        suffix = "cart" if cart else "sph"
        intor = "ECPscalar_%s_%s" % (name, suffix)
        if atom is None:
            ref = mol.intor(intor, comp=3)
        else:
            with mol.with_rinv_at_nucleus(atom):
                ref = mol.intor(intor, comp=3)
        nao = ref.shape[1]
        # PySCF returns early for an atom with no ECP and leaves its output
        # unwritten, so what comes back is whatever the buffer held.  The
        # reference for that case is zero by definition, not by PySCF.
        expect_zero = atom is not None and atom not in set(mol._ecpbas[:, 0])
        ref_name = "%s_REF%d" % (u, k)
        lines += ["", "   ! %s%s" % (intor, "" if atom is None else ", atom %d" % atom)]
        if expect_zero:
            lines.append("   real(dp), parameter :: %s(%d) = 0.0_dp" % (ref_name, ref.size))
        else:
            # (comp, i, j) -> comp slowest, then column-major within a component
            lines += reference_lines(
                ref_name, ["%.17e_dp" % v for v in ref.transpose(0, 2, 1).ravel()])
        what = "%s %s%s" % (label, intor.replace("ECPscalar_", ""),
                            "" if atom is None else " atom %d" % atom)
        calls.append(
            '   call check("%s", %s, %s, %d, %s_NATM, %s_NBAS, %d, &\n'
            '              %s_ATM, %s_BAS, %s_ENV, %s, %s, bad)'
            % (what, ".true." if name == "iprinv" else ".false.",
               ".true." if cart else ".false.", -1 if atom is None else atom,
               u, u, nao, u, u, u, ref_name,
               ".true." if expect_zero else ".false."))
    if sum(1 for n, c, a in integrals if n == "iprinv" and not c) >= mol.natm:
        calls.append(
            '   call check_sum("%s", %s_NATM, %s_NBAS, %s_ATM, %s_BAS, %s_ENV, bad)'
            % (label, u, u, u, u, u))
    return "\n".join(lines), calls


def main():
    blocks, calls = [], []
    for system in SYSTEMS:
        blk, cs = emit_system(*system)
        blocks.append(blk)
        calls += cs
    src = PROGRAM % {
        "pyscf_version": __import__("pyscf").__version__,
        "data": "\n".join(blocks),
        "calls": "\n".join(calls),
        "abi_calls": "\n".join(
            '   call check_layers("%s", %d, %s_NATM, %s_NBAS, %s_ATM, %s_BAS, %s_ENV, bad)'
            % (label, atom, u, u, u, u, u)
            for label, atom, u in (("Cu/H LANL2DZ", 0, "CU"), ("Cu/H LANL2DZ", 1, "CU"),
                                   ("Na/Cl LANL2DZ+f", 1, "NACL"))),
    }
    with open(os.path.join(HERE, "ecp_deriv_check.F90"), "w") as fh:
        fh.write(src)
    print("%d systems, %d comparisons -> ecp_deriv_check.F90" % (len(SYSTEMS), len(calls)))


PROGRAM = '''!
! First derivatives of the scalar ECP integrals against PySCF.
!
! GENERATED by test/gen_ecp_deriv_check.py -- do not edit.
! Reference: PySCF %(pyscf_version)s, mol.intor('ECPscalar_ipnuc_sph', comp=3) and
! its Cartesian and one-atom (iprinv) forms.
!
! The reference is committed because libcint has no ECP code to link; see
! gen_ecp_deriv_check.py for what each system is there to reach.
!
program ecp_deriv_check
   use cint_const, only: dp
   use cint_bas, only: ANG_OF, NCTR_OF, BAS_SLOTS
   use cint_ecp_num, only: AS_RINV_ORIG_ATOM, AS_NECPBAS
   use cint_ecp_deriv, only: ecp_ipnuc_sph, ecp_ipnuc_cart, ecp_iprinv_sph, ecp_iprinv_cart
   use libcint_fortran, only: libcint_ecp_ipnuc_sph, libcint_ecp_ipnuc_cart, &
                              libcint_ecp_iprinv_sph, libcint_ecp_iprinv_cart
   use, intrinsic :: iso_c_binding, only: c_ptr, c_int, c_double, c_loc, c_null_ptr
   implicit none

   ! The C entry points, under the names PySCF's library exports.  Declared
   ! here rather than used from a module so that a rename fails at link.
   interface
      function c_ipnuc_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
         result(ret) bind(C, name="ECPscalar_ipnuc_sph")
         import :: c_ptr, c_int
         type(c_ptr), value :: out, dims, shls, atm, bas, env, opt, cache
         integer(c_int), value :: natm, nbas
         integer(c_int) :: ret
      end function c_ipnuc_sph
      function c_ipnuc_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
         result(ret) bind(C, name="ECPscalar_ipnuc_cart")
         import :: c_ptr, c_int
         type(c_ptr), value :: out, dims, shls, atm, bas, env, opt, cache
         integer(c_int), value :: natm, nbas
         integer(c_int) :: ret
      end function c_ipnuc_cart
      function c_iprinv_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
         result(ret) bind(C, name="ECPscalar_iprinv_sph")
         import :: c_ptr, c_int
         type(c_ptr), value :: out, dims, shls, atm, bas, env, opt, cache
         integer(c_int), value :: natm, nbas
         integer(c_int) :: ret
      end function c_iprinv_sph
      function c_iprinv_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
         result(ret) bind(C, name="ECPscalar_iprinv_cart")
         import :: c_ptr, c_int
         type(c_ptr), value :: out, dims, shls, atm, bas, env, opt, cache
         integer(c_int), value :: natm, nbas
         integer(c_int) :: ret
      end function c_iprinv_cart
   end interface

   !> Relative to the largest element of the matrix, not to each element
   real(dp), parameter :: TOL = 1.0e-11_dp
%(data)s

   integer :: bad

   bad = 0
%(calls)s
%(abi_calls)s
   call check_no_ecp(CU_NATM, CU_NBAS, CU_ATM, CU_BAS, CU_ENV, bad)

   if (bad > 0) then
      write (*, '(a,i0,a)') "ecp_deriv_check: FAILED (", bad, " problem(s))"
      stop 1
   end if
   write (*, '(a)') "ecp_deriv_check: OK"

contains

   !> Functions in a shell, contractions included
   pure function shell_dim(cart, sh, bas) result(n)
      logical, intent(in) :: cart
      integer, intent(in) :: sh, bas(0:)
      integer :: n, l

      l = bas(ANG_OF + sh*BAS_SLOTS)
      if (cart) then
         n = ((l + 1)*(l + 2))/2
      else
         n = 2*l + 1
      end if
      n = n*bas(NCTR_OF + sh*BAS_SLOTS)
   end function shell_dim

   !> One shell pair through the entry point the flags select
   function one_block(one_atom, cart, blk, shls, atm, natm, bas, nbas, env, dims) result(hv)
      logical, intent(in) :: one_atom, cart
      real(dp), intent(inout) :: blk(0:)
      integer, intent(in) :: shls(0:1), atm(0:), natm, bas(0:), nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: hv

      if (one_atom .and. cart) then
         hv = ecp_iprinv_cart(blk, shls, atm, natm, bas, nbas, env, dims)
      else if (one_atom) then
         hv = ecp_iprinv_sph(blk, shls, atm, natm, bas, nbas, env, dims)
      else if (cart) then
         hv = ecp_ipnuc_cart(blk, shls, atm, natm, bas, nbas, env, dims)
      else
         hv = ecp_ipnuc_sph(blk, shls, atm, natm, bas, nbas, env, dims)
      end if
   end function one_block

   !> The whole three-component matrix, assembled shell pair by shell pair
   !>
   !> Every block is computed twice: without `dims` into a poisoned buffer, and
   !> with `dims` straight into a poisoned matrix.  The two have to agree to
   !> the last bit, being one computation placed two ways.
   subroutine build(one_atom, cart, atom, natm, nbas, nao, atm, bas, env_in, mat, &
                    nscreened, nproblem)
      logical, intent(in) :: one_atom, cart
      integer, intent(in) :: atom, natm, nbas, nao
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env_in(0:)
      real(dp), intent(out) :: mat(nao, nao, 3)
      integer, intent(out) :: nscreened, nproblem

      real(dp), allocatable :: env(:), blk(:), placed(:)
      integer :: ish, jsh, di, dj, i0, j0, j, c, hv, hv2

      allocate (env(0:size(env_in) - 1))
      env = env_in
      if (one_atom) env(AS_RINV_ORIG_ATOM) = real(atom, dp)

      nscreened = 0
      nproblem = 0
      mat = 0.0_dp
      allocate (placed(0:3*nao*nao - 1))
      i0 = 0
      do ish = 0, nbas - 1
         di = shell_dim(cart, ish, bas)
         j0 = 0
         do jsh = 0, nbas - 1
            dj = shell_dim(cart, jsh, bas)
            allocate (blk(0:3*di*dj - 1))
            blk = 1.0_dp
            hv = one_block(one_atom, cart, blk, [ish, jsh], atm, natm, bas, nbas, env)
            if (hv == 0) then
               nscreened = nscreened + 1
               if (any(blk /= 0.0_dp)) nproblem = nproblem + 1
            end if

            placed = 1.0_dp
            hv2 = one_block(one_atom, cart, placed, [ish, jsh], atm, natm, bas, nbas, env, &
                            [nao, nao])
            if (hv2 /= hv) nproblem = nproblem + 1
            do c = 0, 2
               do j = 0, dj - 1
                  mat(i0 + 1:i0 + di, j0 + 1 + j, c + 1) = blk(c*di*dj + j*di:c*di*dj + j*di + di - 1)
                  if (any(placed(c*nao*nao + j*nao:c*nao*nao + j*nao + di - 1) &
                          /= blk(c*di*dj + j*di:c*di*dj + j*di + di - 1))) nproblem = nproblem + 1
               end do
            end do
            deallocate (blk)
            j0 = j0 + dj
         end do
         i0 = i0 + di
      end do
   end subroutine build

   subroutine check(label, one_atom, cart, atom, natm, nbas, nao, atm, bas, env, ref, &
                    expect_zero, bad)
      character(len=*), intent(in) :: label
      logical, intent(in) :: one_atom, cart, expect_zero
      integer, intent(in) :: atom, natm, nbas, nao
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:), ref(:)
      integer, intent(inout) :: bad

      real(dp), allocatable :: mat(:, :, :)
      real(dp) :: scale, wabs
      integer :: nscreened, nproblem

      allocate (mat(nao, nao, 3))
      call build(one_atom, cart, atom, natm, nbas, nao, atm, bas, env, mat, nscreened, nproblem)
      if (nproblem > 0) then
         write (*, '(a,a,a,i0,a)') "FAIL ", label, ": ", nproblem, &
            " block(s) unwritten when rejected, or placed differently with dims"
         bad = bad + 1
      end if

      scale = maxval(abs(ref))
      wabs = maxval(abs(reshape(mat, [size(ref)]) - ref))
      write (*, '(a,a)') "  ", label
      write (*, '(a,es11.3,a,es11.3,a,i0)') "     largest |ref| ", scale, &
         "   worst ", wabs, "   screened pairs ", nscreened

      if (expect_zero) then
         ! An atom with no ECP: the reference is identically zero, and every
         ! block has to come back zero and report that nothing was added.
         if (wabs /= 0.0_dp .or. nscreened /= nbas*nbas) then
            write (*, '(a,a,a)') "FAIL ", label, ": an empty potential gave a value"
            bad = bad + 1
         end if
      else if (scale <= 0.0_dp) then
         write (*, '(a,a,a)') "FAIL ", label, ": the reference is identically zero"
         bad = bad + 1
      else if (wabs > TOL*scale) then
         write (*, '(a,a,a)') "FAIL ", label, ": outside tolerance"
         bad = bad + 1
      end if
   end subroutine check

   !> The one-atom integrals add up to the whole
   subroutine check_sum(label, natm, nbas, atm, bas, env, bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: natm, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, intent(inout) :: bad

      real(dp), allocatable :: whole(:, :, :), part(:, :, :), total(:, :, :)
      integer :: nao, sh, atom, nscreened, nproblem
      real(dp) :: d

      nao = 0
      do sh = 0, nbas - 1
         nao = nao + shell_dim(.false., sh, bas)
      end do
      allocate (whole(nao, nao, 3), part(nao, nao, 3), total(nao, nao, 3))
      call build(.false., .false., -1, natm, nbas, nao, atm, bas, env, whole, nscreened, nproblem)
      total = 0.0_dp
      do atom = 0, natm - 1
         call build(.true., .false., atom, natm, nbas, nao, atm, bas, env, part, nscreened, nproblem)
         total = total + part
      end do
      d = maxval(abs(total - whole))/maxval(abs(whole))
      write (*, '(a,a,a,es11.3)') "  ", label, ": one-atom integrals summed against the whole ", d
      if (d > TOL) then
         write (*, '(a,a,a)') "FAIL ", label, ": the one-atom integrals do not add up"
         bad = bad + 1
      end if
   end subroutine check_sum

   !> The two outer layers against the module they wrap, bit for bit
   !>
   !> `check` drives cint_ecp_deriv directly, which leaves the bind(C) entry
   !> points and the libcint_fortran wrappers called by nothing.  Both only
   !> size and place, so they have to reproduce the module exactly -- with and
   !> without `dims` for the C one.
   subroutine check_layers(label, atom, natm, nbas, atm, bas, env_in, bad)
      character(len=*), intent(in) :: label
      integer, intent(in) :: atom, natm, nbas
      integer, intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env_in(0:)
      integer, intent(inout) :: bad

      integer(c_int), allocatable, target :: catm(:), cbas(:)
      integer(c_int), target :: cshls(2), cdims(2)
      real(c_double), allocatable, target :: env(:), cblk(:), cplaced(:)
      real(dp), allocatable :: blk(:), fblk(:)
      logical :: one_atom, cart
      integer :: kind, ish, jsh, di, dj, hv, hvc, hvd, hvf, c, j, nwrong

      allocate (catm(size(atm)), cbas(size(bas)), env(0:size(env_in) - 1))
      catm = int(atm, c_int)
      cbas = int(bas, c_int)
      env = env_in
      env(AS_RINV_ORIG_ATOM) = real(atom, dp)

      nwrong = 0
      do kind = 1, 4
         one_atom = kind > 2
         cart = mod(kind, 2) == 0
         do ish = 0, nbas - 1
            do jsh = 0, nbas - 1
               di = shell_dim(cart, ish, bas)
               dj = shell_dim(cart, jsh, bas)
               allocate (blk(0:3*di*dj - 1), cblk(3*di*dj), fblk(3*di*dj))
               allocate (cplaced(3*(di + 1)*(dj + 2)))
               blk = 1.0_dp; cblk = 2.0_dp; fblk = 3.0_dp; cplaced = 4.0_dp
               cshls = int([ish, jsh], c_int)
               cdims = int([di + 1, dj + 2], c_int)
               hv = one_block(one_atom, cart, blk, [ish, jsh], atm, natm, bas, nbas, env)
               select case (kind)
               case (1)
                  hvc = c_ipnuc_sph(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                    int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                    c_loc(env), c_null_ptr, c_null_ptr)
                  hvd = c_ipnuc_sph(c_loc(cplaced), c_loc(cdims), c_loc(cshls), c_loc(catm), &
                                    int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                    c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_ipnuc_sph(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case (2)
                  hvc = c_ipnuc_cart(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvd = c_ipnuc_cart(c_loc(cplaced), c_loc(cdims), c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_ipnuc_cart(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case (3)
                  hvc = c_iprinv_sph(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvd = c_iprinv_sph(c_loc(cplaced), c_loc(cdims), c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_iprinv_sph(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case default
                  hvc = c_iprinv_cart(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                      int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                      c_loc(env), c_null_ptr, c_null_ptr)
                  hvd = c_iprinv_cart(c_loc(cplaced), c_loc(cdims), c_loc(cshls), c_loc(catm), &
                                      int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                      c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_iprinv_cart(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               end select
               if (hvc /= hv .or. hvd /= hv .or. hvf /= hv) nwrong = nwrong + 1
               if (any(cblk /= blk) .or. any(fblk /= blk)) nwrong = nwrong + 1
               do c = 0, 2
                  do j = 0, dj - 1
                     if (any(cplaced(c*(di + 1)*(dj + 2) + j*(di + 1) + 1: &
                                     c*(di + 1)*(dj + 2) + j*(di + 1) + di) &
                             /= blk(c*di*dj + j*di:c*di*dj + j*di + di - 1))) nwrong = nwrong + 1
                  end do
               end do
               ! Outside the block's corner the caller's matrix is the caller's:
               ! the extra row and the two extra columns keep their poison.
               do c = 0, 2
                  do j = 0, dj + 1
                     if (cplaced(c*(di + 1)*(dj + 2) + j*(di + 1) + di + 1) /= 4.0_dp) &
                        nwrong = nwrong + 1
                     if (j >= dj) then
                        if (any(cplaced(c*(di + 1)*(dj + 2) + j*(di + 1) + 1: &
                                        c*(di + 1)*(dj + 2) + j*(di + 1) + di) /= 4.0_dp)) &
                           nwrong = nwrong + 1
                     end if
                  end do
               end do
               deallocate (blk, cblk, fblk, cplaced)
            end do
         end do
      end do

      write (*, '(a,a,a,i0,a,i0)') "  ", label, ": C ABI and libcint_fortran, atom ", atom, &
         ", blocks differing from the module ", nwrong
      if (nwrong > 0) then
         write (*, '(a,a,a)') "FAIL ", label, ": an outer layer does not reproduce the module"
         bad = bad + 1
      end if
   end subroutine check_layers

   !> A molecule with no ECP rows: every layer zeroes its block and returns 0
   !>
   !> Cu/H with env's row count set to zero, and `bas` cut back to the orbital
   !> shells so that a layer reaching past them is a bounds error, not luck.
   subroutine check_no_ecp(natm, nbas, atm, bas_in, env_in, bad)
      integer, intent(in) :: natm, nbas
      integer, intent(in) :: atm(0:), bas_in(0:)
      real(dp), intent(in) :: env_in(0:)
      integer, intent(inout) :: bad

      integer(c_int), allocatable, target :: catm(:), cbas(:)
      integer(c_int), target :: cshls(2)
      real(c_double), allocatable, target :: env(:), cblk(:)
      integer, allocatable :: bas(:)
      real(dp), allocatable :: blk(:), fblk(:)
      integer :: kind, ish, jsh, di, dj, hv, hvc, hvf, nwrong
      logical :: one_atom, cart

      allocate (bas(0:nbas*BAS_SLOTS - 1), env(0:size(env_in) - 1))
      bas = bas_in(0:nbas*BAS_SLOTS - 1)
      env = env_in
      env(AS_NECPBAS) = 0.0_dp
      env(AS_RINV_ORIG_ATOM) = 0.0_dp
      allocate (catm(size(atm)), cbas(size(bas)))
      catm = int(atm, c_int)
      cbas = int(bas, c_int)

      nwrong = 0
      do kind = 1, 4
         one_atom = kind > 2
         cart = mod(kind, 2) == 0
         do ish = 0, nbas - 1
            do jsh = 0, nbas - 1
               di = shell_dim(cart, ish, bas)
               dj = shell_dim(cart, jsh, bas)
               allocate (blk(0:3*di*dj - 1), cblk(3*di*dj), fblk(3*di*dj))
               blk = 1.0_dp; cblk = 2.0_dp; fblk = 3.0_dp
               cshls = int([ish, jsh], c_int)
               hv = one_block(one_atom, cart, blk, [ish, jsh], atm, natm, bas, nbas, env)
               select case (kind)
               case (1)
                  hvc = c_ipnuc_sph(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                    int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                    c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_ipnuc_sph(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case (2)
                  hvc = c_ipnuc_cart(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_ipnuc_cart(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case (3)
                  hvc = c_iprinv_sph(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                     int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                     c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_iprinv_sph(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               case default
                  hvc = c_iprinv_cart(c_loc(cblk), c_null_ptr, c_loc(cshls), c_loc(catm), &
                                      int(natm, c_int), c_loc(cbas), int(nbas, c_int), &
                                      c_loc(env), c_null_ptr, c_null_ptr)
                  hvf = libcint_ecp_iprinv_cart(fblk, [ish, jsh], atm, natm, bas, nbas, env)
               end select
               if (hv /= 0 .or. hvc /= 0 .or. hvf /= 0) nwrong = nwrong + 1
               if (any(blk /= 0.0_dp) .or. any(cblk /= 0.0_dp) .or. any(fblk /= 0.0_dp)) &
                  nwrong = nwrong + 1
               deallocate (blk, cblk, fblk)
            end do
         end do
      end do

      write (*, '(a,i0)') "  no ECP (necpbas = 0): blocks not zeroed or not reported empty ", nwrong
      if (nwrong > 0) then
         write (*, '(a)') "FAIL no ECP: a layer gave a value for an empty potential"
         bad = bad + 1
      end if
   end subroutine check_no_ecp

end program ecp_deriv_check
'''


if __name__ == "__main__":
    sys.exit(main())
