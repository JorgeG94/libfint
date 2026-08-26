# Fortran port of libcint.
#
# So far: the numeric tables (D2) and the Rys quadrature subsystem (D3).
# The tables are generated from the C by scripts/c_tables_to_fortran.py and the
# closed-form Rys roots by scripts/c_rysroot_to_fortran.py; both outputs are
# committed, so a build needs no Python.

set(cintFortranSrc
  cint_const.f90
  cint_ecp_num.f90
  cint_dd.f90
  cint_bas.f90
  cint_blas.f90
  cint_fmt.F90
  cint_find_roots.f90
  cint_schmidt.F90
  cint_eigh_dstemr.f90
  cint_eigh.F90
  cint_rys_root_n.f90
  cint_tab_roots_x0.f90
  cint_tab_jacobi.f90
  cint_tab_jacobi_ext.f90
  cint_wheeler.F90
  cint_rys_roots.f90
  cint_tab_cart2sph.f90
  cint_cart2sph.f90
  cint_cart2spinor.f90
  cint_ecp_ang.f90
  cint_ecp_rad.f90
  cint_ecp.f90
  cint_ecp_drv.f90
  cint_envs.f90
  cint_screen.f90
  cint_workspace.f90
  cint_g1e.f90
  cint_1e.f90
  cint_1e_spinor.f90
  cint_1e_grids.f90
  cint_g2e_unrolled.f90
  cint_g2e.f90
  cint_2e.f90
  cint_g2e_ops.f90
  cint_3c2e.f90
  cint_2e_spinor.f90
  cint_3c2e_spinor.f90
  cint_1e_grids_spinor.f90
  # The generated catalogue: one Fortran module per output file of
  # scripts/auto_intor.cl, from the same descriptions the C is generated from.
  cint_gen_intor1.f90
  cint_gen_intor2.f90
  cint_gen_intor3.f90
  cint_gen_intor4.f90
  cint_gen_dkb.f90
  cint_gen_grad1.f90
  cint_gen_grad2.f90
  cint_gen_gaunt1.f90
  cint_gen_breit1.f90
  cint_gen_hess.f90
  cint_gen_int3c2e.f90
  cint_gen_deriv3.f90
  cint_gen_deriv4.f90
  cint_gen_lresc.f90
  cint_gen_int1e_grids1.f90
  cint_gen_int3c1e.f90
  cint_gen_a.f90
  cint_3c1e.f90
  cint_opt.f90
  cint_breit_gauge.f90
  # libcint's C entry points under their C names, for consumers that bind
  # them directly instead of using the module.  See the file for which.
  cint_c_abi.f90
  # Always listed.  The F12 three are guarded internally by COMPILE_F12 and
  # compile to nothing without it; they sit in fortran/f12 rather than
  # fortran/src because fpm scans `use` statements lexically -- it cannot see
  # a preprocessor conditional, so a guarded `use cint_tab_stg_roots` still
  # reads to it as a dependency on a table that is not in the repository.  A
  # directory it does not look in is the only thing that hides them.
  #
  # The two polynomial-fit tables are out of fortran/src for a different
  # reason: their array constructors exceed gfortran's 65535 limit and need
  # -fmax-array-constructor, which fpm accepts only through --flag, and
  # --flag REPLACES the profile flags rather than adding to them.  Passing it
  # would cost -fcheck=bounds and -fbacktrace, which is what found the
  # short-circuit bug two commits ago.  They are committed leaf modules that
  # nothing USEs unless the fit is enabled, so their being compiled always
  # costs archive size and nothing else.
  #
  # What all this buys is a source list that does not vary, which is a source
  # list nothing can drift from -- and fpm can then mirror fortran/src with a
  # single directory symlink instead of a hand-maintained file list.  See
  # scripts/check_f90_integral_coverage.py for the bug that argument comes
  # from.
  ../f12/cint_stg_roots.F90
  ../f12/cint_g2e_f12.F90
  ../f12/cint_gen_f12.F90
  # The public API.  It lives in include/ rather than src/ because it is
  # what a consumer reads, but it has to be IN the library: without it,
  # find_package(libfint) hands back a target whose headline module --
  # libcint_fortran, the one every consumer writes `use` for -- was
  # never compiled.
  ../include/libcint_interface.f90
  ../include/libfint.f90)

# The range-separation tables are 91% of all the table data and, like the C
# they come from, sit behind a switch.  WITH_POLYNOMIAL_FIT selects a fast path
# in CINTsr_rys_roots for 3..5 roots that falls through to the general solver
# when absent, so this costs accuracy nowhere -- only speed, and only for
# range-separated functionals.  It is off by default upstream and off here.
#
# The cost is worth knowing: with it on, these tables take ~27 s and ~1.5 GB to
# compile; without, the whole set is ~1.9 s and ~134 MB.
#
# That was the intent; until now it was only the comment.  WITH_POLYNOMIAL_FIT
# was never declared and both tables were listed unconditionally above, so the
# switch described here has never existed and the tables have always been
# compiled.  cint_tab_sr_roots.f90 alone is 29.6 MB of source and 49.5 MB of
# the 65.4 MB archive.
#
# Note what the switch does *not* do yet.  Nothing under src/ or include/ uses
# either module -- cint_rys_roots_sr (src/cint_rys_roots.f90:156) goes to
# Schmidt, Jacobi and Laguerre for every root count, and never consults a
# fitted table.  So the fast path these tables exist to serve is not wired up
# on the Fortran side, and turning this ON today buys compile time and archive
# size and nothing else.  It is declared anyway, and defaulted OFF to match
# upstream, so that the day CINTsr_rys_roots grows its fitted path the tables
# are already gated the way the C gates them.
option(WITH_POLYNOMIAL_FIT
       "Range-separated Rys root tables (29.6 MB of source; no consumer yet)" OFF)
if(WITH_POLYNOMIAL_FIT)
  list(APPEND cintFortranSrc
       ../polyfit/cint_tab_polyfits.f90
       ../polyfit/cint_tab_sr_roots.f90)
endif()
# F12: the Yukawa potential and the Slater-type geminal, off by default here
# because they are off by default upstream.
#
if(WITH_F12)
  # The quadrature table is committed, in f12/, and split across two files
  # because together they are 108 MB of Fortran and GitHub refuses any single
  # file over 100 MB. `cint_tab_stg_roots` is a facade over the two so that
  # `cint_stg_roots` needs one name and does not care about the split.
  #
  # It used to be generated at build time from libcint's roots_xw.dat, which
  # made WITH_F12 unbuildable for anyone who had libfint without libcint --
  # every consumer that fetches this project, since libcint is only ever
  # fetched to validate against in CI. scripts/stg_roots_to_fortran.py still
  # exists and still regenerates these three files; nothing calls it during a
  # build.
  #
  # The cost is real and measured: 106 s and 3.46 GB to compile, taking
  # libfint.a from 13.5 MB to 93.5 MB. That is why the option is off by
  # default rather than merely available.
  list(APPEND cintFortranSrc
    ${PROJECT_SOURCE_DIR}/f12/cint_tab_stg_x.f90
    ${PROJECT_SOURCE_DIR}/f12/cint_tab_stg_w.f90
    ${PROJECT_SOURCE_DIR}/f12/cint_tab_stg_roots.f90)
endif()

# STATIC, not OBJECT: CMake does not propagate one OBJECT library's objects
# through another to the final consumer, and with WITH_FORTRAN_BACKEND the
# libcint_fortran module sits in between.
add_library(fint STATIC ${cintFortranSrc})

# THREAD SAFETY.  libcint is thread-safe -- every call takes its scratch off
# the C stack -- and this library has to be too, because callers drive it from
# inside `!$omp parallel`.  It gets there differently: the workspace is a
# module variable, reused so that a shell quartet costs no allocation, and
# marked `!$omp threadprivate`.
#
# That directive is a COMMENT unless the file is compiled with OpenMP enabled.
# Without this block the library builds, links, and corrupts the heap the
# moment two threads evaluate integrals at once -- which is exactly how it was
# found, as a "double free or corruption" in a caller with no bug in it.
#
# So OpenMP is not optional here, and not merely "found if available": if it is
# missing the threadprivate workspaces are silently shared and the library is
# quietly wrong under threads.
find_package(OpenMP COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(fint PUBLIC OpenMP::OpenMP_Fortran)
else()
  message(WARNING
    "OpenMP not found: the per-thread workspaces become shared, so this build "
    "is only safe to call from one thread at a time.")
endif()
if(WITH_F12)
  target_compile_definitions(fint PRIVATE COMPILE_F12)
endif()
set_target_properties(fint PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
# BUILD_INTERFACE, not a bare path: fint is in the install export
# set, and a bare build-tree path there is a hard CMake error at generate
# time -- it would be baked into cintTargets.cmake and point at a directory
# the consumer does not have.
target_include_directories(fint PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
# cint_fmt/cint_schmidt/cint_wheeler include their algorithm bodies once per
# precision; the .inc files live beside them in src/.
target_include_directories(fint PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# gfortran caps an array constructor at 65535 elements; the largest table here
# has 223,440.  ifx has no such cap.  nvfortran has neither this flag nor
# real128, which is the constraint recorded in PORT_TO_FORTRAN.md 3.1.
# ifx and icx default to a relaxed floating-point model, and it is not a
# small effect: built that way the port disagrees with the C in 75,593 of
# int2e_check's 209,254 values, worst 4.2e-10.  That is reassociation and
# contraction, not a bug in either -- but bit-identity is this port's whole
# acceptance criterion, so the flag is not optional here.  With it, the
# Intel build is bit-identical again.
#
# Nothing equivalent is needed for GNU, which is strict by default.
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  target_compile_options(fint PRIVATE -fp-model=precise)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(fint PRIVATE -fmax-array-constructor=4000000)
  # -fno-tree-vectorize is worth 12% here, which is the opposite of what a
  # flag called that usually does.  The inner loops of the gout kernels index
  # g by three values pulled out of idx -- g(ix)*g(iy)*g(iz) -- and gfortran
  # turns that into AVX gathers.  On this Broadwell a gather is slower than
  # the scalar loads it replaces, so the vectorised form loses; the port is
  # slower at -O3 than at -O2 for exactly this reason.  gcc, given the same
  # source in C, does not vectorise it at all: the C's CINTgout2e is 646
  # instructions of plain mulsd and addsd against the port's 423.
  #
  # It changes no results.  Vectorisation without -ffast-math preserves
  # floating-point semantics, so switching it off preserves them too, and the
  # catalogue sweep is bit-identical either way.
  #
  # Measured alternatives, all rejected.  Release, 5 reps over 518,400
  # quartets, the port's own CPU time; the C is 1.95 s on the same run:
  #
  #   -O3                                   3.199 s   1.64x    <- was
  #   -O3 -fno-tree-vectorize               2.825 s   1.45x    <- is
  #   -O3 -funroll-loops                    3.213 s   1.65x
  #   -O3 -march=native                     3.223 s   1.65x
  #   -O3 -march=native -mno-gather         3.365 s   1.73x
  #   -O3 -march=native -fno-tree-vectorize 2.718 s   1.40x
  #
  # The last one is faster still, and is not taken: -march=native lets
  # gfortran contract FMAs where gcc does not, and the catalogue then differs
  # from the C in 371,752 values.  Bit-identity is the contract, and 3.8% does
  # not buy it out.  -mno-gather was the obvious repair for the gather
  # finding and made things worse, which is why the flag here is the blunt
  # one: it is not the gathers alone, it is the whole vectorised shape.
  # ...and the reason above is half wrong, which matters if anyone retries
  # this.  The 371,752 differing values were FMA contraction, not
  # vectorisation: -ftree-vectorize -ffp-contract=off is bit-identical
  # across the whole suite.  The flag stays anyway, because measured that
  # way gfortran 13's vectoriser is a large regression -- int2e_sph 3,055M
  # instructions scalar against 3,964M vectorised -- including on the very
  # gout kernel gcc vectorises well.  See PERFORMANCE.md 7b.
  target_compile_options(fint PRIVATE -fno-tree-vectorize)
  # ...except where the loops are the shape a vectoriser is actually for.
  # cint_g1e holds prim_to_ctr_0 and _1, whose inner loop is an axpy --
  # gc(nf*j+n) += c0 * gp(n) -- and gcc vectorises the C's freely.  The
  # blanket flag above is right for the gout kernels, which index g by three
  # values out of idx and get gathers; it is wrong here.
  set_source_files_properties(cint_g1e.f90 PROPERTIES
                              COMPILE_OPTIONS "-ftree-vectorize")

  # LTO is worth 15% on the two-electron path, and it is worth it for a
  # reason that has nothing to do with cross-file optimisation in general.
  # The C keeps rys_root1..5 `static` in rys_roots.c, so gcc inlines them
  # into CINTrys_roots.  The port has them in cint_rys_root_n, a different
  # module from cint_rys_roots, and gfortran cannot inline across a module
  # boundary without seeing both at once.  Measured on the Rys layer alone,
  # nroots 1 to 5, the port was 1.5x to 1.8x the C for arithmetic that is
  # transcribed line for line.
  #
  # Submodules do not fix this.  A submodule is a separate translation unit
  # by design -- that is the point of them, breaking the recompilation
  # dependency -- and a caller in the parent module emits an undefined
  # reference to the submodule's procedure, exactly as it would across two
  # ordinary modules.  Only one giant module in one file would let gfortran
  # inline these without LTO, and LTO buys the same thing without asking the
  # port to be unreadable.
  #
  #   -O3 -fno-tree-vectorize             2.825 s   1.45x
  #   -O3 -fno-tree-vectorize -flto       2.395 s   1.23x
  #
  # -ffat-lto-objects so the archive still links under a toolchain with no
  # LTO plugin; it just does not get the inlining.  With GNU's plugin, which
  # is on by default, a consumer picks the benefit up without having to pass
  # -flto itself -- verified by linking the benchmark without it.
  #
  # The cost is link time, not compile time: building this library is 19.3 s
  # either way, and each executable that links it pays about 45 s more.
  # WITH_FORTRAN_LTO=OFF for anyone who would rather have the link back.
  option(WITH_FORTRAN_LTO "Link-time optimisation for the Fortran port" ON)
  if(WITH_FORTRAN_LTO)
    target_compile_options(fint PRIVATE -flto=auto -ffat-lto-objects)
    target_link_options(fint PUBLIC -flto=auto)
  endif()
endif()

# The Wheeler solver diagonalises a tridiagonal matrix with LAPACK's DSTEMR.
# The C ships a hand-translated copy of that routine (src/eigh.c, 1,477 lines)
# for builds with no LAPACK; calling the real one instead is what lets that
# file be dropped.
# LAPACK only when the caller asks for it.  By default the port uses its own
# translation of libcint's eigh.c and links nothing -- see cint_eigh.F90.
option(WITH_EXTERNAL_LAPACK "Use a system dstemr instead of the translated one" OFF)
if(WITH_EXTERNAL_LAPACK)
  target_compile_definitions(fint PRIVATE WITH_EXTERNAL_LAPACK)
  target_link_libraries(fint PUBLIC ${LAPACK_LIBRARIES})
endif()
