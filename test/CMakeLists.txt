find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Everything here compares libfint against libcint, so every target links
# both.  The C comes from extern/libcint, pinned to the commit libfint is
# bit-identical to -- which is what makes that claim checkable rather than
# asserted.

# The public interface is part of the library now, so the tests just link it.
add_library(fint_test_systems OBJECT cint_test_systems.F90)
target_link_libraries(fint_test_systems PUBLIC fint)
target_include_directories(fint_test_systems PRIVATE ${PROJECT_BINARY_DIR}/include)
target_include_directories(fint_test_systems PRIVATE ${PROJECT_BINARY_DIR}/include)
set_target_properties(fint_test_systems PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# catalogue_check.F90 is REGENERATED here, against whichever libcint was
# fetched, rather than used as committed.  That is not a convenience: the
# generator filters on the symbols libcint.so actually exports, so it drops
# any entry point the reference does not have -- libfint carries integrals
# upstream has not merged (int2c2e_ipip2), and a committed catalogue would
# simply fail to link against upstream master.  Regenerating means the
# comparison covers exactly the intersection, and the coverage check is what
# tells you about the difference.
add_custom_command(
  OUTPUT ${PROJECT_BINARY_DIR}/catalogue_check.F90
  COMMAND ${CMAKE_COMMAND} -E env
          ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/test/gen_catalogue_check.py
          ${PROJECT_BINARY_DIR}/extern/libcint
          ${PROJECT_BINARY_DIR}/catalogue_check.F90
  DEPENDS cint ${PROJECT_SOURCE_DIR}/test/gen_catalogue_check.py
  COMMENT "Generating catalogue_check against the fetched libcint")

# catalogue_check.F90 is generated -- see gen_catalogue_check.py, which reads
# each integral's ng array out of the generated modules rather than trusting a
# hand-written interface block.
set(LIBFINT_CHECKS
  backend_check spinor_grids_check opt_check
  int1e_check int2e_check int3c1e_check int3c2e_check
  gen1e_check gen2e_check gen3c_check grids_check rel2e_check
  spinor1e_check spinor2e_check spinor3c_check)

# The ECP checks, kept apart from LIBFINT_CHECKS above.
#
# Everything in that list links `cint` and compares against libcint at run
# time.  These cannot and need not: libcint has no ECP code, so their
# references are committed (generated from PySCF by gen_ecp_check.py and
# gen_ecp_num_check.py) and they link `fint` alone.  Keeping them separate is
# what would let them run in the toolchain lanes that never fetch libcint.
#
# ecp_abi_check is C, and drives ECPscalar_sph/ECPscalar_cart through
# libcint's own C names.  It exists because ecp_check calls the Fortran module
# directly and so cannot reach the bind(C) layer -- which is exactly where the
# two faults fixed in "Fix the radial grid size, and two latent faults in the
# C ABI" lived, invisible to every test at the time.
foreach(chk ecp_num_check ecp_check)
  add_executable(${chk} ${chk}.F90)
  target_link_libraries(${chk} fint m)
  target_include_directories(${chk} PRIVATE ${PROJECT_BINARY_DIR}/include)
  add_test(NAME ${chk} COMMAND $<TARGET_FILE:${chk}>)
endforeach()

add_executable(ecp_abi_check ecp_abi_check.c)
target_link_libraries(ecp_abi_check fint m)
add_test(NAME ecp_abi_check COMMAND $<TARGET_FILE:ecp_abi_check>)

# Double-double carries no libcint reference: it is checked against real128 in
# the same binary, which is exact and needs nothing external. Kept out of
# LIBFINT_CHECKS for that reason -- everything in that list links `cint`.
add_executable(dd_check dd_check.f90)
target_link_libraries(dd_check fint)
add_test(NAME dd_check COMMAND $<TARGET_FILE:dd_check>)

add_executable(dd_ops_check dd_ops_check.f90)
target_link_libraries(dd_ops_check fint)
add_test(NAME dd_ops_check COMMAND $<TARGET_FILE:dd_ops_check>)

add_executable(catalogue_check ${PROJECT_BINARY_DIR}/catalogue_check.F90)
target_link_libraries(catalogue_check fint fint_test_systems cint m)
target_include_directories(catalogue_check PRIVATE ${PROJECT_BINARY_DIR}/include)
add_test(NAME catalogue_check COMMAND $<TARGET_FILE:catalogue_check>)

foreach(chk ${LIBFINT_CHECKS})
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${chk}.F90)
    add_executable(${chk} ${chk}.F90)
    target_link_libraries(${chk} fint fint_test_systems cint m)
    target_include_directories(${chk} PRIVATE ${PROJECT_BINARY_DIR}/include)
    add_test(NAME ${chk} COMMAND $<TARGET_FILE:${chk}>)
  endif()
endforeach()

# The layer checks are two-step: a small C program dumps libcint's answers to
# a binary, then the Fortran reads it back and compares.  A CTest fixture
# sequences the pair, so `ctest -R rys_check` still runs its dump first.
# libcint's private headers, for the C reference programs: cint.h and
# cint_config.h are generated into the submodule's build tree.
set(LIBCINT_PRIVATE_INCLUDES
  ${PROJECT_BINARY_DIR}/extern/libcint/include
  ${PROJECT_BINARY_DIR}/extern/libcint/src
  ${PROJECT_SOURCE_DIR}/extern/libcint/src)

function(fint_layer_check name ref_src chk_src fixture)
  add_executable(${name}_reference ${ref_src})
  target_link_libraries(${name}_reference cint m)
  target_include_directories(${name}_reference PRIVATE ${LIBCINT_PRIVATE_INCLUDES})
  add_executable(${name} ${chk_src})
  target_link_libraries(${name} fint m)
  target_include_directories(${name} PRIVATE ${PROJECT_BINARY_DIR}/include)
  add_test(NAME ${name}_dump COMMAND $<TARGET_FILE:${name}_reference> ${ARGN})
  set_tests_properties(${name}_dump PROPERTIES FIXTURES_SETUP ${fixture})
  add_test(NAME ${name} COMMAND $<TARGET_FILE:${name}> ${ARGN})
  set_tests_properties(${name} PROPERTIES FIXTURES_REQUIRED ${fixture})
endfunction()

fint_layer_check(rys_core_check rys_core_reference.c rys_core_check.f90 fint_rys_core
  ${PROJECT_BINARY_DIR}/rys_roots_c.bin ${PROJECT_BINARY_DIR}/rys_schmidt_c.bin)
fint_layer_check(support_check support_reference.c support_check.f90 fint_support
  ${PROJECT_BINARY_DIR}/c2s_c.bin ${PROJECT_BINARY_DIR}/blas_c.bin)
fint_layer_check(fmt_check fmt_reference.c fmt_check.F90 fint_fmt
  ${PROJECT_BINARY_DIR}/fmt_c.bin)

# rys_check takes a root count the dump does not, so it is not the shape above.
add_executable(rys_reference rys_reference.c)
target_link_libraries(rys_reference cint m)
target_include_directories(rys_reference PRIVATE ${LIBCINT_PRIVATE_INCLUDES})
add_executable(rys_check rys_check.f90)
target_link_libraries(rys_check fint m)
add_test(NAME rys_dump
  COMMAND $<TARGET_FILE:rys_reference> ${PROJECT_BINARY_DIR}/rys_golden.bin 14)
set_tests_properties(rys_dump PROPERTIES FIXTURES_SETUP fint_rys)
add_test(NAME rys_check
  COMMAND $<TARGET_FILE:rys_check> ${PROJECT_BINARY_DIR}/rys_golden.bin)
set_tests_properties(rys_check PROPERTIES FIXTURES_REQUIRED fint_rys)

# L shells, which libcint cannot express and so cannot be compared against.
# The reference is libfint's own split-shell path, which int2e_check holds to
# the last bit against the C -- so this is a comparison with libcint by one
# more step than the others.  Links `fint` alone for that reason.
add_executable(l_shell_check l_shell_check.f90)
target_link_libraries(l_shell_check fint m)
target_include_directories(l_shell_check PRIVATE ${PROJECT_BINARY_DIR}/include)
add_test(NAME l_shell_check COMMAND $<TARGET_FILE:l_shell_check>)

# The same contract for the gradient kernel int2e_ip1, which exercises the
# one thing int2e cannot: contraction buffers with n_comp > 1, where the
# stride of the s/p coefficient pattern carries a factor the plain integral
# never tests.
add_executable(l_shell_grad_check l_shell_grad_check.f90)
target_link_libraries(l_shell_grad_check fint m)
target_include_directories(l_shell_grad_check PRIVATE ${PROJECT_BINARY_DIR}/include)
add_test(NAME l_shell_grad_check COMMAND $<TARGET_FILE:l_shell_grad_check>)

# .f90 rather than .F90, but still compared against the C.
foreach(chk spinor_check)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${chk}.f90)
    add_executable(${chk} ${chk}.f90)
    target_link_libraries(${chk} fint cint m)
    add_test(NAME ${chk} COMMAND $<TARGET_FILE:${chk}>)
  endif()
endforeach()

if(WITH_F12)
  add_executable(f12_check f12_check.F90)
  target_link_libraries(f12_check fint fint_test_systems cint m)
  add_test(NAME f12_check COMMAND $<TARGET_FILE:f12_check>)
endif()

# Source-level, no build needed: does the Fortran generator still know every
# integral the C generator knows?  The C list comes from the submodule, so
# this checks against the pinned upstream rather than a copy.
if(Python3_Interpreter_FOUND)
  add_test(NAME integral_coverage
           COMMAND ${Python3_EXECUTABLE}
                   ${PROJECT_SOURCE_DIR}/scripts/check_f90_integral_coverage.py)
endif()

# The C ABI, driven from C.  Deliberately links `fint` and NOT `cint`: the
# point is that libfint alone resolves libcint's own entry-point names, so a
# missing or misnamed bind(C) fails here at link time instead of at a
# consumer's.  That also means it cannot compare against libcint at runtime --
# both define these symbols -- so its reference values are recorded in the
# source, taken from this same program built against upstream.
add_executable(c_abi_check c_abi_check.c)
target_link_libraries(c_abi_check fint m)
add_test(NAME c_abi_check COMMAND $<TARGET_FILE:c_abi_check>)
