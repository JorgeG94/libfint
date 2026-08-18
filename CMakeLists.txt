cmake_minimum_required(VERSION 3.15)
project(libfint VERSION 6.1.3 LANGUAGES Fortran)

include(GNUInstallDirs)

# A default build type, because CMake has none and every flag this project
# tunes lives in a per-configuration variable.
#
# Without this the command the README gives -- `cmake -S . -B build` -- leaves
# CMAKE_BUILD_TYPE empty, and then no -O reaches the compiler at all.  Not -O0
# by request: -O0 by omission.  On the 6-31G quartet sweep that is 11.71G
# instructions against 4.97G, so the documented way to build the library was
# 2.4x slower than the one every CI job uses.
#
# It survived for exactly that reason.  Every workflow under .github/ passes
# -DCMAKE_BUILD_TYPE=Release explicitly, so no automated build has ever taken
# the path a reader of the README takes.
#
# Guarded on CMAKE_CONFIGURATION_TYPES as well: a multi-config generator picks
# its configuration at build time, where setting this would mean nothing.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Build type; Release unless the caller asks for another" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY
               STRINGS Debug Release RelWithDebInfo MinSizeRel)
  message(STATUS "CMAKE_BUILD_TYPE was not set; defaulting to Release")
endif()

option(LIBFINT_BUILD_TESTS "Build the verification suite (needs libcint)" OFF)
option(WITH_F12 "Slater-type geminal integrals (generates a 100 MB table)" OFF)
option(WITH_EXTERNAL_LAPACK "Use a system dstemr instead of the translated one" OFF)

# No LAPACK, and no BLAS.  libfint carries its own translation of libcint's
# eigh.c, exactly as libcint carries the C of it -- see src/cint_eigh.F90.
# WITH_EXTERNAL_LAPACK swaps in a system dstemr for callers that would rather
# have a vendor implementation; it is faster and not bit-identical.
if(WITH_EXTERNAL_LAPACK)
  find_package(LAPACK REQUIRED)
endif()

add_subdirectory(src)

if(LIBFINT_BUILD_TESTS)
  # libfint contains no C.  The suite compares against libcint at runtime, so
  # it is fetched rather than vendored -- and against upstream master rather
  # than a pin, so that a red run means upstream moved.  See
  # scripts/fetch_libcint.sh.
  if(NOT EXISTS ${PROJECT_SOURCE_DIR}/extern/libcint/CMakeLists.txt)
    message(FATAL_ERROR
      "LIBFINT_BUILD_TESTS needs libcint to compare against.
"
      "  scripts/fetch_libcint.sh          # upstream master
"
      "  scripts/fetch_libcint.sh v6.1.3   # or a specific ref")
  endif()
  enable_language(C)

  # Intel's C compiler defaults to a relaxed floating-point model, and the
  # Fortran side is already built -fp-model=precise (see src/CMakeLists.txt).
  # Comparing the two then compares two different arithmetics: on the Rys
  # Schmidt path that is 6,686 of 7,192 values "differing" with nothing wrong
  # on either side.  So the reference gets the same model the port has.
  #
  # This is set before add_subdirectory so libcint inherits it, and only for
  # the fetched reference -- it does not follow libfint anywhere.
  if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fp-model=precise")
  endif()

  add_subdirectory(extern/libcint EXCLUDE_FROM_ALL)
  enable_testing()
  add_subdirectory(test)
endif()
