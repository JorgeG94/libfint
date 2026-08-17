cmake_minimum_required(VERSION 3.15)
project(libfint VERSION 6.1.3 LANGUAGES Fortran)

include(GNUInstallDirs)

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
  add_subdirectory(extern/libcint EXCLUDE_FROM_ALL)
  enable_testing()
  add_subdirectory(test)
endif()
