################################################################################
#
# \file      charm.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 for
# Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration. All rights in the program are reserved by Triad National
# Security, LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting on its
# behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
# to reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so.
#
# Additionally, redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
#
# 3. Neither the name of Triad National Security, LLC, Los Alamos National
# Laboratory, LANL, the U.S. Government, nor the names of its contributors may be
# used to endorse or promote products derived from this software without specific
# prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# \brief     Function used to setup a Charm++ module.
#
# 2021/07/05 pgrete: adapted from https://github.com/quinoacomputing/cmake-modules
# 2023/02/05 mabruzzo: refactored the addCharmModule so that the primary is now an INTERFACE library rather than a custom target. The benefit of doing this is that the include_directory dependencies are now automatically managed
# 2023/03/19 mabruzzo: refactored addCharmModule so that it now puts the generated headers into custom subdirectories placed in the build-tree and declaring all headers in those subdirectories to be system headers (in order to suppress compiler warnings)
#
################################################################################

if(__addCharmModule)
  return()
endif()
set(__addCharmModule YES)

# Function 'addCharmModule' is used to add custom build commands and dependency
# for a Charm++ module
#
# This function creates an INTERFACE library target called ${MODULE}CharmModule
# that represents the header-files generated by charmc. You should treat
# ${MODULE}CharmModule just like any other header-only library.
#
# Just as with any other header-only library:
# - to denote that a given target depends on the .decl.h and/or .decl.h headers
#   produced by this command (& represented by ${MODULE}CharmModule), users
#   should specify ${MODULE}CharmModule as a dependency in a call to
#   target_link_libraries.
#   * This ensures that the necessary include-directories are automatically
#     used when building that target (or any dependents on that target)
#   * users should specify PUBLIC/PRIVATE/INTERFACE in the call to
#     target_link_libraries based on the visibility of the .decl.h/.def.h files
#     in the target's public header files (if applicable). This handles
#     transitive dependencies correctly.
# - users should specify dependencies of this ${MODULE}CharmModule on another
#   header-only libraries, ${OTHER_HEADER_ONLY}, to properly propagate
#   dependencies. This is accomplished with a call to:
#     target_link_libraries(${MODULE}CharmModule
#       INTERFACE ${OTHER_HEADER_ONLY}
#     )
#   Users should also do this for all other externally defined modules that are
#   referenced within ${MODULE}.ci
#
# Doing all of this ensures any targets that depends (directly or transitively)
# on the generated .decl.h/.def.h files will carry a transitive dependence to
# the underlying .ci file (i.e. a change to that file triggers an appropriate
# rebuild)
#
# Generation Details
# ==================
# All header files generated from *.ci files in a given source directory are
# are placed within the corresponding build-directory inside of the
# _generated_charm_module_headers subdirectory. That subdirectory is created at
# configuration time while the headers are generated later at build time.
#
# This makes it easier to change how charm-generated headers are treated (e.g.
# declaring that this newly created subdirectory includes SYSTEM headers to
# suppress compiler warnings) without affecting the treatment of other
# generated header files
#
# Suppressing Compiler Warnings
# =============================
# Headers generated by certain versions of Charm++ have been found to trigger
# compiler warnings (e.g. gcc's -Wsign-compare and -Wunused-variable). Given
# that applications have no control over these warnings, it's desirable for
# these warnings to be suppressed. This is done by declaring the header files
# generated by Charm++ to be system header files.
#
# If you would prefer to not declare them as system headers, define the
# ``SUPPRESS_CHARM_HEADER_WARNING`` global variable with a value of ``0``.

function(addCharmModule MODULE)
  # Arguments:
  #   MODULE:    Name of the Charm++ module.
  # Add custom command generating .decl.h and .def.h from .ci
  # IMPORTANT: CHARM_PREPROC_DEFS need to include all preprocessor defines that
  #   are used in the *.ci files.

  # determine the path to subdirectory where we should place the header files
  # generated by charm++.
  set(OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/_generated_charm_module_headers")

  # create OUT_DIR (if it doesn't already exist)
  if (NOT IS_ABSOLUTE "${OUT_DIR}")
    # The behavior of IS_DIRECTORY isn't well-defined for relative paths
    MESSAGE(FATAL_ERROR "sanity check failed - something went horribly wrong")
  elseif (NOT IS_DIRECTORY "${OUT_DIR}")
    file(MAKE_DIRECTORY "${OUT_DIR}")
  endif()

  add_custom_command(
    OUTPUT ${OUT_DIR}/${MODULE}.decl.h ${OUT_DIR}/${MODULE}.def.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.ci
    COMMAND ${CHARM_COMPILER} ${CHARM_PREPROC_DEFS} ${ARGN}
            ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.ci
    WORKING_DIRECTORY ${OUT_DIR})

  # Add custom target dependency for Charm++ module
  add_custom_target(${MODULE}CharmModule_custom_target_
    DEPENDS ${OUT_DIR}/${MODULE}.decl.h
            ${OUT_DIR}/${MODULE}.def.h)

  add_library(${MODULE}CharmModule INTERFACE)
  add_dependencies(${MODULE}CharmModule ${MODULE}CharmModule_custom_target_)

  # The following command informs all targets that will be linked against
  # ${MODULE}CharmModule where to find the generated headers
  if(SUPPRESS_CHARM_HEADER_WARNING OR NOT DEFINED SUPPRESS_CHARM_HEADER_WARNING)
    target_include_directories(${MODULE}CharmModule SYSTEM INTERFACE ${OUT_DIR})
  else()
    target_include_directories(${MODULE}CharmModule INTERFACE ${OUT_DIR})
  endif()

endfunction(addCharmModule)
