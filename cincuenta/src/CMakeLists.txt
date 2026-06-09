# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_executable(cincuenta cincuenta.cpp)
target_link_libraries(cincuenta PUBLIC dmrgpp_utils lanczos_globals dmrg_runner kronutil)

add_executable(freeGimp freeGimp.cpp)
target_link_libraries(freeGimp PUBLIC dmrgpp_utils kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# tests below
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

# exact impurity solver
add_test(NAME cincuenta0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta0.ain) # energy 0.456977
# dmrg impurity solver
add_test(NAME cincuenta1 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta1.ain) # energy 0.456977

# Regression: G_imp must equal G_cluster0 at U=0 with FitOptions=particleholesymmetric.
# Physical basis: U=0 => Sigma_imp=0 exactly, so the impurity solver must reproduce
# the non-interacting cluster propagator G_cluster0 to solver precision.
# This was broken by a bath-parameter layout bug fixed in PR #195 (Anderson_symPH).
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_FOUND)
  if(NOT PATH_TO_TESTSUITE)
    set(PATH_TO_TESTSUITE ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite)
  endif()

  set(U0_ED_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0ExactDiag)
  file(MAKE_DIRECTORY ${U0_ED_DIR})

  add_test(
    NAME u0_gImpEqGcluster0_exactdiag
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cincuenta -I ${PATH_TO_INPUTS} -f inputU0ExactDiag.ain
    WORKING_DIRECTORY ${U0_ED_DIR})

  add_test(NAME u0_gImpEqGcluster0_exactdiag_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_gcluster0.py ${U0_ED_DIR}/gimp_exactdiag.txt
                   ${U0_ED_DIR}/gcluster0_exactdiag.txt --tolerance 1e-3)
  set_tests_properties(u0_gImpEqGcluster0_exactdiag_check PROPERTIES DEPENDS u0_gImpEqGcluster0_exactdiag)

  # Same check without FitOptions=particleholesymmetric (unconstrained 2*nBath fit).
  # Verifies that half-filling and G_imp=G_cluster0 hold without the PH constraint.
  set(U0_ED_NOSYM_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0ExactDiagNoSym)
  file(MAKE_DIRECTORY ${U0_ED_NOSYM_DIR})

  add_test(
    NAME u0_gImpEqGcluster0_exactdiag_nosym
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cincuenta -I ${PATH_TO_INPUTS} -f inputU0ExactDiagNoSym.ain
    WORKING_DIRECTORY ${U0_ED_NOSYM_DIR})

  add_test(NAME u0_gImpEqGcluster0_exactdiag_nosym_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_gcluster0.py
                   ${U0_ED_NOSYM_DIR}/gimp_exactdiag.txt ${U0_ED_NOSYM_DIR}/gcluster0_exactdiag.txt --tolerance 1e-3)
  set_tests_properties(u0_gImpEqGcluster0_exactdiag_nosym_check PROPERTIES DEPENDS u0_gImpEqGcluster0_exactdiag_nosym)

  # DMRG U=0 smoke test: nBath=4 (even, no forced zero-energy site), no PH constraint.
  # Even nBath avoids the CorrectionVector ill-conditioning at low omega_n that
  # plagues nBath=5 PH-symmetric at beta=20 (forced eps=0 site, error ~8e-2).
  # nUp=nDown=2 is the nearest Sz=0 half-filling for 5 total sites.
  set(U0_DMRG_4BATH_DIR ${CMAKE_CURRENT_BINARY_DIR}/u0Dmrg4Bath)
  file(MAKE_DIRECTORY ${U0_DMRG_4BATH_DIR})

  add_test(
    NAME u0_gImpEqGcluster0_dmrg_4bath
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cincuenta -I ${PATH_TO_INPUTS} -f inputU0Dmrg4Bath.ain
    WORKING_DIRECTORY ${U0_DMRG_4BATH_DIR})

  add_test(NAME u0_gImpEqGcluster0_dmrg_4bath_check
           COMMAND ${Python3_EXECUTABLE} ${PATH_TO_TESTSUITE}/compare_gimp_gcluster0.py
                   ${U0_DMRG_4BATH_DIR}/gimp_dmrg.txt ${U0_DMRG_4BATH_DIR}/gcluster0_dmrg.txt --tolerance 1e-3)
  set_tests_properties(u0_gImpEqGcluster0_dmrg_4bath_check PROPERTIES DEPENDS u0_gImpEqGcluster0_dmrg_4bath)
endif()
