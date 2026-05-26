# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_executable(cincuenta cincuenta.cpp)
target_link_libraries(cincuenta PUBLIC dmrgpp_utils lanczos_globals dmrg_runner kronutil)

add_executable(freeGimp freeGimp.cpp)
target_link_libraries(freeGimp PUBLIC dmrgpp_utils kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# tests below
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

# General helper: check one captured regex group against an expected value.
# check_value uses the LAST match by default; pass --first as part of pattern
# arg (prepend to pattern string) to stop at the first match instead.
function(add_check name pattern expected tolerance file parent_test)
  add_test(NAME ${name}
    COMMAND $<TARGET_FILE:check_value> "${pattern}" "${expected}" "${tolerance}" "${file}")
  set_tests_properties(${name} PROPERTIES DEPENDS ${parent_test})
endfunction()

# Wrapper: check one captured value from gimp_lanczos.txt.
# gimp_lanczos.txt lists negative frequencies first, then positive.  check_value
# uses the LAST match by default, which for pattern "0\\.0628319\\s+..." picks up
# the positive-ω₀ line (Im[G] < 0, causal) rather than the negative-ω₀ line.
function(add_gimp_check name pattern expected tolerance parent_test)
  add_check(${name} "${pattern}" "${expected}" "${tolerance}"
            ${CMAKE_CURRENT_BINARY_DIR}/gimp_lanczos.txt ${parent_test})
endfunction()

# exact impurity solver
add_test(NAME cincuenta0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta0.ain) # energy 0.456977
# dmrg impurity solver
add_test(NAME cincuenta1 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputCincuenta1.ain) # energy 0.456977
# equilibrium Lanczos solver, NstatesEq=0 (full diag), 1 iteration
# G_imp(iw_0=0.0628) = (0.00551, -2.482) -- 1st-iteration values from zero sigma start
add_test(NAME eqLanczos0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputEqLanczosRef.ain)
# Re[G_imp(+iw_0)]: last match = positive-frequency line
add_gimp_check(eqLanczos0_Re
  "0\\.0628319\\s+([-0-9.eE+]+)\\s+[-0-9.eE+]+"  "0.00551" "1e-4" eqLanczos0)
# Im[G_imp(+iw_0)]: causal (negative), magnitude confirms both sectors contribute
add_gimp_check(eqLanczos0_Im
  "0\\.0628319\\s+[-0-9.eE+]+\\s+([-0-9.eE+]+)"  "-2.48174" "1e-4" eqLanczos0)

# Non-equilibrium DMFT: interaction quench U_i=0 -> U_f=2, 5-bath exact diag,
# beta=16, t_max=5 (N_t=200).  Reference: Tsuji et al. noneq-dmft (IPT, infinite bath).
# Values marked [Tsuji] are analytically exact for the Bethe lattice; others are
# regression anchors for the 5-bath finite-bath discretisation.
add_test(NAME neqTsuji0 COMMAND ./cincuenta -I ${PATH_TO_INPUTS} -f inputNeqTsuji.ain)

# G^R(t=0, t'=0): boundary condition G^R(t,t) = -i  [Tsuji: -1.0000, exact]
add_check(neqTsuji0_GR00_Im
  "0\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)"  "-1.0"  "1e-8"
  ${CMAKE_CURRENT_BINARY_DIR}/green-retarded  neqTsuji0)

# G^M(tau=0) = -n_imp: 5-bath approximation to half-filling  [Tsuji: -0.5000]
# ^ anchor required: "0.0000000000" appears as substring inside "10.0000000000".
add_check(neqTsuji0_GM0_Re
  "^0\\.0000000000 ([-0-9.eE+]+) [-0-9.eE+]+"  "-0.472861"  "1e-4"
  ${CMAKE_CURRENT_BINARY_DIR}/green-matsubara-t  neqTsuji0)

# G^<(t=0, t'=0): initial occupation Im[G^<(0,0)] = n(0)  [Tsuji: 0.5000]
add_check(neqTsuji0_GL00_Im
  "0\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)"  "0.527139"  "1e-4"
  ${CMAKE_CURRENT_BINARY_DIR}/green-lesser  neqTsuji0)

# G^R(t=1, t'=0): retarded GF one time unit after quench — regression anchor
add_check(neqTsuji0_GR10_Im
  "1\\.0000000000 0\\.0000000000 [^ ]+ ([-0-9.eE+]+)"  "-0.294763"  "1e-4"
  ${CMAKE_CURRENT_BINARY_DIR}/green-retarded  neqTsuji0)

# n(t=1) = Im[G^<(1,1)]: occupation after quench — regression anchor
add_check(neqTsuji0_n1
  "1\\.0000000000 1\\.0000000000 [^ ]+ ([-0-9.eE+]+)"  "0.506618"  "1e-4"
  ${CMAKE_CURRENT_BINARY_DIR}/green-lesser  neqTsuji0)

# n(t=5) = Im[G^<(5,5)]: final occupation — catches drift in time propagation
add_check(neqTsuji0_n5
  "5\\.0000000000 5\\.0000000000 [^ ]+ ([-0-9.eE+]+)"  "0.507376"  "1e-4"
  ${CMAKE_CURRENT_BINARY_DIR}/green-lesser  neqTsuji0)

# tDMRG non-equilibrium test: U_i=0 -> U_f=2 quench, compare G^R and G^< against
# exact-diag reference.  Requires Python 3; skipped silently if not available.
if(Python3_FOUND)
  set(PATH_TO_TESTSUITE ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite)

  # Each run writes green-retarded / green-lesser to its own subdirectory so
  # they don't overwrite each other or the neqTsuji0 outputs above.
  set(TDMRG_DIR ${CMAKE_CURRENT_BINARY_DIR}/neqTsujiTdmrg)
  set(ED_DIR    ${CMAKE_CURRENT_BINARY_DIR}/neqTsujiExactDiag)
  file(MAKE_DIRECTORY ${TDMRG_DIR} ${ED_DIR})

  # Run 1: particle + hole tDMRG, produces green-retarded and green-lesser
  add_test(NAME neqTsujiTdmrg0
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqTsujiTdmrg.ain
    WORKING_DIRECTORY ${TDMRG_DIR})

  # Run 2: exact-diag reference (single-step, no DMFT loop)
  add_test(NAME neqTsujiExactDiag0
    COMMAND $<TARGET_FILE:cincuenta> -I ${PATH_TO_INPUTS} -f inputNeqTsujiExactDiag.ain
    WORKING_DIRECTORY ${ED_DIR})

  # Run 3: compare tDMRG vs ED; tolerances set with ~2x headroom over current results
  # (G^R mean |dRe|~0.009, mean |dIm|~0.007 for U_i=0->U_f=2, m=200, dt=0.1, 6 sites)
  add_test(NAME neqTsujiTdmrg0_compare
    COMMAND ${Python3_EXECUTABLE}
      ${PATH_TO_TESTSUITE}/compare_neq_gr.py
      --tdmrg-retarded ${TDMRG_DIR}/green-retarded
      --tdmrg-lesser   ${TDMRG_DIR}/green-lesser
      --ed-retarded    ${ED_DIR}/green-retarded
      --ed-lesser      ${ED_DIR}/green-lesser
      --tolerance-re   0.02
      --tolerance-im   0.02)
  set_tests_properties(neqTsujiTdmrg0_compare PROPERTIES
    DEPENDS "neqTsujiTdmrg0;neqTsujiExactDiag0")
endif()
