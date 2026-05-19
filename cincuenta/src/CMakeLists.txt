# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_executable(cincuenta cincuenta.cpp)
target_link_libraries(cincuenta PUBLIC dmrgpp_utils lanczos_globals dmrg_runner kronutil)

add_executable(freeGimp freeGimp.cpp)
target_link_libraries(freeGimp PUBLIC dmrgpp_utils kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# tests below
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

# Helper: check one captured value from gimp_lanczos.txt after a parent test runs.
# gimp_lanczos.txt lists negative frequencies first, then positive.  check_value
# uses the LAST match by default, which for pattern "0\\.0628319\\s+..." picks up
# the positive-ω₀ line (Im[G] < 0, causal) rather than the negative-ω₀ line.
# Pass --first to explicitly select the first match instead.
function(add_gimp_check name pattern expected tolerance parent_test)
  add_test(NAME ${name}
    COMMAND $<TARGET_FILE:check_value> "${pattern}" "${expected}" "${tolerance}"
            ${CMAKE_CURRENT_BINARY_DIR}/gimp_lanczos.txt)
  set_tests_properties(${name} PROPERTIES DEPENDS ${parent_test})
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
