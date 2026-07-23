# Custom main (not Catch2::Catch2WithMain) because csr_kron_mult's imethod3
# runs a Kokkos kernel and needs Kokkos::ScopeGuard constructed before any
# test body executes.
add_executable(test_csr_kron_mult_imethod3 test_csr_kron_mult_imethod3.cpp)
target_link_libraries(test_csr_kron_mult_imethod3 PRIVATE Catch2::Catch2 kronutil)
catch_discover_tests(test_csr_kron_mult_imethod3)
