add_executable(test_KadanoffBaym test_KadanoffBaym.cpp)
target_link_libraries(test_KadanoffBaym PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_KadanoffBaym PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src)
catch_discover_tests(test_KadanoffBaym)

add_executable(test_spectralWeight test_spectralWeight.cpp)
target_link_libraries(test_spectralWeight PRIVATE Catch2::Catch2WithMain psimaglite)
target_include_directories(test_spectralWeight PRIVATE ${CMAKE_SOURCE_DIR}/cincuenta/src
                                                       ${CMAKE_SOURCE_DIR}/dmrg/Engine)
catch_discover_tests(test_spectralWeight)
