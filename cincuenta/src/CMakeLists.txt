# This is repeated here but should possibly just use the custom
# command and info from dmrgpp/src

add_library(provenance OBJECT Provenance.cpp)
target_link_libraries(provenance PUBLIC dmrgpp psimaglite)

set(cincuenta_src cincuenta.cpp)
add_executable(cincuenta ${cincuenta_src} $<TARGET_OBJECTS:provenance>)
target_link_libraries(cincuenta PUBLIC dmrgpp kronutil)

add_executable(freeGimp freeGimp.cpp $<TARGET_OBJECTS:provenance>)
target_link_libraries(freeGimp PUBLIC dmrgpp kronutil)

install(TARGETS cincuenta freeGimp RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# tests below
set(PATH_TO_INPUTS ${CMAKE_CURRENT_SOURCE_DIR}/../TestSuite/inputs)

add_test(NAME cincuenta0 COMMAND ./cincuenta -f ${PATH_TO_INPUTS}/inputCincuenta0.ain) # energy 0.456977
