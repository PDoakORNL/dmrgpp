add_library(lanczos_globals LanczosGlobals.cpp)
target_include_directories(lanczos_globals PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Engine)
target_link_libraries(lanczos_globals PRIVATE psimaglite)

add_executable(
  lanczos
  lanczos.cpp
  LanczosDriver0.cpp
  LanczosDriver1.cpp
  LanczosDriver2.cpp
  LanczosDriver3.cpp
  LanczosDriver4.cpp
  LanczosDriver5.cpp)
target_include_directories(lanczos PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Engine
                                           ${CMAKE_CURRENT_SOURCE_DIR}/Model)
target_link_libraries(lanczos PRIVATE lanczos_globals)

target_link_libraries(lanczos PUBLIC psimaglite::psimaglite psimaglite)

add_executable(thermal thermal.cpp)
target_include_directories(thermal PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/Engine)
target_link_libraries(thermal psimaglite::psimaglite)

add_executable(lorentzian lorentzian.cpp)

target_link_libraries(lorentzian psimaglite::psimaglite)
