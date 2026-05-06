set(SRC_FILES
    apply_Htarget_vbatch.cpp
    cal_kron_flops.cpp
    dmrg_lapack.cpp
    dmrg_malloc.cpp
    dmrg_vbatch.cpp
    gen_patches_comb.cpp
    setup_matrix.cpp
    setup_vbatch.cpp
    unsetup_vbatch.cpp
    setup_nC.cpp
    setup_sparse_batch.cpp
    unsetup_sparse_batch.cpp
    apply_Htarget_sparse_not_ready.cpp
    estimate_work.cpp
    BatchedGemm.cpp)

add_library(gpuplugin ${SRC_FILES})
target_include_directories(gpuplugin PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(gpuplugin PUBLIC psimaglite::psimaglite)
target_compile_definitions(gpuplugin PUBLIC FpType=double)

add_executable(test_vbatch test_vbatch.cpp)
target_link_libraries(test_vbatch gpuplugin)
add_test(NAME vbatch COMMAND test_vbatch)

add_executable(test_gen_patches_comb test_gen_patches_comb.cpp)
target_link_libraries(test_gen_patches_comb gpuplugin)
add_test(NAME gen_patches_comb COMMAND test_gen_patches_comb)
