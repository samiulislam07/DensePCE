cmake_minimum_required(VERSION 3.16)
project(dense_pce_ablations CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Sources
set(BASE_SRC dense-pce.cpp)
set(MOD_SRC  dense-pce-mod-edge-order.cpp)

# Helper to make a variant target
function(make_variant name src)
  add_executable(${name} ${src})
  target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/EBBkC/src)
endfunction()

# V0 baseline
make_variant(dense-pce-v0 ${BASE_SRC})

# A) EBBkC-only (seeds on, DS off, pruners off)
make_variant(ebbkc_only ${MOD_SRC})
target_compile_definitions(ebbkc_only PRIVATE
  USE_EBBKC_SEEDS=1
  USE_NEW_DATASTRUCTS=0
  APPLY_ORDER_BOUND=0
  APPLY_EDGE_BOUND=0
)
# Try to link EBBkC static lib if available
set(EBBKC_LIB_A ${CMAKE_SOURCE_DIR}/build_integrated/libebbkc_core.a)
set(EBBKC_LIB_A_ALT ${CMAKE_SOURCE_DIR}/build_integrated/libebbkc_core(1).a)
if(EXISTS ${EBBKC_LIB_A})
  target_link_libraries(ebbkc_only PRIVATE ${EBBKC_LIB_A})
elseif(EXISTS ${EBBKC_LIB_A_ALT})
  target_link_libraries(ebbkc_only PRIVATE ${EBBKC_LIB_A_ALT})
endif()

# B) Data-structure-only (seeds off, DS on, pruners off)
make_variant(ds_only ${MOD_SRC})
target_compile_definitions(ds_only PRIVATE
  USE_EBBKC_SEEDS=0
  USE_NEW_DATASTRUCTS=1
  APPLY_ORDER_BOUND=0
  APPLY_EDGE_BOUND=0
)

# (Optional) Full modified (everything on)
make_variant(full_mod ${MOD_SRC})
target_compile_definitions(full_mod PRIVATE
  USE_EBBKC_SEEDS=1
  USE_NEW_DATASTRUCTS=1
  APPLY_ORDER_BOUND=1
  APPLY_EDGE_BOUND=1
)
if(EXISTS ${EBBKC_LIB_A})
  target_link_libraries(full_mod PRIVATE ${EBBKC_LIB_A})
elseif(EXISTS ${EBBKC_LIB_A_ALT})
  target_link_libraries(full_mod PRIVATE ${EBBKC_LIB_A_ALT})
endif()


