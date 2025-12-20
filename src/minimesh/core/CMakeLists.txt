set(TARGET_NAME minimeshcore)

##
## SOURCE FILES
##
set(SOURCE_FILES
  ##
  ## UTIL
  ##
  util/assert.hpp
  util/console_colors.hpp
  util/numbers.hpp
  util/foldertools.hpp
  util/foldertools.cpp
  util/macros.hpp
  util/disable_irrelevant_warnings.hpp
  util/enable_irrelevant_warnings.hpp
  util/bit_vector.hpp
  util/get_basefilename.hpp
  util/get_basefilename.cpp
  ##
  ## Manifold oriented half edge
  ##
  mohe/mesh_connectivity.hpp
  mohe/mesh_connectivity.cpp
  mohe/mesh_io.hpp
  mohe/mesh_io.cpp
  mohe/mesh_modifier_loop_subdivision.hpp
  mohe/mesh_modifier_loop_subdivision.cpp
  mohe/mesh_modifier_edge_collapse.hpp
  mohe/mesh_modifier_edge_collapse.cpp
  mohe/mesh_modifier_arap.hpp
  mohe/mesh_modifier_arap.cpp
  mohe/mesh_analysis.hpp
  mohe/mesh_analysis.cpp
  mohe/fixed_uv_param.hpp
  mohe/fixed_uv_param.cpp
  mohe/lscm_uv_param.hpp
  mohe/lscm_uv_param.cpp
 )

 
# create a library that the viewer can use to link to.
add_library(${TARGET_NAME} ${SOURCE_FILES})

##
## Create a folder structure inside the visual studio project
##
minimesh_folders_lib(${TARGET_NAME}  ${SOURCE_FILES})
