set(TARGET_NAME minimesh_tests)

add_executable(${TARGET_NAME}
  test_mesh_properties.cpp
)

# Link your core lib
target_link_libraries(${TARGET_NAME} PRIVATE minimeshcore)

# Include path for vendored doctest and Eigen
target_include_directories(${TARGET_NAME} PRIVATE
  ${CMAKE_SOURCE_DIR}/third-party/doctest
  ${CMAKE_SOURCE_DIR}/src
)

# Add Eigen include directories (using same approach as main CMakeLists)
target_include_directories(${TARGET_NAME} SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIRS})

# Set C++ standard to match the project
set_target_properties(${TARGET_NAME} PROPERTIES
  CXX_STANDARD 14
  CXX_STANDARD_REQUIRED ON
)

# Register a single CTest (shows up in VS Code too)
# Set working dir so 'mesh/cube.obj' paths resolve.
add_test(NAME minimesh.mesh COMMAND ${TARGET_NAME})
set_tests_properties(minimesh.mesh PROPERTIES
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
