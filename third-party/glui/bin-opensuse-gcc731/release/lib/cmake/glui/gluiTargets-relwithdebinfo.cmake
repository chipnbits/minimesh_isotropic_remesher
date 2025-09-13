#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "glui::glui_static" for configuration "RelWithDebInfo"
set_property(TARGET glui::glui_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(glui::glui_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "CXX"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libglui_static.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui_static "${_IMPORT_PREFIX}/lib/libglui_static.a" )

# Import target "glui::glui" for configuration "RelWithDebInfo"
set_property(TARGET glui::glui APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(glui::glui PROPERTIES
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/libglui.so.2.37"
  IMPORTED_SONAME_RELWITHDEBINFO "libglui.so.2.37"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui "${_IMPORT_PREFIX}/lib/libglui.so.2.37" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
