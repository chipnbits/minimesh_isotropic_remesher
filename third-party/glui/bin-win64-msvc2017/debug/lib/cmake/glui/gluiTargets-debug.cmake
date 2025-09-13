#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "glui::glui_static" for configuration "Debug"
set_property(TARGET glui::glui_static APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(glui::glui_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/glui_staticd.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui_static "${_IMPORT_PREFIX}/lib/glui_staticd.lib" )

# Import target "glui::glui" for configuration "Debug"
set_property(TARGET glui::glui APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(glui::glui PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/gluid.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/gluid.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui "${_IMPORT_PREFIX}/lib/gluid.lib" "${_IMPORT_PREFIX}/bin/gluid.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
