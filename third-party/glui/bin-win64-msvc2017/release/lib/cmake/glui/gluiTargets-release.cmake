#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "glui::glui_static" for configuration "Release"
set_property(TARGET glui::glui_static APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(glui::glui_static PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/glui_static.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui_static )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui_static "${_IMPORT_PREFIX}/lib/glui_static.lib" )

# Import target "glui::glui" for configuration "Release"
set_property(TARGET glui::glui APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(glui::glui PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/glui.lib"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/glui.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS glui::glui )
list(APPEND _IMPORT_CHECK_FILES_FOR_glui::glui "${_IMPORT_PREFIX}/lib/glui.lib" "${_IMPORT_PREFIX}/bin/glui.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
