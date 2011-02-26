set(PROJECT_DESCRIPTION  "Blender is a very fast and versatile 3D modeller/renderer.")
set(PROJECT_COPYRIGHT    "Copyright (C) 2001-2011 Blender Foundation")
set(PROJECT_CONTACT      "foundation@blender.org")
set(PROJECT_VENDOR       "Blender Foundation")
set(ORG_WEBSITE          "www.blender.org")

set(MAJOR_VERSION ${BLENDER_VERSION_MAJOR})
set(MINOR_VERSION ${BLENDER_VERSION_MINOR})
set(PATCH_VERSION ${BLENDER_VERSION_CHAR_INDEX})

set(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME})
set(CPACK_PACKAGE_DESCRIPTION ${PROJECT_DESCRIPTION})
set(CPACK_PACKAGE_VENDOR ${PROJECT_VENDOR})
set(CPACK_PACKAGE_CONTACT ${PROJECT_CONTACT})
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
SET(CPACK_PACKAGE_VERSION_MAJOR "${MAJOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_MINOR "${MINOR_VERSION}")
SET(CPACK_PACKAGE_VERSION_PATCH "${PATCH_VERSION}")


# Get the build revision, note that this can get out-of-sync, so for packaging run cmake first.
include(FindSubversion)
set(MY_WC_REVISION "unknown")
if(EXISTS ${CMAKE_SOURCE_DIR}/.svn/)
	if(Subversion_FOUND)
		Subversion_WC_INFO(${CMAKE_SOURCE_DIR} MY)
	endif()
endif()
set(BUILD_REV ${MY_WC_REVISION})


# Force Package Name
set(CPACK_PACKAGE_FILE_NAME ${PROJECT_NAME}-${BLENDER_VERSION}-r${BUILD_REV}-${CPACK_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})

# RPM packages
include(build_files/cmake/RpmBuild.cmake)
if(RPMBUILD_FOUND AND NOT WIN32)
	set(CPACK_GENERATOR "RPM")
	set(CPACK_SET_DESTDIR TRUE)
	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
	set(CPACK_RPM_PACKAGE_LICENSE "GPLv2")
	set(CPACK_RPM_PACKAGE_GROUP "Amusements/Graphics")
endif()

# Mac Bundle
if(APPLE)
	set(CPACK_GENERATOR "DragNDrop")

	# Libraries are bundled directly
	set(CPACK_COMPONENT_LIBRARIES_HIDDEN TRUE)

	# Bundle Properties
	set(MACOSX_BUNDLE_BUNDLE_NAME blender)
	set(MACOSX_BUNDLE_BUNDLE_VERSION ${BLENDER_VERSION})
	set(MACOSX_BUNDLE_SHORT_VERSION_STRING ${BLENDER_VERSION})
	set(MACOSX_BUNDLE_LONG_VERSION_STRING "Version ${BLENDER_VERSION}-r${BUILD_REV}")
endif(APPLE)

set(CPACK_PACKAGE_EXECUTABLES "blender")
include(CPack)
