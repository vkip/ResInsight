set(SOURCE_GROUP_SOURCE_FILES
    ${CMAKE_CURRENT_LIST_DIR}/RicEclipseWellFeatureImpl.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RicEclipseWellShowFeatures.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RicCreateWellPathFromSimulationWellFeature.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RicCreateApproximateWellPathFeature.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RicApproximateWellPathTools.cpp
)

list(APPEND COMMAND_CODE_SOURCE_FILES ${SOURCE_GROUP_SOURCE_FILES})
