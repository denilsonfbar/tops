FIND_PATH(RMATH_INCLUDE_DIR Rmath.h /usr/include/ /usr/local/include/ /Library/Frameworks/R.framework/Resources/include/)

FIND_LIBRARY(RMATH_LIBRARY Rmath PATHS /usr/lib /usr/local/lib /Library/Frameworks/R.framework/Resources/lib ) 

IF (RMATH_INCLUDE_DIR AND RMATH_LIBRARY)
   SET(RMATH_FOUND TRUE)
ENDIF (RMATH_INCLUDE_DIR AND RMATH_LIBRARY)


IF (RMATH_FOUND)
   IF (NOT RMATH_FIND_QUIETLY)
      MESSAGE(STATUS "Found Rmath: ${RMATH_LIBRARY}")
   ENDIF (NOT RMATH_FIND_QUIETLY)
ELSE (RMATH_FOUND)
   IF (RMATH_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find Rmath")
   ENDIF (RMATH_FIND_REQUIRED)
ENDIF (RMATH_FOUND)
