# CMake script to regression-test --extract-scans: the extracted mzML must
# contain exactly the requested spectra — no more, no fewer.
#
# Regression: the parsed scan numbers used to be stored through a uint32_t*
# view of the long* array returned by string_to_array(), so on LP64 every
# other scan was replaced by a zero ("3,9" extracted only scan 3, "2-25"
# extracted scans 2-13).
#
# Expected variables:
# EXECUTABLE: Path to mscompress executable
# INPUT_FILE: mzML fixture to compress into an msz first
# TEMP_DIR: Path to temporary directory
# SCAN_SPEC: Argument passed to --extract-scans (e.g. "3", "3,9", "2-5", "[1-3,5-6]")
# EXPECTED_SCANS: Semicolon-separated list of scan numbers that must be extracted

if(NOT EXISTS "${EXECUTABLE}")
    message(FATAL_ERROR "Executable not found: ${EXECUTABLE}")
endif()

if(NOT EXISTS "${INPUT_FILE}")
    message(FATAL_ERROR "Input file not found: ${INPUT_FILE}")
endif()

if(NOT DEFINED SCAN_SPEC)
    message(FATAL_ERROR "SCAN_SPEC not defined")
endif()

if(NOT DEFINED EXPECTED_SCANS)
    message(FATAL_ERROR "EXPECTED_SCANS not defined")
endif()

if(NOT EXISTS "${TEMP_DIR}")
    file(MAKE_DIRECTORY "${TEMP_DIR}")
endif()

set(MSZ_FILE "${TEMP_DIR}/input.msz")
set(OUTPUT_FILE "${TEMP_DIR}/output.mzML")

# Compress the fixture into an msz.
execute_process(
    COMMAND ${EXECUTABLE} "${INPUT_FILE}" "${MSZ_FILE}"
    RESULT_VARIABLE COMPRESS_RESULT
    OUTPUT_QUIET
    ERROR_QUIET
)

if(NOT COMPRESS_RESULT EQUAL 0)
    message(FATAL_ERROR "Compression failed with return code: ${COMPRESS_RESULT}")
endif()

if(NOT EXISTS "${MSZ_FILE}")
    message(FATAL_ERROR "msz file was not created: ${MSZ_FILE}")
endif()

# Extract the requested scans.
file(REMOVE "${OUTPUT_FILE}")

message(STATUS "Running: ${EXECUTABLE} --extract --extract-scans ${SCAN_SPEC} ${MSZ_FILE} ${OUTPUT_FILE}")

execute_process(
    COMMAND ${EXECUTABLE} --extract --extract-scans "${SCAN_SPEC}" "${MSZ_FILE}" "${OUTPUT_FILE}"
    RESULT_VARIABLE EXTRACT_RESULT
    OUTPUT_VARIABLE EXTRACT_OUTPUT
    ERROR_VARIABLE EXTRACT_ERROR
)

if(NOT EXTRACT_RESULT EQUAL 0)
    message(STATUS "Extraction output: ${EXTRACT_OUTPUT}")
    message(STATUS "Extraction error: ${EXTRACT_ERROR}")
    message(FATAL_ERROR "Extraction failed with return code: ${EXTRACT_RESULT}")
endif()

if(NOT EXISTS "${OUTPUT_FILE}")
    message(FATAL_ERROR "Output file was not created: ${OUTPUT_FILE}")
endif()

file(READ "${OUTPUT_FILE}" CONTENT)

# Every emitted <spectrum> element must carry a scan= number...
string(REGEX MATCHALL "<spectrum [^>]*scan=[0-9]+" SPECTRUM_TAGS "${CONTENT}")
string(REGEX MATCHALL "<spectrum " SPECTRUM_ELEMENTS "${CONTENT}")

list(LENGTH SPECTRUM_ELEMENTS TOTAL_SPECTRA)
list(LENGTH EXPECTED_SCANS EXPECTED_COUNT)

if(NOT TOTAL_SPECTRA EQUAL EXPECTED_COUNT)
    message(FATAL_ERROR "--extract-scans ${SCAN_SPEC}: expected ${EXPECTED_COUNT} spectra, got ${TOTAL_SPECTRA}")
endif()

# ...and exactly the requested ones. The pattern is anchored to the opening
# <spectrum tag so offsets in a trailing indexList (idRef="...scan=N") never
# match.
set(FOUND_SCANS "")
foreach(TAG IN LISTS SPECTRUM_TAGS)
    string(REGEX REPLACE "^.*scan=" "" SCAN "${TAG}")
    list(APPEND FOUND_SCANS "${SCAN}")
endforeach()

if(NOT FOUND_SCANS STREQUAL "${EXPECTED_SCANS}")
    message(FATAL_ERROR "--extract-scans ${SCAN_SPEC}: expected scans [${EXPECTED_SCANS}], got [${FOUND_SCANS}]")
endif()

# Cleanup
file(REMOVE "${OUTPUT_FILE}")
file(REMOVE "${MSZ_FILE}")

message(STATUS "--extract-scans ${SCAN_SPEC} extracted scans [${FOUND_SCANS}]")
message(STATUS "Test passed!")
