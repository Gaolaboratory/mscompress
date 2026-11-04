# Re-encoded arrays must leave a valid indexedmzML behind.
#
# TEST_MZML is an indexedmzML whose binary arrays are zlib-compressed by a
# deflate other than mscompress's (raw2ms), so the round-trip re-deflates them
# to different lengths. After compress -> decompress the output must still be
# internally consistent:
#   1. every <binaryDataArray encodedLength> equals its <binary> text length,
#   2. every <offset> points at the '<' of a <spectrum element,
#   3. <indexListOffset> points at the '<' of <indexList,
#   4. <fileChecksum> is the SHA-1 of byte 0 .. end of "<fileChecksum>".
# The same holds for --extract of a subset, from the .msz and from the
# .mzML: the index must list exactly the extracted spectra. Then a document
# with uncompressed arrays (BYTE_EXACT_MZML) must still round-trip
# byte-identically, its own index untouched.
#
# Expected variables: EXECUTABLE, TEST_MZML, BYTE_EXACT_MZML, TEMP_DIR

foreach(F "${EXECUTABLE}" "${TEST_MZML}" "${BYTE_EXACT_MZML}")
    if(NOT EXISTS "${F}")
        message(FATAL_ERROR "Missing: ${F}")
    endif()
endforeach()

file(REMOVE_RECURSE "${TEMP_DIR}")
file(MAKE_DIRECTORY "${TEMP_DIR}")

function(run_ok)
    execute_process(COMMAND ${ARGN}
        RESULT_VARIABLE R OUTPUT_VARIABLE O ERROR_VARIABLE E)
    if(NOT R EQUAL 0)
        message(STATUS "out: ${O}\nerr: ${E}")
        message(FATAL_ERROR "command failed (${R}): ${ARGN}")
    endif()
endfunction()

# Asserts OUT is a self-consistent indexedmzML; WANT_SPECTRA (optional) is
# the number of spectra it must index.
function(check_indexed OUT)
    file(READ "${OUT}" DOC)
    # --- 1. encodedLength ---
    set(REST "${DOC}")
    set(N_ARRAYS 0)
    foreach(_unused RANGE 1000000)
        string(FIND "${REST}" "encodedLength=\"" P)
        if(P EQUAL -1)
            break()
        endif()
        math(EXPR P "${P} + 15")
        string(SUBSTRING "${REST}" ${P} -1 REST)
        string(REGEX MATCH "^[0-9]+" WANT "${REST}")
        string(FIND "${REST}" "<binary>" B0)
        string(FIND "${REST}" "</binary>" B1)
        math(EXPR GOT "${B1} - ${B0} - 8")
        if(NOT GOT EQUAL WANT)
            message(FATAL_ERROR "array ${N_ARRAYS}: encodedLength=${WANT} but "
                "<binary> holds ${GOT} characters")
        endif()
        string(SUBSTRING "${REST}" ${B1} -1 REST)
        math(EXPR N_ARRAYS "${N_ARRAYS} + 1")
    endforeach()
    if(N_ARRAYS EQUAL 0)
        message(FATAL_ERROR "no binaryDataArray found in ${OUT}")
    endif()

    # --- 2. spectrum offsets ---
    string(REGEX MATCHALL "<offset idRef=\"[^\"]*\">[0-9]+</offset>" OFFS "${DOC}")
    list(LENGTH OFFS N_OFFS)
    if(N_OFFS EQUAL 0)
        message(FATAL_ERROR "no <offset> entries in ${OUT}")
    endif()
    foreach(E IN LISTS OFFS)
        string(REGEX REPLACE ".*\">([0-9]+)</offset>" "\\1" N "${E}")
        string(SUBSTRING "${DOC}" ${N} 10 AT)
        if(NOT AT STREQUAL "<spectrum ")
            message(FATAL_ERROR "${E} points at '${AT}', not a <spectrum")
        endif()
    endforeach()

    # --- 3. indexListOffset ---
    string(REGEX MATCH "<indexListOffset>([0-9]+)</indexListOffset>" M "${DOC}")
    set(ILO "${CMAKE_MATCH_1}")
    string(SUBSTRING "${DOC}" ${ILO} 10 AT)
    if(NOT AT STREQUAL "<indexList")
        message(FATAL_ERROR "indexListOffset ${ILO} points at '${AT}'")
    endif()

    # --- 4. fileChecksum ---
    string(FIND "${DOC}" "<fileChecksum>" CK)
    math(EXPR CK_END "${CK} + 14")
    string(SUBSTRING "${DOC}" 0 ${CK_END} SPAN)
    string(SHA1 WANT_SHA "${SPAN}")
    string(SUBSTRING "${DOC}" ${CK_END} 40 GOT_SHA)
    if(NOT GOT_SHA STREQUAL WANT_SHA)
        message(FATAL_ERROR "fileChecksum ${GOT_SHA} != SHA-1 ${WANT_SHA}")
    endif()

    list(LENGTH OFFS N_OFFS)
        if(DEFINED WANT_SPECTRA AND NOT N_OFFS EQUAL WANT_SPECTRA)
            message(FATAL_ERROR "${OUT}: ${N_OFFS} index entries, want ${WANT_SPECTRA}")
        endif()
        string(FIND "${DOC}" "><<" DBL)
        if(NOT DBL EQUAL -1)
            message(FATAL_ERROR "${OUT}: malformed \"><<\" at byte ${DBL}")
        endif()
        message(STATUS "${OUT}: ${N_ARRAYS} arrays, ${N_OFFS} offsets, checksum ${GOT_SHA}")

endfunction()

set(MSZ "${TEMP_DIR}/rt.msz")
set(OUT "${TEMP_DIR}/rt.mzML")
run_ok("${EXECUTABLE}" "${TEST_MZML}" "${MSZ}")
run_ok("${EXECUTABLE}" "${MSZ}" "${OUT}")

file(READ "${OUT}" DOC)
string(LENGTH "${DOC}" DOC_LEN)

# The fixture must actually exercise re-encoding: if the output were
# byte-identical to the source, this test would prove nothing.
file(READ "${TEST_MZML}" SRC)
if(SRC STREQUAL DOC)
    message(FATAL_ERROR "fixture round-tripped byte-exact; it no longer "
        "exercises re-encoded arrays")
endif()

check_indexed("${OUT}")

# --- extraction: a subset gets a fresh index of exactly that subset, from
# the .msz (re-encoded arrays) and straight from the .mzML (source bytes) ---
set(WANT_SPECTRA 3)
run_ok("${EXECUTABLE}" --extract --extract-indices "[2,5,7]" "${MSZ}"
    "${TEMP_DIR}/ex_msz.mzML")
check_indexed("${TEMP_DIR}/ex_msz.mzML")
run_ok("${EXECUTABLE}" --extract --extract-indices "[2,5,7]" "${TEST_MZML}"
    "${TEMP_DIR}/ex_mzml.mzML")
check_indexed("${TEMP_DIR}/ex_mzml.mzML")
unset(WANT_SPECTRA)

# --- 5. arrays that round-trip unchanged leave the document byte-exact ---
run_ok("${EXECUTABLE}" "${BYTE_EXACT_MZML}" "${TEMP_DIR}/be.msz")
run_ok("${EXECUTABLE}" "${TEMP_DIR}/be.msz" "${TEMP_DIR}/be.mzML")
file(SHA1 "${BYTE_EXACT_MZML}" H_IN)
file(SHA1 "${TEMP_DIR}/be.mzML" H_OUT)
if(NOT H_IN STREQUAL H_OUT)
    message(FATAL_ERROR "${BYTE_EXACT_MZML} no longer round-trips byte-exact")
endif()
