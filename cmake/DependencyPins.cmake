# fpm.toml owns the dependency revisions for every build entrypoint.
file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/../fpm.toml" _libneo_dependencies
    REGEX "^(fortio|fortnum) = ")
foreach(_libneo_dependency IN ITEMS fortio fortnum)
    set(_libneo_found OFF)
    foreach(_libneo_line IN LISTS _libneo_dependencies)
        if(_libneo_line MATCHES "^${_libneo_dependency} = ")
            string(REGEX MATCH "rev = \"([0-9a-f]+)\"" _libneo_match
                "${_libneo_line}")
            set(_libneo_revision "${CMAKE_MATCH_1}")
            string(LENGTH "${_libneo_revision}" _libneo_length)
            if(NOT _libneo_length EQUAL 40)
                message(FATAL_ERROR "Invalid ${_libneo_dependency} revision in fpm.toml")
            endif()
            string(TOUPPER "${_libneo_dependency}" _libneo_upper)
            set(${_libneo_upper}_REF "${_libneo_revision}" CACHE STRING
                "${_libneo_dependency} commit, tag, or branch used for the build")
            set(_libneo_found ON)
        endif()
    endforeach()
    if(NOT _libneo_found)
        message(FATAL_ERROR "Missing ${_libneo_dependency} revision in fpm.toml")
    endif()
endforeach()
