set(fpotencia_test_SOURCES
    main.cpp)

add_executable(fpotencia_test ${fpotencia_test_SOURCES})
add_test(fpotencia_test fpotencia_test)
target_link_libraries(fpotencia_test
    fpotencia
    ${GTEST_LIBRARIES})
set_property(TARGET fpotencia PROPERTY CXX_STANDARD 11)
