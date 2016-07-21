set(fpotencia_test_SOURCES
    SolverTest.cpp
    JacobiSolverTest.cpp
    ZBusGSSolverTest.cpp
    NRrectSolverTest.cpp
    IwamotoSolverTest.cpp
    NRpolarSolverTest.cpp
    NRcurrentSolverTest.cpp
    main.cpp)

include_directories(${eigen3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_executable(fpotencia_test ${fpotencia_test_SOURCES})
add_test(fpotencia_test fpotencia_test)
target_include_directories(fpotencia_test
    PUBLIC ${fPotencia_SOURCE_DIR}/src)
target_link_libraries(fpotencia_test
    fpotencia
    gtest
    pthread)
set_property(TARGET fpotencia_test PROPERTY CXX_STANDARD 11)
