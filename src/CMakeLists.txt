set(fpotencia_SOURCES
    Solver.cpp
    UnfittingSolverException.cpp
    LineConstructors.cpp
    Solver_ZBusGS.cpp
    Solution.cpp
    Circuit.cpp
    Generator.cpp
    TransformerConstructors.cpp
    Solver_Iwamoto.cpp
    Graph.cpp
    Solver_NRcurrent.cpp
    Shunt.cpp
    Solver_Jacobi.cpp
    Circuit3.cpp
    Load.cpp
    Line.cpp
    Solver_NRrect.cpp
    Bus.cpp
    Solver_NRpolar.cpp
    CommonFunctions.cpp
    Transformer.cpp
    ExternalGrid.cpp)

set(fpotencia_HEADERS
    Solver.h
    UnfittingSolverException.h
    Circuit3.h
    Circuit.h
    TransformerConstructors.h
    Shunt.h
    Transformer.h
    Graph.h
    Line.h
    enumaratons.h
    Solver_ZBusGS.h
    Load.h
    Bus.h
    fpotencia_libs.h
    Solution.h
    Solver_NRrect.h
    Solver_Jacobi.h
    fpotencia.h
    ExternalGrid.h
    Solver_NRpolar.h
    LineConstructors.h
    Solver_Iwamoto.h
    Solver_NRcurrent.h
    CommonFunctions.h
    Generator.h)

include_directories(${eigen3_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
add_library(fpotencia ${fpotencia_SOURCES})
target_link_libraries(fpotencia
    Boost::graph)
set_property(TARGET fpotencia PROPERTY CXX_STANDARD 11)
set_property(TARGET fpotencia PROPERTY SOVERSION ${VERSION})
install(TARGETS fpotencia
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${fpotencia_HEADERS}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/fpotencia)
