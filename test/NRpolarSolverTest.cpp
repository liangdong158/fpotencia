#include <gtest/gtest.h>

#include "fpotencia.h"
#include "SolverTest.h"
#include "NRpolarSolverTest.h"


TEST_F(NRpolarSolverTest, ConvergesOnIeee14BusTest)
{
    auto model = generateIeee14Bus();

    bool estimate_angles = false;
    model.compile(estimate_angles);

    fPotencia::Solver_NRpolar NRcs(model);
    NRcs.maxIterations = 6;
    NRcs.tolerance = 1e-9;

    auto state = NRcs.solve();

    ASSERT_EQ(fPotencia::Solver_State::Converged, state);
}


TEST_F(NRpolarSolverTest, SolvesLynnPowellBusTest)
{
    auto model = generateLynnPowell();
    model.compile(false);

    fPotencia::Solver_NRpolar solver(model);
    solver.tolerance = 1e-9;
    solver.maxIterations = 100;
    auto state = solver.solve();

    ASSERT_EQ(fPotencia::Solver_State::Converged, state);

    for (auto& bus: solver.Model.buses) {
        bus.print();
    }
}
