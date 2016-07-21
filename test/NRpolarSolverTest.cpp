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
