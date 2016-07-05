#include <gtest/gtest.h>

#include "fpotencia.h"
#include "SolverTest.h"
#include "NRcurrentSolverTest.h"


TEST_F(NRcurrentSolverTest, ConvergesOnIeee14BusTest)
{
    auto model = generateIeee14Bus();

    bool estimate_angles = false;
    model.compile(estimate_angles);

    fPotencia::Solver_NRcurrent NRcs(model);
    NRcs.Max_Iter = 6;
    NRcs.EPS = 1e-9;

    auto state = NRcs.solve();

    ASSERT_EQ(fPotencia::Solver_State::Converged, state);
}
