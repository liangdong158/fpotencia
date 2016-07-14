#include <gtest/gtest.h>

#include "fpotencia.h"
#include "SolverTest.h"
#include "IwamotoSolverTest.h"


TEST_F(IwamotoSolverTest, ConvergesOnIeee14BusTest)
{
    auto model = generateIeee14Bus();

    bool estimate_angles = false;
    model.compile(estimate_angles);

    fPotencia::Solver_Iwamoto NRcs(model);
    auto state = NRcs.solve();

    ASSERT_EQ(fPotencia::Solver_State::Converged, state);
}
