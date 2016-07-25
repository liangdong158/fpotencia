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

TEST_F(NRcurrentSolverTest, SolvesLynnPowellBusTest)
{
    auto model = generateLynnPowellWithoutGenerator();

    bool estimate_angles = false;
    model.compile(estimate_angles);

    fPotencia::Solver_NRcurrent NRcs(model);
    NRcs.Max_Iter = 60;
    NRcs.EPS = 1e-9;

    auto state = NRcs.solve();

    ASSERT_EQ(fPotencia::Solver_State::Converged, state);

    ASSERT_NEAR(NRcs.Model.buses.at(1).voltage_pu.imag(), -0.026648, NRcs.EPS);
    ASSERT_NEAR(NRcs.Model.buses.at(1).voltage_pu.real(), 0.974707, NRcs.EPS);
}
