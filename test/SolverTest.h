#ifndef SOLVERTEST_H
#define SOLVERTEST_H


#include <gtest/gtest.h>
#include "Circuit.h"


class SolverTest: public ::testing::Test
{
public:
    SolverTest();
    virtual ~SolverTest() noexcept;
    fPotencia::Circuit generateIeee14Bus() const;
};

#endif // SOLVERTEST_H
