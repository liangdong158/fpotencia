/*
 * File:   Solver_NRpolar.h
 * Author: Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#ifndef SOLVER_IWAMOTO_H
#define	SOLVER_IWAMOTO_H


#include <cmath>
#include "Circuit.h"
#include "Solution.h"

//using namespace arma;
using namespace std;

namespace fPotencia {

    /*
     * This class implements the Iwamoto method in current equations
     * for PQ buses and Power equations for PV buses according to the paper:
     * Improvement of Power Flow Calculation with Optimization Factor Based on 
     * Current Injection Method
     */
    class Solver_Iwamoto {
    public:


        /*! 
         * \brief The default tolerance for the solution
         *
         * Solving the power flow equations is an interative process. For
         * every iteration, the parameters are adjusted in order to reach
         * convergence. The tolerance defines the allowable deviation/error
         * after which the process halts.
         */
        static constexpr const double DEFAULT_SOLUTION_TOLERANCE = 1e-3;


        /*!
         * \brief Default maximum number of iterations after which the solver
         *  declares failure
         */
        static constexpr const unsigned DEFAULT_MAX_ITERATIONS = 10;


        Solver_Iwamoto(Circuit model);

        Solver_Iwamoto(Circuit model, solution sol_);

        virtual ~Solver_Iwamoto();

        /*Properties*/
        Circuit Model;


        /*!
         * \brief Allowable tolerance of the solver instance
         *
         * \sa DEFAULT_SOLUTION_TOLERANCE
         */
        double tolerance;


        /*!
         * \brief Maximum number of iterations
         *
         * \sa DEFAULT_MAX_ITERATIONS
         */
        int maxIterations;

        Solver_State solve(); //Solves the grid

        void update_solution_power_from_circuit();

    private:

        vector<int> PQPV;

        vector<int> PQPVPV;

        solution Sol;

        vec Pesp;

        vec Qesp;

        /*
         * This function returns the calculated increments of x
         */
        vec y(solution sol, uint N, uint npv);
        
        /*
         * This function returns the specified increments of x 
         */
        vec ys(solution sol, uint N, uint npv);

        /*
         * Calculates the Jacobian wich is passed by refference
         */
        mat Jacobian(solution & sol, uint N, uint npv);

        /*Calculates the optimal Iwamoto multiplicator mu
         */
        double compute_mu(vec a, vec b, vec c, uint Nj, double mu0);

        /*
         * Solves a polynomial defined by the coeffients g0, g1, g3 and g3 such
         * that:  d + c*x + b*x^2 + a*x^3 = 0
         * 
         * Provides the real solution using the Newton-Raphson technique
         */
        double solve_poly_deg3(double d, double c, double b, double a, double x);


        bool converged(vec X, uint Nj); //check if the solution converged


        vec pack_solution(solution sol, uint N, uint npv);


        solution unpack_solution(vec x,  uint N, uint npv);


        void calculate_slack_power(); //calculate the slack bus power   


        bool checks(); //check the solvability with this method


        void fill_especifyed_values();
    };

#endif	/* SOLVER_IWAMOTO_H */

}
