/*!
 * \file Solver_NRpolar.h
 * \author Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SOLVER_NRPOLAR_H
#define	SOLVER_NRPOLAR_H


#include <cmath>
#include "Circuit.h"
#include "Solution.h"

//using namespace arma;
using namespace std;

namespace fPotencia {

    /*!
     * \brief This class implements the Nerwton-Raphson method of load flow
     *  analysis using polar coordinates.
     */
    class Solver_NRpolar
    {
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


        /*!
         * \brief Creates a new solver
         *
         * \param[in] model The circuit the solver should perform load flow
         *  analysis on
         */
        Solver_NRpolar(Circuit const& model);


        /*!
         * \brief Constructs a new solver instance that works off an initial
         *  solution
         *
         * \param[in] model The circuit the solver should perform load flow
         *  analysis on
         *
         * \param[in] sol_ The initial solution the solver should start
         *  working with
         */
        Solver_NRpolar(Circuit const& model, solution const& sol_);


        virtual ~Solver_NRpolar() noexcept;


        //!  \brief The circuit model the solver analyses
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
        unsigned maxIterations;


        /*!
         * \brief Solves the grid
         *
         * \return Solver_State::converged if the grid was solved
         *
         * \sa Solver_State
         */
        Solver_State solve();


        void update_solution_power_from_circuit();
        
    private:

        vector<int> BUSES;

        vector<int> PQPV;

        vector<int> LastPQ;
        
        vector<int> LastPV;

        vec Pesp;

        vec Qesp;

        solution Sol;

        void Jacobian(mat &J, vec &V, vec &D, uint npq, uint npv); //calculate the jacobian, J is passed by refference
        
        double mu(mat &J, mat &J2, vec &F, vec &dV, vec &dD, vec & dx, uint npq, uint npv);
        
        double solve_poly_deg3(double d, double c, double b, double a, double x) ;

        void get_power_inc(vec &PQinc, uint npq, uint npv); //PQinc is passed by refference

        void calculate_Q(uint npq, uint npv); //calculate the reative power at the PV buses

        double Q(uint k);

        double P(uint k);

        bool converged(vec PQinc, uint npqpvpq); //check if the solution converged

        void update_solution(vec X, uint npq, uint npv);
        
        void get_increments(vec X, vec &incV, vec &incD, uint npq, uint npv);

        void calculate_slack_power(); //calculate the slack bus power        


        //! \brief Checks whether the solver can work on the given model
        bool checks() const;


        void correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X);
        
    };

#endif	/* SOLVER_NRPOLAR_H */

}
