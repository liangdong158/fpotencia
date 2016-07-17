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


        virtual ~Solver_Iwamoto() noexcept;


        /*!
         * \brief The model we solve
         *
         * \todo This should become a per-call variable, and not remain an
         *  instance variable, to allow for multithreading.
         */
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
         * \brief Computes the partial derivatives of power injection
         *  with regards to the voltage.
         *
         * The following explains the expressions used to form the
         * matrices:
         *
         *  S = diag(V) * conj(ibus) = diag(conj(ibus)) * V
         *
         * Partials of V & ibus w.r.t. voltage magnitudes:
         *
         *  dV/dVm = diag(V / abs(V))
         *  dI/dVm = Ybus * dV/dVm = Ybus * diag(V / abs(V))
         *
         * Partials of V & Ibus w.r.t. voltage angles:
         *
         *  dV/dVa = j * diag(V)
         *  dI/dVa = ybus * dV/dVa = ybus * j * diag(V)
         *
         * Partials of S w.r.t. voltage magnitudes:
         *
         *  dS/dVm = diag(V) * conj(dI/dVm) + diag(conj(ibus)) * dV/dVm
         *      = diag(V) * conj(ybus * diag(V / abs(V)))
         *          + conj(diag(ibus)) * diag(V / abs(V))
         *
         * Partials of S w.r.t. voltage angles:
         *
         *  dS/dVa = diag(V) * conj(dI/dVa) + diag(conj(Ibus)) * dV/dVa
         *      = diag(V) * conj(ybus * j * diag(V))
         *          + conj(diag(ibus)) * j * diag(V)
         *      = -j * diag(V) * conj(ybus * diag(V))
         *          + conj(diag(Ibus)) * j * diag(V)
         *      = j * diag(V) * conj(diag(ibus) - ybus * diag(V))
         *
         * For more details on the derivations behind the derivative code
         * used (in MATPOWER), see:
         * R. D. Zimmerman, "AC Power Flows, Generalized OPF Costs and
         * their Derivatives using Complex Matrix Notation", MATPOWER
         * Technical Note 2, February 2010.
         * http://www.pserc.cornell.edu/matpower/TN2-OPF-Derivatives.pdf
         *
         * \return Two matrices containing the partial derivatives of the
         *  complex bus power injections with regards to the (1) voltage
         *  magnitude and (2) the voltage angle respectively for all buses.
         *  If `ybus` is a sparse matrix, the return values will also be
         *  sparse.
         */
        std::pair<mat, mat> partialPowerInjectionDerivatives(
                mat const& ybus,
                vec const& voltages)
                const;


        /*!
         * \brief Calculates µ, the Iwamoto acceleration parameter.
         *
         * Calculates the Iwamoto acceleration parameter as described in:
         * "A Load Flow Calculation Method for Ill-Conditioned Power Systems"
         * by Iwamoto, S. and Tamura, Y.
         *
         * \param[in] ybus Admittance matrix
         *
         * \param[in] jacobian Jacobian matrix
         *
         * \param[in] mismatches Mismatch vector
         *
         * \param[in] dV Voltage increment in complex form
         *
         * \param[in] solution Current solution vector as calculated by
         *  `solve(jacobian, mismatches)`
         *
         * \param[in] pvpq PQ and PV indices
         *
         * \param[in] pq PQ indices
         *
         * \return The optimal multiplier for ill-conditioned systems
         */
        double mu(
                mat const& ybus,
                mat const& jacobian,
                vec const& mismatches,
                mat const& voltages,
                vec const& solution,
                std::vector<int> const& pvpq,
                std::vector<int> const& pq)
                const;



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


        /*!
         * \brief Checks whether the current solution marks convergence
         *
         * The solver has converged on a solution if all coefficients in
         * the given vector are smaller than the allowed tolerance.
         *
         * \param[in] solution solution vector
         *
         * \return `true` on convergence, `false` otherwise
         */
        bool converged(vec const& solution) const;


        vec pack_solution(solution sol, uint N, uint npv);


        solution unpack_solution(vec x,  uint N, uint npv);


        void calculate_slack_power(); //calculate the slack bus power


        bool checks(); //check the solvability with this method


        void fill_especifyed_values();
    };

#endif	/* SOLVER_IWAMOTO_H */

}
