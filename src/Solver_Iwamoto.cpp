/*
 * File:   Solver_NRpolar.cpp
 * Author: Santiago Peñate Vera
 *
 * Created on 25 of January of 2015, 23:05
 * Copyright (C) 2015 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/Polynomials>

#include "Solver_Iwamoto.h"

//using namespace arma;
using namespace std;


namespace fPotencia {

    /* Constructor
     */
    Solver_Iwamoto::Solver_Iwamoto(Circuit model):
            tolerance(DEFAULT_SOLUTION_TOLERANCE),
            maxIterations(DEFAULT_MAX_ITERATIONS)
    {
        Model = model;

        Sol = Model.get_initial_solution();
        if (!Sol.initialized) {
            Model.compile(false);
            Sol = Model.get_initial_solution();
        }

        PQPV.reserve(2 * Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        PQPV.insert(PQPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPV.insert(PQPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        PQPVPV.reserve(2 * Model.PQ_list.size() + 2 * Model.PV_list.size()); // preallocate memory
        PQPVPV.insert(PQPVPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPVPV.insert(PQPVPV.end(), Model.PV_list.begin(), Model.PV_list.end());
        PQPVPV.insert(PQPVPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        fill_especifyed_values();

        if (!checks())
            throw std::invalid_argument("The circuit failed the solver compatibility test.");
    }

    /* Constructor where the solution is given
     */
    Solver_Iwamoto::Solver_Iwamoto(Circuit model, solution sol_) {
        Model = model;

        Sol = sol_;

        PQPV.reserve(2 * Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        PQPV.insert(PQPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPV.insert(PQPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        PQPVPV.reserve(2 * Model.PQ_list.size() + 2 * Model.PV_list.size()); // preallocate memory
        PQPVPV.insert(PQPVPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPVPV.insert(PQPVPV.end(), Model.PV_list.begin(), Model.PV_list.end());
        PQPVPV.insert(PQPVPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        fill_especifyed_values();

        if (!checks())
            throw std::invalid_argument("The circuit failed the solver compatibility test.");
    }


    Solver_Iwamoto::~Solver_Iwamoto() noexcept
    {
    }


    std::pair<cx_mat, cx_mat>
    Solver_Iwamoto::partialPowerInjectionDerivatives(
            cx_mat const& ybus,
            cx_vec const& voltages)
            const
    {

        // Python code is in the comments.

        // Ibus = Ybus * asmatrix(V).T
        auto ibus = ybus * voltages.transpose();

        // diagV = asmatrix(diag(V))
        auto diagV = voltages.diagonal();

        // diagIbus = asmatrix(diag( asarray(Ibus).flatten()))
        cx_mat diagIbus = cx_mat::Zero(
                ibus.cols() * ibus.rows(),
                ibus.cols() * ibus.rows());
        for (int i = 0, r = 0; r != ibus.rows(); ++r) {
            for (int c = 0; c != ibus.cols(); ++c, i++) {
                diagIbus(i, i) = ibus(r, c);
            }
        }

        // diagVnorm = asmatrix(diag(V / abs(V)))
        mat diagVnorm = voltages;
        for (int i = 0; i != voltages.size(); ++i) {
            diagVnorm[i] /= std::abs(voltages[i]);
        }
        diagVnorm = diagVnorm.diagonal();

        // dS_dVm = diagV * conj(Ybus * diagVnorm)
        //  + conj(diagIbus) * diagVnorm
        auto dS_dVm = diagV * (ybus * diagVnorm).conjugate()
                + diagIbus.conjugate() * diagVnorm;

        // dS_dVa = 1j * diagV * conj(diagIbus - Ybus * diagV)
        mat dS_dVa = 1j * diagV * (diagIbus - ybus * diagV).conjugate();

        // return dS_dVm, dS_dVa
        return std::make_pair(dS_dVm, dS_dVa);
    }


    double Solver_Iwamoto::mu(
            const cx_mat& ybus,
            const cx_mat& jacobian,
            const vec& mismatches,
            const cx_vec& voltages,
            const vec& solution,
            const Indices& pvpq,
            const Indices& pq)
            const
    {
        assert(pvpq.size() == pq.size());

        auto partials = partialPowerInjectionDerivatives(ybus, voltages);
        auto dSdVm = partials.first, dSdVa = partials.second;

        // J11 = dS_dVa[array([pvpq]).T, pvpq].real
        mat j11(pvpq.size());
        for (auto const& i: pvpq) {
            j11 << dSdVa(i, i).real();
        }

        // J12 = dS_dVm[array([pvpq]).T, pq].real     # [row, col]?
        mat j12(pvpq.size());
        for (std::vector<int>::size_type i = 0; i != pvpq.size(); ++i) {
            j12 << dSdVm(pvpq[i], pq[i]).real();
        }

        // J21 = dS_dVa[array([pq]).T, pvpq].imag
        mat j21(pvpq.size());
        for (std::vector<int>::size_type i = 0; i != pvpq.size(); ++i) {
            j21 << dSdVa(pq[i], pvpq[i]).imag();
        }

        // J22 = dS_dVm[array([pq]).T, pq].imag
        mat j22(pvpq.size());
        for (std::vector<int>::size_type i = 0; i != pvpq.size(); ++i) {
            j22 << dSdVa(pq[i], pq[i]).imag();
        }

        // Theoretically this is the second derivative matrix
        // since the Jacobian has been calculated with dV instead of V
        // J2 = vstack([
        //     hstack([J11, J12]),
        //     hstack([J21, J22]) ], format="csr")
        mat j11j12(j11.rows(), j11.cols() + j12.cols());
        j11j12<< j11, j12;
        mat j21j22(j21.rows(), j21.cols() + j22.cols());
        j21j22 << j11, j12;
        mat j2(j11j12.rows() + j21j22.rows(), j11j12.cols());
        j2 << j11j12, j21j22;

        // a = F
        mat a = mismatches;

        // b = J * dx
        mat b = jacobian * solution;

        // c = 0.5 * dx * J2 * dx
        mat c = 0.5 * solution * j2 * solution;

        // g0 = -a.dot(b)
        auto g0 = -a.dot(b);

        // g1 = b.dot(b) + 2 * a.dot(c)
        auto g1 = b.dot(b) + 2 * a.dot(c);

        // g2 = -3.0 * b.dot(c)
        auto g2 = -3.0 * b.dot(c);

        //  g3 = 2.0 * c.dot(c)
        auto g3 = 2.0 * c.dot(c);

        // roots = np.roots([g3, g2, g1, g0])
        cx_vec polynomial(4), roots;
        polynomial << g3, g2, g1, g0;
        roots_to_monicPolynomial(polynomial, roots);

        // Three solutions are provided, the first two are complex;
        // only the real solution is valid:
        // return roots[2].real
        return roots[2].real();
    }


    //check the solvability with this method

    bool Solver_Iwamoto::checks() {
        bool val = true;

        //Only one slack bus allowed
        if (Model.VD_list.size() > 1)
            val = false;

        return val;
    }

    /*Fills the initial power values
     */
    void Solver_Iwamoto::fill_especifyed_values() {
        Pesp = vec(Model.buses.size());
        Qesp = vec(Model.buses.size());

        for (uint k : PQPV) {
            Pesp(k) = Sol.P.coeff(k); //P at PQ and PV buses
            Qesp(k) = Sol.Q.coeff(k); //Q at PQ buses
        }

    }

    vec Solver_Iwamoto::y(solution sol, uint N, uint npv) {

        uint Nj = 2 * N + npv;
        vec x(Nj);
        uint a = 0;
        uint k;
        double val1, val2;
        for (uint idx = 0; idx < N + npv; idx++) { //filas
            k = PQPVPV[idx];
            if (idx < N) {
                //Increment of real current
                val1 = 0;
                val2 = 0;
                for (uint i = 0; i < Model.buses.size(); i++) {
                    val1 += sol.V.coeff(i) * (Model.G(k, i) * cos(sol.D.coeff(k)) - Model.B(k, i) * sin(sol.D.coeff(k)));
                    val2 += sol.V.coeff(i) * (Model.G(k, i) * sin(sol.D.coeff(k)) + Model.B(k, i) * cos(sol.D.coeff(k)));
                }
                x[a] = val1; //1

                //increment of imaginary current
                x[a + 1] = val2; //2

                a += 2;
            } else {
                //Increment of active power
                val1 = 0;
                for (uint i = 0; i < Model.buses.size(); i++) {
                    val1 += sol.V.coeff(i)*(Model.G(k, i) * cos(sol.D.coeff(k) - sol.D.coeff(i)) + Model.B(k, i) * sin(sol.D.coeff(k) - sol.D.coeff(i)));
                }
                x[a] = sol.V.coeff(k) * val1; //3

                a += 1;
            }
        }

        return x;
    }

    vec Solver_Iwamoto::ys(solution sol, uint N, uint npv) {

        uint Nj = 2 * N + npv;
        vec x(Nj);
        uint a = 0;
        uint k;

        for (uint idx = 0; idx < N + npv; idx++) { //filas
            k = PQPVPV[idx];
            if (idx < N) {
                //Increment of real current
                x[a] = ((Pesp.coeff(k) * cos(sol.D.coeff(k)) + Qesp.coeff(k) * sin(sol.D.coeff(k))) / sol.V.coeff(k)); //1

                //increment of imaginary current
                x[a + 1] = ((Pesp.coeff(k) * sin(sol.D.coeff(k)) - Qesp.coeff(k) * cos(sol.D.coeff(k))) / sol.V.coeff(k)); //2

                a += 2;
            } else {
                //Increment of active power
                x[a] = Pesp.coeff(k); //3

                a += 1;
            }
        }

        return x;
    }

    /*
     * Calculates the Jacobian wich is passed by refference
     */
    mat Solver_Iwamoto::Jacobian(solution &sol, uint N, uint npv) {

        /* indices: x, y: conceptual bus index
         *          i, k: real bus indices
         *          a, b: jacobian indices
         *
         */
        uint Nj = 2 * N + npv;
        mat J(Nj, Nj);
        J.setZero(2 * N + npv, 2 * N + npv);

        //W
        uint a = 0;
        uint b, k, i;
        double val;
        for (uint x = 0; x < N; x++) { //rows
            b = 0;
            k = PQPVPV[x];
            for (uint y = 0; y < N; y++) { //cols
                i = PQPVPV[y];

                //H
                if (i == k)//a == b
                    J(a, b) = sol.V.coeff(k) * (Model.G(k, k) * sin(sol.D.coeff(k)) + Model.B(k, k) * cos(sol.D.coeff(k)))
                    + (Qesp.coeff(k) * cos(sol.D.coeff(k)) - Pesp.coeff(k) * sin(sol.D.coeff(k))) / sol.V.coeff(k); //11
                else
                    J(a, b) = sol.V.coeff(i) * (Model.G(k, i) * sin(sol.D.coeff(i)) + Model.B(k, i) * cos(sol.D.coeff(i))); //1

                //N
                if (i == k)//a == b + 1
                    J(a, b + 1) = -(Model.G(k, k) * cos(sol.D.coeff(k)) - Model.B(k, k) * sin(sol.D.coeff(k)))
                    - (Pesp.coeff(k) * cos(sol.D.coeff(k)) + Qesp.coeff(k) * sin(sol.D.coeff(k))) / sol.V.coeff(k) / sol.V.coeff(k); //22
                else
                    J(a, b + 1) = Model.B(k, i) * sin(sol.D.coeff(k)) - Model.G(k, i) * cos(sol.D.coeff(k)); //2

                //J
                if (i == k) //a + 1 == b
                    J(a + 1, b) = -sol.V.coeff(k) * (Model.G(k, k) * cos(sol.D.coeff(k)) - Model.B(k, k) * sin(sol.D.coeff(k)))
                    + (Qesp.coeff(k) * sin(sol.D.coeff(k)) + Pesp.coeff(k) * cos(sol.D.coeff(k))) / sol.V.coeff(k); //33
                else
                    J(a + 1, b) = 0; //3

                //L
                if (i == k) //a + 1 == b + 1
                    J(a + 1, b + 1) = -Model.B(k, k) * cos(sol.D.coeff(k)) - Model.G(k, k) * sin(sol.D.coeff(k))
                    -(Pesp.coeff(k) * sin(sol.D.coeff(k)) - Qesp.coeff(k) * cos(sol.D.coeff(k))) / sol.V.coeff(k) / sol.V.coeff(k); //44
                else
                    J(a + 1, b + 1) = -Model.B(k, i) * cos(sol.D.coeff(k)) - Model.G(k, i) * sin(sol.D.coeff(k)); //4


                //only to debug the J structure
                //J(a, b) = i * 10 + k;
                //J(a + 1, b) = i * 10 + k;
                //J(a, b + 1) = i * 10 + k;
                //J(a + 1, b + 1) = i * 10 + k;

                b += 2;
            }
            a += 2;
        }

        //K
        a = 0;
        for (uint x = 0; x < N; x++) { //rows
            k = PQPVPV[x];
            b = 2 * N;
            for (uint y = 0; y < npv; y++) { //cols
                i = PQPVPV[y];
                if (Model.buses[k].Type == BusType::PV) {

                    J(a, b) = sin(sol.D.coeff(k)) / sol.V.coeff(k); //8

                    J(a + 1, b) = -cos(sol.D.coeff(k)) / sol.V.coeff(k); //81

                    //only to debug the J structure
                    //J(a, b) = i * 10 + k;
                    //J(a + 1, b) = i * 10 + k;
                }
                b += 1;
            }
            a += 2;
        }

        //M
        a = 2 * N;
        for (uint x = 0; x < npv; x++) { //rows
            k = PQPVPV[x];
            b = 0;
            for (uint y = 0; y < N; y++) { //cols
                i = PQPVPV[y];
                if (i == k) { // a == b
                    //Fkk
                    val = 0;
                    for (uint j = 0; j < Model.buses.size(); j++)
                        if (j != k)
                            val += sol.V.coeff(i) * (Model.G(k, i) * sin(sol.D.coeff(k) - sol.D.coeff(i)) - Model.B(k, i) * cos(sol.D.coeff(k) - sol.D.coeff(i)));
                    J(a, b) = sol.V.coeff(k) * val; //9

                    //Rkk
                    J(a, b + 1) = 0; //91
                } else {
                    //Fki
                    J(a, b) = -sol.V.coeff(k) * sol.V.coeff(i) *
                            (Model.G(k, i) * sin(sol.D.coeff(k) - sol.D.coeff(i)) - Model.B(k, i) * cos(sol.D.coeff(k) - sol.D.coeff(i))); //9

                    //Rki
                    J(a, b + 1) = -sol.V.coeff(k) * sol.V.coeff(i) *
                            (Model.G(k, i) * cos(sol.D.coeff(k) - sol.D.coeff(i)) + Model.B(k, i) * sin(sol.D.coeff(k) - sol.D.coeff(i))); //91
                }

                //only to debug the J structure
                //J(a, b) = i * 10 + k;
                //J(a, b + 1) = i * 10 + k;

                b += 2;
            }
            a += 1;
        }

        return J;
    }

    /*
     * Solves a polynomial defined by the coeffients g0, g1, g3 and g3 such that
     * d + c*x + b*x^2 + a*x^3 = 0
     *
     * Provides the real solution using the Newton-Raphson technique
     */
    double Solver_Iwamoto::solve_poly_deg3(double d, double c, double b, double a, double x) {
        double eps = 1e-9;

        double fx = a * x * x * x + b * x * x + c * x + d;
        double fxd = 3 * a * x * x + 2 * b * x + c;
        double incx = fx / fxd;

        while (abs(incx) > eps) {
            x -= incx;
            fx = a * x * x * x + b * x * x + c * x + d;
            fxd = 3 * a * x * x + 2 * b * x + c;
            incx = fx / fxd;
        }
        return x;
    }

    /*Computes the optimal multiplier
     */
    double Solver_Iwamoto::compute_mu(vec a, vec b, vec c, uint Nj, double mu0) {
        double g0 = 0.0;
        double g1 = 0.0;
        double g2 = 0.0;
        double g3 = 0.0;

        for (uint i = 0; i < Nj; i++) {
            g0 += a.coeff(i) * b.coeff(i);
            g1 += b.coeff(i) * b.coeff(i) + 2.0 * a.coeff(i) * c.coeff(i);
            g2 += b.coeff(i) * c.coeff(i);
            g3 += c.coeff(i) * c.coeff(i);
        }
        g2 *= 3.0;
        g3 *= 2.0;

        // calculo del multiplicador optimo
        return solve_poly_deg3(g3, g2, g1, g0, mu0);
    }


    bool Solver_Iwamoto::converged(const vec &solution) const
    {
        return solution.lpNorm<Infinity>()< tolerance;
    }


    /* Generates a vector X from a solution object
     */
    vec Solver_Iwamoto::pack_solution(solution sol, uint N, uint npv) {

        uint Nj = 2 * N + npv;
        vec x(Nj);

        uint a = 0;
        uint k;
        for (uint idx = 0; idx < N + npv; idx++) { //filas
            k = PQPVPV[idx];
            if (idx < N) {
                //Increment of voltage angle
                x[a] = sol.D.coeff(k); //1

                //increment of voltage magnitude
                x[a + 1] = sol.V.coeff(k); //2

                a += 2;
            } else {
                //Increment of reactive power
                x[a] = sol.Q.coeff(k); //3

                a += 1;
            }
        }

        return x;
    }

    /* Generates a solution object from a vector X
     */
    solution Solver_Iwamoto::unpack_solution(vec x, uint N, uint npv) {

        solution sol;
        sol.resize(Sol.Lenght);
        //sol.clear();

        uint a = 0;
        uint k;
        for (uint idx = 0; idx < N + npv; idx++) { //filas
            k = PQPVPV[idx];
            if (idx < N) {
                //Increment of voltage angle
                sol.D(k) = x[a]; //1

                //increment of voltage magnitude
                sol.V(k) = x[a + 1]; //2

                a += 2;
            } else {
                //Increment of reactive power
                sol.Q(k) = x[a]; //3

                a += 1;
            }
        }

        //copy the slack values
        for (uint k : Model.VD_list) {
            sol.V(k) = Sol.V.coeff(k);
            sol.D(k) = Sol.D.coeff(k);
            sol.P(k) = Sol.P.coeff(k);
            sol.Q(k) = Sol.Q.coeff(k);
        }

        return sol;
    }

    /*calculate the slack bus power
     */
    void Solver_Iwamoto::calculate_slack_power() {
        for (uint k : Model.VD_list) {
            cx_double I(0.0, 0.0);
            for (uint j = 0; j < Model.buses.size(); j++) {
                I += Model.Y.coeff(k, j) * Sol.Vcx(j);
            }
            I = Sol.Vcx(k) * conj(I); //now this is the power
            Sol.P(k) = I.real();
            Sol.Q(k) = I.imag();
        }
    }

    /*Solves the grid
     */
    Solver_State Solver_Iwamoto::solve() {
        auto npv = Model.PV_list.size();
        auto npq = Model.PQ_list.size();
        auto N = npq + npv;
        auto Nj = 2 * N + npv;


        solution inc_sol;
        inc_sol.resize(Sol.Lenght);
        vec b, c;
        double mu = 0.0;

        mat J = Jacobian(Sol, N, npv);
        Eigen::FullPivLU<mat>LU(J); //Full pivot LU (only once for Iwamoto)

        vec y_s = ys(Sol, N, npv);

        vec y_xe = y(Sol, N, npv);

        vec a = y_s - y_xe;

        vec inc_x = -1 * LU.solve(a);

        vec x_e = pack_solution(Sol, N, npv);


        for (unsigned i = 0; ! converged(inc_x) && i < maxIterations; ++i) {

            b = -1.0 * a;

            inc_sol = unpack_solution(x_e + inc_x, N, npv);
            c = -1.0 * y(inc_sol, N, npv);

            mu = compute_mu(a, b, c, Nj, mu);

            // update solution vector
            x_e += mu * inc_x;

            // Update calculation values
            inc_sol = unpack_solution(x_e, N, npv); //recycling inc_sol
            //J = Jacobian(inc_sol, N, npv);
            y_xe = y(inc_sol, N, npv);
            a = y_s - y_xe;

            inc_x = -1 * LU.solve(a);
        }

        //Unpacks the solution vector incX to the solution
        Sol = unpack_solution(x_e, N, npv);

        return converged(inc_x)
                ? Solver_State::Converged
                : Solver_State::Not_Converged;
    }

    /*
     */
    void Solver_Iwamoto::update_solution_power_from_circuit() {

    }

}
