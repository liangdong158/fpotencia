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

#include "Solver_NRpolar.h"
#include "Circuit.h"

namespace fPotencia {

    /*
     * constructor
     */
    Solver_NRpolar::Solver_NRpolar(Circuit model) {

        Model = model;

        Sol = Model.get_initial_solution();
        if (!Sol.initialized) {
            Model.compile(false);
            Sol = Model.get_initial_solution();
        }

        BUSES.reserve(Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        BUSES.insert(BUSES.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        BUSES.insert(BUSES.end(), Model.PV_list.begin(), Model.PV_list.end());

        PQPV.reserve(2 * Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        PQPV.insert(PQPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPV.insert(PQPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        LastPQ.reserve(Model.PQ_list.size()); // preallocate memory
        LastPQ.insert(LastPQ.end(), Model.PQ_list.begin(), Model.PQ_list.end());

        LastPV.reserve(Model.PV_list.size()); // preallocate memory
        LastPV.insert(LastPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        Pesp = Sol.P;
        Qesp = Sol.Q;
        
        if (!checks())
            throw std::invalid_argument("The circuit failed the solver compatibility test.");
    }

    /*
     * constructor, given previous solution
     */
    Solver_NRpolar::Solver_NRpolar(Circuit model, solution sol_) {

        Model = model;

        Sol = sol_;

        BUSES.reserve(Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        BUSES.insert(BUSES.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        BUSES.insert(BUSES.end(), Model.PV_list.begin(), Model.PV_list.end());

        PQPV.reserve(2 * Model.PQ_list.size() + Model.PV_list.size()); // preallocate memory
        PQPV.insert(PQPV.end(), Model.PQ_list.begin(), Model.PQ_list.end());
        PQPV.insert(PQPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        LastPQ.reserve(Model.PQ_list.size()); // preallocate memory
        LastPQ.insert(LastPQ.end(), Model.PQ_list.begin(), Model.PQ_list.end());

        LastPV.reserve(Model.PV_list.size()); // preallocate memory
        LastPV.insert(LastPV.end(), Model.PV_list.begin(), Model.PV_list.end());

        Pesp = Sol.P;
        Qesp = Sol.Q;
        
        if (!checks())
            throw std::invalid_argument("The circuit failed the solver compatibility test.");
    }

    /*
     * destructor
     */
    Solver_NRpolar::~Solver_NRpolar() {
    }

    /*//////////////////////////////////////////////////////////////////////////
     * checks wether if the grid can be solved with this method
     */
    bool Solver_NRpolar::checks() {
        bool val = true;

        //Only one slack bus allowed
        if (Model.VD_list.size() > 1)
            val = false;

        return val;
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the slack bus power
     */
    void Solver_NRpolar::calculate_slack_power() {        
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

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the reactive power of the bus k (usefull for PV uses)
     */
    void Solver_NRpolar::calculate_Q(uint npq, uint npv) {
        double val;
        uint k;
        for (uint i = npq - 1; i < npq + npv; i++) {
            k = PQPV[i];
            val = Q(k);
            Sol.Q(k) = val;
        }
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the active power at a bus
     */
    double Solver_NRpolar::P(uint k) {
        double val = 0.0;
        for (uint j = 0; j < Model.buses.size(); j++) {
            val += Sol.V.coeff(j)
                    *(Model.G(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j))
                    + Model.B(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j)));
        }
        return Sol.V.coeff(k) * val;
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the reactive power at a bus
     */
    double Solver_NRpolar::Q(uint k) {
        double val = 0.0;
        for (uint j = 0; j < Model.buses.size(); j++) {
            val += Sol.V.coeff(j)
                    *(Model.G(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j))
                    - Model.B(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j)));
        }
        return Sol.V.coeff(k) * val;
    }

    /*
     * This function corects the PV buses that exeed the reative power limit
     */
    void Solver_NRpolar::correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X) {

        /*Find the PV buses that violate the reative power limit*/
        vector<int> lst;
        uint k;
        for (uint i = npq; i < npq + npv; i++) {
            k = PQPV[i];
            if (Model.buses[k].max_q != 0.0) {
                if (Sol.Q(k) >= Model.buses[k].max_q) {
                    Sol.Q(k) = Model.buses[k].max_q; //truncate Q to the limit
                    lst.push_back(k); //add the PV bus to the list to be treated as a PQ bus
                    std::cout << "PV to PQ: " << Model.buses[k].Name << std::endl;
                }
            } else {
                std::cout << "Probably invalid reactie power value at bus " << Model.buses[k].Name << std::endl;
            }
        }

        /*Change the lists and arrays size to accomodate the new situation*/
        npq += lst.size();
        npv -= lst.size();
        uint npqpvpq = 2 * npq + npv; //size of the arrays

        //Resize the linear system, since if npq and npv vary their size vary as well
        J = mat(npqpvpq, npqpvpq);
        K = vec(npqpvpq);
        X = vec(npqpvpq);

        //add the PV buses indices from lst to the PQ list
        LastPQ.insert(LastPQ.end(), lst.begin(), lst.end());

        //Remove the same lst indices fromthe PV list
        for (uint k : lst)
            LastPV.erase(std::remove(LastPV.begin(), LastPV.end(), k), LastPV.end());

        PQPV.clear();
        PQPV.reserve(LastPQ.size() + LastPV.size()); // preallocate memory
        PQPV.insert(PQPV.end(), LastPQ.begin(), LastPQ.end());
        PQPV.insert(PQPV.end(), LastPV.begin(), LastPV.end());
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the jacobian of the circuit
     */
    void Solver_NRpolar::Jacobian(mat &J, uint npq, uint npv) {
        //matrix(rows, cols)
        uint npqpv = npq + npv;
        double val;
        uint k, j;
        uint da, db;

        J.setZero();

        //J1 
        for (uint a = 0; a < npqpv; a++) { //rows
            k = PQPV[a];
            //diagonal
            J(a, a) = -Q(k) - Model.B(k, k) * Sol.V.coeff(k) * Sol.V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npqpv; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = Sol.V.coeff(k) * Sol.V.coeff(j)
                            *(Model.G(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j))
                            - Model.B(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j)));
                    //if (val != 0.0)
                    J(a, b) = val;
                }
            }
        }

        //J2
        da = 0;
        db = npqpv;
        for (uint a = 0; a < npqpv; a++) { //rows
            k = PQPV[a];
            //diagonal
            //std::cout << "J2D:" << (a + da) << "," << (a + db) << std::endl;
            if (a < npq)
                J(a + da, a + db) = P(k) + Model.G(k, k) * Sol.V.coeff(k) * Sol.V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npq; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = Sol.V.coeff(k) * Sol.V.coeff(j)
                            *(Model.G(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j))
                            + Model.B(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j)));
                    //if (val != 0.0)
                    //std::cout << "J2ij:" << (a + da) << "," << (b + db) << std::endl;
                    J(a + da, b + db) = val;
                }
            }
        }


        //J3
        da = npqpv;
        db = 0;
        for (uint a = 0; a < npq; a++) { //rows
            k = PQPV[a];
            //diagonal
            //std::cout << "J3:" << (a + da) << "," << (a + db) << std::endl;
            J(a + da, a + db) = P(k) - Model.G(k, k) * Sol.V.coeff(k) * Sol.V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npqpv; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = Sol.V.coeff(k) * Sol.V.coeff(j)
                            *(Model.G(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j))
                            + Model.B(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j)));
                    //if (val != 0.0)
                    //std::cout << "J3:" << (a + da) << "," << (b + db) << std::endl;
                    J(a + da, b + db) = -val;
                }
            }
        }

        //J4
        da = npqpv;
        db = npqpv;
        for (uint a = 0; a < npq; a++) { //rows
            k = PQPV[a];
            //diagonal
            //std::cout << "J4:" << (a + da) << "," << (a + db) << std::endl;
            J(a + da, a + db) = Q(k) - Model.B(k, k) * Sol.V.coeff(k) * Sol.V.coeff(k);

            //non diagonal elements
            for (uint b = 0; b < npq; b++) {
                if (b != a) {
                    j = PQPV[b];
                    val = Sol.V.coeff(k) * Sol.V.coeff(j)
                            *(Model.G(k, j) * sin(Sol.D.coeff(k) - Sol.D.coeff(j))
                            - Model.B(k, j) * cos(Sol.D.coeff(k) - Sol.D.coeff(j)));
                    if (val != 0.0) {
                        //std::cout << "J4:" << (a + da) << "," << (b + db) << std::endl;
                        J(a + da, b + db) = val;
                    }
                }
            }
        }


    }

    /*//////////////////////////////////////////////////////////////////////////
     * Calculate the power increments
     */
    void Solver_NRpolar::get_power_inc(vec& PQinc, uint npq, uint npv) {

        uint npqpv = npq + npv;
        uint k;
        PQinc.setZero();

        for (uint a = 0; a < npqpv; a++) {
            //For PQ and PV buses; calculate incP
            k = PQPV[a];
            PQinc(a) = Pesp.coeff(k) - P(k);

            if (a < npq) //only for PQ buses, calculate incQ
                PQinc(a + npqpv) = Qesp.coeff(k) - Q(k);
        }

    }

    /*//////////////////////////////////////////////////////////////////////////
     * Check the convergence
     */
    bool Solver_NRpolar::converged(vec PQinc, uint npqpvpq) {

        for (uint k = 0; k < npqpvpq; k++)
            if (abs(PQinc.coeff(k)) > EPS)
                return false;

        return true;
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Check the convergence
     */
    void Solver_NRpolar::update_solution(vec X, uint npq, uint npv) {

        uint npqpv = npq + npv;
        uint k;

        for (uint a = 0; a < npqpv; a++) {
            k = PQPV[a];
            Sol.D(k) += X.coeff(a);

            if (a < npq)
                Sol.V(k) = Sol.V.coeff(k) *(1.0 + X.coeff(a + npqpv));
        }

        //Correct for PV buses
        for (uint i = npq; i < npq + npv; i++) {
            k = PQPV[i];
            cx_double v = Sol.Vcx(k);
            v *= Model.buses[k].v_set_point / abs(v);
            Sol.V(k) = abs(v);
            Sol.D(k) = arg(v);
        }
    }

    /*//////////////////////////////////////////////////////////////////////////
     * Run the solve process
     */
    Solver_State Solver_NRpolar::solve() {

        uint npq = Model.PQ_list.size();
        uint npv = Model.PV_list.size();
        uint npqpvpq = 2 * npq + npv; //size of the arrays

        //System : J*X = K
        mat J(npqpvpq, npqpvpq);
        vec K(npqpvpq);
        vec X(npqpvpq);

        bool conv;
        Iterations = 0;

        //first calculations
        get_power_inc(K, npq, npv);
        conv = converged(K, npqpvpq);

        //std::cout << "Converged: " << conv << std::endl;

        while (!conv && Iterations < Max_Iter) {
            //std::cout << "-----------------------Iter: " << Iterations << std::endl;
            //std::cout << "K:\n" << K << std::endl;
            //Calculate the jacobian
            Jacobian(J, npq, npv);
            //std::cout << "J:\n" << J << std::endl;

            Eigen::FullPivLU<mat>lu(J); //Full pivot LU
            X = lu.solve(K);
            //std::cout << "X:\n" << X << std::endl;

            //upgrade the solution
            update_solution(X, npq, npv);

            calculate_Q(npq, npv); //Calculate the reactive power for the PV buses

            correct_PVbuses_violating_Q(npq, npv, J, K, X); //pass PV buses to PQ

            //Calculate the increment of power for the new iteration
            get_power_inc(K, npq, npv);

            //Check convergency
            conv = converged(K, npqpvpq);

            //Sol.print("Soltution:");

            Iterations++;
        }

        if (!conv || Iterations == Max_Iter)
            return Solver_State::Not_Converged;
        else {
            //std::cout << "Converged in " << Iterations << " iterations." << std::endl;
            calculate_slack_power();
            Model.set_solution(Sol.get_cx());
            //Sol.print("Final Soltution:");
            return Solver_State::Converged;
        }
    }
    
    
    
    /*This function updates the solver solution object power values using the
     * circuit's own solution power values. this is specially usefull when updating
     * the circuit power values while keeping the previous voltage solution
     */
    void Solver_NRpolar::update_solution_power_from_circuit(){    
        Sol.P = Model.get_initial_solution().P;
        Sol.Q = Model.get_initial_solution().Q;
        Pesp = Sol.P;
        Qesp = Sol.Q;
    }

}