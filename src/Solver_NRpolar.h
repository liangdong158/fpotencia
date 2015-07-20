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
#include <cmath>
#include "Circuit.h"
#include "Solution.h"

//using namespace arma;
using namespace std;

namespace fPotencia {

#ifndef SOLVER_NRPOLAR_H
#define	SOLVER_NRPOLAR_H

    /*
     * This class implements the Nerwton Raphson method
     * in polar coordinates to solve the circuit.
     */
    class Solver_NRpolar {
    public:

        Solver_NRpolar(Circuit model);

        Solver_NRpolar(Circuit model, solution sol_);

        virtual ~Solver_NRpolar();

        /*Properties*/
        Circuit Model;

        double EPS = 1e-6;

        int Iterations = 0;

        int Max_Iter = 100;

        Solver_State solve(); //Solves the grid

        void update_solution_power_from_circuit();
        
    private:

        vector<int> BUSES;

        vector<int> PQPV;

        vector<int> LastPQ;
        
        vector<int> LastPV;

        vec Pesp;

        vec Qesp;

        solution Sol;

        void Jacobian(mat &J, uint npq, uint npv); //calculate the jacobian, J is passed by refference

        void get_power_inc(vec &PQinc, uint npq, uint npv); //PQinc is passed by refference

        void calculate_Q(uint npq, uint npv); //calculate the reative power at the PV buses

        double Q(uint k);

        double P(uint k);

        bool converged(vec PQinc, uint npqpvpq); //check if the solution converged

        void update_solution(vec X, uint npq, uint npv);

        void calculate_slack_power(); //calculate the slack bus power        

        bool checks(); //check the solvability with this method

        void correct_PVbuses_violating_Q(uint &npq, uint &npv, mat &J, vec &K, vec &X);
        
    };

#endif	/* SOLVER_NRPOLAR_H */

}