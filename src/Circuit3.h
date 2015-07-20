/*
 * File:   Circuit3.h
 * Author: Santiago Peñate Vera
 *
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <vector>
#include <math.h>

#include "fpotencia_libs.h"

#include "Bus.h"
#include "Load.h"
#include "Line.h"
#include "Transformer.h"
#include "Generator.h"
#include "Shunt.h"
#include "ExternalGrid.h"
#include "Solution.h"

#ifndef CIRCUIT3_H
#define	CIRCUIT3_H

using namespace Eigen;
//using namespace arma;
using namespace std;


namespace fPotencia {

    class Circuit3 {
    public:
        Circuit3(string name);

        virtual ~Circuit3();

        /*
         * Function to assemble the circuit admittance matrix and initial solutions
         * from the specifyed data.
         */
        void Compile();
        
        /*
         * Function that accesses correctly Y at the bus i,j at the phase a,b
         * 
         * i: bus i
         * j: bus j
         * a: phase {0, 1, 2}
         * b: phase {0, 1, 2}
         */
        cx_double Y(int i, int j, int a, int b);

        /*
         * Function that accesses correctly to the real part of Y
         *  the bus i,j at the phase a,b
         * 
         * i: bus i
         * j: bus j
         * a: phase {0, 1, 2}
         * b: phase {0, 1, 2}
         */
        double G(int i, int j, int a, int b);

        /*
         * Function that accesses correctly to the imaginary part of Y
         *  the bus i,j at the phase a,b
         * 
         * i: bus i
         * j: bus j
         * a: phase {0, 1, 2}
         * b: phase {0, 1, 2}
         */
        double B(int i, int j, int a, int b);


        /*
         * Function that sets a value to Y at the bus i,j at the phase a,b
         * 
         * i: bus i
         * j: bus j
         * a: phase {0, 1, 2}
         * b: phase {0, 1, 2}
         * val: value to set
         */
        void setY(int i, int j, int a, int b, cx_double val);
        
        
        string Name;
        
        std::vector<Bus> buses;

        std::vector<Load> loads;

        std::vector<Generator> generators;

        std::vector<Line> lines;

        std::vector<Transformer> transformers;

        std::vector<Shunt> shunts;

        std::vector<ExternalGrid> external_grids;

        std::vector<unsigned int> PQ_list; //consumption buses indices list

        std::vector<unsigned int> PV_list; //generation buses indices list

        std::vector<unsigned int> VD_list; //Slack buses indices list

    private:

        sp_cx_mat Ymat;
    };

}

#endif	/* CIRCUIT_H */