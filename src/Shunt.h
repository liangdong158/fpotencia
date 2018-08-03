/* 
 * File:   Shunt.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:06
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

//#include "armadillo"
#include "fpotencia_libs.h"
#include "Solution.h"

//using namespace arma;
using namespace std;

namespace fPotencia {
#ifndef SHUNT_H
#define	SHUNT_H

	class Shunt {
		public:                
                Shunt(string name, int bus, double R, double X);
                
		virtual ~Shunt();

		void get_element_Y(int n, sp_cx_mat &Yret);

		void calculate_current(cx_solution sol);

		void print();

		/*************************************************************************
		 * Properties
		 *************************************************************************/
		string Name;

		int bus1 = 0;

		cx_double impedance;

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double current;

		cx_double power;

		cx_double power_losses;

		private:

		/*************************************************************************
		 * Calculated variables: Results
		 *************************************************************************/

		cx_double Y_element;

	};

#endif	/* SHUNT_H */

}