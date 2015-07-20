/* 
 * File:   Load.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

//#include "armadillo"
#include "fpotencia_libs.h"

//using namespace arma;
using namespace std;

namespace fPotencia {
#ifndef LOAD_H
#define	LOAD_H

	class Load {
		public:
		Load(string name, int connection_bus, double P, double Q);
		virtual ~Load();

		/*Properties*/
		string Name;

		int bus = -1;

		cx_double power;


		private:

	};

#endif	/* LOAD_H */

}