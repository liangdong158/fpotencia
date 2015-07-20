/* 
 * File:   ExternalGrid.h
 * Author: Santiago Peñate Vera
 *
 * Created on 8 de agosto de 2014, 14:45
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

#ifndef EXTERNALGRID_H
#define	EXTERNALGRID_H

	class ExternalGrid {
		public:
		ExternalGrid(string name, int connection_bus);
		virtual ~ExternalGrid();

		/*Properties*/
		string Name;

		int bus;

		cx_double power;

		private:

	};

#endif	/* EXTERNALGRID_H */

}