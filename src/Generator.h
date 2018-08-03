/* 
 * File:   Generator.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 10:05
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "fpotencia_libs.h"

using namespace std;

namespace fPotencia {
#ifndef GENERATOR_H
#define	GENERATOR_H

	class Generator {
		public:
		Generator(string name, uint connection_bus, double P, double Q);

		Generator(string name, uint connection_bus, double P, double Vset, double Qmin, double Qmax, bool Vset_per_unit);

		virtual ~Generator();


		/*Properties*/
		string Name;

		uint bus;

		cx_double power = cx_double(0.0, 0.0);

		double min_Q;

		double max_Q;

		double voltage_set_point;

		bool voltage_controlled = false;
                
                bool Vset_in_per_unit = false;

		private:

	};

#endif	/* GENERATOR_H */

}