/* 
 * File:   Bus.h
 * Author: Santiago Peñate Vera
 *
 * Created on 6 de agosto de 2014, 9:54
 * 
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
//#include "armadillo"
//using namespace arma;

#include "fpotencia_libs.h"

using namespace std;

namespace fPotencia {
#ifndef BUS_H
#define	BUS_H

    class Bus {
    public:
        Bus(string name, const BusType& type, double Bus_Nominal_Voltage);
        virtual ~Bus();

        //Properties
        int index = -1;

        string Name;

        BusType Type;

        double nominal_voltage;

        double max_voltage;

        double min_voltage;

        cx_double voltage;

        cx_double voltage_pu;

        cx_double power;

        cx_double connected_power = cx_double(0, 0); //To store the sum of load and generation

        /*-----Only for PV buses--------------------------------------*/

        double min_q = 0; //minimum reactive power per unit (changed in the circuit class)

        double max_q = 0; //maximum reactive power per unit (changed in the circuit class)

        double v_set_point = 1.0; //only used if the bus is PV
        /*------------------------------------------------------------*/

        void print();

    private:

    };

#endif	/* BUS_H */

}