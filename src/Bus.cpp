/* 
 * File:   Bus.cpp
 * Author: Santiago Peñate Vera
 * 
 * Created on 6 de agosto de 2014, 9:55
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Bus.h"

namespace fPotencia {

    /*
     * Bus class constructor
     */
    Bus::Bus(string name, const BusType& type, double Bus_Nominal_Voltage) {
        Name = name;
        Type = type;
        nominal_voltage = Bus_Nominal_Voltage;
        max_voltage = 1.05 * nominal_voltage;
        min_voltage = 0.95 * nominal_voltage;
    }

    /*
     * Bus class destructor
     */
    Bus::~Bus() {
    }

    /*
     * Print the bus results
     */
    void Bus::print() {
        cout << Name + " -> " + BusType_name[Type] << endl;
        cout << "\tPower: " << power << endl;

        cout << "\tVoltage: " << voltage << "\tVoltage p.u.: " << voltage_pu << endl;
    }

}