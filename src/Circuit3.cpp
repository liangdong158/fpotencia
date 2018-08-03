/*
 * File:   Circuit3.cpp
 * Author: Santiago Peñate Vera
 *
 * Copyright (C) 2014 Santiago Peñate Vera
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "Circuit3.h"

namespace fPotencia {

    Circuit3::Circuit3(string name) {
        Name = name;
    }

    Circuit3::~Circuit3() {
    }

    void Circuit3::Compile() {
    }

    cx_double Circuit3::Y(int i,int j, int a, int b) {
        return Ymat.coeff(i * 3 + a, j * 3 + b);
    }

    void Circuit3::setY(int i, int j, int a, int b, cx_double val) {
        Ymat.coeffRef(i * 3 + a, j * 3 + b) = val;
    }

    double Circuit3::G(int i, int j, int a, int b) {
        return Ymat.coeff(i * 3 + a, j * 3 + b).real();
    }

    double Circuit3::B(int i, int j, int a, int b) {
        return Ymat.coeff(i * 3 + a, j * 3 + b).imag();
    }

}