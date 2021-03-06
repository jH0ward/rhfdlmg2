/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef _myelst_h
#define _myelst_h

#include "psi4/libmints/electrostatic.h"

namespace psi {

class ElectrostaticInt;
class BasisSet;
class GaussianShell;
class ObaraSaikaTwoCenterRecursion;
class OneBodyAOInt;
class PotentialInt;
class IntegralFactory;
class SphericalTransform;
class Vector3;

/*! \ingroup MINTS
 *  \class PotentialInt
 *  \brief Computes potential integrals.
 * Use an IntegralFactory to create this object.
 */
class MyElstInt : public ElectrostaticInt
{
    virtual void compute_pair(const GaussianShell&, const GaussianShell&)
    {}

public:
    /// Constructor
    MyElstInt(std::vector<SphericalTransform>&, std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>, int deriv=0);
    //MyElstInt(vector<SphericalTransform>&, BasisSet, BasisSet, int deriv=0);
    ~MyElstInt();

    // Intel C++ 12 thinks we're trying to overload the "void compute_shell(int, int)" and warns us about it.
    // The following line is to shut it up.
    #pragma warning disable 1125
    /// Computes integrals between two shells.
    virtual void compute_shell(int, int, const Vector3&);
    /// Computes integrals between two shells.
    virtual void compute_pair(const GaussianShell&, const GaussianShell&, const Vector3&);

#if __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#endif
    /// Computes integrals and stores in result.
    virtual void compute(SharedMatrix& result, const Vector3&);
#if __clang__
#pragma clang diagnostic pop
#endif
    /// Does the method provide first derivatives?
    bool has_deriv1() { return false; }

    static SharedVector nuclear_contribution(std::shared_ptr<Molecule> mol);
};

}

#endif
