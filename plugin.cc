/*
 *@BEGIN LICENSE
 *
 * rhf_dlmg2 by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libscf_solver/rhf.h>
#include <psi4/libscf_solver/hf.h>
#include <psi4/libciomr/libciomr.h>
#include "rhfdlmg2.h"
#include <iostream>

namespace psi{ namespace scf {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "RHF_DLMG2"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
        options.add_int("eflag",0);
        options.add_int("nx",161);
        options.add_int("ny",161);
        options.add_int("nz",161);
        options.add_double("g_step",0.075);
        options.add_double("sigma",0.5);
        options.add_double("eps",1.0);
        //options.add_double("eps_sol",78.34);
        options.add_double("eps_sol",78.54);
        options.add_double("rho0",0.00055);
        options.add_double("beta",1.6);
        options.add_int("dens_algo",0);
        options.add_int("dens_only",0);
        options.add_int("esp_threaded",1);
        options.add_int("dlmg_tol",7);
    }

    return true;
}

extern "C"
SharedWavefunction rhf_dlmg2(SharedWavefunction ref_wfn, Options& options)
{
    // need to use dynamic cast
    //shared_ptr<psi::scf::RHF> rhf_shptr = reinterpret_cast<*RHF>(ref_wfn.get());
    std::shared_ptr<RHF> rhf_shptr( dynamic_cast<RHF*>(ref_wfn.get()));
    int print = options.get_int("PRINT");
    std::shared_ptr<PSIO> psio;
    std::cout << "Eps is " << options.get_double("eps") << "\n";
    std::cout << "eps_sol is " << options.get_double("eps_sol") << "\n";
   // exit(0);
    //std::cout << ref_wfn->functional() << "\n";

    /* Your code goes here */
    //boost::shared_ptr<Wavefunction> scf(new RHFDLMG2(ref_wfn,options,psio));
    std::shared_ptr<Wavefunction> scf(new RHFDLMG2(ref_wfn,options,psio));
    scf->compute_energy();
    return scf;
}

}} // End namespaces

