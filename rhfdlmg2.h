#ifndef _rhfdlmg2
#define _rhfdlmg2

#include <vector>

#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/oeprop.h"
#include "psi4/libscf_solver/rhf.h"
#include "psi4/libscf_solver/hf.h"


namespace psi{

namespace scf{

class RHFDLMG2: public RHF{
public:
  RHFDLMG2(SharedWavefunction ref_wfn,Options &options, std::shared_ptr<PSIO> psio);
  virtual ~RHFDLMG2();
protected:

  void calc_eps_symm();
  void read_dens();
  std::string dens_file_;
  void read_ca();
  std::string ca_file_;
  void calc_sigsum();
  void compute_esp_ao();
  void compute_esp_old();
  void compute_esp();
  int esp_threaded;
  virtual double compute_E();
  virtual void form_G();
  virtual void form_H();
  virtual double compute_energy();
  double calc_sg_pot();
  double compute_oneE();
  void calc_eps();
	void calc_V_vac();
  void calc_BC_vals();
  void calc_BC_vals_vac();
  std::shared_ptr<OEProp> oe_;
	void calc_V_sol();
	void calc_V_sgao2();
	void calc_V_sgao_P();
	void calc_V_sgao();
	void calc_V_eff();
	void calc_V_eff2();
	void calc_V_eff3();
	void calc_dens_grid();
  void calc_dens_grid_old();
	void new_dens_grid();
	double calc_sg_rep();
	void coarse_rho();
	double num_int();
	double nearest_dist(const double gp[]);
	double nearest_comp(double q,int indx);
  double phiao_contrib(double gp[3]);
	double numint2(int ijk0[3],int ijkf[3],double h[3]);
  SharedMatrix V_eff_;
  SharedMatrix V_sgao_;
  SharedMatrix V_ion_;
  SharedMatrix myDa_;


    /// The number of atoms in the current molecule.
  int natom_;
  void calc_V_bc(std::vector<double> &pot_vac_,double eps);
int nx_,ny_,nz_;
int dens_algo;
int dens_only;
std::vector <std::vector<double> > cart_grid_;//(nx_*ny_*nz_,std::vector<double>(3));
std::vector<double> sigsum_;
double g_step_,sigma_;

	// Vector to hold density at grid points
	std::vector<double> gmat_;
	std::vector<double> nmat_;
	std::vector<double> emat_;
	double centerx_,centery_,centerz_;
	int npts;
	double rho0_;//=0.00078;
  int test_elec_flag_;
    //double beta2_=2.6;
  double beta2_;
  int dlmg_tol;
   double eps_inf_;
	std::vector<double> pot_sol_;
	std::vector<double> pot_vac_;
  std::vector<double> bound_;
	std::vector<double> residual_;
	std::vector<double> error_;
	std::vector<double> eps_half_;
	std::vector<double> eps_full_;
	std::vector<double> eps_half_vac;
	std::vector<double> eps_full_vac;
  std::shared_ptr<BasisSet> basis_;
  int narray[3];
  double darray[3];
	
};

}} // Namespaces
#endif
