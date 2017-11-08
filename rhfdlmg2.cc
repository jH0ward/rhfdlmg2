#include <typeinfo>
#include "rhfdlmg2.h"
#include <psi4/psi4-dec.h>
#include <psi4/liboptions/liboptions.h>
#include <fstream>
#include <iostream>
#include <psi4/libpsio/psio.hpp>
#include <psi4/libciomr/libciomr.h>
#include <psi4/libscf_solver/hf.h>
#include <psi4/libscf_solver/rhf.h>
#include <psi4/libmints/oeprop.h>
#include <psi4/psifiles.h>
#include <psi4/libpsi4util/PsiOutStream.h>
#include <psi4/libmints/matrix.h>
#include <psi4/libpsi4util/process.h>
#include <psi4/libmints/molecule.h>
#include <psi4/libmints/basisset.h>
#include <psi4/libmints/mintshelper.h>
#include <psi4/libmints/petitelist.h>
#include <psi4/libmints/psimath.h>
#include <psi4/libqt/qt.h>
#include <psi4/libmints/integral.h>
#include <omp.h>
#include <chrono>
#include <stdio.h>
#include <iomanip>
#include "myelst.h"


namespace psi{ namespace scf {

//#define DL_MG_INIT dl_mg_wrap_mp_dl_mg_init_wrap_
//#define DL_MG_SOLVER dl_mg_wrap_mp_dl_mg_solver_wrap_
//#define DL_MG_FREE dl_mg_mp_dl_mg_free_
//
#define DL_MG_INIT __dl_mg_wrap_MOD_dl_mg_init_wrap
#define DL_MG_SOLVER __dl_mg_wrap_MOD_dl_mg_solver_wrap
#define DL_MG_FREE __dl_mg_MOD_dl_mg_free

extern "C"
{
void DL_MG_INIT(int narray[3], double darray[3], int bc[3], int gstart[3],
                int gend[3], int * maxcycles, int * report_unit,
                const char * filename, int * filenamelen,
                int * pbe_mode, int * is_damping, double * error_tol,
                double * defect_tol, int * max_iters,
                int * fullagg_size, int * fullagg_level, int * minblocksize,
                int * block_size);


void DL_MG_SOLVER(int sizes[3], double * pot, double * rho, int * order, double * bound, int * report_unit, double * eps_full = NULL, double * eps_half = NULL);

void DL_MG_FREE();
}

    //std::shared_ptr<RHF> rhf_shptr( dynamic_cast<RHF*>(ref_wfn.get()));

//RHFDLMG2::RHFDLMG2(ref_wfn, Options
RHFDLMG2::RHFDLMG2(SharedWavefunction ref_wfn, Options
                  &options,std::shared_ptr<PSIO> psio): RHF(ref_wfn,
                  dynamic_cast<RHF*>(ref_wfn.get())->functional(), options,psio)
                  //&options,std::shared_ptr<PSIO> psio): RHF(ref_wfn,options,psio)
{
  test_elec_flag_=options.get_int("eflag");
  nx_=options.get_int("nx");
  ny_=options.get_int("ny");
  nz_=options.get_int("nz");
  dlmg_tol=options.get_int("dlmg_tol");
  dens_algo=options.get_int("dens_algo");
  dens_only=options.get_int("dens_only");
  esp_threaded=options.get_int("esp_threaded");
  std::cout << "Check 1\n";
  narray[0]=nx_;
  narray[1]=ny_;
  narray[2]=nz_;
  rho0_=options.get_double("rho0");
  beta2_=2*options.get_double("beta");
  std::cout << "beta2 = " << beta2_ << "\n";
  std::cout << "Nx = " << nx_ << "\n";
  std::cout << "Ny = " << ny_ << "\n";
  std::cout << "Nz = " << nz_ << "\n";
  g_step_=options.get_double("g_step");
  darray[0]=g_step_;
  darray[1]=g_step_;
  darray[2]=g_step_;
  std::cout << "Check 2\n";
  sigma_=options.get_double("sigma");
  eps_inf_=options.get_double("eps");
  outfile->Printf("eps = %15.10f\n",eps_inf_);
  basis_=ref_wfn->basisset();	
  V_eff_ = SharedMatrix (new Matrix("V eff",nirrep_,nsopi_,nsopi_));
  V_sgao_=SharedMatrix(new Matrix("V sgao",nirrep_,nsopi_,nsopi_));

  //std::cout << " about to guess\n";
  std::cout << "instead of guessing, I'm going to set ca and da\n";
  Ca_->copy(ref_wfn->Ca());
  Da_->copy(ref_wfn->Da());
  H_->copy(ref_wfn->H());
  std::cout << "I Just set Da and Ca\n";
  Ca_->print_out();
  //RHF::guess();
  //std::cout << " just guessed\n";
  centerx_=g_step_*(nx_-1)/2.0;
  centery_=g_step_*(ny_-1)/2.0;
  centerz_=g_step_*(nz_-1)/2.0;
	pot_sol_.resize(nx_*ny_*nz_);
	bound_.resize(nx_*ny_*nz_);
  gmat_.resize(nx_*ny_*nz_);
  nmat_.resize(nx_*ny_*nz_);
  // We don't have to pass eps_half if it's vacuum so
  // don't allocate it
  if (eps_inf_>1.0){
    eps_half_.resize(nx_*ny_*nz_*3);
    eps_full_.resize(nx_*ny_*nz_);
    calc_eps_symm();
  }
  else{
    eps_half_.resize(nx_*ny_*nz_*3);
    eps_full_.resize(nx_*ny_*nz_);
    for (int i=0;i<nx_*ny_*nz_*3;i++){
      eps_half_[i]=1.0;
    }
    for (int i=0;i<nx_*ny_*nz_;i++){
      eps_full_[i]=1.0;
    }
  }

  std::cout << "Box is centered at " << centerx_ << ", " << centery_ << ", "
	          << centerz_ << std::endl;
  //std::cout << "printing density before reading attempt\n";
  Da_->print_out();
  //std::cout << "Going to create a mintshelper to calculate integrals\n";
  //std::shared_ptr<MintsHelper> mints (new MintsHelper(basisset_, options_, 0));
  //mints->one_electron_integrals();
  // Look for density file
  //std::string dens_file_;
  //exit(0);
}
RHFDLMG2::~RHFDLMG2()
{
  std::cout << "Entering destructor\n";
  //std::vector<double>().swap(pot_sol_);
  //std::vector<double>().swap(bound_);
  //std::vector<double>().swap(gmat_);
  //std::vector<double>().swap(nmat_);
  //std::vector<double>().swap(eps_half_);
  //std::vector<double>().swap(eps_full_);


	//pot_sol_.resize(0);
	//bound_.resize(0);
  //gmat_.resize(0);
  //nmat_.resize(0);
  //eps_half_.resize(0);
  //eps_full_.resize(0);


}

void
RHFDLMG2::read_dens()
{
  std::ifstream densin(dens_file_);
  std::string myline;
  double read_val;
  for (int i=0;i<nirrep_;i++){
    std::getline(densin,myline);
    int last_col=0;
    while(last_col<nsopi_[i]){
      for (int j=0;j<3;j++){
        std::getline(densin,myline);
      }
      // Read and throw away row_id
      // Loop over nrows in irrep
      for (int r=0;r<nsopi_[i];r++){
        densin >> read_val;
        // Read in 5 values
        for (int j=last_col;j<std::min(last_col+5,nsopi_[i]);j++){
          densin >> read_val;
          Da_->set(i,r,j,read_val);
        }
        std::getline(densin,myline);
      }
      last_col+=5;
    }
    for (int j=0;j<3;j++){
      std::getline(densin,myline);
    }
  }
}


void
RHFDLMG2::read_ca()
{
  std::ifstream casin(ca_file_);
  std::string myline;
  double read_val;
  for (int i=0;i<nirrep_;i++){
    std::getline(casin,myline);
        //if (i>0){
          //std::cout << "myline = " << myline << "\n";
        //}
    int last_col=0;
    //while(last_col<nsopi_[i]+1){
    while(last_col<nsopi_[i]){
      //std::cout << "last_col is " << last_col << "\n";
      for (int j=0;j<3;j++){
        std::getline(casin,myline);
        //if (i>0){
          //std::cout << "in j loop, myline  = " << myline << "\n";
       // }
      }
      // Read and throw away row_id
      // Loop over nrows in irrep
      for (int r=0;r<nsopi_[i];r++){
        casin >> read_val;
        //if (i>0){
        //  std::cout << "read_val from r loop = " << read_val << "\n";
        //}
        // Read in 5 values
        for (int j=last_col;j<std::min(last_col+5,nsopi_[i]);j++){
          casin >> read_val;
          //if (i>0){
          //}
          Ca_->set(i,r,j,read_val);
        }
        //std::cout << "last read val = " << read_val << "\n";
        std::getline(casin,myline);
      }
      last_col+=5;
    }
    //std::cout << "Entering throwaway 3 lines on i = " << i << "\n";
    for (int j=0;j<3;j++){
      std::getline(casin,myline);
     // std::cout << myline << "\n";
    }
  }
}

double
RHFDLMG2::compute_energy()
{
  std::cout << "about to compute energy\n";
  double energy =  RHF::compute_energy();
  std::cout << " done with compute energy\n";
  std::cout << " swapping vectors\n";
  std::vector<double>().swap(pot_sol_);
  std::vector<double>().swap(bound_);
  std::vector<double>().swap(gmat_);
  std::vector<double>().swap(nmat_);
  std::vector<double>().swap(eps_half_);
  std::vector<double>().swap(eps_full_);
  return energy;
}

double
RHFDLMG2::compute_E()
{
  outfile->Printf("I'm in compute_E; gonna call RHF compute_E to start\n");
	double Etot=RHF::compute_E();
  outfile->Printf("Etot from RHF: %24.16f\n",Etot);
  
  double sg_pot = calc_sg_pot();
  Etot+=sg_pot;
  outfile->Printf("Interaction of cores with Veff = %24.16f\n",sg_pot); 

  double sgrep=calc_sg_rep();
  Etot+=sgrep;
  outfile->Printf("Returning total E = %24.16f\n",Etot);

	return Etot;
}

void RHFDLMG2::form_H()
{
  outfile->Printf("I'm in my own form_H\n");
  //RHF::form_H();
  outfile->Printf("Just called RHF::form_H()\nNow printing H_");
  H_->print_out();
  calc_V_sgao_P();
  H_->add(V_sgao_);
  printf("Leaving my own form_H\n");
}

void RHFDLMG2::form_G()
{
    outfile->Printf("I'm in my own form_G\n");
    outfile->Printf("Calling calc_V_sol from form_G\n");
    std::cout << "Iteration = " << iteration_ << "\n";
    if (iteration_==0){
      if (eps_inf_>1.0){
        ca_file_="sol_cas";
      }
      else{
        ca_file_="vac_cas";
      }
      std::ifstream casin(ca_file_.c_str());
      if (casin){
        std::cout << "cas file exists\n";
        // call read_density function
        read_ca();
        std::cout << "Printing out cas after read\n";
        Ca_->print_out();
        //exit(0);
      }
      else{
        std::cout << "cas file does not exist\n";
      }
      // don't, ok do
      // repeat for density
      if (eps_inf_>1.0){
        dens_file_="sol_dens";
      }
      else{
        dens_file_="vac_dens";
      }
      std::ifstream densin(dens_file_.c_str());
      if (densin){
        std::cout << "dens file exists\n";
        // call read_density function
        read_dens();
        std::cout << "Printing out density after read\n";
        Da_->print_out();
      }
      else{
        std::cout << "density file does not exist\n";
      }
    }
    RHF::form_G();
    G_->subtract(J_);
    //G_->print_out();
    calc_V_sol();
    outfile->Printf("Calling calc_V_eff from form_G\n");
    calc_V_eff3();
    G_->add(V_eff_);
    outfile->Printf("Leaving my form_G\n");
    outfile->Printf("But I'm going to print out the coefficients first\n");
    Ca_->print_out();
}

void
RHFDLMG2::calc_dens_grid()
{
  std::shared_ptr<Molecule> molecule = Process::environment.molecule();
  Matrix geom=molecule->geometry();
  int natom=molecule->natom();
  molecule->print_in_bohr();
  Da_->print_out();
  int npts=nx_*ny_*nz_;
  std::cout << "about to set gmat and nmat\n";
  for (int i=0;i<npts;i++){
    gmat_[i]=0.;
    nmat_[i]=0.;
  }
  std::cout << "getting nao\n";
  int nao = basisset_->nao();
  std::cout << "counting nso\n";
  int nso=0;
  for (int h=0;h<nirrep_;h++) nso+=nsopi_[h];

	double tot=0.;
  std::cout << "nso = " << nso << "\n";
  double ndens=0;
  double rrmin=100;
  int savei,savej,savek;
  std::cout << "Beginning loop over x\n";
  double bfprod,xt,yt,zt;
  double Z[natom];
  double Zreal[natom];

  for (int a=0;a<natom;a++){
    Z[a]=molecule->Z(a);
    Zreal[a]=molecule->Z(a);
  }


  std::chrono::time_point<std::chrono::system_clock> dens_start, dens_end;
  dens_start = std::chrono::system_clock::now();
  double phiao[nao];
  double phiso[nso];
  int *col_offset = new int[nirrep_];
  MintsHelper helper(basis_, RHF::options_, 0);
  SharedMatrix aotoso = helper.petite_list(true)->aotoso();
  col_offset[0] = 0;
  for(int h=1; h < nirrep_; h++){
    col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);
  }
  double **u = block_matrix(nao, nso);
  for(int h=0; h < nirrep_; h++){
    for(int j=0; j < aotoso->coldim(h); j++){
      for(int ii=0; ii < nao; ii++){
        u[ii][j+col_offset[h]] = aotoso->get(h, ii, j);
      }
    }
  }
  #pragma omp parallel for reduction(+:tot,ndens) private(phiao,phiso)
  for (int i=0;i<nx_;i++){
    double xt=g_step_*i-centerx_;
    for (int j=0;j<ny_;j++){
      double yt=g_step_*j-centery_;
    for (int k=0;k<nz_;k++){
        double zt=g_step_*k-centerz_;
        int g=i+nx_*j+nx_*ny_*k;
        basisset_->compute_phi(phiao,xt,yt,zt);
        C_DGEMV('t',nao, nso, 1.0, &(u[0][0]), nso, &(phiao[0]), 
                1, 0.0, &(phiso[0]), 1);
        for (int h=0;h<nirrep_;h++){
          for(int u=0; u < nsopi_[h]; u++){
            for(int v=0; v < u; v++){
            bfprod=(Da_->get(h,u,v))*4*phiso[u+col_offset[h]]*
                    phiso[v+col_offset[h]];
            gmat_[g]+=bfprod;
            }
          gmat_[g]+=(Da_->get(h,u,u))*2*phiso[u+col_offset[h]]*
                    phiso[u+col_offset[h]];
          }
        }
        tot+=fabs(gmat_[g]);
        double rr;
        for (int l=0;l<natom;l++){
          double dx=xt-geom.get(l,0);
          double dy=yt-geom.get(l,1);
          double dz=zt-geom.get(l,2);
          rr=dx*dx+dy*dy+dz*dz;
          nmat_[g]-=1.0*Zreal[l]/pow(sigma_,3)*pow(M_PI,-1.5)*
                        exp(-1*rr/pow(sigma_,2));
          gmat_[g]-=1.0*Z[l]/pow(sigma_,3)*pow(M_PI,-1.5)*
                        exp(-1*rr/pow(sigma_,2));
          ndens+=Z[l]/pow(sigma_,3)*pow(M_PI,-1.5)*
                        exp(-1*rr/pow(sigma_,2));
        }
      }
    }
  }
	std::cout << "total charge = " << std::setprecision(14) << tot*pow(g_step_,3) << std::endl;
	std::cout << "total nuc charge  = " << ndens*pow(g_step_,3) << std::endl;
  dens_end=std::chrono::system_clock::now();
  std::chrono::duration<double> dens_seconds=dens_end-dens_start;
  std::cout << "Calc_dens_grid took " << dens_seconds.count() << " seconds\n";
}

void RHFDLMG2::calc_eps_symm()
{
  std::chrono::time_point<std::chrono::system_clock> eps_start, eps_end;
  eps_start = std::chrono::system_clock::now();
	std::cout << "I'm in eps" << "\n";
	double max_eps=0.;
	double tot=0;
	double min_eps=100.;
  std::cout << "reserving eps_half\n";
  std::cout << "just reserved eps_half\n";
	for (int i=0;i<nx_*ny_*nz_*3;i++){
	  eps_half_[i]=0.;
	}
  std::cout << "reserving eps_full\n";
  for (int i=0;i<nx_*ny_*nz_;i++){
    eps_full_[i]=0.;
  }
  int nao = basisset_->nao();
  int nso=0;
  std::cout << "just got nao\n";
  std::cout << "counting nso  in eps symm\n";
  for (int h=0;h<nirrep_;h++) nso+=nsopi_[h];
  std::cout << "nso = " << nso << "\n";
  double phiao[nao];
  double phiso[nso];
  int *col_offset = new int[nirrep_];
  MintsHelper helper(basis_, RHF::options_, 0);
  SharedMatrix aotoso = helper.petite_list(true)->aotoso();
  col_offset[0] = 0;
  for(int h=1; h < nirrep_; h++){
    col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);
  }
  double **u = block_matrix(nao, nso);
  for(int h=0; h < nirrep_; h++){
    for(int j=0; j < aotoso->coldim(h); j++){
      for(int ii=0; ii < nao; ii++){
        u[ii][j+col_offset[h]] = aotoso->get(h, ii, j);
      }
    }
  }

  #pragma omp parallel for private(phiao,phiso)
	for (int i=0;i<nx_;i++){
    double xt=g_step_*i-centerx_;
    for (int j=0;j<ny_;j++){
      double yt=g_step_*j-centery_;
      for (int k=0;k<nz_;k++){
        double zt=g_step_*k-centerz_;
        int findx=i+nx_*j+nx_*ny_*k;
        basisset_->compute_phi(phiao,g_step_*i-centerx_,g_step_*j-centery_,
        g_step_*k-centerz_);
        C_DGEMV('t',nao, nso, 1.0, &(u[0][0]), nso, &(phiao[0]), 1, 0.0,
                &(phiso[0]), 1);
        double this_rho=0.;
        for (int h=0;h<nirrep_;h++){
          for(int u=0; u < nsopi_[h]; u++){
            for(int v=0; v < u; v++){
              double bfprod=(Da_->get(h,u,v))*4*phiso[u+col_offset[h]]*
                                              phiso[v+col_offset[h]];
              this_rho+=bfprod;
            }
          this_rho+=(Da_->get(h,u,u))*2*phiso[u+col_offset[h]]*phiso[u+col_offset[h]];
          }
        }
        eps_full_[findx] = 1.+(eps_inf_-1.)/2.0*(1.+(1.-pow((this_rho/rho0_),beta2_))/
             (1.+pow((this_rho/rho0_),beta2_)));
        this_rho=0.;
        basisset_->compute_phi(phiao,g_step_*(i+0.5)-centerx_,g_step_*j-centery_,
                              g_step_*k-centerz_);
        C_DGEMV('t',nao, nso, 1.0, &(u[0][0]), nso, &(phiao[0]), 1, 0.0, &(phiso[0]), 1);
        for (int h=0;h<nirrep_;h++){
          for(int u=0; u < nsopi_[h]; u++){
            for(int v=0; v < u; v++){
              double bfprod=(Da_->get(h,u,v))*4*phiso[u+col_offset[h]]*
                                              phiso[v+col_offset[h]];
              this_rho+=bfprod;
            }
          this_rho+=(Da_->get(h,u,u))*2*phiso[u+col_offset[h]]*phiso[u+col_offset[h]];
          }
        }


        int indx=i+nx_*j+nx_*ny_*k+nx_*ny_*nz_*0;
        eps_half_[indx] = 1.+(eps_inf_-1.)/2.0*(1.+(1.-pow((this_rho/rho0_),beta2_))/
                       (1.+pow((this_rho/rho0_),beta2_)));


        this_rho=0.;
        basisset_->compute_phi(phiao,g_step_*i-centerx_,g_step_*(j+0.5)-centery_,
                              g_step_*k-centerz_);
        indx=i+nx_*j+nx_*ny_*k+nx_*ny_*nz_*1;
        C_DGEMV('t',nao, nso, 1.0, &(u[0][0]), nso, &(phiao[0]), 1, 0.0, &(phiso[0]), 1);
        for (int h=0;h<nirrep_;h++){
          for(int u=0; u < nsopi_[h]; u++){
            for(int v=0; v < u; v++){
              double bfprod=(Da_->get(h,u,v))*4*phiso[u+col_offset[h]]*
                                              phiso[v+col_offset[h]];
              this_rho+=bfprod;
            }
          this_rho+=(Da_->get(h,u,u))*2*phiso[u+col_offset[h]]*phiso[u+col_offset[h]];
          }
        }

        eps_half_[indx] = 1.+(eps_inf_-1.)/2.0*(1.+(1.-pow((this_rho/rho0_),beta2_))/
                       (1.+pow((this_rho/rho0_),beta2_)));

        this_rho=0.;
        basisset_->compute_phi(phiao,g_step_*i-centerx_,g_step_*j-centery_,
                            g_step_*(k+0.5)-centerz_);
        indx=i+nx_*j+nx_*ny_*k+nx_*ny_*nz_*2;
        C_DGEMV('t',nao, nso, 1.0, &(u[0][0]), nso, &(phiao[0]), 1, 0.0, &(phiso[0]), 1);
        for (int h=0;h<nirrep_;h++){
          for(int u=0; u < nsopi_[h]; u++){
            for(int v=0; v < u; v++){
              double bfprod=(Da_->get(h,u,v))*4*phiso[u+col_offset[h]]*
                                              phiso[v+col_offset[h]];
              this_rho+=bfprod;
            }
          this_rho+=(Da_->get(h,u,u))*2*phiso[u+col_offset[h]]*phiso[u+col_offset[h]];
          }
        }

        eps_half_[indx] = 1.+(eps_inf_-1.)/2.0*(1.+(1.-pow((this_rho/rho0_),beta2_))/
                       (1.+pow((this_rho/rho0_),beta2_)));
      }
    }
	}
	std::cout << "done looping over grid pts \n" << std::endl;
	for (int i=0;i<nx_*ny_*nz_;i++){
    if (eps_full_[i]>max_eps) max_eps=eps_full_[i];
    if (eps_full_[i]<min_eps) min_eps=eps_full_[i];
	}
	std::cout << "Max eps = " << max_eps << ", ";
	std::cout << "Min eps = " << min_eps << std::endl;
  eps_end=std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds=eps_end-eps_start;
  std::cout << "Calc_eps took " << elapsed_seconds.count() << " seconds\n";
}

void RHFDLMG2::calc_V_sol()
{
  std::cout << "made it to v_sol\n";
  int bc[3] = { 2, 2, 2};
  int gstart[3] = {1,1,1};
  int gend[3] = {nx_,ny_,nz_};
  int maxcyc = 40;
  int rptu = 160; 
  std::cout << "bout to get tmpdir\n";
  const char* pathc;
  std::string path;
  if (pathc = std::getenv("TMPDIR")){
    path=pathc;
  }
  else {
    pathc = std::getenv("PWD");
    path=pathc;
  }
  std::cout << "done got tmpdir\n";
  std::string logname = "/dlmg_log.vac";
  if (eps_inf_ > 2.0){
    logname="/dlmg_log.sol";
  }
  std::string movename = "/dlmg_log.old";
  path.append(logname);
  const char * finalname = path.c_str();
  int plen = strlen(finalname);
  std::cout << "filename = " << finalname << "\n";
  std::cout << "Plen = " << plen << "\n";
	for (int i=0;i<nx_*ny_*nz_;i++){
    bound_[i]=0.;
	}
  std::cout << "Bout to call my own calc_bc_vals\n";
  std::cout << "calling espao" << std::endl;
  std::cout << "The espao is gonna be " << std::getenv("TMPDIR")+std::string("/espao");
  std::ifstream look4thisfile(std::getenv("TMPDIR")+std::string("/espao")
       ,std::ios::binary);
  if (!look4thisfile){
  compute_esp_ao();
  }
  compute_esp();
  //compute_esp_old();
  //exit(0);
  std::cout << "just got back from bc vals\n";
  calc_dens_grid();
  int is_damping = 1;
  int pbe_mode = 0;
  double error_tol=pow(10.,-1*dlmg_tol);
  double defect_tol=pow(10.,-1*dlmg_tol);
  int max_iters=100;
  int block_size=-1;
  int min_block_size=3;
  int full_agg_level=0;
  int full_agg_size=0;
  std::cout << "bound_0[0] = " << bound_[0] << "\n";
  DL_MG_INIT(narray,darray,bc,gstart,gend,&maxcyc,&rptu,
            finalname,&plen,&pbe_mode,&is_damping,&error_tol,
            &defect_tol,&max_iters,&full_agg_size,&full_agg_level,
            &min_block_size,&block_size);
  double alpha=-4.0*M_PI;
  double tol=1e-5;
  int ierror=0;
  int npt[3] = {nx_,ny_,nz_};
  std::cout << "Before dl_mg_solver" << std::endl;
	int wrtflag = -1;
  int order=12; 
  std::cout << "Calling solver now\n";
	for (int i=0;i<nx_*ny_*nz_;i++){
    pot_sol_[i]=0.;
	}
  //if (eps_inf_>1.0){
    DL_MG_SOLVER(narray,&pot_sol_[0], &gmat_[0],&order,&bound_[0],&rptu,
        &eps_full_[0],&eps_half_[0]);
  //}
  //else{
  //  DL_MG_SOLVER(narray,&pot_sol_[0], &gmat_[0],&order,&bound_[0],&rptu);
  //}
  std::cout << "Done calling solver now\n";
  
  DL_MG_FREE();
}


double RHFDLMG2::calc_sg_pot()
{
  double tot=0.;
  for (int i=0;i<nx_;i++){
    for (int j=0;j<ny_;j++){
      for (int k=0;k<nz_;k++){
        int indx=i+nx_*j+nx_*ny_*k;
        tot+=nmat_[indx]*pot_sol_[indx];
      }
    }
  }
  tot/=2;
  tot*=(g_step_*g_step_*g_step_);
  return tot;
}

double RHFDLMG2::calc_sg_rep()
{
  double tot=0.;
  std::shared_ptr<Molecule> molecule=Process::environment.molecule();
  Matrix geom=molecule->geometry();
  int natom=molecule->natom();
  double Z[natom];
  double Zreal[natom];
  for (int a=0;a<natom;a++){
    Z[a]=molecule->Z(a);
    Zreal[a]=molecule->Z(a);
  }
  double f1=-1/sqrt(2*M_PI);
  double t1=0.;
	for (int l=0;l<natom;l++){
    t1+=Zreal[l]*Z[l];
	}
  t1=t1*f1*1/sigma_;
  std::cout << "T1 = " << t1 << std::endl;
  double t2=0;
  for (int l=0;l<natom;l++){
    for (int m=0;m<l;m++){
      double dx=geom.get(l,0)-geom.get(m,0);
      double dy=geom.get(l,1)-geom.get(m,1);
      double dz=geom.get(l,2)-geom.get(m,2);
      double r2=dx*dx+dy*dy+dz*dz;
      double r=sqrt(r2);
      double signorm=sqrt(2*sigma_*sigma_);
      t2+=Zreal[l]*Z[m]/r*erf(r/signorm);
    }
      std::cout << "T2 now = " << t2 << std::endl;
  }
  t2*=-1;
  return t1+t2;
}

void RHFDLMG2::calc_V_eff3()
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  std::cout << "entering v eff calc" << std::endl;
	V_eff_->zero();
  //int nthread=omp_get_max_threads();
  int nthread=1;
  std::cout << "num threads = " << nthread << std::endl;
  std::vector< SharedMatrix > thread_mats(nthread);
  for (int i=0;i<nthread;i++){
    thread_mats[i]=SharedMatrix(new Matrix(nirrep_,nsopi_,nsopi_));
    thread_mats[i]->zero();
  }

  MintsHelper helper(basis_,RHF::options_,0);
  SharedMatrix aotoso = helper.petite_list(true)->aotoso();
  int *col_offset = new int [nirrep_];
  col_offset[0]=0;
  for (int h=1;h<nirrep_;h++){
    col_offset[h]=col_offset[h-1]+aotoso->coldim(h-1);
  }

  int nao=basisset_->nao();
  int nso=nso_;
  double **u = block_matrix(nao,nso);
  for(int h=0; h < nirrep_; h++){
    for(int j=0; j < aotoso->coldim(h); j++){
      for(int i=0; i < nao; i++){
          u[i][j+col_offset[h]] = aotoso->get(h, i, j);
      }
    }
  }

      #pragma omp parallel for schedule(guided)
      for (int i=0;i<nx_;i++){
        int myid=omp_get_thread_num();
        double * phiao = new double [nao];
        double * phiso = new double [nso];
        double xt=g_step_*i-centerx_;
        for (int j=0;j<ny_;j++){
          double yt=g_step_*j-centery_;
          for (int k=0;k<nz_;k++){
            int indx=i+nx_*j+nx_*ny_*k;
            double zt=g_step_*k-centerz_;
            basisset_->compute_phi(phiao,xt,yt,zt);
            C_DGEMV('t',nao,nso,1.0,&(u[0][0]),nso,&(phiao[0]),1,0.0,
                    &(phiso[0]),1);
            for (int h=0;h<nirrep_;h++){
              for (int u=0;u<nsopi_[h];u++){
                for (int v=0;v<=u;v++){
                  double val=phiso[u+col_offset[h]]*phiso[v+col_offset[h]]
                            *pot_sol_[indx];
                  thread_mats[myid]->add(h,u,v,val);
                }
              }
            }
          }
        }
      delete [] phiao;
      delete [] phiso;
      }
  // Fill in rest of symmetric matrix
  for (int h=0;h<nirrep_;h++){
    for (int u=0;u<nsopi_[h];u++){
      for (int v=0;v<u;v++){
        for (int t=0;t<nthread;t++){
          V_eff_->add(h,u,v,thread_mats[t]->get(h,u,v));
          V_eff_->add(h,v,u,thread_mats[t]->get(h,u,v));
        }
      }
    }
  }
  for (int h=0;h<nirrep_;h++){
    for (int u=0;u<nsopi_[h];u++){
      for (int t=0;t<nthread;t++){
        V_eff_->add(h,u,u,thread_mats[t]->get(h,u,u));
      }
    }
  }

	V_eff_->scale(pow(g_step_,3));
	V_eff_->print_out();
  std::cout << "leaving v eff calc" << std::endl;
  end=std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds=end-start;
  std::cout << "Calc_V_eff took " << elapsed_seconds.count() << " seconds\n";
}

void RHFDLMG2::compute_esp()
{
    std::chrono::time_point<std::chrono::system_clock> esp_start, esp_end,z_end,y_end,x_end,write_end;
    esp_start = std::chrono::system_clock::now();



    int nbf = basisset_->nbf();

    FILE *gridout = fopen("grid_esp.dat", "w");
    FILE *egridout = fopen("grid_espe.dat", "w");
    FILE *ngridout = fopen("grid_espn.dat", "w");
    if(!gridout)
        throw PSIEXCEPTION("Unable to write to grid_esp.dat");

    //int nthread=omp_get_max_threads();
    int nthread=1;
    if (esp_threaded==0){
      nthread=1;
    }

    SharedMatrix Dtot = shared_from_this()->D_subset_helper(shared_from_this()->Da(), shared_from_this()->Ca(), "AO");
    Dtot->scale(2.0);
    std::ifstream ifile(std::getenv("TMPDIR")+std::string("/espao"),
        std::ios::binary);
    SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
    ints->zero();
    Vector3 origin;
    double * read_value = new double[1];
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
      std::shared_ptr<Molecule> mol = basisset_->molecule();
      origin[0]=xt;
      int myid=omp_get_thread_num();
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        origin[1]=yt;
        for (int k=0;k<nz_;k=k+nz_-1){
          double zt=g_step_*k-centerz_;
          origin[2]=zt;
          int g=i+nx_*j+nx_*ny_*k;
          // Read in AO integrals
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              ifile.read((char *)read_value,sizeof(double));
              //if (g==0){
              //  std::cout << "ints[" << row << "]["<<col<<"] = " << 
              //    *read_value << "\n";
              //}
              ints->set(row,col,*read_value);
              ints->set(col,row,*read_value);
            }
          }
          // Contract density matrix with AOxAO integrals matrix
          double Velec=Dtot->vector_dot(ints);
          double Vnuc = 0.0;
          int natom = mol->natom();
          for(int a=0; a < natom; a++) {
              Vector3 dR = origin - mol->xyz(a);
              double r = dR.norm();
              if(r > 1.0E-8) Vnuc += mol->Z(a)/r;
          }
          bound_[g]=-1.*Velec-Vnuc;
          if (g==0){
            std::cout << "Coleman it's g=0, bound = " << bound_[g]  << "\n";
            //ints->print_out();
          }
        }
      }
    }
    z_end=std::chrono::system_clock::now();
    std::chrono::duration<double> z_seconds=z_end-esp_start;
    std::cout << "Z faces took " << z_seconds.count() << " seconds\n";
    ints->zero();
    std::shared_ptr<Molecule> mol = basisset_->molecule();
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
      origin[0]=xt;
      for (int k=0;k<nz_;k++){
        double zt=g_step_*k-centerz_;
        origin[2]=zt;
        for (int j=0;j<ny_;j=j+ny_-1){
          double yt=g_step_*j-centery_;
          origin[1]=yt;
          int g=i+nx_*j+nx_*ny_*k;
          // Read in AO integrals
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              ifile.read((char *)read_value,sizeof(double));
              ints->set(row,col,*read_value);
              ints->set(col,row,*read_value);
            }
          }
          double Velec=Dtot->vector_dot(ints);
          double Vnuc = 0.0;
          int natom = mol->natom();
          for(int a=0; a < natom; a++) {
              Vector3 dR = origin - mol->xyz(a);
              double r = dR.norm();
              if(r > 1.0E-8)
                  Vnuc += mol->Z(a)/r;
          }
          bound_[g]=-1.*Velec-Vnuc;
          if (g==0){
            std::cout << "Coleman it's g=0, bound = " << bound_[g]  << "\n";
            //ints->print_out();
          }
        }
      }
    }
    y_end=std::chrono::system_clock::now();
    std::chrono::duration<double> y_seconds=y_end-z_end;
    std::cout << "y faces took " << y_seconds.count() << " seconds\n";
    for (int j=0;j<ny_;j++){
      double yt=g_step_*j-centery_;
      origin[1]=yt;
      for (int k=0;k<nz_;k++){
        double zt=g_step_*k-centerz_;
        origin[2]=zt;
        for (int i=0;i<nx_;i=i+nx_-1){
          int g=i+nx_*j+nx_*ny_*k;
          double xt=g_step_*i-centerx_;
          origin[0]=xt;
          // read in ao integrals
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              ifile.read((char *)read_value,sizeof(double));
              ints->set(row,col,*read_value);
              ints->set(col,row,*read_value);
            }
          }
          double Velec=Dtot->vector_dot(ints);
          double Vnuc = 0.0;
          int natom = mol->natom();
          for(int a=0; a < natom; a++) {
              Vector3 dR = origin - mol->xyz(a);
              double r = dR.norm();
              if(r > 1.0E-8)
                  Vnuc += mol->Z(a)/r;
          }
          bound_[g]=-1.*Velec-Vnuc;
          if (g==0){
            std::cout << "Coleman it's g=0, bound = " << bound_[g]  << "\n";
            //ints->print_out();
          }
      }
    }
    }
    x_end=std::chrono::system_clock::now();
    std::chrono::duration<double> x_seconds=x_end-y_end;
    std::cout << "x faces took " << x_seconds.count() << " seconds\n";
    for (int g=0;g<nx_*ny_*nz_;g++){
      bound_[g]/=eps_inf_;
    }
   for (int i=0;i<nx_;i++){
    for (int j=0;j<ny_;j++){
      for (int k=0;k<nz_;k=k+nz_-1){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

   for (int i=0;i<nx_;i++){
    for (int j=0;j<ny_;j=j+ny_-1){
      for (int k=0;k<nz_;k++){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

   for (int i=0;i<nx_;i=i+nx_-1){
    for (int j=0;j<ny_;j++){
      for (int k=0;k<nz_;k++){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

    write_end=std::chrono::system_clock::now();
    std::chrono::duration<double> write_seconds=write_end-x_end;
    std::cout << "write faces took " << write_seconds.count() << " seconds\n";
    fclose(gridout);
    fclose(egridout);
    fclose(ngridout);
    esp_end=std::chrono::system_clock::now();
    std::chrono::duration<double> esp_seconds=esp_end-esp_start;
    std::cout << "Calc_esp_vals took " << esp_seconds.count() << " seconds\n";
}

void RHFDLMG2::compute_esp_ao()
{
    std::chrono::time_point<std::chrono::system_clock> esp_start, esp_end,z_end,y_end,x_end,write_end;
    esp_start = std::chrono::system_clock::now();


    outfile->Printf( "\n Electrostatic potential computed on the grid and written to grid_esp.dat\n");


    int nbf = basisset_->nbf();

   // int nthread=omp_get_max_threads();
   int nthread=1;
    if (esp_threaded==0){
      nthread=1;
    }
    std::vector< SharedMatrix> thread_mats(nthread);
    std::cout << "Creating myelstin shared_ptr\n";
    std::cout << "getting bs1\n";
    std::shared_ptr<BasisSet> mybs1 = integral_->basis1();

    std::vector< std::shared_ptr<MyElstInt>> epot_vec(nthread);
    std::cout << "entering nthread loop with nthread = " << nthread << "\n";
    for (int i=0;i<nthread;i++){
    epot_vec[i]=std::shared_ptr<MyElstInt>((dynamic_cast<MyElstInt*>(new MyElstInt(integral_->spherical_transform(),mybs1,mybs1,0))));
    //epot_vec[i]=std::shared_ptr<MyElstInt>((dynamic_cast<MyElstInt*>(new MyElstInt(integral_->electrostatic()))))
    }
    std::cout << "Ended nthread loop\n";

    SharedMatrix Dtot = shared_from_this()->D_subset_helper(shared_from_this()->Da(), shared_from_this()->Ca(), "AO");
    std::cout << "Dtot just got got\n";
    Dtot->scale(2.0);
    std::cout << "Opening binary file\n";
    std::ofstream ofile(std::getenv("TMPDIR")+std::string("/espao"),
          std::ios::binary);
    std::cout << "about to get in loop" << std::endl;
   //#pragma omp parallel for ordered schedule(static,1)
    for (int i=0;i<nx_;i++){
      std::cout << "i = " << i << "\n";
      double xt=g_step_*i-centerx_;
      double * buf4=new double [ny_*(nbf*(nbf+1))];
      SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
      std::shared_ptr<Molecule> mol = basisset_->molecule();
      int natom=mol->natom();
      Vector3 origin;
      int myid=omp_get_thread_num();
      int maxid=0;
      double Vnuc=0.;
      double Velec=0.;
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        for (int k=0;k<nz_;k=k+nz_-1){
          double zt=g_step_*k-centerz_;
          int g=i+nx_*j+nx_*ny_*k;
          origin[0]=xt;
          origin[1]=yt;
          origin[2]=zt;
          ints->zero();
          int k2=k;
          if (k2>0){
            k2=1;
          }
          epot_vec[myid]->compute(ints,origin);
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              //int indx=j*2*nbf*nbf+k2*nbf*nbf+row*nbf+col;
              int indx=j*2*(nbf*(nbf+1))/2+k2*(nbf*(nbf+1))/2+(row*(row+1))/2+col;
              if (indx>maxid) maxid=indx;
              //std::cout << "indx = \n" << indx;
              buf4[indx]=ints->get(row,col);
            }
          }
          double Velec=Dtot->vector_dot(ints);
          for(int a=0; a < natom; a++) {
            Vector3 dR = origin - mol->xyz(a);
            double r = dR.norm();
            if(r > 1.0E-8) Vnuc += mol->Z(a)/r;
            }
          //std::cout << -1.*Velec-Vnuc << "\n";
          //exit(0);
        }
      }
      #pragma omp ordered
      for (int j=0;j<ny_;j++){
        for (int k=0;k<2;k++){
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              //int indx=j*2*nbf*nbf+k*nbf*nbf+row*nbf+col;
              int indx=j*2*(nbf*(nbf+1))/2+k*(nbf*(nbf+1))/2+(row*(row+1))/2+col;
              ofile.write((char*) &buf4[indx],sizeof(double)); 
            }
          }
        }
      }
      delete [] buf4;
    }
    z_end=std::chrono::system_clock::now();
    std::chrono::duration<double> z_seconds=z_end-esp_start;
    std::cout << "Z faces took " << z_seconds.count() << " seconds\n";
    //exit(0);
    std::cout << "Entering y faces \n";
   //#pragma omp parallel for ordered schedule(static,1)
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
      double * buf4=new double [nz_*nbf*(nbf+1)];
      std::shared_ptr<Molecule> mol = basisset_->molecule();
      SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
      int myid=omp_get_thread_num();
      Vector3 origin;
      for (int k=0;k<nz_;k++){
        double zt=g_step_*k-centerz_;
        for (int j=0;j<ny_;j=j+ny_-1){
          double yt=g_step_*j-centery_;
          origin[0]=xt;
          origin[1]=yt;
          origin[2]=zt;
          ints->zero();
          int j2=j;
          if (j2>0){
            j2=1;
          }
          epot_vec[myid]->compute(ints,origin);
          for (int row=0;row<nbf;row++){
            for (int col=0;col<=row;col++){
              int indx=k*2*(nbf*(nbf+1))/2+j2*(nbf*(nbf+1))/2+(row*(row+1))/
                  2+col;
              buf4[indx]=ints->get(row,col);
            }
          }
        }
      }
      //#pragma omp ordered
        for (int k=0;k<nz_;k++){
          for (int j=0;j<2;j++){
            for (int row=0;row<nbf;row++){
              for (int col=0;col<=row;col++){
                int indx=k*2*(nbf*(nbf+1))/2+j*(nbf*(nbf+1))/2+
                      (row*(row+1))/2+col;
                ofile.write((char*) &buf4[indx],sizeof(double));
              }
            }
          }
        }
      delete [] buf4;
    }
    y_end=std::chrono::system_clock::now();
    std::chrono::duration<double> y_seconds=y_end-z_end;
    std::cout << "y faces took " << y_seconds.count() << " seconds\n";
    std::cout << "entering x faces\n";
    //#pragma omp parallel for ordered schedule(static,1)
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        double * buf4=new double[nz_*(nbf*(nbf+1))];
        std::shared_ptr<Molecule> mol = basisset_->molecule();
        SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
        int myid=omp_get_thread_num();
        Vector3 origin;
        for (int k=0;k<nz_;k++){
          double zt=g_step_*k-centerz_;
          for (int i=0;i<nx_;i=i+nx_-1){
            double xt=g_step_*i-centerx_;
            origin[0]=xt;
            origin[1]=yt;
            origin[2]=zt;
            ints->zero();
            int i2=i;
            if (i2>0){
              i2=1;
            }
            epot_vec[myid]->compute(ints,origin);
            for (int row=0;row<nbf;row++){
              for (int col=0;col<=row;col++){
                int indx=k*2*(nbf*(nbf+1))/2+i2*(nbf*(nbf+1))/2+(row*(row+1))/2
                  +col;
                buf4[indx]=ints->get(row,col);
              }
            }
          }
        }
       // #pragma omp ordered
          for (int k=0;k<nz_;k++){
            for (int i=0;i<2;i++){
              for (int row=0;row<nbf;row++){
                for (int col=0;col<=row;col++){
                  int indx=k*2*(nbf*(nbf+1))/2+i*(nbf*(nbf+1))/2+(row*(row+1))/2
                          +col;
                  ofile.write((char*)&buf4[indx],sizeof(double));
                }
              }
            }
          }
        delete [] buf4;
      }
    x_end=std::chrono::system_clock::now();
    std::chrono::duration<double> x_seconds=x_end-y_end;
    std::cout << "x faces took " << x_seconds.count() << " seconds\n";
    esp_end=std::chrono::system_clock::now();
    std::chrono::duration<double> esp_seconds=esp_end-esp_start;
    std::cout << "Calc_esp_aos took " << esp_seconds.count() << " seconds\n";
}

void RHFDLMG2::compute_esp_old()
{
    std::chrono::time_point<std::chrono::system_clock> esp_start, esp_end,z_end,y_end,x_end,write_end;
    esp_start = std::chrono::system_clock::now();


    outfile->Printf( "\n Electrostatic potential computed on the grid and written to grid_esp.dat\n");


    int nbf = basisset_->nbf();

    FILE *gridout = fopen("grid_esp.dat", "w");
    FILE *egridout = fopen("grid_espe.dat", "w");
    FILE *ngridout = fopen("grid_espn.dat", "w");
    if(!gridout)
        throw PSIEXCEPTION("Unable to write to grid_esp.dat");

    int nthread=omp_get_max_threads();
    if (esp_threaded==0){
      nthread=1;
    }
    std::vector< SharedMatrix> thread_mats(nthread);
    std::cout << "Creating myelstin shared_ptr\n";
    std::cout << "getting bs1\n";
    std::shared_ptr<BasisSet> mybs1 = integral_->basis1();
    std::cout << mybs1 << "\n";
    std::shared_ptr<MyElstInt> epot(dynamic_cast<MyElstInt*>(new MyElstInt(integral_->spherical_transform(),mybs1,mybs1,0)));
    std::cout << epot << "\n";
    std::cout << "did we survive dynamic cast\n";

    std::shared_ptr<ElectrostaticInt> epot_orig(dynamic_cast<ElectrostaticInt*>(integral_->electrostatic()));

    std::vector< std::shared_ptr<MyElstInt>> epot_vec(nthread);
    std::cout << "entering nthread loop with nthread = " << nthread << "\n";
    for (int i=0;i<nthread;i++){
    epot_vec[i]=std::shared_ptr<MyElstInt>((dynamic_cast<MyElstInt*>(new MyElstInt(integral_->spherical_transform(),mybs1,mybs1,0))));
    }

    SharedMatrix Dtot = shared_from_this()->D_subset_helper(shared_from_this()->Da(), shared_from_this()->Ca(), "AO");
    Dtot->scale(2.0);
      #pragma omp parallel for schedule(auto)
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
      SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
      std::shared_ptr<Molecule> mol = basisset_->molecule();
      Vector3 origin;
      int myid=omp_get_thread_num();
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        for (int k=0;k<nz_;k=k+nz_-1){
          double zt=g_step_*k-centerz_;
          int g=i+nx_*j+nx_*ny_*k;
          origin[0]=xt;
          origin[1]=yt;
          origin[2]=zt;
          ints->zero();
          epot_vec[myid]->compute(ints,origin);
          double Velec=Dtot->vector_dot(ints);
          double Vnuc = 0.0;
          int natom = mol->natom();
          for(int a=0; a < natom; a++) {
              Vector3 dR = origin - mol->xyz(a);
              double r = dR.norm();
              if(r > 1.0E-8) Vnuc += mol->Z(a)/r;
          }
          bound_[g]=-1.*Velec-Vnuc;
          std::cout << "bound_g[0] = " << bound_[0] << "\n";
          //exit(0);
        }
      }
    }
    z_end=std::chrono::system_clock::now();
    std::chrono::duration<double> z_seconds=z_end-esp_start;
    std::cout << "Z faces took " << z_seconds.count() << " seconds\n";
      #pragma omp parallel for schedule(auto)
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
        std::shared_ptr<Molecule> mol = basisset_->molecule();
        SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
        int myid=omp_get_thread_num();
        Vector3 origin;
      for (int j=0;j<ny_;j=j+ny_-1){
        double yt=g_step_*j-centery_;
        for (int k=0;k<nz_;k++){
          double zt=g_step_*k-centerz_;
          int g=i+nx_*j+nx_*ny_*k;
          origin[0]=xt;
          origin[1]=yt;
          origin[2]=zt;
          ints->zero();
          epot_vec[myid]->compute(ints,origin);
          double Velec=Dtot->vector_dot(ints);
          double Vnuc = 0.0;
          int natom = mol->natom();
          for(int a=0; a < natom; a++) {
              Vector3 dR = origin - mol->xyz(a);
              double r = dR.norm();
              if(r > 1.0E-8)
                  Vnuc += mol->Z(a)/r;
          }
          bound_[g]=-1.*Velec-Vnuc;
        }
      }
    }
    y_end=std::chrono::system_clock::now();
    std::chrono::duration<double> y_seconds=y_end-z_end;
    std::cout << "y faces took " << y_seconds.count() << " seconds\n";
      #pragma omp parallel for schedule(auto)
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        std::shared_ptr<Molecule> mol = basisset_->molecule();
        SharedMatrix ints(new Matrix("Ex integrals", nbf, nbf));
        int myid=omp_get_thread_num();
        Vector3 origin;
        for (int k=0;k<nz_;k++){
          double zt=g_step_*k-centerz_;
          for (int i=0;i<nx_;i=i+nx_-1){
            int g=i+nx_*j+nx_*ny_*k;
            double xt=g_step_*i-centerx_;
            origin[0]=xt;
            origin[1]=yt;
            origin[2]=zt;
            ints->zero();
            epot_vec[myid]->compute(ints,origin);
            double Velec=Dtot->vector_dot(ints);
            double Vnuc = 0.0;
            int natom = mol->natom();
            for(int a=0; a < natom; a++) {
                Vector3 dR = origin - mol->xyz(a);
                double r = dR.norm();
                if(r > 1.0E-8)
                    Vnuc += mol->Z(a)/r;
            }
            bound_[g]=-1.*Velec-Vnuc;
        }
      }
    }
    x_end=std::chrono::system_clock::now();
    std::chrono::duration<double> x_seconds=x_end-y_end;
    std::cout << "x faces took " << x_seconds.count() << " seconds\n";
    for (int g=0;g<nx_*ny_*nz_;g++){
      bound_[g]/=eps_inf_;
    }
   for (int i=0;i<nx_;i++){
    for (int j=0;j<ny_;j++){
      for (int k=0;k<nz_;k=k+nz_-1){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

   for (int i=0;i<nx_;i++){
    for (int j=0;j<ny_;j=j+ny_-1){
      for (int k=0;k<nz_;k++){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

   for (int i=0;i<nx_;i=i+nx_-1){
    for (int j=0;j<ny_;j++){
      for (int k=0;k<nz_;k++){
        int indx = i+nx_*j+nx_*ny_*k;
        fprintf(gridout, "%16.10f\n", bound_[indx]);
      }
    }
   }

    write_end=std::chrono::system_clock::now();
    std::chrono::duration<double> write_seconds=write_end-x_end;
    std::cout << "write faces took " << write_seconds.count() << " seconds\n";
    fclose(gridout);
    fclose(egridout);
    fclose(ngridout);
    esp_end=std::chrono::system_clock::now();
    std::chrono::duration<double> esp_seconds=esp_end-esp_start;
    std::cout << "Calc_esp_vals took " << esp_seconds.count() << " seconds\n";

}
void RHFDLMG2::calc_V_sgao_P()
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  std::cout << "entering sgVaoP calc" << std::endl;
	V_sgao_->zero();
  std::shared_ptr<Molecule> molecule=Process::environment.molecule();
  Matrix geom=molecule->geometry();
  int natom=molecule->natom();
  double Z[natom];
  for (int a=0;a<natom;a++){
    Z[a]=molecule->Z(a);
  }

  //const int nthread=omp_get_max_threads();
  const int nthread=1;
  std::cout << "num threads = " << nthread << std::endl;
  std::vector< SharedMatrix > thread_mats(nthread);
  for (int i=0;i<nthread;i++){
    thread_mats[i]=SharedMatrix(new Matrix(nirrep_,nsopi_,nsopi_));
    thread_mats[i]->zero();
  }
  std::cout << "After thread_mats created\n"; 
  MintsHelper helper(basis_,RHF::options_,0);
  SharedMatrix aotoso = helper.petite_list(true)->aotoso();
  int *col_offset = new int [nirrep_];
  col_offset[0]=0;
  for (int h=1;h<nirrep_;h++){
    col_offset[h]=col_offset[h-1]+aotoso->coldim(h-1);
  }
  std::cout << "After col_offset filled in\n";
  const int nao=basisset_->nao();
  const int nso=nso_;
  double **u = block_matrix(nao,nso);
  for(int h=0; h < nirrep_; h++){
    for(int j=0; j < aotoso->coldim(h); j++){
      for(int i=0; i < nao; i++){
        u[i][j+col_offset[h]] = aotoso->get(h, i, j);
      }
    }
  }
  std::cout << "After transformation matrix filled in\n";
  std::cout << "Beginning omp loop\n";
  #pragma omp parallel for schedule(auto)
    for (int i=0;i<nx_;i++){
      double xt=g_step_*i-centerx_;
      int myid=omp_get_thread_num();
      for (int j=0;j<ny_;j++){
        double yt=g_step_*j-centery_;
        for (int k=0;k<nz_;k++){
          double zt=g_step_*k-centerz_;
          //int indx=i+nx_*j+nx_*ny_*k;
          double * phiao = new double [nao];
          double * phiso = new double [nso];
          basisset_->compute_phi(phiao,xt,yt,zt);

          C_DGEMV('t',nao,nso,1.0,&(u[0][0]),nso,&(phiao[0]),1,0.0,
                  &(phiso[0]),1);
         
          // calc sigsum here first 
          double sigsum=0.;
          for (int l=0;l<natom;l++){
            double dx=xt-geom.get(l,0);
            double dy=yt-geom.get(l,1);
            double dz=zt-geom.get(l,2);
            double r=sqrt(dx*dx+dy*dy+dz*dz);
            sigsum+=(Z[l]/r*erf(r/sigma_));  
          }

          for (int h=0;h<nirrep_;h++){
            for (int u=0;u<nsopi_[h];u++){
              for (int v=0;v<=u;v++){
                double val=phiso[u+col_offset[h]]*phiso[v+col_offset[h]]*sigsum;
                thread_mats[myid]->add(h,u,v,val);
              }
            }
          }
        delete [] phiao;
        delete [] phiso;
        }
      }
    }
  // Fill in rest of symmetric matrix
  for (int h=0;h<nirrep_;h++){
    for (int u=0;u<nsopi_[h];u++){
      for (int v=0;v<u;v++){
        for (int t=0;t<nthread;t++){
          V_sgao_->add(h,u,v,thread_mats[t]->get(h,u,v));
          V_sgao_->add(h,v,u,thread_mats[t]->get(h,u,v));
        }
      }
    }
  }

  // fill in diagonal part
  for (int h=0;h<nirrep_;h++){
    for (int u=0;u<nsopi_[h];u++){
      for (int t=0;t<nthread;t++){
        V_sgao_->add(h,u,u,thread_mats[t]->get(h,u,u));
      }
    }
  }
	V_sgao_->scale(pow(g_step_,3));
	V_sgao_->print_out();
  std::cout << "leaving v sgao calc" << std::endl;
  end=std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds=end-start;
  std::cout << "Calc_V_sgao took " << elapsed_seconds.count() << " seconds\n";
}


}} // End namespaces

