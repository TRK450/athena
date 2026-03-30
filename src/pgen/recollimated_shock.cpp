// C++ headers
#include <cmath>      // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../units/units.hpp"

namespace {
Real rho_amb, T_amb, bz_amb;
Real power_jet, v_jet, bz_jet;
Real power_iso, v_iso, r_source, theta , theta_rad;
Real gm1;
// 预计算的转换值 (Code Units)
Real d_amb, p_amb, bz_amb_ism, bz_jet_ism;
Real power_iso_code, power_jet_code;
Real M_dot_jet, M_dot_iso;
Real code_mass_cgs, code_time_cgs, code_velocity_cgs, code_energydensity_cgs, code_magneticfield_cgs;
} // namespace

void Source(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Units units(pin); 
  code_mass_cgs = units.code_mass_cgs;
  code_time_cgs = units.code_time_cgs;
  code_velocity_cgs  = units.code_velocity_cgs;
  code_energydensity_cgs = units.code_energydensity_cgs;
  code_magneticfield_cgs = units.code_magneticfield_cgs;
  // initialize global variables
  rho_amb  = pin->GetReal("problem", "rho_amb");//number density of ISM in 1/cm^3
  T_amb  = pin->GetReal("problem", "T_amb");//temperature of ISM in K
  if (MAGNETIC_FIELDS_ENABLED) {
    bz_amb = pin->GetReal("problem", "bz_amb");//magnet field in Gauss in ISM and wind
  }
  power_jet  = pin->GetReal("problem", "power_jet");//jet power in erg/s
  v_jet = pin->GetReal("problem", "v_jet");//speed of jet in km/s
  if (MAGNETIC_FIELDS_ENABLED) {
    bz_jet = pin->GetReal("problem", "bz_jet");//magnetic field in Gauss at jet base
  }
  r_source = pin->GetReal("problem", "r_source");//radius of source in pc
  power_iso = pin->GetReal("problem","power_iso");//wind power in erg/s
  v_iso = pin->GetReal("problem", "v_iso");//speed of wind in km/s
  theta=pin->GetReal("problem","theta");//half opening angle in degree
  
  theta_rad=theta*PI/180.0;
  power_iso_code=power_iso*code_time_cgs/(code_mass_cgs*code_velocity_cgs*code_velocity_cgs);
  M_dot_iso=2.0*power_iso_code/SQR(v_iso);//mass lossing rate in code unit
  power_jet_code=power_jet*code_time_cgs/(code_mass_cgs*code_velocity_cgs*code_velocity_cgs);
  M_dot_jet=2.0*power_jet_code/SQR(v_jet);//mass lossing rate in code unit
  
  EnrollUserExplicitSourceFunction(Source);
  
  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;
  d_amb=rho_amb/1.4;//mass density of ISM in ism units(1.4m_H/cm^3).
  p_amb=rho_amb*k_boltzmann_cgs*T_amb/code_energydensity_cgs;
  bz_amb_ism=bz_amb/code_magneticfield_cgs;
  bz_jet_ism=bz_jet/code_magneticfield_cgs;
  

  // initialize conserved variables
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = d_amb;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = p_amb/gm1;
        }
      }
    }
  }

  // initialize interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz_amb_ism;
        }
      }
    }
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real R = pcoord->x1f(i);
          Real z = pcoord->x3f(k);
          Real r = std::sqrt(SQR(R)+SQR(z));
          Real mu = z/r;
          if (r < 1.0*r_source && r> 0.25*r_source) {
            Real angle=std::acos(z/r);
            Real theta_rad=theta*PI/180.0;
            Real ftheta_rad=1.0/(1.0+std::exp(10.0*(angle-theta_rad)/theta_rad))
             +1.0/(1.0+std::exp(10.0*(PI-angle-theta_rad)/theta_rad));
            Real btemp=0.45*r_source;
            Real fr=1.0/(1.0+std::exp(-40.0*(r-btemp)/btemp))
             *1.0/(1.0+std::exp(40.0*(r-2.0*btemp)/btemp));
            pfield->b.x3f(k,j,i) += bz_jet_ism*ftheta_rad*fr;
          }
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += 0.5*SQR(bz_amb_ism);
          }
        }
      }
    }
  }

  return;
}

void Source(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar){
  Real g = pmb->peos->GetGamma();
  //Real temp_goal = 10.0;
  //Real tau = 0.01;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real R = pmb->pcoord->x1f(i);
        Real z = pmb->pcoord->x3f(k);
        Real r = std::sqrt(SQR(R)+SQR(z));
        Real mu = z/r;
        if (r < 1.0*r_source && r> 0.25*r_source) {
           Real cos_val = std::fmin(1.0, std::fmax(-1.0, z/r));
           Real angle = std::acos(cos_val);
           
           Real ftheta_rad=1.0/(1.0+std::exp(10.0*(angle-theta_rad)/theta_rad))
             +1.0/(1.0+std::exp(10.0*(PI-angle-theta_rad)/theta_rad));
           Real btemp=0.45*r_source;
           Real fr=1.0/(1.0+std::exp(-40.0*(r-btemp)/btemp))
             *1.0/(1.0+std::exp(40.0*(r-2.0*btemp)/btemp));
           
           Real d_iso=M_dot_iso/(4*PI*r*r*v_iso);
           
           Real d_jet=M_dot_jet/(2*PI*r*r*(1-std::cos(theta_rad))*v_jet);

           Real dens_current=cons(IDN,k,j,i);
           Real p1_current=cons(IM1,k,j,i);
           Real p2_current=cons(IM2,k,j,i);
           Real p3_current=cons(IM3,k,j,i);
           Real energy_current=cons(IEN,k,j,i);

           Real dens_target=d_jet*ftheta_rad+d_iso*(1-ftheta_rad);
           Real p1_target=d_jet*v_jet*R/r*ftheta_rad+d_iso*v_iso*R/r*(1-ftheta_rad);
           Real p2_target=0.0;
           Real p3_target=d_jet*v_jet*z/r*ftheta_rad+d_iso*v_iso*z/r*(1-ftheta_rad);
           Real energy_target=0.5*(p1_target*p1_target+p2_target*p2_target+p3_target*p3_target)/dens_target
             +0.5*(bcc(IB1,k,j,i)*bcc(IB1,k,j,i)+bcc(IB2,k,j,i)*bcc(IB2,k,j,i)+bcc(IB3,k,j,i)*bcc(IB3,k,j,i))
             +p_amb/(g-1.0);

           cons(IDN,k,j,i)=dens_target*fr+dens_current*(1-fr);
           cons(IM1,k,j,i)=p1_target*fr+p1_current*(1-fr);
           cons(IM2,k,j,i)=p2_target*fr+p2_current*(1-fr);
           cons(IM3,k,j,i)=p3_target*fr+p3_current*(1-fr);
           cons(IEN,k,j,i)=energy_target*fr+energy_current*(1-fr);
        }
      }
    }
  }
  return;
}
