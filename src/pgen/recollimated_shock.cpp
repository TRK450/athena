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
Real theta;
Real gm1;
} // namespace

void Source(MeshBlock *pmb, const Real time, const Real dt,
              const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
              const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
              AthenaArray<Real> &cons_scalar);

void Mesh::InitUserMeshData(ParameterInput *pin) {
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
  //Real g = pmb->peos->GetGamma();
  //Real temp_goal = 10.0;
  //Real tau = 0.01;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real R = pmb->pcoord->x1v(i);
        Real z = pmb->pcoord->x3v(k);
        Real r = std::sqrt(SQR(R)+SQR(z));
        Real mu = z/r;
        if (r < 1.25*r_source && r> 0.25*r_source) {
           // 比如注入动量或能量

           Real angle=std::acos(z/r);
           Real theta_rad=theta*PI/180.0;
           Real ftheta_rad=1.0/(1.0+std::exp(10.0*(angle-theta_rad)/theta_rad))
             +1.0/(1.0+std::exp(10.0*(PI-angle-theta_rad)/theta_rad));
           Real btemp=0.5*r_source;
           Real fr=1.0/(1.0+std::exp(-40.0*(r-btemp)/btemp))
             *1.0/(1.0+std::exp(40.0*(r-2.0*btemp)/btemp));
           Real ftheta_rad_r=ftheta_rad*fr;
           Real one_minus_ftheta_rad_r=1.0-ftheta_rad_r;
           Real power_iso_code=power_iso*code_time_cgs/(code_mass_cgs*code_velocity_cgs*code_velocity_cgs);
           Real M_dot_iso=2.0*power_iso_code/SQR(v_iso);//mass lossing rate in code unit
           Real d_iso=M_dot_iso/(4*PI*r*r*v_iso);
           Real power_jet_code=power_jet*code_time_cgs/(code_mass_cgs*code_velocity_cgs*code_velocity_cgs);
           Real M_dot_jet=2.0*power_jet_code/SQR(v_jet);//mass lossing rate in code unit
           Real d_jet=M_dot_jet/(2*PI*r*r*(1-std::cos(theta_rad))*v_jet);
           cons(IDN,k,j,i)=d_jet*ftheta_rad_r+d_iso*one_minus_ftheta_rad_r;
           cons(IM1,k,j,i)=d_jet*v_jet*R/r*ftheta_rad_r+d_iso*v_iso*R/r*one_minus_ftheta_rad_r;
           cons(IM2,k,j,i)=0.0;
           cons(IM3,k,j,i)=d_jet*v_jet*z/r*ftheta_rad_r+d_iso*v_iso*z/r*one_minus_ftheta_rad_r;
           cons(IEN,k,j,i)=0.5*d_jet*v_jet*v_jet*ftheta_rad_r+0.5*d_iso*v_iso*v_iso*one_minus_ftheta_rad_r;
           bcc(IB1,k,j,i)=0.0;
           bcc(IB2,k,j,i)=0.0;
           bcc(IB3,k,j,i)=bz_jet_ism*ftheta_rad_r+bz_amb_ism*one_minus_ftheta_rad_r;
        }
      }
    }
  }
  return;
}
