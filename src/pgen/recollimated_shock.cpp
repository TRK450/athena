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

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  rho_amb  = pin->GetReal("problem", "rho_amb");//number density of ISM in 1/cm^3
  T_amb  = pin->GetReal("problem", "T_amb");//temperature of ISM in K
  if (MAGNETIC_FIELDS_ENABLED) {
    bz_amb = pin->GetReal("problem", "bz_amb");//magnet field in Gauss in ISM and wind
  }
  power_jet  = pin->GetReal("problem", "power_jet");//log10 of jet power in erg/s
  v_jet = pin->GetReal("problem", "v_jet");//speed of jet in km/s
  if (MAGNETIC_FIELDS_ENABLED) {
    bz_jet = pin->GetReal("problem", "bz_jet");//magnetic field in Gauss at jet base
  }
  r_source = pin->GetReal("problem", "r_source");//radius of source in pc
  power_iso = pin->GetReal("problem","power_iso");//log10 of wind power in erg/s
  theta=pin->GetReal("problem","theta");//half opening angle in degree

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;
  d_amb=rho_amb/1.4;//mass density of ISM in ism units(1.4m_H/cm^3).
  p_amb=rho_amb*k_boltzmann_cgs*T_amb/code_energydensity_cgs;
  bz_amb_ism=bz_amb/code_magneticfield_cgs;
  

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

