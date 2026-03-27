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
