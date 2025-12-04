//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jet.cpp
//! \brief Sets up a nonrelativistic jet introduced through L-x1 boundary (left edge)
//========================================================================================

// C headers

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

// BCs on L-x1 (left edge) of grid with jet inflow conditions
void ShockOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh);

namespace {
// Make shock variables global so they can be accessed by BC functions
Real d_sh, p_sh, vx_sh, vy_sh, vz_sh, bx_sh, by_sh, bz_sh;
Real gm1;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // initialize global variables
  d_sh  = pin->GetReal("problem", "d");
  p_sh  = pin->GetReal("problem", "p");
  vx_sh = pin->GetReal("problem", "vx");
  vy_sh = pin->GetReal("problem", "vy");
  vz_sh = pin->GetReal("problem", "vz");
  if (MAGNETIC_FIELDS_ENABLED) {
    bx_sh = pin->GetReal("problem", "bx");
    by_sh = pin->GetReal("problem", "by");
    bz_sh = pin->GetReal("problem", "bz");
  }
  // enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, ShockOuterX1);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Jet problem
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gm1 = peos->GetGamma() - 1.0;

  // initialize conserved variables
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = d_sh;
        phydro->u(IM1,k,j,i) = d_sh*vx_sh;
        phydro->u(IM2,k,j,i) = d_sh*vy_sh;
        phydro->u(IM3,k,j,i) = d_sh*vz_sh;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = p_sh/gm1
                                 + 0.5*d_sh*(SQR(vx_sh)+SQR(vy_sh)+SQR(vz_sh));
        }
      }
    }
  }

  // initialize interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) = bx_sh;
        }
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = by_sh;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = bz_sh;
        }
      }
    }
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            phydro->u(IEN,k,j,i) += 0.5*(SQR(bx_sh) + SQR(by_sh) + SQR(bz_sh));
          }
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ShockOuterX1()
//  \brief Sets boundary condition on right X boundary (iib) for shock problem

void ShockOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                Real time, Real dt,
                int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables in inlet ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
          prim(IDN,k,j,iu+i) = d_sh;
          prim(IVX,k,j,iu+i) = vx_sh;
          prim(IVY,k,j,iu+i) = vy_sh;
          prim(IVZ,k,j,iu+i) = vz_sh;
          prim(IPR,k,j,iu+i) = p_sh;
        
      }
    }
  }

  // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
            b.x1f(k,j,iu+i) = bx_sh;
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
            b.x2f(k,j,iu+i) = by_sh;
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
            b.x3f(k,j,iu+i) = bz_sh;
        }
      }
    }
  }
}
