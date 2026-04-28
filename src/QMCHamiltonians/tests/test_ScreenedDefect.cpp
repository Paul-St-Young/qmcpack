#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
#include "Particle/LongRange/Screen2DKernel.h"
#include "QMCHamiltonians/ScreenedDefect.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Utilities/RuntimeOptions.h"

#include <cmath>
#include <sstream>
#include <vector>

namespace qmcplusplus
{
namespace
{
using Real = QMCTraits::RealType;
using PosType = ParticleSet::PosType;

constexpr Real pi = 3.14159265358979323846;

ParticleSet makeElectrons(const SimulationCell& simulation_cell, const std::vector<PosType>& positions)
{
  ParticleSet elec(simulation_cell);
  elec.setName("e");
  elec.create(std::vector<int>{static_cast<int>(positions.size())});

  for (size_t i = 0; i < positions.size(); ++i)
    elec.R[i] = positions[i];

  SpeciesSet& species         = elec.getSpeciesSet();
  const int up_idx            = species.addSpecies("u");
  const int charge_idx        = species.addAttribute("charge");
  const int mass_idx          = species.addAttribute("mass");
  species(charge_idx, up_idx) = -1.0;
  species(mass_idx, up_idx)   = 1.0;

  elec.resetGroups();
  elec.update();
  return elec;
}

ParticleSet makeDefects(const SimulationCell& simulation_cell, const std::vector<PosType>& positions, Real charge)
{
  ParticleSet defects(simulation_cell);
  defects.setName("ion0");
  defects.create(std::vector<int>{static_cast<int>(positions.size())});

  for (size_t i = 0; i < positions.size(); ++i)
    defects.R[i] = positions[i];

  SpeciesSet& species                    = defects.getSpeciesSet();
  const int defect_idx                   = species.addSpecies("D");
  const int defect_charge_idx            = species.addAttribute("charge");
  species(defect_charge_idx, defect_idx) = charge;

  defects.update();
  return defects;
}

// Brute-force lattice image sum reference of the screened potential.
// Loops far enough that the spline cutoff in the operator is reproduced.
Real bruteLatticeSum(const PosType& disp,
                     const PosType& a1,
                     const PosType& a2,
                     int nrange,
                     Real dgate,
                     int mimg,
                     Real rcut)
{
  Real sum = 0.0;
  for (int ix = -nrange; ix <= nrange; ++ix)
    for (int iy = -nrange; iy <= nrange; ++iy)
    {
      PosType R = ix * a1 + iy * a2;
      PosType d = disp + R;
      Real r    = std::sqrt(dot(d, d));
      if (r >= rcut) continue;
      sum += screen2DKernel<Real>(r, dgate, mimg);
    }
  return sum;
}
} // namespace

TEST_CASE("ScreenedDefect open boundary single electron", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = false;
  lattice.R         = 0.0;
  lattice.R(0, 0)   = 4.0;
  lattice.R(1, 1)   = 5.0;
  lattice.R(2, 2)   = 6.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet defects = makeDefects(simulation_cell, {{0.0, 0.0, 0.0}}, 2.0);
  ParticleSet elec    = makeElectrons(simulation_cell, {{1.0, 0.0, 0.0}});

  elec.addTable(defects);
  elec.update();

  ScreenedDefect screened_defect(elec, defects);

  Libxml2Document doc;
  const char* xml = R"(<pairpot type="screened" source="ion0" target="e" dgate="0.5" max_image="10"/>)";
  REQUIRE(doc.parseFromString(xml));
  REQUIRE(screened_defect.put(doc.getRoot()));

  const Real dgate          = 0.5;
  const int mimg            = 10;
  const Real expected_v     = screen2DKernel<Real>(1.0, dgate, mimg);
  const Real area           = 4.0 * 5.0;
  const Real q_e            = -1.0;
  const Real q_d            = 2.0;
  const Real expected_const = -1.0 * 1 * q_e * q_d * 2.0 * pi * dgate / area;
  const Real expected_value = expected_const + q_e * q_d * expected_v;

  CHECK(screened_defect.evaluate(elec) == Approx(expected_value));

  std::ostringstream output;
  CHECK(screened_defect.get(output));
  CHECK(output.str().find("Screened defect potential") != std::string::npos);
  CHECK(output.str().find("dgate   = 0.5") != std::string::npos);
  CHECK(output.str().find("mimg    = 10") != std::string::npos);
  CHECK(output.str().find("mlat    = 0") != std::string::npos);
}

TEST_CASE("ScreenedDefect open boundary multi-electron and clone", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds = false;
  lattice.R         = 0.0;
  lattice.R(0, 0)   = 6.0;
  lattice.R(1, 1)   = 3.0;
  lattice.R(2, 2)   = 8.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet defects = makeDefects(simulation_cell, {{0.0, 0.0, 0.0}}, 3.0);
  ParticleSet elec    = makeElectrons(simulation_cell, {{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}});

  elec.addTable(defects);
  elec.update();

  ScreenedDefect screened_defect(elec, defects);

  Libxml2Document doc;
  const char* xml = R"(<pairpot type="screened" source="ion0" target="e" dgate="1.5" max_image="10"/>)";
  REQUIRE(doc.parseFromString(xml));
  REQUIRE(screened_defect.put(doc.getRoot()));

  const Real dgate          = 1.5;
  const int mimg            = 10;
  const Real v1             = screen2DKernel<Real>(1.0, dgate, mimg);
  const Real v2             = screen2DKernel<Real>(2.0, dgate, mimg);
  const Real area           = 6.0 * 3.0;
  const Real q_e            = -1.0;
  const Real q_d            = 3.0;
  const Real expected_const = -2.0 * 1 * q_e * q_d * 2.0 * pi * dgate / area;
  const Real expected_value = expected_const + q_e * q_d * (v1 + v2);

  CHECK(screened_defect.evaluate(elec) == Approx(expected_value));

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  std::unique_ptr<OperatorBase> cloned = screened_defect.makeClone(elec, psi);
  REQUIRE(cloned);
  CHECK(cloned->evaluate(elec) == Approx(expected_value));
}

TEST_CASE("ScreenedDefect 2D periodic lattice image sum", "[hamiltonian]")
{
  CrystalLattice<OHMMS_PRECISION, OHMMS_DIM> lattice;
  lattice.BoxBConds    = 0;
  lattice.BoxBConds[0] = 1;
  lattice.BoxBConds[1] = 1;
  lattice.R            = 0.0;
  lattice.R(0, 0)      = 4.0;
  lattice.R(1, 1)      = 4.0;
  lattice.R(2, 2)      = 8.0;
  lattice.reset();

  const SimulationCell simulation_cell(lattice);
  ParticleSet defects = makeDefects(simulation_cell, {{0.0, 0.0, 0.0}}, 1.0);
  ParticleSet elec    = makeElectrons(simulation_cell, {{0.5, 0.0, 0.0}});

  elec.addTable(defects);
  elec.update();

  ScreenedDefect screened_defect(elec, defects);

  Libxml2Document doc;
  const char* xml = R"(<pairpot type="screened" source="ion0" target="e" dgate="1.0" max_image="20"/>)";
  REQUIRE(doc.parseFromString(xml));
  REQUIRE(screened_defect.put(doc.getRoot()));

  const Real dgate     = 1.0;
  const int mimg       = 20;
  const Real rmax      = 15.0 * dgate;
  const PosType a1{4.0, 0.0, 0.0};
  const PosType a2{0.0, 4.0, 0.0};
  const PosType disp{0.5, 0.0, 0.0};
  const Real expected_sum   = bruteLatticeSum(disp, a1, a2, 8, dgate, mimg, rmax);
  const Real area           = 4.0 * 4.0;
  const Real q_e            = -1.0;
  const Real q_d            = 1.0;
  const Real expected_const = -1.0 * 1 * q_e * q_d * 2.0 * pi * dgate / area;
  const Real expected_value = expected_const + q_e * q_d * expected_sum;

  CHECK(screened_defect.evaluate(elec) == Approx(expected_value).epsilon(1e-4));

  std::ostringstream output;
  CHECK(screened_defect.get(output));
  // for dgate=1, rmax=15, a=4 -> mlat = ceil((15+4)/4) = 5
  CHECK(output.str().find("mlat    = 5") != std::string::npos);
}

} // namespace qmcplusplus
