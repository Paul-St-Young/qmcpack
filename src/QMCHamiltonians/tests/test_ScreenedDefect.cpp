#include "catch.hpp"

#include "OhmmsData/Libxml2Doc.h"
#include "Particle/ParticleSet.h"
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

  SpeciesSet& species        = elec.getSpeciesSet();
  const int up_idx           = species.addSpecies("u");
  const int charge_idx       = species.addAttribute("charge");
  const int mass_idx         = species.addAttribute("mass");
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

  SpeciesSet& species            = defects.getSpeciesSet();
  const int defect_idx           = species.addSpecies("D");
  const int defect_charge_idx    = species.addAttribute("charge");
  species(defect_charge_idx, defect_idx) = charge;

  defects.update();
  return defects;
}
} // namespace

TEST_CASE("ScreenedDefect evaluates configured image sum", "[hamiltonian]")
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
  const char* xml = R"(<pairpot type="screened" source="ion0" target="e" dgate="0.5" max_image="1"/>)";
  REQUIRE(doc.parseFromString(xml));
  REQUIRE(screened_defect.put(doc.getRoot()));

  const Real expected_sum    = 1.0 - 2.0 / std::sqrt(2.0);
  const Real expected_vconst = -pi / 40.0;
  const Real expected_value  = expected_vconst - 2.0 * expected_sum;

  CHECK(screened_defect.evaluate(elec) == Approx(expected_value));

  std::ostringstream output;
  CHECK(screened_defect.get(output));
  CHECK(output.str().find("Screened defect potential") != std::string::npos);
  CHECK(output.str().find("dgate   = 0.5") != std::string::npos);
  CHECK(output.str().find("mimg    = 1") != std::string::npos);
}

TEST_CASE("ScreenedDefect sums over electrons and preserves state in clone", "[hamiltonian]")
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
  const char* xml = R"(<pairpot type="screened" source="ion0" target="e" dgate="1.5" max_image="0"/>)";
  REQUIRE(doc.parseFromString(xml));
  REQUIRE(screened_defect.put(doc.getRoot()));

  const Real expected_sum    = 1.0 + 0.5;
  const Real expected_vconst = -2.0 * pi * 1.5 / 18.0;
  const Real expected_value  = expected_vconst - 3.0 * expected_sum;

  CHECK(screened_defect.evaluate(elec) == Approx(expected_value));

  RuntimeOptions runtime_options;
  TrialWaveFunction psi(runtime_options);
  std::unique_ptr<OperatorBase> cloned = screened_defect.makeClone(elec, psi);
  REQUIRE(cloned);
  CHECK(cloned->evaluate(elec) == Approx(expected_value));
}

} // namespace qmcplusplus
