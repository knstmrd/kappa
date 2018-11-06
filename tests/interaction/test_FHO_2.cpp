/*
 * \file test_FHO.cpp
 * \brief Test for FHO model for k_VT and VT probability
 */
 
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>        
#include <stdlib.h>

#include "kappa.hpp"
 
int main(int argc, char** argv) {
   
  std::cout << "Start test: computation of thermal diffusion coefficients" << std::endl;
  
  std::string m_source = std::getenv("KAPPA_DATA_DIRECTORY");
  std::cout << "KAPPA_DATA_DIRECTORY is: " << m_source << '\n';
  std::string particle_source    = m_source + "particles.yaml";
  std::string interaction_source = m_source + "interaction.yaml";

  std::cout << "Loading particles data" << std::endl;

  kappa::Molecule mol_N2("N2", false, true, particle_source);
  kappa::Molecule mol_NO("NO", false, true, particle_source);
  kappa::Atom at_N("N", particle_source);
  kappa::Atom at_O("O", particle_source);

  std::cout << "Finished loading particles data" << std::endl;

  std::cout << "N2 vibrational levels: " << mol_N2.num_vibr_levels[0] << std::endl;
  std::cout << "NO vibrational levels: " << mol_NO.num_vibr_levels[0] << std::endl;

  kappa::Approximation appr{};
  kappa::Interaction inter_N2_N(mol_N2, at_N, interaction_source);
  kappa::Interaction inter_N2_O(mol_N2, at_O, interaction_source);
  kappa::Interaction inter_NO_N(at_N, mol_NO, interaction_source);
  kappa::Interaction inter_NO_O(at_O, mol_NO, interaction_source);

  std::vector<double> T_vals = {2000.0, 15000.};
  int i1 = 10, i2 = 9;

  for (auto const& T: T_vals) {
    printf("T=%f\n", T);
    printf("FHO N2(%d)+N -> N2(%d)+N: %.5e\n", i1, i2, appr.k_VT(T, mol_N2, inter_N2_N, i1, (i2-i1)));
    printf("FHO N2(%d)+O -> N2(%d)+O: %.5e\n", i1, i2, appr.k_VT(T, mol_N2, inter_N2_O, i1, (i2-i1)));

    printf("FHO NO(%d)+N -> NO(%d)+O: %.5e\n", i1, i2, appr.k_VT(T, mol_NO, inter_NO_N, i1, (i2-i1)));
    printf("FHO NO(%d)+O -> NO(%d)+O: %.5e\n", i1, i2, appr.k_VT(T, mol_NO, inter_NO_O, i1, (i2-i1)));
  }

  // Output from Python:
  // N2 N 2000.0 2.2566808616748333e-19
  // N2 O 2000.0 3.2893309887433872e-18
  // NO N 2000.0 1.7340272009044324e-18
  // NO O 2000.0 1.881905690986804e-16

  // N2 N 15000.0 1.2551839907129305e-16
  // N2 O 15000.0 1.542073027452873e-16
  // NO N 15000.0 1.4648337853806354e-16
  // NO O 15000.0 2.338709618494171e-16



   return 0;
}
