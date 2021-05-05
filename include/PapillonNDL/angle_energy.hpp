/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#ifndef PAPILLON_NDL_ANGLE_ENERGY_H
#define PAPILLON_NDL_ANGLE_ENERGY_H

#include <PapillonNDL/tabulated_1d.hpp>
#include <functional>
#include <memory>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief A struct to hold a sampled angle and energy.
 */
struct AngleEnergyPacket {
  double cosine_angle; /**< Sampled cosine of scattering angle */
  double energy;       /**< Sampled outgoing energy in MeV */
};

/**
 * @brief Interface to represent any secondary angle-energy distribution.
 */
class AngleEnergy {
 public:
  virtual ~AngleEnergy() = default;

  /**
   * @brief Samples an angle and energy from the distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Randum number generation function.
   * @return Sampled cosine of the scattering angle and energy in an
   *         AngleEnergyPacket.
   */
  virtual AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const = 0;
};

}  // namespace pndl

#endif
