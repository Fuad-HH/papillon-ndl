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
#ifndef PAPILLON_NDL_UNCORRELATED_H
#define PAPILLON_NDL_UNCORRELATED_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <memory>

namespace pndl {

/**
 * @brief AngleDistribution - EnergyLaw pair for secondary distributions
 *        where the angle and energy are not correlated.
 */
class Uncorrelated : public AngleEnergy {
 public:
  /**
   * @param angle Angle distribution for all incoming energies.
   * @param energy Shared pointer to the energy distribution.
   */
  Uncorrelated(std::shared_ptr<AngleDistribution> angle,
               std::shared_ptr<EnergyLaw> energy);
  ~Uncorrelated() = default;

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  /**
   * @brief Returns a pointer to the angular distribution.
   */
  std::shared_ptr<AngleDistribution> angle() const;

  /**
   * @brief Returns a pointer to the energy distribution.
   */
  std::shared_ptr<EnergyLaw> energy() const;

 private:
  std::shared_ptr<AngleDistribution> angle_;
  std::shared_ptr<EnergyLaw> energy_;
};

}  // namespace pndl

#endif
