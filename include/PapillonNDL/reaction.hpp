/*
 * Copyright 2020, Hunter Belanger
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
#ifndef PAPILLON_NDL_REACTION_H
#define PAPILLON_NDL_REACTION_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/frame.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <memory>

namespace pndl {

class Reaction {
 public:
  Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid);
  ~Reaction() = default;

  uint32_t MT() const;
  double Q() const;
  double yield(double E) const;
  double threshold() const;
  Frame frame() const;
  double xs(double E) const;
  double xs(double E, size_t i) const;
  AngleEnergyPacket sample_angle_energy(double E_in,
                                        std::function<double()> rng) const;

  const CrossSection& cross_section() const;
  const AngleEnergy& angle_energy() const;
  const Function1D& yield() const;

 private:
  uint32_t mt_;
  double q_;
  double awr_;
  double threshold_;
  Frame frame_;
  CrossSection xs_;
  std::shared_ptr<AngleEnergy> angle_energy_;
  std::shared_ptr<Function1D> yield_;
};

}  // namespace pndl

#endif
