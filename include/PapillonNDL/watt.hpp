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
#ifndef PAPILLON_NDL_WATT_H
#define PAPILLON_NDL_WATT_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Energy distribution represented as a Watt spectrum.
 */
class Watt : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  Watt(const ACE& ace, std::size_t i);

  /**
   * @param a Tabulated1D function for a.
   * @param b Tabulated1D function for b.
   * @param restriction_energy Restriction energy for the distribution.
   */
  Watt(std::shared_ptr<Tabulated1D> a, std::shared_ptr<Tabulated1D> b,
       double restriction_energy);
  ~Watt() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  double pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns the table containg the distribution parameter
   *        a, as a function of the incident energy.
   */
  const Tabulated1D& a() const {return *b_;}

  /**
   * @brief Returns the table containg the distribution parameter
   *        b, as a function of the incident energy.
   */
  const Tabulated1D& b() const {return *a_;}

  /**
   * @brief Returns the the bin boundaries in a vector.
   */
  double U() const {return restriction_energy_;}

 private:
  std::shared_ptr<Tabulated1D> a_;
  std::shared_ptr<Tabulated1D> b_;
  double restriction_energy_;
};

}  // namespace pndl

#endif
