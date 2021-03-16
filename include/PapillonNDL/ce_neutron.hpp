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
#ifndef PAPILLON_NDL_CE_NEUTRON_H
#define PAPILLON_NDL_CE_NEUTRON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/fission_data.hpp>
#include <PapillonNDL/reaction.hpp>
#include <unordered_map>

namespace pndl {

/**
 * @brief Holds all continuous energy neutron data for a single nuclide
 *        and at a single temperature.
 */
class CENeutron {
 public:
  /**
   * @param ace ACE file from which to construct the data.
   */
  CENeutron(const ACE& ace);

  /**
   * @param ace ACE file from which to take the new cross sections.
   * @param nuclide CENeutron containing another instance of the desired
   *                nuclide. Secondary distributions and fission data
   *                will be shared between the two data sets.
   */
  CENeutron(const ACE& ace, const CENeutron& nuclide);

  ~CENeutron() = default;

  /**
   * @brief Returns the nuclide ZAID.
   */
  uint32_t zaid() const { return zaid_; }

  /**
   * @brief Returns the nuclide Atomic Weight Ratio.
   */
  double awr() const { return awr_; }

  /**
   * @brief Returns the temperature at which the data has been prepared.
   */
  double temperature() const { return temperature_; }

  /**
   * @brief Returns true if the nuclide is fissile, and false otherwise.
   */
  bool fissile() const { return fissile_; }

  /**
   * @brief Returns the energy grid for the nuclide.
   */
  const EnergyGrid& energy_grid() const;

  /**
   * @brief Returns a pointer to the total CrossSection for the nuclide.
   */
  std::shared_ptr<CrossSection> total_cross_section() const;

  /**
   * @brief Returns a pointer to the elastic scattering CrossSection for the nuclide.
   */
   std::shared_ptr<CrossSection> elastic_cross_section() const;

  /**
   * @brief Returns a pointer to the disappearance CrossSection for the nuclide.
   */
  std::shared_ptr<CrossSection> disappearance_cross_section() const;

  /**
   * @brief Returns a pointer to the photon production CrossSection for the nuclide.
   */
  std::shared_ptr<CrossSection> photon_production_cross_section() const;

  /**
   * @brief Returns a pointer to the AngleDistribution for elastic scattering.
   */
  std::shared_ptr<AngleDistribution> elastic_angle_distribution() const;

  /**
   * @brief Retrieves the index in the energy grid for an energy.
   * @param E Energy to find in the energy grid.
   */
  size_t energy_grid_index(double E) const {
    return energy_grid_.get_lower_index(E);
  }

  /**
   * @brief Evaluates the total cross section at E using bisection search.
   * @param E Energy.
   */
  double total_xs(double E) const {
    return total_xs_->evaluate(E);
  }

  /**
   * @brief Evaluates the total cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double total_xs(double E, size_t i) const {
    return total_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the elastic scattering cross section at E using bisection
   * search.
   * @param E Energy.
   */
  double elastic_xs(double E) const {
    return elastic_xs_->evaluate(E);
  }

  /**
   * @brief Evaluates the elastic scattering cross section at energy E and index
   * i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double elastic_xs(double E, size_t i) const {
    return elastic_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the disappearance cross section at E using bisection search.
   * @param E Energy.
   */
  double disappearance_xs(double E) const {
    return disappearance_xs_->evaluate(E);
  }

  /**
   * @brief Evaluates the disappearance cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double disappearance_xs(double E, size_t i) const {
    return disappearance_xs_->evaluate(E, i);
  }

  /**
   * @brief Evaluates the photon production cross section at E using bisection search.
   * @param E Energy.
   */
  double photon_production_xs(double E) const {
    // Need to check photon production XS exists, as this one is not
    // necessarily present.
    if(photon_production_xs_) return photon_production_xs_->evaluate(E);
    return 0.;
  }

  /**
   * @brief Evaluates the photon production cross section at energy E and index i.
   * @param E Energy.
   * @param i Index to the energy grid.
   */
  double photon_production_xs(double E, size_t i) const {
    // Need to check photon production XS exists, as this one is not
    // necessarily present.
    if(photon_production_xs_) return photon_production_xs_->evaluate(E, i);
    return 0.;
  }

  /**
   * @brief Samples a scattering angle from the elastic scattering angular
   *        distribution.
   * @param E Incident energy.
   * @param rng Random number generation function.
   */
  double sample_elastic_angle(double E, std::function<double()> rng) const {
    return elastic_angle_->sample_angle(E, rng);
  }

  /**
   * @brief Checks to see if a nucldie has a given reaction.
   * @param mt MT reaction to search for.
   */
  bool has_reaction(uint32_t mt) const {
    if (reactions_.find(mt) == reactions_.end()) return false;
    return true;
  }

  /**
   * @brief Retrieved a given MT reaction.
   * @param mt MT reaction to return.
   */
  const Reaction& reaction(uint32_t mt) const {
    return reactions_.find(mt)->second;
  }

  /**
   * @brief Returns the cross section for a perscriped reaction at a
   *        provided energy. Uses bisection search.
   * @param mt MT value of the reaction.
   * @param E Energy to evaluate cross section at.
   */
  double reaction_xs(uint32_t mt, double E) const {
    if (!has_reaction(mt)) return 0.;

    return reactions_.find(mt)->second.xs(E);
  }

  /**
   * @brief Returns the cross section for a perscriped reaction at a
   *        provided energy, and energy grid index.
   * @param mt MT value of the reaction.
   * @param E Energy to evaluate cross section at.
   * @param i Index to the energy grid for energy E.
   */
  double reaction_xs(uint32_t mt, double E, size_t i) const {
    if (!has_reaction(mt)) return 0.;

    return reactions_.find(mt)->second.xs(E, i);
  }

  /**
   * @brief Returns the fission data for the nuclide.
   */
  const FissionData& fission_data() const { return fission_data_; }

 private:
  uint32_t zaid_;
  double awr_;
  double temperature_;
  bool fissile_;

  EnergyGrid energy_grid_;

  std::shared_ptr<CrossSection> total_xs_;
  std::shared_ptr<CrossSection> disappearance_xs_;
  std::shared_ptr<CrossSection> elastic_xs_;
  std::shared_ptr<CrossSection> photon_production_xs_;

  std::shared_ptr<AngleDistribution> elastic_angle_;

  FissionData fission_data_;

  std::unordered_map<uint32_t, Reaction> reactions_;
};

}  // namespace pndl

#endif