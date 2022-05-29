/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef PAPILLON_NDL_URRPTABLES_H
#define PAPILLON_NDL_URRPTABLES_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/reaction.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <vector>

namespace pndl {

class URRPTables {
  public:
    /**
     * @brief A struct to hold a probability table for a single incident
     *        energy.
     */
    struct PTable {
      /**
       * @brief A struct to hold the cross section values for a single
       *        probability band.
       */
      struct XSBand {
        double total; /**< Total cross section */
        double elastic; /**< Elastic cross section (MT 2) */
        double fission; /**< Fission cross section (MT 18) */
        double capture; /**< Radiative capture cross section (MT 102) */
        double heating; /**< Heating number */
      };

      std::vector<double> cdf; /**< Probability CDF for cross section bands */
      std::vector<XSBand> xs_bands; /**< Cross section bands */
    };
    
    /**
     * @brief A struct to hold all of the microscopic cross sections which
     *        are calculated from the sampled PTables.
     */
    struct MicroXS {
      double total; /**< Total cross section */
      double elastic; /**< Elastic cross section (MT 2) */
      double inelastic; /**< Inelastic cross section */
      double absorption; /**< Absorption cross section */
      double capture; /**< Capture cross section (MT 102) */
      double fission; /**< Fission cross section (MT 18) */
      double heating; /**< Heating number */
    };

  public:
    URRPTables(const ACE& ace,
               const std::shared_ptr<CrossSection>& elastic,
               const std::shared_ptr<CrossSection>& capture,
               const std::shared_ptr<CrossSection>& fission,
               const std::shared_ptr<CrossSection>& heating,
               const std::vector<STReaction>& reactions);
    
    /**
     * @brief Returns true if the PTables are present, and false if not.
     */
    bool is_valid() const {
      return energy_->size() > 2; 
    }
    
    /**
     * @brief Determines the cross section band from the incident energy
     *        and the random variable.
     * @param E Incident energy (MeV).
     * @param xi Random variable in the interval [0,1).
     */ 
    std::size_t sample_xs_band(double E, double xi) const {
      // Find the energy index
      std::size_t i = 0;
      auto Eit = std::lower_bound(energy_->begin(), energy_->end(), E);
      if (Eit == energy_->begin()) {
        i = 0;  
      } else if (Eit == energy_->end()) {
        i = energy_->size() - 1; 
      } else {
        i = std::distance(energy_->begin(), Eit) - 1; 
      }

      // Get reference to the PTable
      const PTable& ptable = (*ptables_)[i];

      // Figure out which band we have sampled
      for (std::size_t b = 0; b < ptable.cdf.size(); b++) {
        if (xi >= ptable.cdf[b]) return b;
      }

      // SHOULD NEVER GET HERE
      return this->n_xs_bands()-1;
    }
    
    /**
     * @brief Calculates the cross section for a given incident energy and
     *        cross section band.
     * @param E Incident energy (MeV).
     * @param i Index of E in the global energy grid, for evaluating cross
     *          sections.
     * @param b Index of the sampled cross section band.
     */
    MicroXS evaluate_xs_band(double E, std::size_t i, std::size_t b) const {
      // Find the energy index and interpolation factor
      std::size_t j = 0;
      double f = 0.;
      std::vector<double>& energy = *energy_;
      auto Eit = std::lower_bound(energy.begin(), energy.end(), E);
      if (Eit == energy.begin()) {
        j = 0;
        f = 0.; 
      } else if (Eit == energy.end()) {
        j = energy.size() - 2; 
        f = 1.;
      } else {
        j = std::distance(energy.begin(), Eit) - 1; 
        if (interp_ == Interpolation::LinLin) {
          f = (E - energy[j])/(energy[j+1] - energy[j]);
        } else {
          f = std::log(E/energy[j]) / std::log(energy[j+1]/energy[j]); 
        }
      }
      
      // XSBand struct which will contain the returned cross sections
      MicroXS xsout{0., 0., 0., 0., 0., 0., 0.};
      std::vector<PTable>& ptabs = *ptables_;

      // Evaluate the cross sections depending on interpolation
      if (interp_ == Interpolation::LinLin) {
        xsout.elastic = ptabs[j].xs_bands[b].elastic + f*(ptabs[j+1].xs_bands[b].elastic - ptabs[j].xs_bands[b].elastic);
        xsout.capture = ptabs[j].xs_bands[b].capture + f*(ptabs[j+1].xs_bands[b].capture - ptabs[j].xs_bands[b].capture);
        xsout.fission = ptabs[j].xs_bands[b].fission + f*(ptabs[j+1].xs_bands[b].fission - ptabs[j].xs_bands[b].fission);
        xsout.heating = ptabs[j].xs_bands[b].heating + f*(ptabs[j+1].xs_bands[b].heating - ptabs[j].xs_bands[b].heating);
      } else {
        if (ptabs[j].xs_bands[b].elastic > 0. &&
            ptabs[j+1].xs_bands[b].elastic > 0.) {
          xsout.elastic = std::exp(std::log(ptabs[j].xs_bands[b].elastic) + f*std::log(ptabs[j+1].xs_bands[b].elastic/ptabs[j].xs_bands[b].elastic));
        } else {
          xsout.elastic = 0.; 
        }

        if (ptabs[j].xs_bands[b].capture > 0. &&
            ptabs[j+1].xs_bands[b].capture > 0.) {
          xsout.capture = std::exp(std::log(ptabs[j].xs_bands[b].capture) + f*std::log(ptabs[j+1].xs_bands[b].capture/ptabs[j].xs_bands[b].capture));
        } else {
          xsout.capture = 0.; 
        }
        
        if (ptabs[j].xs_bands[b].fission > 0. &&
            ptabs[j+1].xs_bands[b].fission > 0.) {
          xsout.fission = std::exp(std::log(ptabs[j].xs_bands[b].fission) + f*std::log(ptabs[j+1].xs_bands[b].fission/ptabs[j].xs_bands[b].fission));
        } else {
          xsout.fission = 0.; 
        }

        if (ptabs[j].xs_bands[b].heating > 0. &&
            ptabs[j+1].xs_bands[b].heating > 0.) {
          xsout.heating = std::exp(std::log(ptabs[j].xs_bands[b].heating) + f*std::log(ptabs[j+1].xs_bands[b].heating/ptabs[j].xs_bands[b].heating));
        } else {
          xsout.heating = 0.; 
        }
      }

      // Check if these are factors. If so, we mulitply by smooth cross sections.
      if (factors_) {
        xsout.elastic *= elastic_->evaluate(E, i); 
        xsout.capture *= capture_->evaluate(E, i);
        xsout.fission *= fission_->evaluate(E, i);
        // I don't know if we are supposed to do this with heating too ??
        // Will need to check this one day.
        xsout.heating *= heating_->evaluate(E, i);
      }

      // Set any negatives to zero.
      if (xsout.elastic < 0.) xsout.elastic = 0.;
      if (xsout.capture < 0.) xsout.capture = 0.;
      if (xsout.fission < 0.) xsout.fission = 0.;
      if (xsout.heating < 0.) xsout.heating = 0.;

      // Now get the inelastic portion (if there is any)
      if (inelastic_) {
        xsout.inelastic = inelastic_->evaluate(E, i); 
      }

      // Now get other absorption portion (if there is any)
      double other_absorption = 0.;
      if (absorption_) {
        other_absorption = absorption_->evaluate(E, i); 
      }
      
      // Calculate full absorption
      xsout.absorption = xsout.capture + xsout.fission + other_absorption;
      
      // Calculate total cross section
      xsout.total = xsout.elastic + xsout.inelastic + xsout.absorption;

      return xsout;
    }
    
    /**
     * @brief Returns the minimum energy of the URR probability tables.
     */
    double min_energy() const {
      if (energy_->size() == 0) return -1.;
      return energy_->front();
    }
    
    /**
     * @brief Returns the maximum energy of the URR probability tables.
     */
    double max_energy() const {
      if (energy_->size() == 0) return -1.;
      return energy_->back();
    }
    
    /**
     * @brief Returns true if provided energy is in the URR energy range.
     * @param E Energy to check, in MeV.
     */
    bool energy_in_range(double E) const {
      if (energy_->size() < 2) return false;
      return this->min_energy() < E && E , this->max_energy(); 
    }
    
    /**
     * @brief Energies for which a PTable is given.
     */
    const std::vector<double>& energy() const { return *energy_; }
    
    /**
     * @brief All PTables for the nuclide.
     */
    const std::vector<PTable>& ptables() const { return *ptables_; }
    
    /**
     * @brief Number of cross section bands in each PTable.
     */
    std::size_t n_xs_bands() const {
      if (ptables_->size() == 0) return 0;
      return ptables_->front().xs_bands.size();
    }
    
  private:
    Interpolation interp_; 
    bool factors_;
    std::shared_ptr<CrossSection> elastic_; // MT 2
    std::shared_ptr<CrossSection> capture_; // MT 102
    std::shared_ptr<CrossSection> fission_; // MT 18
    std::shared_ptr<CrossSection> heating_;
    std::shared_ptr<CrossSection> inelastic_;
    std::shared_ptr<CrossSection> absorption_;
    std::shared_ptr<std::vector<double>> energy_;
    std::shared_ptr<std::vector<PTable>> ptables_;
};

}  // namespace pndl

#endif
