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
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_table.hpp>
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/isotropic.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

AngleDistribution::AngleDistribution(const ACE& ace, int locb)
    : energy_grid_(), laws_() {
  // Locb must be >= 0! If locb == -1, it means that there is
  // no angular distribution for the reaction (must use product distribution)
  if (locb < 0) {
    std::string mssg =
        "AngleDistribution::AngleDistribution: Must have locb >= 0.\n";
    mssg += "Was provided with locb = " + std::to_string(locb) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (locb > 0) {
    // Set index
    size_t i = ace.AND() + locb - 1;

    // Get number of energies
    uint32_t NE = ace.xss<uint32_t>(i);

    // Read in energies
    energy_grid_ = ace.xss(i + 1, NE);

    // Get each table
    for (uint32_t j = 0; j < NE; j++) {
      int l = ace.xss<int>(i + 1 + NE + j);
      uint32_t loc = ace.AND() + std::abs(l) - 1;

      try {
        if (l > 0) {
          laws_.push_back(std::make_shared<EquiprobableAngleBins>(ace, loc));
        } else if (l < 0) {
          laws_.push_back(std::make_shared<AngleTable>(ace, loc));
        } else {
          laws_.push_back(std::make_shared<Isotropic>());
        }
      } catch (PNDLException& err) {
        std::string mssg =
            "AngleDistribution::AngleDistribution: Could not construct angular "
            "distribution\n";
        mssg += "for energy index = " + std::to_string(j) +
                ", energy = " + std::to_string(energy_grid_[j]);
        mssg += " MeV.\nOccurred at loc = " + std::to_string(loc) +
                ", locb = " + std::to_string(locb) + ".";
        err.add_to_exception(mssg, __FILE__, __LINE__);
        throw err;
      }
    }
  } else if (locb == 0) {
    energy_grid_.push_back(1.E-5);
    laws_.push_back(std::make_shared<Isotropic>());
  }
}

AnglePacket AngleDistribution::sample_angle(double E_in,
                                            std::function<double()> rng) const {
  auto E_it = std::lower_bound(energy_grid_.begin(), energy_grid_.end(), E_in);
  if (E_it == energy_grid_.begin())
    return laws_.front()->sample_mu(rng());
  else if (E_it == energy_grid_.end())
    return laws_.back()->sample_mu(rng());
  E_it--;

  // Get index of low energy
  size_t l = std::distance(energy_grid_.begin(), E_it);
  double f = (E_in - energy_grid_[l]) / (energy_grid_[l + 1] - energy_grid_[l]);

  if (rng() > f)
    return laws_[l]->sample_mu(rng());
  else
    return laws_[l + 1]->sample_mu(rng());
}

size_t AngleDistribution::size() const { return energy_grid_.size(); }

const std::vector<double>& AngleDistribution::energy() const {
  return energy_grid_;
}

double AngleDistribution::energy(size_t i) const { return energy_grid_[i]; }

std::shared_ptr<AngleLaw> AngleDistribution::law(size_t i) const {
  return laws_[i];
}

}  // namespace pndl
