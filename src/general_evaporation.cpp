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
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <cmath>

namespace pndl {

GeneralEvaporation::GeneralEvaporation(const ACE& ace, std::size_t i)
    : temperature_(), bin_bounds_() {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT;
  std::vector<Interpolation> INT;

  if (NR == 0) {
    NBT = {NE};
    INT = {Interpolation::LinLin};
  } else {
    NBT = ace.xss<uint32_t>(i + 1, NR);
    INT = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);

  std::vector<double> temperature = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Get number of bins
  uint32_t NX = ace.xss<uint32_t>(i + 2 + 2 * NR + 2 * NE);

  // Get bins
  bin_bounds_ = ace.xss(i + 2 + 2 * NR + 2 * NE, NX);

  if (!std::is_sorted(bin_bounds_.begin(), bin_bounds_.end())) {
    std::string mssg =
        "GeneralEvaporation::GeneralEvaporation: Bin bounds for X are not "
        "sorted. Index in the XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Create Function1D pointer
  try {
    if (NBT.size() == 1) {
      temperature_ = std::make_shared<Region1D>(energy, temperature, INT[0]);
    } else {
      temperature_ =
          std::make_shared<MultiRegion1D>(NBT, INT, energy, temperature);
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "GeneralEvaporation::GeneralEvaporation: Could not construct Tabular1D "
        "for the effective nuclear temperature. Index in the XSS block is i "
        "= " +
        std::to_string(i) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }
}

GeneralEvaporation::GeneralEvaporation(std::shared_ptr<Tabulated1D> temperature,
                                       const std::vector<double>& bounds)
    : temperature_(temperature), bin_bounds_(bounds) {
  if (!std::is_sorted(bin_bounds_.begin(), bin_bounds_.end())) {
    std::string mssg =
        "GeneralEvaporation::GeneralEvaporation: Bin bounds for X are not "
        "sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

double GeneralEvaporation::sample_energy(double E_in,
                                         std::function<double()> rng) const {
  double T = (*temperature_)(E_in);
  double xi1 = rng();
  std::size_t bin =
      static_cast<std::size_t>(std::floor(bin_bounds_.size() * xi1));
  double xi2 = rng();
  double Chi =
      (bin_bounds_[bin + 1] - bin_bounds_[bin]) * xi2 + bin_bounds_[bin];
  return Chi * T;
}

double GeneralEvaporation::pdf(double E_in, double E_out) const {
  double T = (*temperature_)(E_in);
  double Chi = E_out / T;

  // Go find Chi in bins
  if (Chi < bin_bounds_.front() || Chi > bin_bounds_.back()) return 0.;

  std::size_t bin = 0;
  for (std::size_t i = 0; i < bin_bounds_.size() - 1; i++) {
    if (bin_bounds_[i] <= Chi && bin_bounds_[i + 1] >= Chi) {
      bin = i;
      break;
    }
  }

  double Chi_low = bin_bounds_[bin];
  double Chi_hi = bin_bounds_[bin + 1];
  double nbins = static_cast<double>(bin_bounds_.size() - 1);
  double prob_per_bin = 1. / nbins;

  return prob_per_bin / ((Chi_hi - Chi_low) * T);
}

}  // namespace pndl
