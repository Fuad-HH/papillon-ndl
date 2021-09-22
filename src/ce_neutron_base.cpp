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
#include <PapillonNDL/ce_neutron_base.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/difference_1d.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/sum_1d.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <memory>
#include <system_error>
#include <vector>

namespace pndl {

CENeutronBase::CENeutronBase(const ACE& ace)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      fissile_(ace.fissile()),
      elastic_angle_(nullptr),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      delayed_groups_(),
      mt_list_(),
      reaction_indices_() {
  // Make elastic AngleDistribution
  elastic_angle_ =
      std::make_shared<AngleDistribution>(ace, ace.xss<int>(ace.LAND()));

  // Read all reaction MTs
  uint32_t NMT = ace.nxs(3);
  mt_list_.resize(NMT, 0);
  reaction_indices_.fill(-1);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    mt_list_[indx] = ace.xss<uint32_t>(ace.MTR() + indx);
  }

  if (fissile_) {
    read_fission_data(ace);
  } else {
    nu_total_ = std::make_shared<Constant>(0.);
    nu_prompt_ = std::make_shared<Constant>(0.);
    nu_delayed_ = std::make_shared<Constant>(0.);
  }
}

void CENeutronBase::read_fission_data(const ACE& ace) {
  // These will temporarily hold the nu functions untill all are made
  // and we can make the Sum1D or Difference1D instances from them.
  std::shared_ptr<Function1D> total(nullptr);
  std::shared_ptr<Function1D> prompt(nullptr);
  std::shared_ptr<Function1D> delayed(nullptr);

  // If prompt and or total neutrons are given
  if (ace.jxs(1) > 0) {
    if (ace.xss(ace.NU()) > 0.) {
      // Either prompt or total given, but not both
      if (ace.jxs(23) > 0) {
        // Prompt is provided, as delayed is present
        prompt = read_nu(ace, ace.NU());
      } else {
        // No delayed, so this must be total
        total = read_nu(ace, ace.NU());
      }
    } else {
      // Both prompt and total given
      uint32_t KNU_prmpt = ace.NU() + 1;
      uint32_t KNU_tot = ace.NU() + std::abs(ace.xss<int32_t>(ace.NU())) + 1;

      total = read_nu(ace, KNU_tot);
      prompt = read_nu(ace, KNU_prmpt);
    }
  }

  // Read delayed nu if given
  if (ace.DNU() > 0) {
    delayed = read_nu(ace, ace.DNU());
  }

  // First we make nu_total_
  if (total) {
    nu_total_ = total;
  } else if (prompt && delayed) {
    nu_total_ = std::make_shared<Sum1D>(prompt, delayed);
  } else if (delayed) {
    nu_total_ = delayed;
  } else {
    // Should never get here, but just in case
    nu_total_ = std::make_shared<Constant>(0.);
  }

  // Now we make nu_prompt_
  if (prompt) {
    nu_prompt_ = prompt;
  } else {
    nu_prompt_ = nu_total_;
  }

  // And finally, delayed
  if (delayed) {
    nu_delayed_ = delayed;
  } else if (total && prompt) {
    nu_delayed_ = std::make_shared<Difference1D>(total, prompt);
  } else {
    nu_delayed_ = std::make_shared<Constant>(0.);
  }

  // Read all delayed group data
  if (ace.BDD() > 0) {
    uint32_t NGRPS = ace.nxs(7);
    std::size_t g = 1;
    std::size_t i = ace.BDD();
    while (g <= NGRPS) {
      delayed_groups_.push_back(DelayedGroup(ace, i, g));
      uint32_t NR = ace.xss<uint32_t>(i + 1);
      uint32_t NE = ace.xss<uint32_t>(i + 2 + 2 * NR);
      i += 3 + 2 * (NR + NE);
      g++;
    }
  }
}

std::shared_ptr<Function1D> CENeutronBase::read_nu(const ACE& ace,
                                                   std::size_t i) {
  uint32_t LNU = ace.xss<uint32_t>(i);

  if (LNU == 1) {  // Polynomial
    return read_polynomial_nu(ace, ++i);
  } else {  // Tabular
    return read_tabular_nu(ace, ++i);
  }
}

std::shared_ptr<Function1D> CENeutronBase::read_polynomial_nu(const ACE& ace,
                                                              std::size_t i) {
  uint32_t NC = ace.xss<uint32_t>(i);
  std::vector<double> coeffs = ace.xss(i + 1, NC);
  return std::make_shared<Polynomial1D>(coeffs);
}

std::shared_ptr<Function1D> CENeutronBase::read_tabular_nu(const ACE& ace,
                                                           std::size_t i) {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    return std::make_shared<Region1D>(energy, y, interp);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);
    return std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
  }
}

}  // namespace pndl
