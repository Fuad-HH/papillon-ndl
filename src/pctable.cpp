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
#include <PapillonNDL/pctable.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

PCTable::PCTable(const ACE& ace, size_t i, double normalization)
    : values_(), pdf_(), cdf_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "PCTable::PCTable: Invalid interpolation of ";
    mssg += std::to_string(static_cast<int>(interp_)) + ".";
    mssg += "\nIndex of PCTable in XSS block is " + std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
  uint32_t NP = ace.xss<uint32_t>(i + 1);
  values_ = ace.xss(i + 2, NP);
  // Apply normalization to values
  for (auto& v : values_) v *= normalization;

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);

  if (!std::is_sorted(values_.begin(), values_.end())) {
    std::string mssg = "PCTable::PCTable: Values are not sorted.";
    mssg += "\nIndex of PCTable in XSS block is " + std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg = "PCTable::PCTable: CDF is not sorted.";
    mssg += "\nIndex of PCTable in XSS block is " + std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (cdf_[cdf_.size() - 1] != 1.) {
    // If last element is close to 1, just set it to exactly 1
    if (std::abs(cdf_[cdf_.size() - 1] - 1.) < 1.E-7) {
      cdf_[cdf_.size() - 1] = 1.;
    } else {
      std::string mssg = "PCTable::PCTable: Last CDF entry is not 1, but ";
      mssg += std::to_string(cdf_[cdf_.size() - 1]) + ".";
      mssg += "\nIndex of PCTable in XSS block is " + std::to_string(i) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  for (const auto& p : pdf_) {
    if (p < 0.) {
      std::string mssg = "PCTable::PCTable: Negative value found in PDF.";
      mssg += "\nIndex of PCTable in XSS block is " + std::to_string(i) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

PCTable::PCTable(const std::vector<double>& values,
                 const std::vector<double>& pdf, const std::vector<double>& cdf,
                 Interpolation interp)
    : values_(values), pdf_(pdf), cdf_(cdf), interp_(interp) {
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    std::string mssg = "PCTable::PCTable: Invalid interpolation of ";
    mssg += std::to_string(static_cast<int>(interp_)) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if ((values_.size() != pdf_.size()) || (pdf_.size() != cdf_.size())) {
    std::string mssg = "PCTable::PCTable: Values, PDF, and CDF must";
    mssg += " have the same length.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(values_.begin(), values_.end())) {
    std::string mssg = "PCTable::PCTable: Values are not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    std::string mssg = "PCTable::PCTable: CDF is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (cdf_[cdf_.size() - 1] != 1.) {
    // If last element is close to 1, just set it to exactly 1
    if (std::abs(cdf_[cdf_.size() - 1] - 1.) < 1.E-7) {
      cdf_[cdf_.size() - 1] = 1.;
    } else {
      std::string mssg = "PCTable::PCTable: Last CDF entry is not 1, but ";
      mssg += std::to_string(cdf_[cdf_.size() - 1]) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }

  for (const auto& p : pdf_) {
    if (p < 0.) {
      std::string mssg = "PCTable::PCTable: Negative value found in PDF.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

}  // namespace pndl
