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
#include <PapillonNDL/kalbach_table.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

KalbachTable::KalbachTable(const ACE& ace, size_t i)
    : energy_(), pdf_(), cdf_(), R_(), A_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) && (interp_ != Interpolation::LinLin)) {
    std::string mssg = "KalbachTable::KalbackTable: Invalid interpolation of ";
    mssg += std::to_string(static_cast<int>(interp_)) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
  uint32_t NP = ace.xss<uint32_t>(i + 1);
  energy_ = ace.xss(i + 2, NP);

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);
  R_ = ace.xss(i + 2 + NP + NP + NP, NP);
  A_ = ace.xss(i + 2 + NP + NP + NP + NP, NP);

  if (!std::is_sorted(energy_.begin(), energy_.end())) {
    throw PNDLException("KalbachTable::KalbackTable: Energies are not sorted.", __FILE__, __LINE__);
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    throw PNDLException("KalbachTable::KalbackTable: CDF is not sorted.", __FILE__, __LINE__);
  }
}

}  // namespace pndl
