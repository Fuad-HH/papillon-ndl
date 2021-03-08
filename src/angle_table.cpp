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
#include <PapillonNDL/angle_table.hpp>

namespace pndl {

AngleTable::AngleTable(const ACE& ace, size_t i) : distribution_(ace, i) {}

AnglePacket AngleTable::sample_mu(double xi) const {
  double mu = distribution_.sample_value(xi);
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
  double pdf = distribution_.pdf(mu);
  return {mu, pdf};
}

size_t AngleTable::size() const { return distribution_.size(); }

const std::vector<double>& AngleTable::cosines() const {
  return distribution_.values();
}

const std::vector<double>& AngleTable::pdf() const {
  return distribution_.pdf();
}

const std::vector<double>& AngleTable::cdf() const {
  return distribution_.cdf();
}

Interpolation AngleTable::interpolation() const {
  return distribution_.interpolation();
}

}  // namespace pndl
