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
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

LevelInelasticScatter::LevelInelasticScatter(const ACE& ace, size_t i)
    : C1_(), C2_() {
  C1_ = ace.xss(i);
  C2_ = ace.xss(i + 1);
}

LevelInelasticScatter::LevelInelasticScatter(double Q, double AWR)
    : C1_(), C2_() {
  if (AWR <= 0.) {
    std::string mssg =
        "LevelInelasticScattter::LevelInelasticScatter: AWR must be greater "
        "than zero.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  C1_ = (AWR + 1.) * std::abs(Q) / AWR;
  double tmp = AWR / (AWR + 1.);
  C2_ = tmp * tmp;
}

double LevelInelasticScatter::sample_energy(double E_in,
                                            std::function<double()>) const {
  return C2_ * (E_in - C1_);
}

double LevelInelasticScatter::pdf(double E_in, double E_out) const {
  auto rng = []() { return 0.5; };
  double smp_E_out = this->sample_energy(E_in, rng);

  // Make sure that the provided E_out is close to the actual
  // outgoing energy for the discrete distribution. If so, return 1,
  // otherwise return 0.
  if (std::abs(E_out - smp_E_out) > 1.E-15) return 0.;
  return 1.;
}

double LevelInelasticScatter::C1() const { return C1_; }

double LevelInelasticScatter::C2() const { return C2_; }

}  // namespace pndl
