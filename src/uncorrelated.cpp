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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/uncorrelated.hpp>

namespace pndl {

Uncorrelated::Uncorrelated(const AngleDistribution& angle,
                           std::shared_ptr<EnergyLaw> energy)
    : angle_(angle), energy_(energy) {
  if (!energy_) {
    std::string mssg =
        "Uncorrelated::Uncorrelated: Provided energy distribution is nullptr.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

AngleEnergyPacket Uncorrelated::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  double mu = angle_.sample_angle(E_in, rng);
  double E_out = energy_->sample_energy(E_in, rng);
  return {mu, E_out};
}

std::optional<double> Uncorrelated::angle_pdf(double E_in, double mu) const {
  return angle_.pdf(E_in, mu);
}

std::optional<double> Uncorrelated::pdf(double E_in, double mu,
                                        double E_out) const {
  std::optional<double> energy_pdf = energy_->pdf(E_in, E_out);

  if (!energy_pdf) {
    return std::nullopt;
  }

  double angle_pdf = angle_.pdf(E_in, mu);

  energy_pdf.value() *= angle_pdf;

  return energy_pdf;
}

}  // namespace pndl
