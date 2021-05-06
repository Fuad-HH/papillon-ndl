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
#ifndef PAPILLON_NDL_TABULAR_ENERGY_H
#define PAPILLON_NDL_TABULAR_ENERGY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/pctable.hpp>

namespace pndl {

/**
 * @brief Energy distribution represented as a tabulated PDF and CDF.
 */
class TabularEnergy : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param JED Index offset to find the PDF and CDF tables.
   */
  TabularEnergy(const ACE& ace, std::size_t i, std::size_t JED);

  /**
   * @param incoming_energy Incoming energy grid.
   * @param tables PCTables for each point in the incoming energy grid.
   */
  TabularEnergy(const std::vector<double>& incoming_energy,
                const std::vector<PCTable>& tables);

  ~TabularEnergy() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  double pdf(double E_in, double E_out) const override final;

  /**
   * @brief Reterns the incoming energy points in MeV for which a PCTable
   *        is stored.
   */
  const std::vector<double>& incoming_energy() const {return incoming_energy_;}

  /**
   * @brief Returns the ith PDF/CDF, corresponding to the ith incoming energy.
   * @param i Index in the incoming energy grid.
   */
  const PCTable& table(std::size_t i) const {return tables_[i];}

  /**
   * @brief Returns the number of incoming energy points.
   */
  std::size_t size() const {return incoming_energy_.size();}

 private:
  std::vector<double> incoming_energy_;
  std::vector<PCTable> tables_;
};

}  // namespace pndl

#endif
