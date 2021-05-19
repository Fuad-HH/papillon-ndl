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
#ifndef PAPILLON_NDL_MULTI_REGION_1D_H
#define PAPILLON_NDL_MULTI_REGION_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/region_1d.hpp>

namespace pndl {

/**
 * @brief Implementation of a Tabulated1D which is comprised of several
 *        interpolation regions.
 */
class MultiRegion1D : public Tabulated1D {
 public:
  /**
   * @param regions Vector of Region1D objects, one for each interpolation
   *                region.
   */
  MultiRegion1D(const std::vector<Region1D>& regions);

  /**
   * @param NBT Vector of the breakpoint location.
   * @param INT Vector of the interpolations for each segment.
   * @param x   Vector containing the x grid.
   * @param y   Vector contianing the y grid.
   */
  MultiRegion1D(const std::vector<uint32_t>& NBT,
                const std::vector<Interpolation>& INT,
                const std::vector<double>& x, const std::vector<double>& y);
  ~MultiRegion1D() = default;

  // Methods required by Function1D
  double operator()(double x) const override final;
  double integrate(double x_low, double x_hi) const override final;

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final;
  std::vector<Interpolation> interpolation() const override final;
  std::vector<double> x() const override final;
  std::vector<double> y() const override final;

  // Extra methods
  /**
   * @brief Returns reference to one of the interpolation segments,
   *        which are stored as Region1D objects.
   * @param i Index to the interpolation segment.
   */
  const Region1D& operator[](std::size_t i) const { return regions_[i]; }

  /**
   * @brief Returns the number of interpolation regions.
   */
  std::size_t size() const { return regions_.size(); }

  /**
   * @brief Returns the lowest x value.
   */
  double min_x() const { return regions_.front().min_x(); }

  /**
   * @brief Returns the highest x value.
   */
  double max_x() const { return regions_.back().max_x(); }

 private:
  std::vector<Region1D> regions_;
};

}  // namespace pndl

#endif
