/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */

#include <PapillonNDL/elastic_svt.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <functional>

#include "constants.hpp"

namespace pndl {

ElasticSVT::ElasticSVT(const AngleDistribution& angle, double awr,
                       double temperature, double tar_threshold)
    : angle_(angle),
      awr_(awr),
      kT_(temperature * K_TO_EV * EV_TO_MEV),
      tar_threshold_(tar_threshold) {
  if (awr_ <= 0.) {
    std::string mssg = "Atomic weight ratio must be greater than zero.";
    throw PNDLException(mssg);
  }

  if (kT_ < 0.) {
    std::string mssg = "Temperature must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (tar_threshold_ < 0.) {
    std::string mssg =
        "Target At Rest threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }
}

double ElasticSVT::temperature() const { return kT_ * MEV_TO_EV * EV_TO_K; }

inline double ElasticSVT::Vector::magnitude() const {
  return std::sqrt(x * x + y * y + z * z);
}

inline ElasticSVT::Vector ElasticSVT::Vector::rotate(double mu,
                                                     double phi) const {
  double xo, yo, zo;
  const double c = std::cos(phi);
  const double s = std::sin(phi);
  const double C = std::sqrt(1. - mu * mu);

  if (std::abs(1. - z * z) > 1.E-10) {
    const double denom = std::sqrt(1. - z * z);

    xo = x * mu + C * (c * x * z - s * y) / denom;
    yo = y * mu + C * (c * y * z + s * x) / denom;
    zo = z * mu - c * C * denom;
  } else {
    const double denom = std::sqrt(1. - y * y);

    xo = x * mu + C * (c * x * y + s * z) / denom;
    yo = y * mu - c * C * denom;
    zo = z * mu + C * (c * y * z - s * x) / denom;
  }

  return {xo, yo, zo};
}

AngleEnergyPacket ElasticSVT::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // Direction in
  const Vector u_n(0., 0., 1.);

  // Get the "velocity" of the incident neutron in the lab frame
  const Vector v_n = u_n * std::sqrt(E_in);

  // Get the "velocity" of the target nuclide
  const Vector v_t = sample_target_velocity(E_in, rng);

  // Calculate the "velocity" of the center of mass.
  const Vector v_cm = (v_n + (v_t * awr_)) / (awr_ + 1.);

  // Calculate the "velocity" of the neutron in the CM frame
  const Vector V_n = v_n - v_cm;

  // Calculate the "speed" in the CM frame
  const double S_n = V_n.magnitude();

  // Calculate direction in the CM frame
  const Vector U_n = V_n / S_n;

  // Sample the scattering angle in the CM frame
  const double mu_cm = angle_.sample_angle(E_in, rng);

  // Get the outgoing velocity in the CM frame
  const Vector V_n_out = U_n.rotate(mu_cm, 2. * PI * rng()) * S_n;

  // Get the outgoing velocity in the LAB frame
  const Vector v_n_out = V_n_out + v_cm;

  // Calculate "speed" in the LAB frame
  const double s_n_out = v_n_out.magnitude();

  // Calculate direction in the LAB frrame
  const Vector u_n_out = v_n_out / s_n_out;

  // Calculate outgoing energy in LAB frame
  const double E_out = s_n_out * s_n_out;

  // Calculate the cosine of the scattering angle in the LAB frame
  const double mu_lab = u_n.dot(u_n_out);

  return {mu_lab, E_out};
}

ElasticSVT::Vector ElasticSVT::sample_target_velocity(
    double Ein, std::function<double()> rng) const {
  if (Ein >= tar_threshold_ * kT_ && awr_ > 1.) {
    return {0., 0., 0.};
  }

  const double y = std::sqrt(awr_ * Ein / kT_);
  double x_sqrd = 0.;
  double mu = 0.;

  const double P_C49 = 2. / (std::sqrt(PI) * y + 2.);

  bool sample_velocity = true;
  while (sample_velocity) {
    if (rng() < P_C49) {
      // Sample x from the distribution C49 in MC sampler
      x_sqrd = -std::log(rng() * rng());
    } else {
      // Sample x from the distribution C61 in MC sampler
      const double c = std::cos(PI / 2.0 * rng());
      x_sqrd = -std::log(rng()) - std::log(rng()) * c * c;
    }

    const double x = std::sqrt(x_sqrd);
    mu = 2. * rng() - 1.;
    const double P_accept =
        std::sqrt(y * y + x_sqrd - 2. * y * x * mu) / (x + y);

    if (rng() < P_accept) sample_velocity = false;
  }

  // Get speed of target
  double s_t = std::sqrt(x_sqrd * kT_ / awr_);

  // Use mu to get the direction vector of the target. We know in the sample
  // method that we always assume the same incident neutron vector:
  const Vector u_n{0., 0., 1.};
  const Vector u_t = u_n.rotate(mu, 2. * PI * rng());

  return u_t * s_t;
}

}  // namespace pndl

/*
 * REFERENCES
 *
 * [1] R. R. Coveyou, R. R. Bate, and R. K. Osborn, “Effect of moderator
 * temperature upon neutron flux in infinite, capturing medium,” J Nucl Energy
 * 1954, vol. 2, no. 3–4, pp. 153–167, 1956, doi: 10.1016/0891-3919(55)90030-9.
 *
 * [2] P. K. Romano and J. A. Walsh, “An improved target velocity sampling
 * algorithm for free gas elastic scattering,” Ann Nucl Energy, vol. 114, no.
 * Ann.  Nucl. Energy 36 2009, pp. 318–324, 2018,
 * doi: 10.1016/j.anucene.2017.12.044.
 */
