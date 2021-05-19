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
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/ce_neutron.hpp>

namespace py = pybind11;

using namespace pndl;

void init_CENeutron(py::module& m) {
  py::class_<CENeutron>(m, "CENeutron")
      .def(py::init<const ACE&>())
      .def(py::init<const ACE&, const CENeutron&>())
      .def("zaid", &CENeutron::zaid)
      .def("awr", &CENeutron::awr)
      .def("temperature", &CENeutron::temperature)
      .def("fissile", &CENeutron::fissile)
      .def("energy_grid", &CENeutron::energy_grid)
      .def("total_xs", &CENeutron::total_xs)
      .def("elastic_xs", &CENeutron::elastic_xs)
      .def("fission_xs", &CENeutron::fission_xs)
      .def("disappearance_xs", &CENeutron::disappearance_xs)
      .def("photon_production_xs", &CENeutron::photon_production_xs)
      .def("elastic_angle_distribution", &CENeutron::elastic_angle_distribution)
      .def("nu_total", &CENeutron::nu_total,
           py::return_value_policy::reference_internal)
      .def("nu_prompt", &CENeutron::nu_prompt,
           py::return_value_policy::reference_internal)
      .def("nu_delayed", &CENeutron::nu_delayed,
           py::return_value_policy::reference_internal)
      .def("n_delayed_groups", &CENeutron::n_delayed_groups)
      .def("delayed_group", &CENeutron::delayed_group)
      .def("mt_list", &CENeutron::mt_list)
      .def("has_reaction", &CENeutron::has_reaction)
      .def("reaction", &CENeutron::reaction);
}