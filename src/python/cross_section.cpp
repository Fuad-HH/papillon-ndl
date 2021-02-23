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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <PapillonNDL/cross_section.hpp>

namespace py = pybind11;

using namespace pndl;

void init_CrossSection(py::module& m) {
  py::class_<CrossSection>(m, "CrossSection")
      .def(py::init<>())
      .def(py::init<const ACE&, size_t, const EnergyGrid&, bool>())
      .def("__getitem__", &CrossSection::operator[])
      .def("__call__",
           py::overload_cast<double>(&CrossSection::operator(), py::const_))
      .def("__call__", py::overload_cast<double, size_t>(
                           &CrossSection::operator(), py::const_))
      .def("size", &CrossSection::size)
      .def("index", &CrossSection::index)
      .def("xs", py::overload_cast<size_t>(&CrossSection::xs, py::const_))
      .def("xs", py::overload_cast<>(&CrossSection::xs, py::const_))
      .def("energy",
           py::overload_cast<size_t>(&CrossSection::energy, py::const_))
      .def("energy", py::overload_cast<>(&CrossSection::energy, py::const_));
}
