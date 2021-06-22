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

#include <PapillonNDL/frame.hpp>
#include <tuple>

namespace py = pybind11;

using namespace pndl;

void init_Frame(py::module& m) {
  py::enum_<Frame>(m, "Frame").value("Lab", Frame::Lab).value("CM", Frame::CM);

  // Because fundamental types (i.e. float) are imutable in Python, we need
  // to bind a lambda function which returns a tuple for the first transform
  // overload.
  py::class_<CMToLab>(m, "CMToLab")
      .def("transform",
           [](double Ein, double A, double mu, double Eout) {
             CMToLab::transform(Ein, A, mu, Eout);
             return std::make_tuple(mu, Eout);
           })
      .def("transform", py::overload_cast<double, double, AngleEnergyPacket&>(
                            &CMToLab::transform))
      .def("angle_jacobian", py::overload_cast<double, double, double, double>(
                                 &CMToLab::angle_jacobian))
      .def("angle_jacobian",
           py::overload_cast<double, double, double, double, double>(
               &CMToLab::angle_jacobian))
      .def("angle_jacobian",
           py::overload_cast<double, double, AngleEnergyPacket>(
               &CMToLab::angle_jacobian))
      .def("jacobian", &CMToLab::jacobian);

  py::class_<LabToCM>(m, "LabToCM")
      .def("transform",
           [](double Ein, double A, double mu, double Eout) {
             LabToCM::transform(Ein, A, mu, Eout);
             return std::make_tuple(mu, Eout);
           })
      .def("transform", py::overload_cast<double, double, AngleEnergyPacket&>(
                            &LabToCM::transform))
      .def("angle", &LabToCM::angle);
}
