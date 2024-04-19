/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
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
#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>

#include <highfive/highfive.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include "constants.hpp"

// Define the compound datatype for the izaw array
typedef struct {
    int32_t first;
    double second;
} PairIntDouble;

HighFive::CompoundType create_compound_PairIntDouble() {
    return {{"first", HighFive::create_datatype<int32_t>()},
            {"second", HighFive::create_datatype<double>()}};
}
HIGHFIVE_REGISTER_TYPE(PairIntDouble, create_compound_PairIntDouble)

namespace pndl {

// Forward declaration of split_line function
static std::vector<std::string> split_line(std::string line);

ACE::ACE(std::string fname, Type type)
    : zaid_(0, 0),
      temperature_(),
      awr_(),
      fissile_(),
      fname_(fname),
      zaid_txt(10, ' '),
      date_(10, ' '),
      comment_(70, ' '),
      mat_(10, ' '),
      izaw_(),
      nxs_(),
      jxs_(),
      xss_() {
  // Make sure file exists
  if (!std::filesystem::exists(fname)) {
    std::string mssg = "File \"" + fname + "\" does not exist.";
    throw PNDLException(mssg);
  }

  // Open ACE file
  std::ios_base::openmode mode = std::ios_base::in;
  if (type == Type::BINARY) mode = std::ios_base::binary;
  std::ifstream file(fname, mode);

  switch (type) {
    case Type::ASCII:
      read_ascii(file);
      break;

    case Type::BINARY:
      read_binary(file);
      break;
    case Type::HDF5:
      read_hdf5(fname); // because the highfive needs the file name
      break;
    case Type::ADIOS2:
      read_adios2(fname); // adios2 also nees the file name
      break;
  }

  file.close();
}

void ACE::read_ascii(std::ifstream& file) {
  // Check first line to determine header type
  bool legacy_header = true;

  char c1 = static_cast<char>(file.get());
  char c2 = static_cast<char>(file.peek());
  if (c1 == '2' && c2 == '.') legacy_header = false;
  file.unget();

  // Parse header
  std::string awr_txt(12, ' ');
  std::string temp_txt(12, ' ');
  if (legacy_header) {
    file.read(zaid_txt.data(), 10);
    file.read(awr_txt.data(), 12);
    file.read(temp_txt.data(), 12);
    awr_ = std::stod(awr_txt);
    temperature_ = std::stod(temp_txt) * MEV_TO_EV * EV_TO_K;

    // Skip blank char
    file.ignore(1);

    // Read date
    file.read(date_.data(), 10);

    // Ignore the newline chars
    if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);
    if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);

    // Read comment
    file.read(comment_.data(), 70);

    // Read mat id
    file.read(mat_.data(), 10);
  } else {
    std::string line;
    std::getline(file, line);
    // Read next line
    std::getline(file, line);
    std::vector<std::string> split = split_line(line);
    awr_ = std::stod(split[0]);
    temperature_ = std::stod(split[1]) * MEV_TO_EV * EV_TO_K;
    int n_skip = std::stoi(split[3]);

    if (n_skip == 2) {
      // These are the legacy header. Read them
      // Read zaid text
      file.read(zaid_txt.data(), 10);

      // Skip duplicate awr and temp and space
      file.ignore(25);

      // Read date
      file.read(date_.data(), 10);

      // Ignore the newline chars
      if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);
      if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);

      // Read comment
      file.read(comment_.data(), 70);

      // Read mat id
      file.read(mat_.data(), 10);
    } else {
      // Skip comment lines
      for (int i = 0; i < n_skip; i++) std::getline(file, line);
    }
  }

  // Parse IZAW
  for (std::size_t i = 0; i < 16; i++) {
    int32_t i_zaid;
    double i_awr;
    file >> i_zaid;
    file >> i_awr;
    izaw_[i] = {i_zaid, i_awr};
  }

  // Parse NXS
  for (std::size_t i = 0; i < 16; i++) {
    file >> nxs_[i];
  }

  // Parse JXS
  for (std::size_t i = 0; i < 32; i++) {
    file >> jxs_[i];
  }

  // Parse XSS
  xss_.resize(static_cast<std::size_t>(nxs_[0]));
  int i = 0;
  while (!file.eof() && i < nxs_[0]) {
    file >> xss_[static_cast<std::size_t>(i)];
    i++;
  }

  if (i != nxs_[0]) {
    std::string mssg =
        "Found incorrect number of entries in XSS array while reading the \"" +
        fname_ +
        "\" ACE file. This is likely due to a numerical entry which is missing "
        "the \"E\". Please correct the ACE file.";
    throw PNDLException(mssg);
  }

  uint32_t zaid_int = static_cast<uint32_t>(nxs_[1]);
  uint8_t Z_ = static_cast<uint8_t>(zaid_int / 1000);
  uint32_t A_ = zaid_int - (Z_ * 1000);
  zaid_ = ZAID(Z_, A_);

  if (jxs_[1] > 0) fissile_ = true;
}

void ACE::read_binary(std::ifstream& file) {
  // Skip first record length
  file.ignore(4);

  // Skip zaid
  file.read(zaid_txt.data(), 10);

  // Read the AWR
  file.read(reinterpret_cast<char*>(&awr_), sizeof(double));

  // Read the temperatuer
  file.read(reinterpret_cast<char*>(&temperature_), sizeof(double));
  temperature_ *= MEV_TO_EV * EV_TO_K;

  // Read date
  file.read(date_.data(), 10);

  // Read comment
  file.read(comment_.data(), 70);

  // Read mat
  file.read(mat_.data(), 10);

  // Parse IZAW
  for (std::size_t i = 0; i < 16; i++) {
    int32_t i_zaid;
    double i_awr;

    file.read(reinterpret_cast<char*>(&i_zaid), sizeof(int32_t));
    file.read(reinterpret_cast<char*>(&i_awr), sizeof(double));
    izaw_[i] = {i_zaid, i_awr};
  }

  // Parse NXS
  for (std::size_t i = 0; i < 16; i++) {
    file.read(reinterpret_cast<char*>(&nxs_[i]), sizeof(int32_t));
  }

  // Parse JXS
  for (std::size_t i = 0; i < 32; i++) {
    file.read(reinterpret_cast<char*>(&jxs_[i]), sizeof(int32_t));
  }

  // Skip end record length
  file.ignore(4);

  // Parse XSS
  xss_.resize(static_cast<std::size_t>(nxs_[0]));
  uint32_t rlen;
  std::size_t i = 0;
  while (!file.eof() && i < xss_.size()) {
    // Get the record entry length
    file.read(reinterpret_cast<char*>(&rlen), 4);

    file.read(reinterpret_cast<char*>(&xss_[i]), rlen);
    i += rlen / sizeof(double);

    // Skip last len entry
    file.ignore(4);
  }

  uint32_t zaid_int = static_cast<uint32_t>(nxs_[1]);
  uint8_t Z_ = static_cast<uint8_t>(zaid_int / 1000);
  uint32_t A_ = zaid_int - (Z_ * 1000);
  zaid_ = ZAID(Z_, A_);

  if (jxs_[1] > 0) fissile_ = true;
}


void ACE::read_hdf5(std::string& fname) {
  HighFive::File hdf5_file(fname, HighFive::File::ReadOnly);
  // Read all the datasets
  HighFive::DataSet xss = hdf5_file.getDataSet("xss");
  xss.read(xss_);

  HighFive::DataSet nxs = hdf5_file.getDataSet("nxs");
  nxs.read(nxs_);

  HighFive::DataSet jxs = hdf5_file.getDataSet("jxs");
  jxs.read(jxs_);

  HighFive::DataSet izaw = hdf5_file.getDataSet("izaw");
  std::vector<PairIntDouble> izaw_converted;
  izaw.read(izaw_converted);
  for (size_t i = 0; i < izaw_converted.size(); ++i) {
    izaw_[i] = {izaw_converted[i].first, izaw_converted[i].second};
  }

  HighFive::DataSet zaid_ds = hdf5_file.getDataSet("zaid");
  std::vector<uint32_t> zaid;
  zaid_ds.read(zaid);
  zaid_ = ZAID((uint8_t)zaid[0], (uint8_t)zaid[1]); // casting back to uint8_t

  HighFive::DataSet awr = hdf5_file.getDataSet("awr");
  awr.read(awr_);

  HighFive::DataSet temperature = hdf5_file.getDataSet("temperature");
  temperature.read(temperature_);

  HighFive::DataSet fissile = hdf5_file.getDataSet("fissile");
  fissile.read(fissile_);

  HighFive::DataSet mat = hdf5_file.getDataSet("mat");
  mat.read(mat_);

  HighFive::DataSet date = hdf5_file.getDataSet("date");
  date.read(date_);

  HighFive::DataSet comment = hdf5_file.getDataSet("comment");
  comment.read(comment_);

  HighFive::DataSet flname = hdf5_file.getDataSet("fname");
  flname.read(fname_);
}

void ACE::read_adios2(std::string& fname) {
  adios2::ADIOS adios;
  adios2::IO io = adios.DeclareIO("ACE-serial-read");
  adios2::Engine bpReader = io.Open(fname, adios2::Mode::Read);

  bpReader.BeginStep();
  // Read the variables
  adios2::Variable<double> bpxss = io.InquireVariable<double>("xss");
  adios2::Variable<int32_t> bpnxs = io.InquireVariable<int32_t>("nxs");
  adios2::Variable<int32_t> bpjxs = io.InquireVariable<int32_t>("jxs");
  //adios2::Variable<std::pair<int32_t, double>> bpizaw = io.InquireVariable<std::pair<int32_t, double>>("izaw");
  //adios2::Variable<uint32_t> bpzaid = io.InquireVariable<uint32_t>("zaid");

  // Read the attributes
  awr_ = io.InquireAttribute<double>("awr").Data()[0];
  temperature_ = io.InquireAttribute<double>("temperature").Data()[0];
  int fissile_flag = io.InquireAttribute<int>("fissile").Data()[0];
  fissile_ = (bool)fissile_flag; // adios2 does not support bool
  mat_ = io.InquireAttribute<std::string>("mat").Data()[0];
  date_ = io.InquireAttribute<std::string>("date").Data()[0];
  comment_ = io.InquireAttribute<std::string>("comment").Data()[0];
  fname_ = io.InquireAttribute<std::string>("fname").Data()[0];
  zaid_txt = io.InquireAttribute<std::string>("zaid_txt").Data()[0];

  // Read the data
  std::vector<double> xss_data;
  bpReader.Get(bpxss, xss_data);
  xss_ = xss_data; //TODO consider std::move

  std::vector<int32_t> nxs_data;
  bpReader.Get(bpnxs, nxs_data);
  std::copy(nxs_data.begin(), nxs_data.end(), nxs_.begin());

  std::vector<int32_t> jxs_data;
  bpReader.Get(bpjxs, jxs_data);
  std::copy(jxs_data.begin(), jxs_data.end(), jxs_.begin());

  // close the file
  bpReader.EndStep();
  bpReader.Close();

}


void ACE::save_binary(std::string& fname) {
  std::ofstream file(fname, std::ios_base::binary);

  // Write first record length which is size of all
  // of the ACE header
  uint32_t rlen = 100 + 64 * sizeof(int32_t) + 18 * sizeof(double);
  file.write(reinterpret_cast<char*>(&rlen), 4);

  // Write zaid
  file.write(zaid_txt.data(), 10);

  // Write the AWR
  file.write(reinterpret_cast<char*>(&awr_), sizeof(double));

  // Write the temperatuer
  temperature_ /= MEV_TO_EV * EV_TO_K;
  file.write(reinterpret_cast<char*>(&temperature_), sizeof(double));

  // Write date, comment, and mat
  file.write(date_.data(), 10);
  file.write(comment_.data(), 70);
  file.write(mat_.data(), 10);

  // Write IZAW
  for (std::size_t i = 0; i < 16; i++) {
    file.write(reinterpret_cast<char*>(&izaw_[i].first), sizeof(int32_t));
    file.write(reinterpret_cast<char*>(&izaw_[i].second), sizeof(double));
  }

  // Write NXS
  for (std::size_t i = 0; i < 16; i++) {
    file.write(reinterpret_cast<char*>(&nxs_[i]), sizeof(int32_t));
  }

  // Write JXS
  for (std::size_t i = 0; i < 32; i++) {
    file.write(reinterpret_cast<char*>(&jxs_[i]), sizeof(int32_t));
  }

  // Write end record length
  file.write(reinterpret_cast<char*>(&rlen), 4);

  // Write XSS
  const uint32_t ner = 512;
  std::size_t ll = 0;
  std::size_t nn = xss_.size();
  while (nn > 0) {
    std::size_t n = nn;

    if (n > ner) n = ner;

    rlen = static_cast<uint32_t>(n * sizeof(double));

    // Write first len header
    file.write(reinterpret_cast<char*>(&rlen), 4);

    for (std::size_t j = ll; j < ll + n; j++) {
      file.write(reinterpret_cast<char*>(&xss_[j]), sizeof(double));
    }
    ll += n;
    nn -= n;

    // Write second len header
    file.write(reinterpret_cast<char*>(&rlen), 4);
  }

  file.close();
}

void ACE::save_hdf5(HighFive::File& file, std::string gname) {
  // create a new hdf5 file
  //HighFive::File file(fname, HighFive::File::ReadWrite | HighFive::File::Create |
  //                               HighFive::File::Truncate);
  
  // ************* all the vectors ****************//
  // write xss
  HighFive::DataSet xss = file.createDataSet<double>(gname+"/xss", HighFive::DataSpace::From(xss_));
  xss.write(xss_);

  // write nxs
  HighFive::DataSet nxs = file.createDataSet<int32_t>(gname+"/nxs", HighFive::DataSpace::From(nxs_));
  nxs.write(nxs_);

  // write jxs
  HighFive::DataSet jxs = file.createDataSet<int32_t>(gname+"/jxs", HighFive::DataSpace::From(jxs_));
  jxs.write(jxs_);
/*
  // Convert izaw_ to a vector of PairIntDouble
  std::vector<PairIntDouble> izaw_converted;
  for (const auto& pair : izaw_) {
      izaw_converted.push_back({pair.first, pair.second});
  }
  // Write izaw_converted to the HDF5 file
  auto t1 = create_compound_PairIntDouble();
  t1.commit(file, "PairIntDouble");
  auto dataset = file.createDataSet(gname+"/izaw", izaw_converted);
*/

  // write zaid (its just two integers: uint32_t Z and uint32_t A)
  std::vector<uint32_t> zaid = {zaid_.Z(), zaid_.A()};
  HighFive::DataSet zaid_ds = file.createDataSet<uint32_t>(gname+"/zaid", HighFive::DataSpace::From(zaid));
  zaid_ds.write(zaid);

  // *************** all the scalars ****************
  // write awr
  HighFive::DataSet awr = file.createDataSet<double>(gname+"/awr", HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar));
  awr.write(awr_);

  // write temperature
  HighFive::DataSet temperature = file.createDataSet<double>(gname+"/temperature", HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar));
  temperature.write(temperature_);

  // write if fissile
  HighFive::DataSet fissile = file.createDataSet<bool>(gname+"/fissile", HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar));
  fissile.write(fissile_);

  // ************************ strings *******************

  // mat
  auto scalar_dataspace = HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar);
  auto variable_stringtype = HighFive::VariableLengthStringType();
  file.createDataSet(gname+"/mat", scalar_dataspace, variable_stringtype).write(mat_);

  // date
  file.createDataSet(gname+"/date", scalar_dataspace, variable_stringtype).write(date_);

  // comment
  file.createDataSet(gname+"/comment", scalar_dataspace, variable_stringtype).write(comment_);

  // fname
  file.createDataSet(gname+"/fname", scalar_dataspace, variable_stringtype).write(fname_);

  // write zaid_text
  file.createDataSet(gname+"/zaid_txt", scalar_dataspace, variable_stringtype).write(zaid_txt);

}

// save_adios2
void ACE::save_adios2(adios2::IO io, adios2::Engine bpWriter, std::string atom, std::string id)
{
  /* set up adios*/
  //adios2::ADIOS adios;
  //adios2::IO io = adios.DeclareIO("ACE-serial");
  //adios2::Engine bpWriter = io.Open(fname, adios2::Mode::Write);

  // atom is the group name
  // id is the sub group name and the following data goes into this sub group
  // ***************** check if the group exists ****************
  std::string gname = atom + "/" + id + "/";

  // sizes of the arrays
  const std::size_t xss_size = xss_.size();
  const std::size_t nxs_size = nxs_.size();
  const std::size_t jxs_size = jxs_.size();
  //const std::size_t izaw_size = izaw_.size();

  // Define the variables
  adios2::Variable<double> bpxss = io.DefineVariable<double>(gname+"xss", {xss_size}, {0}, {xss_size}, adios2::ConstantDims);
  adios2::Variable<int32_t> bpnxs = io.DefineVariable<int32_t>(gname+"nxs", {nxs_size}, {0}, {nxs_size}, adios2::ConstantDims);
  adios2::Variable<int32_t> bpjxs = io.DefineVariable<int32_t>(gname+"jxs", {jxs_size}, {0}, {jxs_size}, adios2::ConstantDims);
  //adios2::Variable<std::pair<int32_t, double>> bpizaw = io.DefineVariable<std::pair<int32_t, double>>("izaw", {izaw_size}, {0}, {izaw_size}, adios2::ConstantDims);

  adios2::Variable<uint32_t> bpzaid = io.DefineVariable<uint32_t>(gname+"zaid", {2}, {0}, {2}, adios2::ConstantDims);

  // string variables as attributes
  io.DefineAttribute<std::string>(gname+"mat", mat_);
  io.DefineAttribute<std::string>(gname+"date", date_);
  io.DefineAttribute<std::string>(gname+"comment", comment_);
  io.DefineAttribute<std::string>(gname+"fname", fname_);
  io.DefineAttribute<std::string>(gname+"zaid_txt", zaid_txt);

  // scalar variables as attributes
  io.DefineAttribute<double>(gname+"awr", awr_);
  io.DefineAttribute<double>(gname+"temperature", temperature_);
  int fissile_flag = (int)fissile_; //looks like adios2 does not support bool
  io.DefineAttribute<int>(gname+"fissile", fissile_flag);

  // Write the data
  //bpWriter.BeginStep();

  bpWriter.Put(bpxss, xss_.data());
  bpWriter.Put(bpnxs, nxs_.data());
  bpWriter.Put(bpjxs, jxs_.data());
  //bpWriter.Put(bpizaw, izaw_.data());
  // zaid array creation
  std::vector<uint32_t> zaid_arr = {zaid_.Z(), zaid_.A()};
  bpWriter.Put(bpzaid, zaid_arr.data());

  bpWriter.PerformPuts();

  //bpWriter.EndStep();
  //bpWriter.Close();
}

std::vector<std::pair<int32_t, double>> ACE::izaw(std::size_t i,
                                                  std::size_t len) const {
  return {izaw_.begin() + static_cast<std::ptrdiff_t>(i),
          izaw_.begin() + static_cast<std::ptrdiff_t>(i) +
              static_cast<std::ptrdiff_t>(len)};
}

std::vector<int32_t> ACE::nxs(std::size_t i, std::size_t len) const {
  return {nxs_.begin() + static_cast<std::ptrdiff_t>(i),
          nxs_.begin() + static_cast<std::ptrdiff_t>(i) +
              static_cast<std::ptrdiff_t>(len)};
}

std::vector<int32_t> ACE::jxs(std::size_t i, std::size_t len) const {
  return {jxs_.begin() + static_cast<std::ptrdiff_t>(i),
          jxs_.begin() + static_cast<std::ptrdiff_t>(i) +
              static_cast<std::ptrdiff_t>(len)};
}

std::vector<double> ACE::xss(std::size_t i, std::size_t len) const {
  return {xss_.begin() + static_cast<std::ptrdiff_t>(i),
          xss_.begin() + static_cast<std::ptrdiff_t>(i) +
              static_cast<std::ptrdiff_t>(len)};
}

const double* ACE::xss_data() const { return xss_.data(); }

static std::vector<std::string> split_line(std::string line) {
  std::vector<std::string> out;

  std::string tmp = "";
  for (std::size_t i = 0; i < line.size(); i++) {
    if (line[i] != ' ')
      tmp += line[i];
    else {
      if (tmp.size() > 0) {
        out.push_back(tmp);
        tmp = "";
      }
    }
  }

  if (tmp.size() > 0) out.push_back(tmp);

  return out;
}

}  // namespace pndl
