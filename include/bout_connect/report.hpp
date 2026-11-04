#pragma once

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bout_connect
{

// The tracked data/ folder of the package (set by CMake); test outputs go here.
inline std::string data_path(const std::string& name)
{
  return std::string(BC_DATA_DIR) + "/" + name;
}

inline FILE* g_rep = nullptr;

template <class... A>
void say(const char* fmt, A... a)
{
  std::printf(fmt, a...);
  if (g_rep)
    std::fprintf(g_rep, fmt, a...);
  std::fflush(stdout);
}

inline double secs(std::chrono::steady_clock::time_point t0)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
    .count();
}

inline bool is_int(double x)
{
  return std::fabs(x - std::round(x)) < 1e-12;
}

inline std::vector<std::vector<double>> read_csv(const std::string& path,
                                                 size_t ncols = 1)
{
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(in, line); // header
  std::vector<std::vector<double>> rows;
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<double> r;
    while (std::getline(ss, cell, ','))
      r.push_back(std::strtod(cell.c_str(), nullptr));
    if (r.size() >= ncols)
      rows.push_back(std::move(r));
  }
  return rows;
}

} // namespace bout_connect
