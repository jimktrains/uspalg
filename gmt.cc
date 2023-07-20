#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <utility>

int main () {
  std::ifstream f("github.com/evansiroky/timezone-boundary-builder/releases/download/2023b/combined-shapefile-with-oceans.gmt");
  std::string line;

  int state = 0;
  int lines = 0;
  std::vector<std::string> names;
  std::vector<std::vector<std::pair<double, double>>> polys;
  std::vector<int> poly_to_name;
  while (std::getline(f, line)) {
    lines++;
    if (line[0] == '#') {
      if (line[2] == '@') {
        if (line[3] == 'D') {
          if (state == 0) {
            names.push_back(line.substr(4));
          }
        } else if (line[3] == 'P') {
        }
      }
      state = 1;
    } else if (line[0] == '>') {
      state = 0;
      polys.emplace_back();
      poly_to_name.emplace_back(names.size() - 1);
    } else {
      double x,y;
      sscanf(line.c_str(), "%lf %lf", &x, &y);
      polys.back().emplace_back(x,y);
    }
  }


  std::ofstream of("polys", std::ios::binary);
  uint32_t n = polys.size();
  of.write(reinterpret_cast<char*>( &n ), sizeof(n));
  for (auto i : polys) {
    n = i.size();
    of.write(reinterpret_cast<char*>( &n ), sizeof(n));
    for (auto j : i) {
      of.write(reinterpret_cast<char*>( &j.first ), sizeof(j.first));
      of.write(reinterpret_cast<char*>( &j.second), sizeof(j.second));
    }
  }

  std::cout << "Lines:    " << lines << std::endl
            << "Objects:  " << names.size() << std::endl
            << "Polygons: " << polys.size() << std::endl;
}
