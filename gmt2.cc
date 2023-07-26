#include "rtree.cc"
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

int main() {

  std::ifstream f("github.com/evansiroky/timezone-boundary-builder/releases/"
                  "download/2023b/combined-shapefile-with-oceans.gmt");
  std::string line;

  int state = 0;
  int lines = 0;
  std::vector<std::string> names;
  std::vector<std::vector<Point>> polys;
  std::vector<int> poly_to_name;
  while (std::getline(f, line)) {
    lines++;
    if (line[0] == '#') {
      if (line[2] == '@') {
        if (line[3] == 'D') {
          names.push_back(line.substr(4));
        } else if (line[3] == 'P') {
        }
      }
      state = 1;
    } else if (line[0] == '>') {
      polys.emplace_back();
      state = 0;
    } else {
      if (state != 2) {
        poly_to_name.emplace_back(names.size() - 1);
        state = 2;
      }
      double x, y;
      sscanf(line.c_str(), "%lf %lf", &x, &y);
      polys.back().emplace_back(Qs10d21(x), Qs10d21(y));
    }
  }

  std::ofstream of("polys", std::ios::binary);
  std::ofstream ofidx("polysidx", std::ios::binary);
  uint32_t n = polys.size();
  unsigned int offset = 0;

  of.write(reinterpret_cast<char *>(&n), sizeof(n));
  offset += sizeof(n);
  for (auto i : polys) {
    ofidx.write(reinterpret_cast<char *>(&offset), sizeof(offset));
    n = i.size();
    of.write(reinterpret_cast<char *>(&n), sizeof(n));
    offset += sizeof(n);
    for (auto j : i) {
      of.write(reinterpret_cast<char *>(&j.x), sizeof(j.x));
      offset += sizeof(j.x);
      of.write(reinterpret_cast<char *>(&j.y), sizeof(j.y));
      offset += sizeof(j.y);
    }
  }

  std::cout << "Lines:    " << lines << std::endl
            << "Objects:  " << names.size() << std::endl
            << "Polygons: " << polys.size() << std::endl;

  auto pittsburgh = Point(Qs10d21(-80.0), Qs10d21(40.44));
  auto chicago = Point(Qs10d21(-87.65), Qs10d21(41.85));
  auto portland = Point(Qs10d21(-122.66), Qs10d21(45.51));

  std::ofstream ofbb("bbox", std::ios::binary);
  RTree rtree;
  std::vector<Entry> entries;
  entries.reserve(polys.size());

  auto p = Polygon{&polys[0][0], polys[0].size()};
  auto bb = p.boundingBox();
  BoundingBox maximal = bb;
  for (size_t i = 0; i < polys.size(); i++) {
    auto p = Polygon{&polys[i][0], polys[i].size()};
    auto bb = p.boundingBox();
    // rtree.insert(bb, i);
    bb.center();
    entries.emplace_back(bb, i);
    maximal = maximal + bb;

    ofbb.write(reinterpret_cast<char *>(&bb.upperleft.x),
               sizeof(bb.upperleft.x));
    ofbb.write(reinterpret_cast<char *>(&bb.upperleft.y),
               sizeof(bb.upperleft.y));
    ofbb.write(reinterpret_cast<char *>(&bb.lowerright.x),
               sizeof(bb.lowerright.x));
    ofbb.write(reinterpret_cast<char *>(&bb.lowerright.y),
               sizeof(bb.lowerright.y));

    if (p && pittsburgh) {
      std::cout << "Found pittsburgh in " << i << " " << names[poly_to_name[i]]
                << std::endl;
    }
    if (bb && pittsburgh) {
      std::cout << "Found pittsburgh in bb for " << i << " "
                << names[poly_to_name[i]] << std::endl;
    }
    if (p && chicago) {
      std::cout << "Found chicago in " << names[poly_to_name[i]] << std::endl;
    }
    if (bb && chicago) {
      std::cout << "Found chicago in bb for " << names[poly_to_name[i]]
                << std::endl;
    }
    if (p && portland) {
      std::cout << "Found portland in " << names[poly_to_name[i]] << std::endl;
    }
    if (bb && portland) {
      std::cout << "Found portland in bb for " << names[poly_to_name[i]]
                << std::endl;
    }
  }

  rtree.SRTLoad(&entries, maximal);

  std::cout << std::endl;
  std::cout << "Node Count: " << nodes.size() << std::endl;
  std::cout << std::endl;

  for (auto r : rtree.find(pittsburgh)) {
    auto p = Polygon{&polys[r][0], polys[r].size()};
    auto b = p.boundingBox();
    std::cout << r << " " << names[poly_to_name[r]] << std::endl;
    std::cout << "[" << (double)b.upperleft.x << " " << (double)b.upperleft.y
              << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y
              << "]" << std::endl;
  }
}
