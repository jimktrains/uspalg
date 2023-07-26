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

  std::cout << "Lines:    " << lines << std::endl
            << "Objects:  " << names.size() << std::endl
            << "Polygons: " << polys.size() << std::endl;

  auto pittsburgh = Point(Qs10d21(-80.0), Qs10d21(40.44));
  auto chicago = Point(Qs10d21(-87.65), Qs10d21(41.85));
  auto portland = Point(Qs10d21(-122.66), Qs10d21(45.51));

  RTree rtree;
  std::vector<Entry> entries;
  entries.reserve(polys.size());

  BoundingBox maximal;
  for (size_t i = 0; i < polys.size(); i++) {
    auto p = Polygon{polys[i]};
    auto bb = p.boundingBox();
    // rtree.insert(bb, i);
    bb.center();
    entries.emplace_back(bb, i);
    maximal = maximal + bb;

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
    auto p = Polygon{polys[r]};
    auto b = p.boundingBox();
    std::cout << r << " " << names[poly_to_name[r]] << std::endl;
    std::cout << "[" << (double)b.upperleft.x << " " << (double)b.upperleft.y
              << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y
              << "]" << std::endl;
  }
}
