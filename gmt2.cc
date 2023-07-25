#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include "rtree.cc"

int main () {


  std::ifstream f("github.com/evansiroky/timezone-boundary-builder/releases/download/2023b/combined-shapefile-with-oceans.gmt");
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
      double x,y;
      sscanf(line.c_str(), "%lf %lf", &x, &y);
      polys.back().emplace_back(x,y);
    }
  }


  std::ofstream of("polys", std::ios::binary);
  std::ofstream ofidx("polysidx", std::ios::binary);
  uint32_t n = polys.size();
  unsigned int offset = 0;

  of.write(reinterpret_cast<char*>( &n ), sizeof(n));
  offset += sizeof(n);
  for (auto i : polys) {
    ofidx.write(reinterpret_cast<char*>( &offset ), sizeof(offset));
    n = i.size();
    of.write(reinterpret_cast<char*>( &n ), sizeof(n));
    offset += sizeof(n);
    for (auto j : i) {
      of.write(reinterpret_cast<char*>( &j.x ), sizeof(j.x));
      offset += sizeof(j.x);
      of.write(reinterpret_cast<char*>( &j.y), sizeof(j.y));
      offset += sizeof(j.y);
    }
  }

  std::cout << "Lines:    " << lines << std::endl
            << "Objects:  " << names.size() << std::endl
            << "Polygons: " << polys.size() << std::endl;

  auto pittsburgh = Point(-80.0, 40.44);
  auto chicago = Point(-87.65, 41.85);
  auto portland = Point(-122.66, 45.51);


  std::ofstream ofbb("bbox", std::ios::binary);
  RTree rtree;
  for (int i = 0; i < polys.size(); i++) {
    auto p = Polygon{&polys[i][0], polys[i].size()};
    auto bb = p.boundingBox();
    rtree.insert(bb, i);

    ofbb.write(reinterpret_cast<char*>( &bb.upperleft.x ), sizeof(bb.upperleft.x));
    ofbb.write(reinterpret_cast<char*>( &bb.upperleft.y ), sizeof(bb.upperleft.y));
    ofbb.write(reinterpret_cast<char*>( &bb.lowerright.x ), sizeof(bb.lowerright.x));
    ofbb.write(reinterpret_cast<char*>( &bb.lowerright.y ), sizeof(bb.lowerright.y));

    if (p && pittsburgh) {
      std::cout << "Found pittsburgh in " << i << " " << names[poly_to_name[i]] << std::endl;
    }
    if (bb && pittsburgh) {
      std::cout << "Found pittsburgh in bb for " << i << " " << names[poly_to_name[i]] << std::endl;
    }
    if (p && chicago) {
      std::cout << "Found chicago in " << names[poly_to_name[i]] << std::endl;
    }
    if (bb && chicago) {
      std::cout << "Found chicago in bb for " << names[poly_to_name[i]] << std::endl;
    }
    if (p && portland) {
      std::cout << "Found portland in " << names[poly_to_name[i]] << std::endl;
    }
    if (bb && portland) {
      std::cout << "Found portland in bb for " << names[poly_to_name[i]] << std::endl;
    }
  }
  std::cout << std::endl;
  std::cout << "Node Count: " << nodes.size() << std::endl;
  std::cout << std::endl;

  for (auto r : rtree.find(pittsburgh)) {
    auto p = Polygon{&polys[r][0], polys[r].size()};
    auto b = p.boundingBox();
    std::cout << r << " " << names[poly_to_name[r]] << std::endl;
    std::cout << "[" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
  }

}
