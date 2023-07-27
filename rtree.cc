/*
 *  uspal - micro spatial algorithms
 *  Copyright (C) 2023  James Keener <jim@jimkeener.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <math.h>

#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

#include "qs10d21.cc"

struct Point {
  Qs10d21 x;
  Qs10d21 y;

  Point(Qs10d21 xx, Qs10d21 yy) : x{xx}, y{yy} {}
};

struct Range {
  Qs10d21 start;
  Qs10d21 end;

  Qs10d21 length() { return end - start; }

  bool overlaps(const Range &other) {
    return start <= other.end && end >= other.start;
  }

  bool operator&&(const Range &other) { return overlaps(other); }
};

struct BoundingBox {
  // This bounding box is flipped such that nothing would be contained
  // in it.
  BoundingBox()
      : upperleft{Point{Qs10d21(MAX_Qs10d21), Qs10d21(-MAX_Qs10d21)}},
        lowerright{Point{Qs10d21(-MAX_Qs10d21), Qs10d21(MAX_Qs10d21)}} {};

  BoundingBox(Point ul, Point lr) : upperleft{ul}, lowerright{lr} {};

  Point upperleft;
  Point lowerright;

  Point center() const {
    return Point(upperleft.x + ((lowerright.x - upperleft.x) / Qs10d21(2)),
                 lowerright.y + ((upperleft.y - lowerright.y) / Qs10d21(2)));
  }

  BoundingBox operator+(const BoundingBox &other) const {
    return BoundingBox{Point{std::min(upperleft.x, other.upperleft.x),
                             std::max(upperleft.y, other.upperleft.y)},
                       Point{std::max(lowerright.x, other.lowerright.x),
                             std::min(lowerright.y, other.lowerright.y)}};
  }

  BoundingBox operator+(const Point &other) const {
    return BoundingBox{
        Point{std::min(upperleft.x, other.x), std::max(upperleft.y, other.y)},
        Point{std::max(lowerright.x, other.x),
              std::min(lowerright.y, other.y)}};
  }

  Range xrange() const { return Range{upperleft.x, lowerright.x}; }

  Range yrange() const { return Range{lowerright.y, upperleft.y}; }

  bool overlaps(const BoundingBox &other) const {
    return (xrange() && other.xrange()) && (yrange() && other.yrange());
  }

  bool operator&&(const BoundingBox &other) const { return overlaps(other); }

  bool operator&&(const Point &other) const {
    return (upperleft.x <= other.x && lowerright.x >= other.x) &&
           (upperleft.y >= other.y && lowerright.y <= other.y);
  }

  double area() const {
    return (static_cast<double>(lowerright.x - upperleft.x)) *
           (static_cast<double>(upperleft.y - lowerright.y));
  }

  double wastedArea(const BoundingBox &other) const {
    return (*this + other).area() - this->area() - other.area();
  }
};

struct Entry {
  Entry() {}
  Entry(BoundingBox bbox, uint64_t tuple_id) : bbox{bbox}, tuple_id{tuple_id} {}

  BoundingBox bbox;
  uint64_t tuple_id;

  // Sorting like this for the STR algorithm.
  bool operator<(const Entry &other) const {
    if (bbox.center().x == other.bbox.center().x) {
      return bbox.center().y < other.bbox.center().y;
    }
    return bbox.center().x < other.bbox.center().x;
  }
};

struct Polygon {
  std::vector<Point> points;

  BoundingBox boundingBox() {
    return std::accumulate(points.begin(), points.end(), BoundingBox());
  }

  /*
   * > The Method
   * > I run a semi-infinite ray horizontally (increasing x, fixed y) out from
   * > the test point, and count how many edges it crosses. At each crossing,
   * the > ray switches between inside and outside. This is called the Jordan
   * curve > theorem.
   * >
   * > The case of the ray going thru a vertex is handled correctly
   * > via a careful selection of inequalities. Don't mess with this code unless
   * > you're familiar with the idea of Simulation of Simplicity. This pretends
   * > to shift the ray infinitesimally down so that it either clearly
   * > intersects, or clearly doesn't touch. Since this is merely a conceptual,
   * > infinitesimal, shift, it never creates an intersection that didn't exist
   * > before, and never destroys an intersection that clearly existed before.
   * >
   * > The ray is tested against each edge thus:
   * >
   * >     Is the point in the half-plane to the left of the extended edge? and
   * >     Is the point's Y coordinate within the edge's Y-range?
   * >
   * > Handling endpoints here is tricky.
   * https://wrfranklin.org/Research/Short_Notes/pnpoly.html
   */
  bool contains(const Point &p) const {
    bool c = false;
    for (size_t i = 0, j = points.size() - 1; i < points.size(); j = i++) {
      if (((points[i].y > p.y) != (points[j].y > p.y)) &&
          (p.x < (points[j].x - points[i].x) * (p.y - points[i].y) /
                         (points[j].y - points[i].y) +
                     points[i].x))
        c = !c;
    }
    return c;
  }

  bool operator&&(const Point &p) const { return contains(p); }
};

constexpr static int RTREE_MAX_CHILDREN_COUNT = 15;
constexpr static int RTREE_MIN_CHILDREN_COUNT = RTREE_MAX_CHILDREN_COUNT / 2;

/*
 * Guttman, Antomn. "R-Trees - A Dynamic Index Structure for Spatial
 *   Searching." ACM SIGMOD Record, vol. 14, no. 2, June 1984, pp. 47–57.,
 *   https://doi.org/10.1145/971697.602266.
 *
 * Leutenegger, Scott T., et al. STR:  A Simple and Efficient Algorithm for
 *   R-Tree Packing. Institute for Computer Applications in Science and
 *   Engineering NASA Langley Research Center, Feb. 1997.
 *   NASA Contract No. NAS1-19480
 */
struct RTreeNode {
  RTreeNode(const Entry &e, bool il, bool hp, uint64_t pid) {
    isLeaf = il;
    hasParent = hp;
    parentid = pid;
    children_count = 1;
    children[0] = e;
  }

  uint8_t children_count;
  uint64_t myid;
  bool isLeaf;
  bool hasParent;
  uint64_t parentid;
  Entry children[RTREE_MAX_CHILDREN_COUNT];

  uint64_t chooseLeaf(const BoundingBox &bbox) {
    if (isLeaf) {
      return myid;
    }

    double min_waste = children[0].bbox.wastedArea(bbox);
    double min_waste_area = children[0].bbox.area();
    uint8_t min_waste_i = 0;
    for (uint8_t i = 1; i < children_count; i++) {
      auto wasted = children[i].bbox.wastedArea(bbox);
      auto area = children[i].bbox.area();
      if ((wasted < min_waste) ||
          (wasted == min_waste && area < min_waste_area)) {
        min_waste = wasted;
        min_waste_area = area;
        min_waste_i = i;
      }
    }

    return children[min_waste_i].tuple_id;
  }

  int64_t split(std::vector<std::unique_ptr<RTreeNode>> *nodes) {
    double max_wasted_space = 0;
    // uint8_t max_wasted_i = 0;
    uint8_t max_wasted_j = 0;
    uint8_t i, j;
    for (i = 0; i < (children_count - 1); i++) {
      for (j = i + 1; j < children_count; j++) {
        auto wasted_area = children[i].bbox.wastedArea(children[j].bbox);
        if (wasted_area > max_wasted_space) {
          max_wasted_space = wasted_area;
          // max_wasted_i = i;
          max_wasted_j = j;
        }
      }
    }

    auto new_nodeb = std::unique_ptr<RTreeNode>(
        new RTreeNode(children[max_wasted_j], isLeaf, true, parentid));
    // Not thread-safe.
    new_nodeb->myid = nodes->size();
    auto new_node = new_nodeb->myid;
    nodes->push_back(std::move(new_nodeb));

    // auto i_bb = children[max_wasted_i].bbox;
    auto j_bb = children[max_wasted_j].bbox;

    std::copy(std::begin(children) + max_wasted_j + 1, std::end(children),
              std::begin(children) + max_wasted_j);
    children_count--;

    for (int k = 0; k < children_count; k++) {
      if (children_count < RTREE_MIN_CHILDREN_COUNT) {
        break;
      }
      if (nodes->at(new_node)->children_count > RTREE_MIN_CHILDREN_COUNT) {
        break;
      }
      double min_wasted_j = std::numeric_limits<double>::max();
      uint8_t min_m = 0;
      for (uint8_t m = 0; m < children_count; m++) {
        auto wasted_j = j_bb.wastedArea(children[m].bbox);

        if (min_wasted_j > wasted_j) {
          min_m = m;
        }
      }

      nodes->at(new_node)->insert(children[min_m]);
      std::copy(std::begin(children) + min_m + 1, std::end(children),
                std::begin(children) + min_m);
      children_count--;
    }

    return new_node;
  }

  bool insert(const Entry &e) {
    if (children_count > (RTREE_MAX_CHILDREN_COUNT - 1)) {
      return false;
    }

    children[children_count] = e;
    children_count++;

    return true;
  }

  void updateChildBoundingBox(uint64_t cid, BoundingBox bb) {
    auto it = std::find_if(std::begin(children), std::end(children),
                           [cid](auto e) { return e.tuple_id == cid; });
    it->bbox = bb;
  }

  BoundingBox computeBoundingBox() {
    return std::accumulate(std::begin(children), std::end(children),
                           BoundingBox(),
                           [](auto a, auto b) { return a + b.bbox; });
  }
};

// This is written to only ever access one preëxisting node plus
// potentially one new node at a time to think about how this
// will work on a microcontroller where I may only have access
// to one node at a time when searching.
struct RTree {
  int32_t current_node = -1;
  int32_t root = -1;
  std::vector<std::unique_ptr<RTreeNode>> nodes;

  void chooseLeaf(const BoundingBox &bbox) {
    uint64_t next = nodes[root]->myid;

    do {
      loadNode(next);
      next = nodes[current_node]->chooseLeaf(bbox);
    } while (next != nodes[current_node]->myid);
  }

  void loadNode(uint64_t i) { current_node = i; }

  void updateParentBoundingBox() {
    if (nodes[current_node]->hasParent) {
      auto cid = nodes[current_node]->myid;
      auto cbb = nodes[current_node]->computeBoundingBox();
      loadNode(nodes[current_node]->parentid);
      nodes[current_node]->updateChildBoundingBox(cid, cbb);
      updateParentBoundingBox();
    }
  }

  std::vector<uint64_t> find(const Point &p) { return find(BoundingBox{p, p}); }

  std::vector<uint64_t> find(const BoundingBox &bbox) {
    std::vector<uint64_t> tocheck;
    std::vector<uint64_t> tupleids;
    int nodes_checked = 0;
    tocheck.push_back(nodes[root]->myid);
    while (!tocheck.empty()) {
      loadNode(tocheck.back());
      nodes_checked++;
      tocheck.pop_back();
      for (int i = 0; i < nodes[current_node]->children_count; i++) {
        if (nodes[current_node]->children[i].bbox && bbox) {
          if (nodes[current_node]->isLeaf) {
            tupleids.push_back(nodes[current_node]->children[i].tuple_id);
          } else {
            tocheck.push_back(nodes[current_node]->children[i].tuple_id);
          }
        }
      }
    }
    return tupleids;
  }

  void insert(const Entry &e) {
    if (root == -1) {
      auto new_root =
          std::unique_ptr<RTreeNode>(new RTreeNode(e, true, false, 0));
      new_root->myid = nodes.size();
      root = new_root->myid;
      current_node = root;
      nodes.push_back(std::move(new_root));
      return;
    }

    chooseLeaf(e.bbox);
    insertPhase2(e);
  }

  void insertPhase2(const Entry &e) {
    if (!nodes[current_node]->insert(e)) {
      auto new_node = nodes[current_node]->split(&nodes);
      if (e.bbox.wastedArea(nodes[new_node]->computeBoundingBox()) <
          e.bbox.wastedArea(nodes[current_node]->computeBoundingBox())) {
        nodes[new_node]->insert(e);
      } else {
        nodes[current_node]->insert(e);
      }
      auto cid = nodes[current_node]->myid;
      auto cbb = nodes[current_node]->computeBoundingBox();

      if (!nodes[new_node]->isLeaf) {
        auto nnid = nodes[new_node]->myid;
        for (int i = 0; i < nodes[new_node]->children_count; i++) {
          loadNode(nodes[new_node]->children[i].tuple_id);
          nodes[current_node]->parentid = nnid;
        }
        loadNode(cid);
      }

      if (nodes[current_node]->hasParent) {
        auto pid = nodes[current_node]->parentid;
        updateParentBoundingBox();
        loadNode(pid);
        insertPhase2(Entry(nodes[new_node]->computeBoundingBox(),
                           nodes[new_node]->myid));
      } else {
        auto new_root = std::unique_ptr<RTreeNode>(
            new RTreeNode(Entry(cbb, cid), false, false, 0));
        nodes[current_node]->hasParent = true;
        nodes[current_node]->parentid = nodes.size();
        nodes[new_node]->parentid = nodes.size();
        new_root->myid = nodes.size();
        new_root->insert(Entry(nodes[new_node]->computeBoundingBox(),
                               nodes[new_node]->myid));
        root = new_root->myid;
        nodes.push_back(std::move(new_root));
      }
    } else {
      updateParentBoundingBox();
    }
  }

  // serial insertion:
  //   Node Count: 222
  //   Nodes checked: 23
  //
  // SRTLoad:
  //  Node Count: 178
  //  Nodes checked: 7
  void SRTLoad(std::vector<Entry> *entries, const BoundingBox &maximal) {
    std::sort(entries->begin(), entries->end());

    nodes.reserve(ceil(1.25 * entries->size() / RTREE_MAX_CHILDREN_COUNT));

    auto num_tiles = entries->size() / RTREE_MAX_CHILDREN_COUNT;
    auto x_slice = Qs10d21(static_cast<double>(maximal.xrange().length()) /
                           ceil(sqrt(num_tiles)));
    auto next_x = MIN_Qs10d21;

    for (auto e : *entries) {
      bool need_new_node = false;
      if (e.bbox.upperleft.x > next_x) {
        need_new_node = true;
        next_x = next_x + x_slice;
      } else {
        need_new_node = !nodes[current_node]->insert(e);
      }
      if (need_new_node) {
        auto current_nodeb =
            std::unique_ptr<RTreeNode>(new RTreeNode(e, true, false, 0));
        current_nodeb->myid = nodes.size();
        current_node = current_nodeb->myid;
        nodes.push_back(std::move(current_nodeb));
      }
    }
    uint16_t start = 0, end = nodes.size() + 1;
    int32_t next_node = -1;

    do {
      for (uint16_t i = start; i < end; i++) {
        auto bb = nodes[i]->computeBoundingBox();
        auto tuple_id = nodes[i]->myid;

        if (static_cast<int32_t>(i) > next_node) {
          auto current_nodeb = std::unique_ptr<RTreeNode>(
              new RTreeNode(Entry(bb, tuple_id), false, false, 0));
          current_nodeb->myid = nodes.size();
          current_node = current_nodeb->myid;
          nodes.push_back(std::move(current_nodeb));
          next_node = (i - 1) + RTREE_MAX_CHILDREN_COUNT;
        } else {
          nodes[current_node]->insert(Entry(bb, tuple_id));
        }
        nodes[tuple_id]->parentid = nodes[current_node]->myid;
        nodes[tuple_id]->hasParent = true;
      }
      start = end;
      end = nodes.size();
      next_node = -1;
    } while ((end - start) > 1);
    root = current_node;
  }
};
