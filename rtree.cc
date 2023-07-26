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
  Entry(BoundingBox &bb, uint64_t &tuple_id) : bb{bb}, tuple_id{tuple_id} {}

  BoundingBox bb;
  uint64_t tuple_id;

  // Sorting like this for the STR algorithm.
  bool operator<(const Entry &other) const {
    if (bb.center().x < other.bb.center().x) {
      return true;
    }
    if (bb.center().x > other.bb.center().x) {
      return false;
    }
    return bb.center().y < other.bb.center().y;
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

struct RTreeNode;

std::vector<RTreeNode *> nodes;

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
  RTreeNode(const BoundingBox &bb, uint64_t i, bool isLeaf, bool hasParent,
            uint64_t pid)
      : isLeaf{isLeaf}, hasParent{hasParent}, parentid{pid} {
    children_count = 1;
    children_bbox[0] = bb;
    child_tuple_or_node_id[0] = i;
  }

  uint8_t children_count;
  uint64_t myid;
  bool isLeaf;
  bool hasParent;
  uint64_t parentid;
  /*
   * If we store the boudning box inside the node, then we can load
   * this block from the sd card and figure out which block to descend
   * into without loading anything else.
   */
  BoundingBox children_bbox[RTREE_MAX_CHILDREN_COUNT];
  uint64_t child_tuple_or_node_id[RTREE_MAX_CHILDREN_COUNT];

  uint64_t chooseLeaf(const BoundingBox &bbox) {
    if (isLeaf) {
      return myid;
    }

    double min_waste = children_bbox[0].wastedArea(bbox);
    double min_waste_area = children_bbox[0].area();
    uint8_t min_waste_i = 0;
    for (uint8_t i = 1; i < children_count; i++) {
      auto wasted = children_bbox[i].wastedArea(bbox);
      auto area = children_bbox[i].area();
      if ((wasted < min_waste) ||
          (wasted == min_waste && area < min_waste_area)) {
        min_waste = wasted;
        min_waste_area = area;
        min_waste_i = i;
      }
    }

    return child_tuple_or_node_id[min_waste_i];
  }

  RTreeNode *split() {
    double max_wasted_space = 0;
    // uint8_t max_wasted_i = 0;
    uint8_t max_wasted_j = 0;
    uint8_t i, j;
    for (i = 0; i < (children_count - 1); i++) {
      for (j = i + 1; j < children_count; j++) {
        auto wasted_area = children_bbox[i].wastedArea(children_bbox[j]);
        if (wasted_area > max_wasted_space) {
          max_wasted_space = wasted_area;
          // max_wasted_i = i;
          max_wasted_j = j;
        }
      }
    }

    auto new_node = new RTreeNode(children_bbox[max_wasted_j],
                                  child_tuple_or_node_id[max_wasted_j], isLeaf,
                                  true, parentid);
    // Not thread-safe.
    new_node->myid = nodes.size();
    nodes.push_back(new_node);

    // auto i_bb = children_bbox[max_wasted_i];
    auto j_bb = children_bbox[max_wasted_j];

    for (int k = max_wasted_j + 1; k < children_count; k++) {
      children_bbox[k - 1] = children_bbox[k];
      child_tuple_or_node_id[k - 1] = child_tuple_or_node_id[k];
    }
    children_count--;

    for (int k = 0; k < children_count; k++) {
      if (children_count < RTREE_MIN_CHILDREN_COUNT) {
        break;
      }
      if (new_node->children_count > RTREE_MIN_CHILDREN_COUNT) {
        break;
      }
      double min_wasted_j = std::numeric_limits<double>::max();
      uint8_t min_m = 0;
      for (uint8_t m = 0; m < children_count; m++) {
        auto wasted_j = j_bb.wastedArea(children_bbox[m]);

        if (min_wasted_j > wasted_j) {
          min_m = m;
        }
      }

      new_node->insert(children_bbox[min_m], child_tuple_or_node_id[min_m]);
      for (uint8_t l = min_m + 1; l < children_count; l++) {
        children_bbox[l - 1] = children_bbox[l];
        child_tuple_or_node_id[l - 1] = child_tuple_or_node_id[l];
      }
      children_count--;
    }

    return new_node;
  }

  bool insert(BoundingBox bbox, uint64_t tuple_id) {
    if (children_count > (RTREE_MAX_CHILDREN_COUNT - 1)) {
      return false;
    }
    children_bbox[children_count] = bbox;
    child_tuple_or_node_id[children_count] = tuple_id;

    children_count++;

    return true;
  }

  void updateChildBoundingBox(uint64_t cid, BoundingBox bb) {
    for (uint8_t i = 0; i < children_count; i++) {
      if (child_tuple_or_node_id[i] == cid) {
        children_bbox[i] = bb;
        break;
      }
    }
  }

  BoundingBox computeBoundingBox() {
    BoundingBox bb = children_bbox[0];
    for (int16_t i = 1; i < children_count; i++) {
      bb = bb + children_bbox[i];
    }
    return bb;
  }
};

// This is written to only ever access one preëxisting node plus
// potentially one new node at a time to think about how this
// will work on a microcontroller where I may only have access
// to one node at a time when searching.
struct RTree {
  RTreeNode *current_node = nullptr;
  RTreeNode *root = nullptr;

  void chooseLeaf(const BoundingBox &bbox) {
    uint64_t next = root->myid;

    do {
      loadNode(next);
      next = current_node->chooseLeaf(bbox);
    } while (next != current_node->myid);
  }

  void loadNode(uint64_t i) { current_node = nodes[i]; }

  void updateParentBoundingBox() {
    if (current_node->hasParent) {
      auto cid = current_node->myid;
      auto cbb = current_node->computeBoundingBox();
      loadNode(current_node->parentid);
      current_node->updateChildBoundingBox(cid, cbb);
      updateParentBoundingBox();
    }
  }

  std::vector<uint64_t> find(const Point &p) { return find(BoundingBox{p, p}); }

  std::vector<uint64_t> find(const BoundingBox &bbox) {
    std::vector<uint64_t> tocheck;
    std::vector<uint64_t> tupleids;
    int nodes_checked = 0;
    tocheck.push_back(root->myid);
    while (!tocheck.empty()) {
      loadNode(tocheck.back());
      nodes_checked++;
      tocheck.pop_back();
      for (int i = 0; i < current_node->children_count; i++) {
        if (current_node->children_bbox[i] && bbox) {
          if (current_node->isLeaf) {
            tupleids.push_back(current_node->child_tuple_or_node_id[i]);
          } else {
            tocheck.push_back(current_node->child_tuple_or_node_id[i]);
          }
        }
      }
    }
    std::cout << "Nodes checked: " << nodes_checked << std::endl;
    return tupleids;
  }

  void insert(const BoundingBox &bbox, uint64_t tuple_id) {
    if (root == nullptr) {
      auto new_root = new RTreeNode(bbox, tuple_id, true, false, 0);
      new_root->myid = nodes.size();
      nodes.push_back(new_root);
      root = new_root;
      current_node = new_root;
      return;
    }

    chooseLeaf(bbox);
    insertPhase2(bbox, tuple_id);
  }

  void insertPhase2(BoundingBox bbox, uint64_t tuple_id) {
    if (!current_node->insert(bbox, tuple_id)) {
      auto new_node = current_node->split();
      if (bbox.wastedArea(new_node->computeBoundingBox()) <
          bbox.wastedArea(current_node->computeBoundingBox())) {
        new_node->insert(bbox, tuple_id);
      } else {
        current_node->insert(bbox, tuple_id);
      }
      auto cid = current_node->myid;
      auto cbb = current_node->computeBoundingBox();

      if (!new_node->isLeaf) {
        auto nnid = new_node->myid;
        for (int i = 0; i < new_node->children_count; i++) {
          loadNode(new_node->child_tuple_or_node_id[i]);
          current_node->parentid = nnid;
        }
        loadNode(cid);
      }

      if (current_node->hasParent) {
        auto pid = current_node->parentid;
        updateParentBoundingBox();
        loadNode(pid);
        insertPhase2(new_node->computeBoundingBox(), new_node->myid);
      } else {
        auto new_root = new RTreeNode(cbb, cid, false, false, 0);
        current_node->hasParent = true;
        current_node->parentid = nodes.size();
        new_node->parentid = nodes.size();
        new_root->myid = nodes.size();
        nodes.push_back(new_root);

        new_root->insert(new_node->computeBoundingBox(), new_node->myid);
        root = new_root;
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

    RTreeNode *current_node;
    for (auto e : *entries) {
      bool need_new_node = false;
      if (e.bb.upperleft.x > next_x) {
        need_new_node = true;
        next_x = next_x + x_slice;
      } else {
        need_new_node = !current_node->insert(e.bb, e.tuple_id);
      }
      if (need_new_node) {
        current_node = new RTreeNode(e.bb, e.tuple_id, true, false, 0);
        current_node->myid = nodes.size();
        nodes.push_back(current_node);
      }
    }
    uint16_t start = 0, end = nodes.size() + 1;
    int32_t next_node = -1;

    do {
      for (uint16_t i = start; i < end; i++) {
        auto bb = nodes[i]->computeBoundingBox();
        auto tuple_id = nodes[i]->myid;

        if (static_cast<int32_t>(i) > next_node) {
          current_node = new RTreeNode(bb, tuple_id, false, false, 0);
          current_node->myid = nodes.size();
          nodes.push_back(current_node);
          next_node = (i - 1) + RTREE_MAX_CHILDREN_COUNT;
        } else {
          current_node->insert(bb, tuple_id);
        }
        nodes[tuple_id]->parentid = current_node->myid;
        nodes[tuple_id]->hasParent = true;
      }
      start = end;
      end = nodes.size();
      next_node = -1;
    } while ((end - start) > 1);
    root = current_node;
  }
};
