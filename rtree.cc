#include <math.h>
#include <limits>
#include <numeric>
#include "fixpoint.cc"

enum Quadrant {
  NE,
  NW,
  SW,
  SE,
};

struct Point {
  Qs10d21 x;
  Qs10d21 y;

  Point(Qs10d21 xx, Qs10d21 yy) : x{xx}, y{yy} {}

  Quadrant operator%(const Point& other) const {
    if (x < other.x) {
      if (y < other.y) {
        return SW;
      }
      return NW;
    } else {
      if (y < other.y) {
        return SE;
      }
      return NE;
    }
  };
};

struct Range {
  Qs10d21 start;
  Qs10d21 end;

  bool overlaps(const Range& other) {
    return start <= other.end and end >= other.start;
  };

  bool operator&&(const Range& other) {
    return overlaps(other);
  };
};

struct BoundingBox {
  BoundingBox()
    : upperleft{Point{Qs10d21(0.), Qs10d21(0.)}}
    , lowerright{Point{Qs10d21(0.), Qs10d21(0.)}}
  {};

  BoundingBox(Point ul, Point lr)
    : upperleft{ul}
    , lowerright{lr}
  {};

  Point upperleft;
  Point lowerright;

  BoundingBox operator+(const BoundingBox& other) const {
    return BoundingBox{
      Point{std::min(upperleft.x, other.upperleft.x), std::max(upperleft.y, other.upperleft.y)},
      Point{std::max(lowerright.x, other.lowerright.x), std::min(lowerright.y, other.lowerright.y)}
    };
  };

  Range xrange() const {
    return Range{upperleft.x, lowerright.x};
  };

  Range yrange() const {
    return Range{lowerright.y, upperleft.y};
  };

  bool overlaps(const BoundingBox& other) const {
    return (xrange() && other.xrange()) && (yrange() && other.yrange());
  };

  bool operator&&(const BoundingBox& other) const {
    return overlaps(other);
  };

  bool operator&&(const Point& other) const {
    return (upperleft.x <= other.x && lowerright.x >= other.x)
        && (upperleft.y >= other.y && lowerright.y <= other.y);
  };

  double area() const {
    return ((double)(lowerright.x - upperleft.x)) * ((double)(upperleft.y - lowerright.y));
  };

  double wastedArea(const BoundingBox& other) const {
    return (*this + other).area() - this->area() - other.area();
  };
};

struct Polygon {
  Point* points;
  size_t npoints;

  BoundingBox boundingBox() {
    Qs10d21 mx = MAX_Qs10d21;
    BoundingBox bb{
      Point{mx, -mx},
      Point{-mx, mx}
    };

    for (int i = 0; i < npoints; i++) {
      bb.upperleft  = Point{std::min(bb.upperleft.x, points[i].x), std::max(bb.upperleft.y, points[i].y)};
      bb.lowerright = Point{std::max(bb.lowerright.x, points[i].x), std::min(bb.lowerright.y, points[i].y)};
    }

    return bb;
  }

  // > The Method
  // > I run a semi-infinite ray horizontally (increasing x, fixed y) out from
  // > the test point, and count how many edges it crosses. At each crossing, the
  // > ray switches between inside and outside. This is called the Jordan curve
  // > theorem.
  // >
  // > The case of the ray going thru a vertex is handled correctly
  // > via a careful selection of inequalities. Don't mess with this code unless
  // > you're familiar with the idea of Simulation of Simplicity. This pretends
  // > to shift the ray infinitesimally down so that it either clearly
  // > intersects, or clearly doesn't touch. Since this is merely a conceptual,
  // > infinitesimal, shift, it never creates an intersection that didn't exist
  // > before, and never destroys an intersection that clearly existed before.
  // >
  // > The ray is tested against each edge thus:
  // >
  // >     Is the point in the half-plane to the left of the extended edge? and
  // >     Is the point's Y coordinate within the edge's Y-range?
  // >
  // > Handling endpoints here is tricky.
  // https://wrfranklin.org/Research/Short_Notes/pnpoly.html
  bool contains(const Point& p) const {
    bool c = false;
    for (int i = 0, j = npoints-1; i < npoints; j = i++) {
      if ( ((points[i].y>p.y) != (points[j].y>p.y)) &&
          (p.x < (points[j].x-points[i].x) * (p.y-points[i].y) / (points[j].y-points[i].y) + points[i].x) )
        c = !c;
    }
    return c;
  };

  bool operator&&(const Point& p) const {
    return contains(p);
  };
};

constexpr static int RTREE_MAX_CHILDREN_COUNT = 15;
constexpr static int RTREE_MIN_CHILDREN_COUNT =  RTREE_MAX_CHILDREN_COUNT / 2;

struct RTreeNode;

#include <vector>
std::vector<RTreeNode*> nodes;

// Guttman, Antomn. "R-Trees - A Dynamic Index Structure for Spatial
//   Searching." ACM SIGMOD Record, vol. 14, no. 2, June 1984, pp. 47–57.,
//   https://doi.org/10.1145/971697.602266.
struct RTreeNode {
  RTreeNode(BoundingBox& bb, uint64_t i, bool isLeaf, bool hasParent, uint64_t pid)
    : isLeaf{isLeaf}
    , hasParent{hasParent}
    , parentid{pid}
  {
    children_count = 1;
    children_bbox[0] = bb;
    child_tuple_or_node_id[0] = i;
  }

  uint8_t children_count;
  uint64_t parentid;
  uint64_t myid;
  bool isLeaf;
  bool hasParent;
  // If we store the boudning box inside the node, then we can load
  // this block from the sd card and figure out which block to descend
  // into without loading anything else.
  BoundingBox children_bbox[RTREE_MAX_CHILDREN_COUNT];
  uint64_t child_tuple_or_node_id[RTREE_MAX_CHILDREN_COUNT];

  uint64_t chooseLeaf(BoundingBox& bbox) {
    if (isLeaf) {
      return myid;
    }

    double max_waste = children_bbox[0].wastedArea(bbox);
    double max_waste_area = children_bbox[0].area();
    size_t max_waste_i = 0;
    for (size_t i = 1; i < children_count; i++) {
      auto wasted = children_bbox[i].wastedArea(bbox);
      auto area   = children_bbox[i].area();
      if ((wasted > max_waste) ||
          (wasted == max_waste && area < max_waste_area)) {
        max_waste = wasted;
        max_waste_area = area;
        max_waste_i = i;
      }
    }

    return child_tuple_or_node_id[max_waste_i];
  };

  RTreeNode* split() {
    // TODO: Do the actual computation to do the split correctly.
    children_count--;
    RTreeNode *new_node = new RTreeNode(children_bbox[children_count], child_tuple_or_node_id[children_count], isLeaf, true, parentid);
    for (children_count--; children_count > RTREE_MIN_CHILDREN_COUNT; children_count--) {
      new_node->insert(children_bbox[children_count], child_tuple_or_node_id[children_count]);
      if (child_tuple_or_node_id[children_count] == 283) {
        auto b = children_bbox[children_count];
            std::cout << "Split " << (isLeaf ? "L" : "I") << " [" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
      }
    }

    // Not thread-safe.
    new_node->myid = nodes.size();
    nodes.push_back(new_node);

    return new_node;
  };

  bool insert(BoundingBox bbox, long tuple_id) {
      if (tuple_id == 283) {
        auto b = bbox;
            std::cout << "Attempted Inser [" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
      }
    if (children_count > (RTREE_MAX_CHILDREN_COUNT - 1)) {
      if (tuple_id == 283) {
      std::cout << "  Failed!" << std::endl;
      }
      return false;
    }
    children_bbox[children_count] = bbox;
    child_tuple_or_node_id[children_count] = tuple_id;

      if (child_tuple_or_node_id[children_count] == 283) {
        auto b = children_bbox[children_count];
            std::cout << "  Insert [" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
      }

    children_count++;

    return true;
  };

  void updateBoundingBox(uint64_t cid, BoundingBox bb) {
    for (int i = 0; i < children_count; i++) {
      if (child_tuple_or_node_id[i] == cid) {
        children_bbox[i] = bb;
        break;
      }
    }
  };

  BoundingBox computeBoundingBox() {
    BoundingBox bb = children_bbox[0];
    for (auto cbb : children_bbox) {
      bb = bb + cbb;
    }
    return bb;
  };
};


struct RTree {
  RTreeNode* current_node = nullptr;
  RTreeNode* root = nullptr;

  void chooseLeaf(BoundingBox& bbox) {
    uint64_t next = root->myid;

    do {
      loadNode(next);
      next = current_node->chooseLeaf(bbox);
    } while(next != current_node->myid );
  };

  void loadNode(uint64_t i) {
    current_node = nodes[i];
  };

  void updateParentBoundingBox() {
    if (current_node->hasParent) {
      auto cid = current_node->myid;
      auto cbb = current_node->computeBoundingBox();
      loadNode(current_node->parentid);
      current_node->updateBoundingBox(cid, cbb);
      updateParentBoundingBox();
    }
  };

  std::vector<uint64_t> find(BoundingBox& bbox) {
    std::vector<uint64_t> tocheck;
    std::vector<uint64_t> tupleids;
    tocheck.push_back(root->myid);
    while (!tocheck.empty()) {
      loadNode(tocheck.back());
      std::cout << "Looking at node " << current_node->myid << std::endl;
      tocheck.pop_back();
      for(int i = 0; i < current_node->children_count; i++) {
        if (current_node->children_bbox[i] && bbox) {
          if (current_node->isLeaf) {
            tupleids.push_back(current_node->child_tuple_or_node_id[i]);
            std::cout << "  Found a leaf " << current_node->child_tuple_or_node_id[i]<< " " << i << std::endl;
            auto b = current_node->children_bbox[i];
            std::cout << "[" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
          } else {
            tocheck.push_back(current_node->child_tuple_or_node_id[i]);
            std::cout << "  Will explore " << current_node->child_tuple_or_node_id[i]<< std::endl;
          }
        }
      }
    }
    return tupleids;
  };

  void insert(BoundingBox bbox, long tuple_id) {
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

  void insertPhase2(BoundingBox bbox, long tuple_id) {
      if (tuple_id == 283) {
        auto b = bbox;
            std::cout << "RTree Attempted Insert " << (current_node->isLeaf ? "L" : "I")  << " [" << (double)b.upperleft.x << " " << (double)b.upperleft.y << ", " << (double)b.lowerright.x << " " << (double)b.lowerright.y << "]" << std::endl;
      }
    if (!current_node->insert(bbox, tuple_id)) {
      if (tuple_id == 283) {
        std::cout << "Splitting" << std::endl;
      }
      auto new_node = current_node->split();
      current_node->insert(bbox, tuple_id);
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
        loadNode(current_node->parentid);
        current_node->updateBoundingBox(cid, cbb);
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
  };
};

