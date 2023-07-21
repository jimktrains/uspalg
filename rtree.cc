#include <math.h>

enum Quadrant {
  NE,
  NW,
  SW,
  SE,
};

struct Point {
  double x;
  double y;

  Point(double xx, double yy) : x{xx}, y{yy} {}

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
  double start;
  double end;

  bool overlaps(const Range& other) {
    return start <= other.end and end >= other.start;
  };

  bool operator&&(const Range& other) {
    return overlaps(other);
  };
};

struct BoundingBox {
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
    return (lowerright.x - upperleft.x) * (upperleft.y - lowerright.y);
  };

  double wastedArea(const BoundingBox& other) const {
    return (*this + other).area() - this->area() - other.area();
  };
};

struct Polygon {
  Point* points;
  size_t npoints;

  BoundingBox boundingBox() {
    double mx = std::numeric_limits<double>::max();
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

// 12 children allows for a 504-byte RTree object, 8 bytes short of a
// 512-byte block.
constexpr static int RTREE_MAX_CHILDREN_COUNT = 12;
constexpr static int RTREE_MIN_CHILDREN_COUNT =  6;

// Guttman, Antomn. "R-Trees - A Dynamic Index Structure for Spatial
//   Searching." ACM SIGMOD Record, vol. 14, no. 2, June 1984, pp. 47–57.,
//   https://doi.org/10.1145/971697.602266.
struct RTree {
  long tuple_id;
  size_t child_count;
  RTree* parent;
  RTree* children[RTREE_MAX_CHILDREN_COUNT];
  // If we store the boudning box inside the node, then we can load
  // this block from the sd card and figure out which block to descend
  // into without loading anything else.
  BoundingBox children_bbox[RTREE_MAX_CHILDREN_COUNT];
};

