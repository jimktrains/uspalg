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

struct BoundingBox {
  Point upperleft;
  Point lowerright;

  BoundingBox operator+(const BoundingBox& other) const {
    return BoundingBox{
      Point{std::min(upperleft.x, other.upperleft.x), std::max(upperleft.y, other.upperleft.y)},
      Point{std::max(lowerright.x, other.lowerright.x), std::min(lowerright.y, other.lowerright.y)}
    };
  };

  bool overlaps(const BoundingBox& other) const {
     if ((upperleft % other.upperleft) == Quadrant::NW) {
       return (lowerright % other.upperleft) == Quadrant::SE;
     }
     if ((upperleft % other.lowerright) == Quadrant::SE) {
       return (lowerright % other.lowerright) == Quadrant::NW;
     }
     return false;
  }

  bool operator&&(const BoundingBox& other) const {
    return overlaps(other);
  }

  bool operator&&(const Point& other) const {
    return (upperleft.x <= other.x && lowerright.x >= other.x)
        && (upperleft.y >= other.y && lowerright.y <= other.y);
  }
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

  //  https://wrfranklin.org/Research/Short_Notes/pnpoly.html
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

constexpr static int RTREE_MAX_CHILDREN_COUNT = 11;
constexpr static int RTREE_MIN_CHILDREN_COUNT = 5;

// Guttman, Antomn. "R-Trees - A Dynamic Index Structure for Spatial
// Searching." ACM SIGMOD Record, vol. 14, no. 2, June 1984, pp. 47–57.,
// https://doi.org/10.1145/971697.602266.
struct RTree {
  BoundingBox bbox;
  long tuple_id;
  size_t child_count;
  RTree* children[RTREE_MAX_CHILDREN_COUNT];
  BoundingBox children_bbox[RTREE_MAX_CHILDREN_COUNT];
};

