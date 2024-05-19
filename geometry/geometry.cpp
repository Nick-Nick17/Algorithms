#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using PointType = long double;

static const int kSetPrecision = 6;
static const PointType kCmp = 1e-6;
//  static const PointType kInf = 1e9;

int GiveSign(PointType value) {
  if (value < -kCmp) {
    return -1;
  }
  if (value > kCmp) {
    return 1;
  }
  return 0;
}

struct Point {
  PointType x;
  PointType y;

  Point() : x(0), y(0) {}

  Point(PointType first, PointType second) : x(first), y(second) {}

  Point(const Point& begin, const Point& end)
      : x(end.x - begin.x), y(end.y - begin.y) {}

  Point(const Point& point) : x(point.x), y(point.y) {}

  Point& operator=(const Point& point) {
    x = point.x;
    y = point.y;
    return *this;
  }

  Point& operator+=(const Point& value) {
    x += value.x;
    y += value.y;
    return *this;
  }

  Point& operator-=(const Point& value) {
    x -= value.x;
    y -= value.y;
    return *this;
  }

  Point& operator*=(const PointType& value) {
    x *= value;
    y *= value;
    return *this;
  }

  Point& operator/=(const PointType& value) {
    x /= value;
    y /= value;
    return *this;
  }

  PointType Length() const { return std::sqrt(x * x + y * y); }
  PointType SquareLinth() const { return x * x + y * y; }
};

bool operator==(const Point& first, const Point& second) {
  return std::abs(first.x - second.x) < kCmp &&
         std::abs(first.y - second.y) < kCmp;
}

bool operator!=(const Point& first, const Point& second) {
  return !(first == second);
}

Point operator+(Point first, const Point& second) {
  first.x += second.x;
  first.y += second.y;
  return first;
}

Point operator-(Point first, const Point& second) {
  first.x -= second.x;
  first.y -= second.y;
  return first;
}

PointType operator^(const Point& first, const Point& second) {
  return first.x * second.y - first.y * second.x;
}  //  векторное произведение

PointType operator*(const Point& first, const Point& second) {
  return first.x * second.x + first.y * second.y;
}  //  скалярное произведение

PointType TripleSquare(const Point& first, const Point& second) {
  return std::abs(first ^ second) / 2;
}

PointType Angle(const Point& first, const Point& second) {
  return std::acos((first * second) / first.Length() / second.Length());
}

PointType FindPointsDist(Point first, Point second) {
  return (first - second).Length();
}

bool IsPointInTriangle(Point aa, Point bb, Point cc, Point point) {
  PointType s1 = std::abs((bb - aa) ^ (cc - aa));
  PointType s2 = std::abs((aa - point) ^ (bb - point));
  s2 += std::abs((bb - point) ^ (cc - point));
  s2 += std::abs((cc - point) ^ (aa - point));
  return GiveSign(s1 - s2) == 0;
}

struct Line {
  PointType a;
  PointType b;
  PointType c;

  Point line_point_first;
  Point line_point_second;

  Line() {}

  Line(const Point& first, const Point& second)
      : a(first.y - second.y),
        b(second.x - first.x),
        c((-a) * first.x - b * first.y),
        line_point_first(first),
        line_point_second(second) {}

  bool IsPointOnLine(const Point& point) const {
    return std::abs(a * point.x + b * point.y + c) < kCmp;
  }

  Point FindNormal() const { return Point(a, b); }

  Point FindProjection(Point point) const {
    if (IsPointOnLine(point)) {
      return point;
    }
    PointType dist = FindDist(point);
    Point new_normal = FindNormal();
    new_normal /= FindNormal().Length();
    new_normal *= dist;

    Point first = point + new_normal;
    Point second = point - new_normal;
    if (FindDist(first) < FindDist(second)) {
      return first;
    }
    return second;
  }

  PointType FindDist(Point point) const {
    PointType sqr = std::sqrt(a * a + b * b);
    return std::abs(a * point.x + b * point.y + c) / sqr;
  }

  bool IsOnSegment(Point projection) const {
    if (!IsPointOnLine(projection)) {
      return false;
    }
    Point ap(line_point_first, projection);
    Point ab(line_point_first, line_point_second);

    Point bp(line_point_second, projection);
    Point ba(line_point_second, line_point_first);

    PointType m1 = std::max(ab.Length(), ap.Length());
    PointType m2 = std::max(ba.Length(), bp.Length());
    return GiveSign((ap + ab).Length() - m1) >= 0 &&
           GiveSign((bp + ba).Length() - m2) >= 0;
  }

  bool IsOnRay(Point projection) const {
    if (!IsPointOnLine(projection)) {
      return false;
    }
    Point ap(line_point_first, projection);
    Point ab(line_point_first, line_point_second);
    PointType m1 = std::max(ab.Length(), ap.Length());
    return GiveSign((ap + ab).Length() - m1) >= 0;
  }

  PointType FindDistToSegmet(Point point) const {
    Point projection = FindProjection(point);
    if (IsOnSegment(projection)) {
      return FindDist(point);
    }
    Point first(point, line_point_first);
    Point second(point, line_point_second);
    return std::min(first.Length(), second.Length());
  }

  PointType FindDistToRay(Point point) const {
    Point projection = FindProjection(point);
    if (IsOnRay(projection)) {
      return FindDist(point);
    }
    Point first(point, line_point_first);
    return first.Length();
  }
};

struct Polygon {
  std::vector<Point> points;

  Polygon() {}

  Polygon(std::vector<Point>& array) : points(array) {}

  void ReorginizeIfNeed() {
    if (GiveSign((points[2] - points[1]) ^ (points[0] - points[1])) >= 0) {
      return;
    }
    std::reverse(points.begin(), points.end());
  }

  size_t Size() const { return points.size(); }

  bool IsConv() {  //  выпуклость
    int nn = static_cast<int>(points.size());
    int res = 0;
    for (int i = 0; i < nn; ++i) {
      int previous = (i - 1 + nn) % nn;
      int next = (i + 1) % nn;
      Point first(points[previous], points[i]);
      Point second(points[i], points[next]);
      if ((first ^ second) < -kCmp) {
        res |= 1;
      }
      if ((first ^ second) > kCmp) {
        res |= 2;
      }
      if (res == 3) {
        return false;
      }
    }
    return true;
  }

  void OrderPoints() {}
};

bool IsPointsOnDiffSemiplanes(Line line, Point first, Point second) {
  Point ab(line.line_point_first, line.line_point_second);
  Point a1(line.line_point_first, first);
  Point a2(line.line_point_first, second);
  return GiveSign(ab ^ a1) != GiveSign(ab ^ a2);
}

bool IsSegmentsAreCrossing(Line fir, Line sec) {
  return IsPointsOnDiffSemiplanes(fir, sec.line_point_first,
                                  sec.line_point_second) &&
         IsPointsOnDiffSemiplanes(sec, fir.line_point_first,
                                  fir.line_point_second);
}

PointType SegsDistHelper(Line line, Point point) {
  Point proj = line.FindProjection(point);
  if (line.IsOnSegment(proj)) {
    return line.FindDist(point);
  }
  PointType ans1 = (line.line_point_first - point).Length();
  PointType ans2 = (line.line_point_second - point).Length();
  return std::min(ans1, ans2);
}

PointType FindSegmentsDist(Point first_start, Point first_end,
                           Point second_start, Point second_end) {
  Line first(first_start, first_end);
  Line second(second_start, second_end);
  if (second.IsPointOnLine(first_start)) {
    return 0;
  }
  if (second.IsPointOnLine(first_end)) {
    return 0;
  }
  if (first.IsPointOnLine(second_start)) {
    return 0;
  }
  if (first.IsPointOnLine(second_end)) {
    return 0;
  }
  if (IsSegmentsAreCrossing(first, second)) {
    return 0;
  }
  PointType min = SegsDistHelper(first, second_start);
  min = std::min(min, SegsDistHelper(first, second_end));
  min = std::min(min, SegsDistHelper(second, first_start));
  min = std::min(min, SegsDistHelper(second, first_end));
  return min;
}

void CinPoint(Point& point) { std::cin >> point.x >> point.y; }

bool CmpY(Point first, Point second) {
  return (GiveSign(first.y - second.y) < 0) ||
         (GiveSign(first.y - second.y) == 0 &&
          GiveSign(first.x - second.x) < 0);
}

bool CmpX(Point first, Point second) {
  return (GiveSign(first.x - second.x) < 0) ||
         (GiveSign(first.x - second.x) == 0 &&
          GiveSign(first.y - second.y) < 0);
}

void ReorderForMinkPolygYY(Polygon& pol) {
  size_t pos = 0;
  for (size_t i = 1; i < pol.Size(); ++i) {
    if (CmpY(pol.points[i], pol.points[pos])) {
      pos = i;
    }
  }
  std::rotate(pol.points.begin(), pol.points.begin() + pos, pol.points.end());
}

void ReorderForMinkPolygXX(Polygon& pol) {
  size_t pos = 0;
  for (size_t i = 1; i < pol.Size(); ++i) {
    if (CmpX(pol.points[i], pol.points[pos])) {
      pos = i;
    }
  }
  std::rotate(pol.points.begin(), pol.points.begin() + pos, pol.points.end());
}

Polygon GiveMinkovskiiSumm(Polygon first, Polygon second) {
  ReorderForMinkPolygYY(first);
  ReorderForMinkPolygYY(second);
  first.points.push_back(first.points[0]);
  first.points.push_back(first.points[1]);
  second.points.push_back(second.points[0]);
  second.points.push_back(second.points[1]);
  Polygon res;
  size_t ii = 0;
  size_t jj = 0;
  while (ii + 2 < first.Size() || jj + 2 < second.Size()) {
    res.points.push_back(first.points[ii] + second.points[jj]);
    auto cross = ((first.points[ii + 1] - first.points[ii]) ^
                  (second.points[jj + 1] - second.points[jj]));
    if (GiveSign(cross) >= 0 && ii + 2 < first.Size()) {
      ++ii;
    }
    if (GiveSign(cross) <= 0 && jj + 2 < second.Size()) {
      ++jj;
    }
  }
  return res;
}

struct FastCheckPointInPolygon {
  Polygon pol;
  std::vector<Point> seq;
  Point translation;

  FastCheckPointInPolygon(const Polygon& pp) : pol(pp) {
    ReorderForMinkPolygXX(pol);
    seq.resize(pol.Size() - 1);
    for (size_t i = 0; i < seq.size(); ++i) {
      seq[i] = pol.points[i + 1] - pol.points[0];
    }
    translation = pol.points[0];
  }

  bool IsPointInPolygon(Point point) {
    point -= translation;

    if (GiveSign(seq[0] ^ point) != 0 &&
        GiveSign(seq[0] ^ point) != GiveSign(seq[0] ^ seq[seq.size() - 1])) {
      return false;
    }
    if (GiveSign(seq[seq.size() - 1] ^ point) != 0 &&
        GiveSign(seq[seq.size() - 1] ^ point) !=
            GiveSign(seq[seq.size() - 1] ^ seq[0])) {
      return false;
    }
    if (GiveSign(seq[0] ^ point) == 0) {
      return GiveSign(seq[0].Length() - point.Length()) >= 0;
    }

    int left = 0;
    int right = static_cast<int>(seq.size()) - 1;
    while (right - left > 1) {
      int mid = (left + right) / 2;
      int pos = mid;
      if (GiveSign(seq[pos] ^ point) >= 0) {
        left = mid;
      } else {
        right = mid;
      }
    }
    int pos = left;
    Point st(0, 0);
    return IsPointInTriangle(seq[pos], seq[pos + 1], st, point);
  }
};

int main() {
  int nn;
  std::cin >> nn;

  Polygon first;
  first.points.resize(nn);
  for (int i = 0; i < nn; ++i) {
    CinPoint(first.points[i]);
  }

  std::cin >> nn;
  Polygon second;
  second.points.resize(nn);
  for (int i = 0; i < nn; ++i) {
    CinPoint(second.points[i]);
  }

  std::cin >> nn;
  Polygon third;
  third.points.resize(nn);
  for (int i = 0; i < nn; ++i) {
    CinPoint(third.points[i]);
  }

  Polygon tmp = GiveMinkovskiiSumm(first, second);
  Polygon mink = GiveMinkovskiiSumm(tmp, third);
  mink.ReorginizeIfNeed();
  FastCheckPointInPolygon is_point_in(mink);

  std::cout << std::fixed << std::setprecision(kSetPrecision);
  int mm;
  std::cin >> mm;
  for (int i = 0; i < mm; ++i) {
    Point point;
    std::cin >> point.x >> point.y;
    point.x *= 3;
    point.y *= 3;
    if (is_point_in.IsPointInPolygon(point)) {
      std::cout << "YES" << '\n';
    } else {
      std::cout << "NO" << '\n';
    }
  }
  return 0;
}
