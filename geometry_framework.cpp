/*
 * Computational Geometry & Mesh Generation Framework
 * High-performance C++ implementation for industrial applications
 * Optimized for speed and memory efficiency
 */

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include <chrono>
#include <set>
#include <map>
#include <queue>
#include <numeric>
#include <limits>

using namespace std;
using namespace std::chrono;

// ============================================================================
// CORE DATA STRUCTURES
// ============================================================================

struct Point {
    double x, y;
    
    Point() : x(0), y(0) {}
    Point(double x_, double y_) : x(x_), y(y_) {}
    
    double distance_to(const Point& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        return sqrt(dx * dx + dy * dy);
    }
    
    bool operator<(const Point& other) const {
        return (x < other.x) || (x == other.x && y < other.y);
    }
    
    bool operator==(const Point& other) const {
        return abs(x - other.x) < 1e-9 && abs(y - other.y) < 1e-9;
    }
};

struct Edge {
    Point p1, p2;
    
    Edge(const Point& a, const Point& b) : p1(a), p2(b) {}
    
    bool operator==(const Edge& other) const {
        return (p1 == other.p1 && p2 == other.p2) ||
               (p1 == other.p2 && p2 == other.p1);
    }
};

struct Triangle {
    Point p1, p2, p3;
    
    Triangle(const Point& a, const Point& b, const Point& c) 
        : p1(a), p2(b), p3(c) {}
    
    Point circumcenter() const {
        double ax = p1.x, ay = p1.y;
        double bx = p2.x, by = p2.y;
        double cx = p3.x, cy = p3.y;
        
        double d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        if (abs(d) < 1e-9) {
            return Point((ax + bx + cx) / 3, (ay + by + cy) / 3);
        }
        
        double ux = ((ax*ax + ay*ay) * (by - cy) + 
                     (bx*bx + by*by) * (cy - ay) + 
                     (cx*cx + cy*cy) * (ay - by)) / d;
        double uy = ((ax*ax + ay*ay) * (cx - bx) + 
                     (bx*bx + by*by) * (ax - cx) + 
                     (cx*cx + cy*cy) * (bx - ax)) / d;
        
        return Point(ux, uy);
    }
    
    double circumradius() const {
        Point center = circumcenter();
        return center.distance_to(p1);
    }
    
    bool circumcircle_contains(const Point& p) const {
        Point center = circumcenter();
        double radius = circumradius();
        return center.distance_to(p) < radius + 1e-9;
    }
    
    bool contains_vertex(const Point& p) const {
        return p == p1 || p == p2 || p == p3;
    }
};

// ============================================================================
// CONVEX HULL ALGORITHMS
// ============================================================================

class ConvexHullSolver {
public:
    static double cross_product(const Point& O, const Point& A, const Point& B) {
        return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
    }
    
    // Graham Scan - O(n log n)
    static vector<Point> graham_scan(vector<Point> points) {
        int n = points.size();
        if (n < 3) return points;
        
        // Find lowest point
        int lowest = 0;
        for (int i = 1; i < n; i++) {
            if (points[i].y < points[lowest].y || 
                (points[i].y == points[lowest].y && points[i].x < points[lowest].x)) {
                lowest = i;
            }
        }
        swap(points[0], points[lowest]);
        Point pivot = points[0];
        
        // Sort by polar angle
        sort(points.begin() + 1, points.end(), 
             [&pivot](const Point& a, const Point& b) {
                 double cp = cross_product(pivot, a, b);
                 if (abs(cp) < 1e-9) {
                     return pivot.distance_to(a) < pivot.distance_to(b);
                 }
                 return cp > 0;
             });
        
        // Build hull
        vector<Point> hull;
        hull.push_back(points[0]);
        
        for (int i = 1; i < n; i++) {
            while (hull.size() > 1 && 
                   cross_product(hull[hull.size()-2], hull[hull.size()-1], points[i]) <= 0) {
                hull.pop_back();
            }
            hull.push_back(points[i]);
        }
        
        return hull;
    }
    
    // Quickhull - O(n log n) average case
    static vector<Point> quickhull(const vector<Point>& points) {
        if (points.size() < 3) return points;
        
        vector<Point> hull;
        
        // Find extreme points
        auto minmax_x = minmax_element(points.begin(), points.end(),
                                       [](const Point& a, const Point& b) { return a.x < b.x; });
        Point left = *minmax_x.first;
        Point right = *minmax_x.second;
        
        hull.push_back(left);
        
        // Split points into upper and lower sets
        vector<Point> upper, lower;
        for (const auto& p : points) {
            if (p == left || p == right) continue;
            if (cross_product(left, right, p) > 0) {
                upper.push_back(p);
            } else {
                lower.push_back(p);
            }
        }
        
        find_hull(hull, left, right, upper);
        hull.push_back(right);
        find_hull(hull, right, left, lower);
        
        return hull;
    }

private:
    static void find_hull(vector<Point>& hull, const Point& p1, const Point& p2, 
                         const vector<Point>& points) {
        if (points.empty()) return;
        
        // Find furthest point
        double max_dist = 0;
        int furthest = 0;
        for (size_t i = 0; i < points.size(); i++) {
            double dist = abs(cross_product(p1, p2, points[i]));
            if (dist > max_dist) {
                max_dist = dist;
                furthest = i;
            }
        }
        
        Point p_far = points[furthest];
        
        // Split remaining points
        vector<Point> left, right;
        for (const auto& p : points) {
            if (p == p_far) continue;
            if (cross_product(p1, p_far, p) > 0) {
                left.push_back(p);
            } else if (cross_product(p_far, p2, p) > 0) {
                right.push_back(p);
            }
        }
        
        find_hull(hull, p1, p_far, left);
        hull.push_back(p_far);
        find_hull(hull, p_far, p2, right);
    }
};

// ============================================================================
// DELAUNAY TRIANGULATION
// ============================================================================

class DelaunayTriangulator {
private:
    vector<Point> points;
    vector<Triangle> triangles;
    
public:
    DelaunayTriangulator(const vector<Point>& pts) : points(pts) {}
    
    vector<Triangle> triangulate() {
        if (points.size() < 3) return {};
        
        // Create super-triangle
        double min_x = points[0].x, max_x = points[0].x;
        double min_y = points[0].y, max_y = points[0].y;
        
        for (const auto& p : points) {
            min_x = min(min_x, p.x);
            max_x = max(max_x, p.x);
            min_y = min(min_y, p.y);
            max_y = max(max_y, p.y);
        }
        
        double dx = max_x - min_x;
        double dy = max_y - min_y;
        double delta_max = max(dx, dy);
        double mid_x = (min_x + max_x) / 2;
        double mid_y = (min_y + max_y) / 2;
        
        Point p1(mid_x - 20 * delta_max, mid_y - delta_max);
        Point p2(mid_x, mid_y + 20 * delta_max);
        Point p3(mid_x + 20 * delta_max, mid_y - delta_max);
        
        triangles.clear();
        triangles.push_back(Triangle(p1, p2, p3));
        
        // Add points incrementally
        for (const auto& point : points) {
            vector<Triangle> bad_triangles;
            
            // Find bad triangles
            for (const auto& tri : triangles) {
                if (tri.circumcircle_contains(point)) {
                    bad_triangles.push_back(tri);
                }
            }
            
            // Find polygon boundary
            vector<Edge> polygon;
            for (const auto& tri : bad_triangles) {
                vector<Edge> edges = {
                    Edge(tri.p1, tri.p2),
                    Edge(tri.p2, tri.p3),
                    Edge(tri.p3, tri.p1)
                };
                
                for (const auto& edge : edges) {
                    bool shared = false;
                    for (const auto& other_tri : bad_triangles) {
                        if (tri.p1 == other_tri.p1 && tri.p2 == other_tri.p2 && 
                            tri.p3 == other_tri.p3) continue;
                        
                        vector<Edge> other_edges = {
                            Edge(other_tri.p1, other_tri.p2),
                            Edge(other_tri.p2, other_tri.p3),
                            Edge(other_tri.p3, other_tri.p1)
                        };
                        
                        for (const auto& other_edge : other_edges) {
                            if (edge == other_edge) {
                                shared = true;
                                break;
                            }
                        }
                        if (shared) break;
                    }
                    
                    if (!shared) {
                        polygon.push_back(edge);
                    }
                }
            }
            
            // Remove bad triangles
            triangles.erase(
                remove_if(triangles.begin(), triangles.end(),
                    [&bad_triangles](const Triangle& tri) {
                        return find_if(bad_triangles.begin(), bad_triangles.end(),
                            [&tri](const Triangle& bt) {
                                return tri.p1 == bt.p1 && tri.p2 == bt.p2 && tri.p3 == bt.p3;
                            }) != bad_triangles.end();
                    }),
                triangles.end()
            );
            
            // Re-triangulate
            for (const auto& edge : polygon) {
                triangles.push_back(Triangle(edge.p1, edge.p2, point));
            }
        }
        
        // Remove triangles with super-triangle vertices
        triangles.erase(
            remove_if(triangles.begin(), triangles.end(),
                [&p1, &p2, &p3](const Triangle& tri) {
                    return tri.contains_vertex(p1) || 
                           tri.contains_vertex(p2) || 
                           tri.contains_vertex(p3);
                }),
            triangles.end()
        );
        
        return triangles;
    }
    
    vector<Triangle> get_triangles() const { return triangles; }
};

// ============================================================================
// SPATIAL INDEXING - QUADTREE
// ============================================================================

class QuadTree {
private:
    struct Boundary {
        double min_x, min_y, max_x, max_y;
        
        bool contains(const Point& p) const {
            return p.x >= min_x && p.x <= max_x && 
                   p.y >= min_y && p.y <= max_y;
        }
        
        bool intersects(const Boundary& other) const {
            return !(max_x < other.min_x || other.max_x < min_x ||
                     max_y < other.min_y || other.max_y < min_y);
        }
    };
    
    Boundary boundary;
    int capacity;
    vector<Point> points;
    bool divided;
    unique_ptr<QuadTree> nw, ne, sw, se;
    
    void subdivide() {
        double mid_x = (boundary.min_x + boundary.max_x) / 2;
        double mid_y = (boundary.min_y + boundary.max_y) / 2;
        
        nw = make_unique<QuadTree>(Boundary{boundary.min_x, mid_y, mid_x, boundary.max_y}, capacity);
        ne = make_unique<QuadTree>(Boundary{mid_x, mid_y, boundary.max_x, boundary.max_y}, capacity);
        sw = make_unique<QuadTree>(Boundary{boundary.min_x, boundary.min_y, mid_x, mid_y}, capacity);
        se = make_unique<QuadTree>(Boundary{mid_x, boundary.min_y, boundary.max_x, mid_y}, capacity);
        
        divided = true;
        
        for (const auto& p : points) {
            nw->insert(p) || ne->insert(p) || sw->insert(p) || se->insert(p);
        }
        points.clear();
    }

public:
    QuadTree(Boundary b, int cap = 4) : boundary(b), capacity(cap), divided(false) {}
    
    bool insert(const Point& p) {
        if (!boundary.contains(p)) return false;
        
        if (points.size() < capacity && !divided) {
            points.push_back(p);
            return true;
        }
        
        if (!divided) subdivide();
        
        return nw->insert(p) || ne->insert(p) || sw->insert(p) || se->insert(p);
    }
    
    vector<Point> query(const Boundary& range) const {
        vector<Point> found;
        
        if (!boundary.intersects(range)) return found;
        
        for (const auto& p : points) {
            if (range.contains(p)) {
                found.push_back(p);
            }
        }
        
        if (divided) {
            auto nw_found = nw->query(range);
            auto ne_found = ne->query(range);
            auto sw_found = sw->query(range);
            auto se_found = se->query(range);
            
            found.insert(found.end(), nw_found.begin(), nw_found.end());
            found.insert(found.end(), ne_found.begin(), ne_found.end());
            found.insert(found.end(), sw_found.begin(), sw_found.end());
            found.insert(found.end(), se_found.begin(), se_found.end());
        }
        
        return found;
    }
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

class GeometryUtils {
public:
    // Point-in-polygon test - O(n)
    static bool point_in_polygon(const Point& p, const vector<Point>& polygon) {
        int n = polygon.size();
        bool inside = false;
        
        for (int i = 0, j = n - 1; i < n; j = i++) {
            double xi = polygon[i].x, yi = polygon[i].y;
            double xj = polygon[j].x, yj = polygon[j].y;
            
            if (((yi > p.y) != (yj > p.y)) && 
                (p.x < (xj - xi) * (p.y - yi) / (yj - yi) + xi)) {
                inside = !inside;
            }
        }
        
        return inside;
    }
    
    // Line segment intersection - O(1)
    static bool line_intersection(const Point& p1, const Point& p2,
                                  const Point& p3, const Point& p4,
                                  Point* intersection = nullptr) {
        double x1 = p1.x, y1 = p1.y;
        double x2 = p2.x, y2 = p2.y;
        double x3 = p3.x, y3 = p3.y;
        double x4 = p4.x, y4 = p4.y;
        
        double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (abs(denom) < 1e-9) return false;
        
        double t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        double u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
        
        if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
            if (intersection) {
                intersection->x = x1 + t * (x2 - x1);
                intersection->y = y1 + t * (y2 - y1);
            }
            return true;
        }
        
        return false;
    }
    
    // Polygon area - O(n)
    static double polygon_area(const vector<Point>& polygon) {
        double area = 0.0;
        int n = polygon.size();
        
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            area += polygon[i].x * polygon[j].y;
            area -= polygon[j].x * polygon[i].y;
        }
        
        return abs(area) / 2.0;
    }
};

// ============================================================================
// BENCHMARK & DEMO
// ============================================================================

void benchmark() {
    cout << "=== COMPUTATIONAL GEOMETRY FRAMEWORK BENCHMARK ===" << endl << endl;
    
    vector<int> sizes = {100, 1000, 10000, 50000};
    
    for (int n : sizes) {
        vector<Point> points;
        for (int i = 0; i < n; i++) {
            points.push_back(Point(rand() % 10000, rand() % 10000));
        }
        
        cout << "Testing with " << n << " points:" << endl;
        
        // Convex Hull
        auto start = high_resolution_clock::now();
        auto hull = ConvexHullSolver::graham_scan(points);
        auto end = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(end - start);
        cout << "  Graham Scan: " << duration.count() / 1000.0 << "ms (" 
             << hull.size() << " hull points)" << endl;
        
        // Delaunay (only for smaller sets)
        if (n <= 1000) {
            start = high_resolution_clock::now();
            DelaunayTriangulator dt(vector<Point>(points.begin(), points.begin() + min(n, 500)));
            auto triangles = dt.triangulate();
            end = high_resolution_clock::now();
            duration = duration_cast<microseconds>(end - start);
            cout << "  Delaunay: " << duration.count() / 1000.0 << "ms ("
                 << triangles.size() << " triangles)" << endl;
        }
        
        cout << endl;
    }
}

int main() {
    cout << "Computational Geometry & Mesh Generation Framework" << endl;
    cout << "===================================================" << endl << endl;
    
    // Demo with sample data
    vector<Point> points = {
        Point(0, 0), Point(10, 0), Point(10, 10), Point(0, 10),
        Point(5, 5), Point(2, 8), Point(8, 2), Point(7, 7)
    };
    
    // Convex Hull
    auto hull = ConvexHullSolver::graham_scan(points);
    cout << "Convex Hull: " << hull.size() << " vertices" << endl;
    
    // Delaunay Triangulation
    DelaunayTriangulator dt(points);
    auto triangles = dt.triangulate();
    cout << "Delaunay Triangulation: " << triangles.size() << " triangles" << endl;
    
    // Point-in-polygon test
    Point test_point(5, 5);
    bool inside = GeometryUtils::point_in_polygon(test_point, hull);
    cout << "Point (5,5) inside convex hull: " << (inside ? "YES" : "NO") << endl;
    
    cout << "\nFramework initialized successfully!" << endl;
    cout << "Ready for high-performance geometric computations." << endl << endl;
    
    // Run benchmark
    benchmark();
    
    return 0;
}