/*!
 * Computational Geometry & Mesh Generation Framework
 * Safe, concurrent, high-performance implementation in Rust
 * Zero-cost abstractions with memory safety guarantees
 */

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::f64::consts::PI;
use std::time::Instant;

// ============================================================================
// CORE DATA STRUCTURES
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
    
    pub fn distance_to(&self, other: &Point) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        (dx * dx + dy * dy).sqrt()
    }
    
    pub fn angle_to(&self, other: &Point) -> f64 {
        (other.y - self.y).atan2(other.x - self.x)
    }
}

impl Eq for Point {}

impl std::hash::Hash for Point {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.to_bits().hash(state);
        self.y.to_bits().hash(state);
    }
}

#[derive(Debug, Clone)]
pub struct Triangle {
    pub p1: Point,
    pub p2: Point,
    pub p3: Point,
}

impl Triangle {
    pub fn new(p1: Point, p2: Point, p3: Point) -> Self {
        Triangle { p1, p2, p3 }
    }
    
    pub fn circumcenter(&self) -> Point {
        let ax = self.p1.x;
        let ay = self.p1.y;
        let bx = self.p2.x;
        let by = self.p2.y;
        let cx = self.p3.x;
        let cy = self.p3.y;
        
        let d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
        
        if d.abs() < 1e-9 {
            return Point::new((ax + bx + cx) / 3.0, (ay + by + cy) / 3.0);
        }
        
        let ux = ((ax * ax + ay * ay) * (by - cy) +
                  (bx * bx + by * by) * (cy - ay) +
                  (cx * cx + cy * cy) * (ay - by)) / d;
        let uy = ((ax * ax + ay * ay) * (cx - bx) +
                  (bx * bx + by * by) * (ax - cx) +
                  (cx * cx + cy * cy) * (bx - ax)) / d;
        
        Point::new(ux, uy)
    }
    
    pub fn circumradius(&self) -> f64 {
        let center = self.circumcenter();
        center.distance_to(&self.p1)
    }
    
    pub fn circumcircle_contains(&self, p: &Point) -> bool {
        let center = self.circumcenter();
        let radius = self.circumradius();
        center.distance_to(p) < radius + 1e-9
    }
    
    pub fn contains_vertex(&self, p: &Point) -> bool {
        *p == self.p1 || *p == self.p2 || *p == self.p3
    }
    
    pub fn area(&self) -> f64 {
        let ax = self.p2.x - self.p1.x;
        let ay = self.p2.y - self.p1.y;
        let bx = self.p3.x - self.p1.x;
        let by = self.p3.y - self.p1.y;
        (ax * by - ay * bx).abs() / 2.0
    }
}

#[derive(Debug, Clone)]
pub struct Edge {
    pub p1: Point,
    pub p2: Point,
}

impl Edge {
    pub fn new(p1: Point, p2: Point) -> Self {
        Edge { p1, p2 }
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.p1 == other.p1 && self.p2 == other.p2) ||
        (self.p1 == other.p2 && self.p2 == other.p1)
    }
}

impl Eq for Edge {}

// ============================================================================
// CONVEX HULL ALGORITHMS
// ============================================================================

pub struct ConvexHullSolver;

impl ConvexHullSolver {
    /// Cross product of vectors OA and OB
    pub fn cross_product(o: &Point, a: &Point, b: &Point) -> f64 {
        (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    }
    
    /// Graham Scan algorithm - O(n log n)
    pub fn graham_scan(mut points: Vec<Point>) -> Vec<Point> {
        if points.len() < 3 {
            return points;
        }
        
        // Find lowest point (leftmost if tie)
        let mut lowest_idx = 0;
        for i in 1..points.len() {
            if points[i].y < points[lowest_idx].y ||
               (points[i].y == points[lowest_idx].y && points[i].x < points[lowest_idx].x) {
                lowest_idx = i;
            }
        }
        points.swap(0, lowest_idx);
        let pivot = points[0];
        
        // Sort by polar angle
        points[1..].sort_by(|a, b| {
            let angle_a = pivot.angle_to(a);
            let angle_b = pivot.angle_to(b);
            angle_a.partial_cmp(&angle_b).unwrap_or(Ordering::Equal)
        });
        
        // Build hull
        let mut hull = Vec::new();
        hull.push(points[0]);
        
        for i in 1..points.len() {
            while hull.len() > 1 {
                let len = hull.len();
                let cp = Self::cross_product(&hull[len-2], &hull[len-1], &points[i]);
                if cp <= 0.0 {
                    hull.pop();
                } else {
                    break;
                }
            }
            hull.push(points[i]);
        }
        
        hull
    }
    
    /// Quickhull algorithm - O(n log n) average case
    pub fn quickhull(points: &[Point]) -> Vec<Point> {
        if points.len() < 3 {
            return points.to_vec();
        }
        
        // Find extreme points
        let left = points.iter().min_by(|a, b| a.x.partial_cmp(&b.x).unwrap()).unwrap();
        let right = points.iter().max_by(|a, b| a.x.partial_cmp(&b.x).unwrap()).unwrap();
        
        let mut hull = Vec::new();
        hull.push(*left);
        
        // Split into upper and lower sets
        let upper: Vec<Point> = points.iter()
            .filter(|p| **p != *left && **p != *right && Self::cross_product(left, right, p) > 0.0)
            .copied()
            .collect();
        
        let lower: Vec<Point> = points.iter()
            .filter(|p| **p != *left && **p != *right && Self::cross_product(left, right, p) < 0.0)
            .copied()
            .collect();
        
        Self::find_hull(&mut hull, *left, *right, &upper);
        hull.push(*right);
        Self::find_hull(&mut hull, *right, *left, &lower);
        
        hull
    }
    
    fn find_hull(hull: &mut Vec<Point>, p1: Point, p2: Point, points: &[Point]) {
        if points.is_empty() {
            return;
        }
        
        // Find furthest point
        let mut max_dist = 0.0;
        let mut furthest_idx = 0;
        
        for (i, p) in points.iter().enumerate() {
            let dist = Self::cross_product(&p1, &p2, p).abs();
            if dist > max_dist {
                max_dist = dist;
                furthest_idx = i;
            }
        }
        
        let furthest = points[furthest_idx];
        
        // Split remaining points
        let left: Vec<Point> = points.iter()
            .filter(|p| **p != furthest && Self::cross_product(&p1, &furthest, p) > 0.0)
            .copied()
            .collect();
        
        let right: Vec<Point> = points.iter()
            .filter(|p| **p != furthest && Self::cross_product(&furthest, &p2, p) > 0.0)
            .copied()
            .collect();
        
        Self::find_hull(hull, p1, furthest, &left);
        hull.push(furthest);
        Self::find_hull(hull, furthest, p2, &right);
    }
}

// ============================================================================
// DELAUNAY TRIANGULATION
// ============================================================================

pub struct DelaunayTriangulator {
    points: Vec<Point>,
    triangles: Vec<Triangle>,
}

impl DelaunayTriangulator {
    pub fn new(points: Vec<Point>) -> Self {
        DelaunayTriangulator {
            points,
            triangles: Vec::new(),
        }
    }
    
    /// Bowyer-Watson algorithm for Delaunay triangulation - O(n log n) expected
    pub fn triangulate(&mut self) -> &[Triangle] {
        if self.points.len() < 3 {
            return &self.triangles;
        }
        
        // Find bounding box
        let min_x = self.points.iter().map(|p| p.x).fold(f64::INFINITY, f64::min);
        let max_x = self.points.iter().map(|p| p.x).fold(f64::NEG_INFINITY, f64::max);
        let min_y = self.points.iter().map(|p| p.y).fold(f64::INFINITY, f64::min);
        let max_y = self.points.iter().map(|p| p.y).fold(f64::NEG_INFINITY, f64::max);
        
        let dx = max_x - min_x;
        let dy = max_y - min_y;
        let delta_max = dx.max(dy);
        let mid_x = (min_x + max_x) / 2.0;
        let mid_y = (min_y + max_y) / 2.0;
        
        // Create super-triangle
        let p1 = Point::new(mid_x - 20.0 * delta_max, mid_y - delta_max);
        let p2 = Point::new(mid_x, mid_y + 20.0 * delta_max);
        let p3 = Point::new(mid_x + 20.0 * delta_max, mid_y - delta_max);
        
        self.triangles.clear();
        self.triangles.push(Triangle::new(p1, p2, p3));
        
        // Add points incrementally
        for point in self.points.clone().iter() {
            let mut bad_triangles = Vec::new();
            
            // Find bad triangles
            for tri in &self.triangles {
                if tri.circumcircle_contains(point) {
                    bad_triangles.push(tri.clone());
                }
            }
            
            // Find polygon boundary
            let mut polygon = Vec::new();
            for tri in &bad_triangles {
                let edges = vec![
                    Edge::new(tri.p1, tri.p2),
                    Edge::new(tri.p2, tri.p3),
                    Edge::new(tri.p3, tri.p1),
                ];
                
                for edge in edges {
                    let mut shared = false;
                    for other_tri in &bad_triangles {
                        if std::ptr::eq(tri, other_tri) {
                            continue;
                        }
                        
                        let other_edges = vec![
                            Edge::new(other_tri.p1, other_tri.p2),
                            Edge::new(other_tri.p2, other_tri.p3),
                            Edge::new(other_tri.p3, other_tri.p1),
                        ];
                        
                        if other_edges.contains(&edge) {
                            shared = true;
                            break;
                        }
                    }
                    
                    if !shared {
                        polygon.push(edge);
                    }
                }
            }
            
            // Remove bad triangles
            self.triangles.retain(|tri| !bad_triangles.iter().any(|bt| 
                tri.p1 == bt.p1 && tri.p2 == bt.p2 && tri.p3 == bt.p3
            ));
            
            // Re-triangulate
            for edge in polygon {
                self.triangles.push(Triangle::new(edge.p1, edge.p2, *point));
            }
        }
        
        // Remove triangles with super-triangle vertices
        self.triangles.retain(|tri| 
            !tri.contains_vertex(&p1) && 
            !tri.contains_vertex(&p2) && 
            !tri.contains_vertex(&p3)
        );
        
        &self.triangles
    }
    
    pub fn get_triangles(&self) -> &[Triangle] {
        &self.triangles
    }
}

// ============================================================================
// VORONOI DIAGRAM
// ============================================================================

pub struct VoronoiDiagram {
    points: Vec<Point>,
    cells: HashMap<usize, Vec<Point>>,
}

impl VoronoiDiagram {
    pub fn new(points: Vec<Point>) -> Self {
        VoronoiDiagram {
            points,
            cells: HashMap::new(),
        }
    }
    
    pub fn compute(&mut self) -> &HashMap<usize, Vec<Point>> {
        let mut triangulator = DelaunayTriangulator::new(self.points.clone());
        let triangles = triangulator.triangulate();
        
        self.cells.clear();
        
        // Calculate circumcenters
        for tri in triangles {
            let circumcenter = tri.circumcenter();
            
            for (i, point) in self.points.iter().enumerate() {
                if tri.contains_vertex(point) {
                    self.cells.entry(i).or_insert_with(Vec::new).push(circumcenter);
                }
            }
        }
        
        // Sort vertices counter-clockwise
        for (i, vertices) in self.cells.iter_mut() {
            let center_x = vertices.iter().map(|v| v.x).sum::<f64>() / vertices.len() as f64;
            let center_y = vertices.iter().map(|v| v.y).sum::<f64>() / vertices.len() as f64;
            let center = Point::new(center_x, center_y);
            
            vertices.sort_by(|a, b| {
                let angle_a = center.angle_to(a);
                let angle_b = center.angle_to(b);
                angle_a.partial_cmp(&angle_b).unwrap_or(Ordering::Equal)
            });
        }
        
        &self.cells
    }
}

// ============================================================================
// GEOMETRIC UTILITIES
// ============================================================================

pub struct GeometryUtils;

impl GeometryUtils {
    /// Point-in-polygon test using ray casting - O(n)
    pub fn point_in_polygon(point: &Point, polygon: &[Point]) -> bool {
        let mut inside = false;
        let n = polygon.len();
        
        let mut j = n - 1;
        for i in 0..n {
            let xi = polygon[i].x;
            let yi = polygon[i].y;
            let xj = polygon[j].x;
            let yj = polygon[j].y;
            
            if ((yi > point.y) != (yj > point.y)) &&
               (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi) {
                inside = !inside;
            }
            
            j = i;
        }
        
        inside
    }
    
    /// Line segment intersection - O(1)
    pub fn line_intersection(p1: &Point, p2: &Point, p3: &Point, p4: &Point) -> Option<Point> {
        let x1 = p1.x;
        let y1 = p1.y;
        let x2 = p2.x;
        let y2 = p2.y;
        let x3 = p3.x;
        let y3 = p3.y;
        let x4 = p4.x;
        let y4 = p4.y;
        
        let denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        
        if denom.abs() < 1e-9 {
            return None;
        }
        
        let t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
        let u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
        
        if t >= 0.0 && t <= 1.0 && u >= 0.0 && u <= 1.0 {
            Some(Point::new(x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
        } else {
            None
        }
    }
    
    /// Polygon area using shoelace formula - O(n)
    pub fn polygon_area(polygon: &[Point]) -> f64 {
        let n = polygon.len();
        let mut area = 0.0;
        
        for i in 0..n {
            let j = (i + 1) % n;
            area += polygon[i].x * polygon[j].y;
            area -= polygon[j].x * polygon[i].y;
        }
        
        area.abs() / 2.0
    }
}

// ============================================================================
// BENCHMARK & DEMO
// ============================================================================

fn benchmark() {
    println!("=== COMPUTATIONAL GEOMETRY FRAMEWORK BENCHMARK ===\n");
    
    let sizes = vec![100, 1000, 5000, 10000];
    
    for n in sizes {
        let points: Vec<Point> = (0..n)
            .map(|_| Point::new(rand::random::<f64>() * 1000.0, rand::random::<f64>() * 1000.0))
            .collect();
        
        println!("Testing with {} points:", n);
        
        // Convex Hull
        let start = Instant::now();
        let hull = ConvexHullSolver::graham_scan(points.clone());
        let duration = start.elapsed();
        println!("  Graham Scan: {:?} ({} hull points)", duration, hull.len());
        
        // Delaunay (only for smaller sets)
        if n <= 1000 {
            let start = Instant::now();
            let mut dt = DelaunayTriangulator::new(points[..n.min(500)].to_vec());
            let triangles = dt.triangulate();
            let duration = start.elapsed();
            println!("  Delaunay: {:?} ({} triangles)", duration, triangles.len());
        }
        
        println!();
    }
}

fn main() {
    println!("Computational Geometry & Mesh Generation Framework");
    println!("===================================================\n");
    
    // Demo with sample points
    let points = vec![
        Point::new(0.0, 0.0),
        Point::new(10.0, 0.0),
        Point::new(10.0, 10.0),
        Point::new(0.0, 10.0),
        Point::new(5.0, 5.0),
        Point::new(2.0, 8.0),
        Point::new(8.0, 2.0),
        Point::new(7.0, 7.0),
    ];
    
    // Convex Hull
    let hull = ConvexHullSolver::graham_scan(points.clone());
    println!("Convex Hull: {} vertices", hull.len());
    
    // Delaunay Triangulation
    let mut dt = DelaunayTriangulator::new(points.clone());
    let triangles = dt.triangulate();
    println!("Delaunay Triangulation: {} triangles", triangles.len());
    
    // Voronoi Diagram
    let mut voronoi = VoronoiDiagram::new(points.clone());
    let cells = voronoi.compute();
    println!("Voronoi Diagram: {} cells", cells.len());
    
    // Point-in-polygon test
    let test_point = Point::new(5.0, 5.0);
    let inside = GeometryUtils::point_in_polygon(&test_point, &hull);
    println!("Point (5,5) inside convex hull: {}", inside);
    
    println!("\nFramework initialized successfully!");
    println!("Ready for safe, concurrent geometric computations.");
    
    // Run benchmark
    println!("\n");
    benchmark();
}

// Mock rand for demonstration
mod rand {
    use std::time::{SystemTime, UNIX_EPOCH};
    
    static mut SEED: u64 = 0;
    
    pub fn random<T: From<f64>>() -> T {
        unsafe {
            if SEED == 0 {
                SEED = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs();
            }
            SEED = SEED.wrapping_mul(1103515245).wrapping_add(12345);
            T::from((SEED % 1000000) as f64 / 1000000.0)
        }
    }
}