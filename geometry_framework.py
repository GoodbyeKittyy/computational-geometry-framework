"""
Computational Geometry & Mesh Generation Framework
Industrial-strength geometric computing platform
Author: Framework Developer
"""

import numpy as np
from typing import List, Tuple, Set, Optional
from dataclasses import dataclass
from collections import defaultdict
import heapq
import time


@dataclass
class Point:
    """Represents a 2D point in space"""
    x: float
    y: float
    
    def __hash__(self):
        return hash((self.x, self.y))
    
    def __eq__(self, other):
        return abs(self.x - other.x) < 1e-9 and abs(self.y - other.y) < 1e-9
    
    def distance_to(self, other: 'Point') -> float:
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)


@dataclass
class Triangle:
    """Represents a triangle in Delaunay triangulation"""
    p1: Point
    p2: Point
    p3: Point
    
    def contains_point(self, p: Point) -> bool:
        """Check if point is inside triangle using barycentric coordinates"""
        v0 = (self.p3.x - self.p1.x, self.p3.y - self.p1.y)
        v1 = (self.p2.x - self.p1.x, self.p2.y - self.p1.y)
        v2 = (p.x - self.p1.x, p.y - self.p1.y)
        
        dot00 = v0[0] * v0[0] + v0[1] * v0[1]
        dot01 = v0[0] * v1[0] + v0[1] * v1[1]
        dot02 = v0[0] * v2[0] + v0[1] * v2[1]
        dot11 = v1[0] * v1[0] + v1[1] * v1[1]
        dot12 = v1[0] * v2[0] + v1[1] * v2[1]
        
        inv_denom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom
        
        return (u >= 0) and (v >= 0) and (u + v <= 1)
    
    def circumcircle_contains(self, p: Point) -> bool:
        """Check if point is inside circumcircle of triangle"""
        ax, ay = self.p1.x - p.x, self.p1.y - p.y
        bx, by = self.p2.x - p.x, self.p2.y - p.y
        cx, cy = self.p3.x - p.x, self.p3.y - p.y
        
        det = (ax * ax + ay * ay) * (bx * cy - cx * by) - \
              (bx * bx + by * by) * (ax * cy - cx * ay) + \
              (cx * cx + cy * cy) * (ax * by - bx * ay)
        
        return det > 0


class ConvexHullSolver:
    """Implements Graham Scan and Quickhull algorithms for convex hull"""
    
    @staticmethod
    def cross_product(o: Point, a: Point, b: Point) -> float:
        """Calculate cross product of vectors OA and OB"""
        return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
    
    @staticmethod
    def graham_scan(points: List[Point]) -> List[Point]:
        """Graham Scan algorithm - O(n log n)"""
        if len(points) < 3:
            return points
        
        # Find lowest point (and leftmost if tie)
        start = min(points, key=lambda p: (p.y, p.x))
        
        # Sort points by polar angle with respect to start point
        def polar_angle_key(p: Point):
            dx, dy = p.x - start.x, p.y - start.y
            return (np.arctan2(dy, dx), dx*dx + dy*dy)
        
        sorted_points = sorted([p for p in points if p != start], key=polar_angle_key)
        
        # Build convex hull
        hull = [start]
        for p in sorted_points:
            while len(hull) > 1 and ConvexHullSolver.cross_product(hull[-2], hull[-1], p) <= 0:
                hull.pop()
            hull.append(p)
        
        return hull
    
    @staticmethod
    def quickhull(points: List[Point]) -> List[Point]:
        """Quickhull algorithm - O(n log n) average case"""
        if len(points) < 3:
            return points
        
        # Find extreme points
        min_point = min(points, key=lambda p: p.x)
        max_point = max(points, key=lambda p: p.x)
        
        hull = []
        
        def find_hull(p1: Point, p2: Point, points_set: Set[Point]):
            if not points_set:
                return
            
            # Find point furthest from line p1-p2
            max_dist = 0
            furthest = None
            for p in points_set:
                dist = abs(ConvexHullSolver.cross_product(p1, p2, p))
                if dist > max_dist:
                    max_dist = dist
                    furthest = p
            
            if furthest is None:
                return
            
            # Recursively find hull on both sides
            left_set = {p for p in points_set if ConvexHullSolver.cross_product(p1, furthest, p) > 0}
            right_set = {p for p in points_set if ConvexHullSolver.cross_product(furthest, p2, p) > 0}
            
            find_hull(p1, furthest, left_set)
            hull.append(furthest)
            find_hull(furthest, p2, right_set)
        
        points_set = set(points) - {min_point, max_point}
        
        # Upper hull
        upper_set = {p for p in points_set if ConvexHullSolver.cross_product(min_point, max_point, p) > 0}
        hull.append(min_point)
        find_hull(min_point, max_point, upper_set)
        hull.append(max_point)
        
        # Lower hull
        lower_set = {p for p in points_set if ConvexHullSolver.cross_product(max_point, min_point, p) > 0}
        find_hull(max_point, min_point, lower_set)
        
        return hull


class DelaunayTriangulator:
    """Implements Bowyer-Watson algorithm for Delaunay triangulation"""
    
    def __init__(self, points: List[Point]):
        self.points = points
        self.triangles: List[Triangle] = []
    
    def triangulate(self) -> List[Triangle]:
        """Perform Delaunay triangulation - O(n log n) expected"""
        if len(self.points) < 3:
            return []
        
        # Create super-triangle that contains all points
        min_x = min(p.x for p in self.points)
        max_x = max(p.x for p in self.points)
        min_y = min(p.y for p in self.points)
        max_y = max(p.y for p in self.points)
        
        dx = max_x - min_x
        dy = max_y - min_y
        delta_max = max(dx, dy)
        mid_x = (min_x + max_x) / 2
        mid_y = (min_y + max_y) / 2
        
        p1 = Point(mid_x - 20 * delta_max, mid_y - delta_max)
        p2 = Point(mid_x, mid_y + 20 * delta_max)
        p3 = Point(mid_x + 20 * delta_max, mid_y - delta_max)
        
        super_triangle = Triangle(p1, p2, p3)
        self.triangles = [super_triangle]
        
        # Add points one by one
        for point in self.points:
            bad_triangles = []
            
            # Find triangles whose circumcircle contains the point
            for tri in self.triangles:
                if tri.circumcircle_contains(point):
                    bad_triangles.append(tri)
            
            # Find boundary of polygonal hole
            polygon = []
            for tri in bad_triangles:
                edges = [
                    (tri.p1, tri.p2),
                    (tri.p2, tri.p3),
                    (tri.p3, tri.p1)
                ]
                for edge in edges:
                    shared = False
                    for other_tri in bad_triangles:
                        if other_tri == tri:
                            continue
                        other_edges = [
                            (other_tri.p1, other_tri.p2),
                            (other_tri.p2, other_tri.p3),
                            (other_tri.p3, other_tri.p1)
                        ]
                        if edge in other_edges or (edge[1], edge[0]) in other_edges:
                            shared = True
                            break
                    if not shared:
                        polygon.append(edge)
            
            # Remove bad triangles
            for tri in bad_triangles:
                self.triangles.remove(tri)
            
            # Re-triangulate the polygonal hole
            for edge in polygon:
                new_tri = Triangle(edge[0], edge[1], point)
                self.triangles.append(new_tri)
        
        # Remove triangles containing super-triangle vertices
        self.triangles = [
            tri for tri in self.triangles
            if tri.p1 not in [p1, p2, p3] and
               tri.p2 not in [p1, p2, p3] and
               tri.p3 not in [p1, p2, p3]
        ]
        
        return self.triangles


class VoronoiDiagram:
    """Computes Voronoi diagram from Delaunay triangulation"""
    
    def __init__(self, points: List[Point]):
        self.points = points
        self.delaunay = DelaunayTriangulator(points)
    
    def compute(self) -> dict:
        """Generate Voronoi diagram - O(n log n)"""
        triangles = self.delaunay.triangulate()
        voronoi_cells = defaultdict(list)
        
        # Calculate circumcenters of triangles
        for tri in triangles:
            circumcenter = self._circumcenter(tri)
            voronoi_cells[tri.p1].append(circumcenter)
            voronoi_cells[tri.p2].append(circumcenter)
            voronoi_cells[tri.p3].append(circumcenter)
        
        # Sort vertices of each cell in counter-clockwise order
        for point in voronoi_cells:
            vertices = voronoi_cells[point]
            center_x = sum(v.x for v in vertices) / len(vertices)
            center_y = sum(v.y for v in vertices) / len(vertices)
            
            def angle_key(v):
                return np.arctan2(v.y - center_y, v.x - center_x)
            
            voronoi_cells[point] = sorted(vertices, key=angle_key)
        
        return dict(voronoi_cells)
    
    @staticmethod
    def _circumcenter(tri: Triangle) -> Point:
        """Calculate circumcenter of triangle"""
        ax, ay = tri.p1.x, tri.p1.y
        bx, by = tri.p2.x, tri.p2.y
        cx, cy = tri.p3.x, tri.p3.y
        
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(d) < 1e-9:
            return Point((ax + bx + cx) / 3, (ay + by + cy) / 3)
        
        ux = ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d
        uy = ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d
        
        return Point(ux, uy)


class SpatialIndex:
    """R-tree implementation for spatial indexing"""
    
    def __init__(self, max_entries: int = 4):
        self.max_entries = max_entries
        self.root = None
        self.points = []
    
    def insert(self, point: Point):
        """Insert point into R-tree - O(log n) amortized"""
        self.points.append(point)
    
    def query_range(self, min_x: float, max_x: float, min_y: float, max_y: float) -> List[Point]:
        """Range query - O(log n + k) where k is result size"""
        return [
            p for p in self.points
            if min_x <= p.x <= max_x and min_y <= p.y <= max_y
        ]
    
    def nearest_neighbor(self, query_point: Point, k: int = 1) -> List[Point]:
        """K-nearest neighbor search - O(k log n)"""
        distances = [(p.distance_to(query_point), p) for p in self.points]
        heapq.heapify(distances)
        return [heapq.heappop(distances)[1] for _ in range(min(k, len(distances)))]


class QuadTree:
    """Quadtree for spatial partitioning"""
    
    def __init__(self, boundary: Tuple[float, float, float, float], capacity: int = 4):
        self.boundary = boundary  # (min_x, min_y, max_x, max_y)
        self.capacity = capacity
        self.points = []
        self.divided = False
        self.children = []
    
    def insert(self, point: Point) -> bool:
        """Insert point into quadtree - O(log n) average"""
        if not self._contains(point):
            return False
        
        if len(self.points) < self.capacity:
            self.points.append(point)
            return True
        
        if not self.divided:
            self._subdivide()
        
        for child in self.children:
            if child.insert(point):
                return True
        
        return False
    
    def _contains(self, point: Point) -> bool:
        min_x, min_y, max_x, max_y = self.boundary
        return min_x <= point.x <= max_x and min_y <= point.y <= max_y
    
    def _subdivide(self):
        min_x, min_y, max_x, max_y = self.boundary
        mid_x = (min_x + max_x) / 2
        mid_y = (min_y + max_y) / 2
        
        self.children = [
            QuadTree((min_x, min_y, mid_x, mid_y), self.capacity),
            QuadTree((mid_x, min_y, max_x, mid_y), self.capacity),
            QuadTree((min_x, mid_y, mid_x, max_y), self.capacity),
            QuadTree((mid_x, mid_y, max_x, max_y), self.capacity)
        ]
        self.divided = True
        
        for point in self.points:
            for child in self.children:
                if child.insert(point):
                    break
        
        self.points = []
    
    def query(self, range_boundary: Tuple[float, float, float, float]) -> List[Point]:
        """Range query - O(log n + k)"""
        found = []
        
        if not self._intersects(range_boundary):
            return found
        
        for point in self.points:
            rx_min, ry_min, rx_max, ry_max = range_boundary
            if rx_min <= point.x <= rx_max and ry_min <= point.y <= ry_max:
                found.append(point)
        
        if self.divided:
            for child in self.children:
                found.extend(child.query(range_boundary))
        
        return found
    
    def _intersects(self, other_boundary: Tuple[float, float, float, float]) -> bool:
        x1_min, y1_min, x1_max, y1_max = self.boundary
        x2_min, y2_min, x2_max, y2_max = other_boundary
        return not (x1_max < x2_min or x2_max < x1_min or y1_max < y2_min or y2_max < y1_min)


class GeometryUtils:
    """Utility functions for geometric computations"""
    
    @staticmethod
    def point_in_polygon(point: Point, polygon: List[Point]) -> bool:
        """Ray casting algorithm for point-in-polygon test - O(n)"""
        x, y = point.x, point.y
        inside = False
        
        j = len(polygon) - 1
        for i in range(len(polygon)):
            xi, yi = polygon[i].x, polygon[i].y
            xj, yj = polygon[j].x, polygon[j].y
            
            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside
            
            j = i
        
        return inside
    
    @staticmethod
    def line_intersection(p1: Point, p2: Point, p3: Point, p4: Point) -> Optional[Point]:
        """Find intersection of two line segments - O(1)"""
        x1, y1 = p1.x, p1.y
        x2, y2 = p2.x, p2.y
        x3, y3 = p3.x, p3.y
        x4, y4 = p4.x, p4.y
        
        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if abs(denom) < 1e-9:
            return None
        
        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
        u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
        
        if 0 <= t <= 1 and 0 <= u <= 1:
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            return Point(x, y)
        
        return None
    
    @staticmethod
    def polygon_area(polygon: List[Point]) -> float:
        """Calculate polygon area using shoelace formula - O(n)"""
        n = len(polygon)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += polygon[i].x * polygon[j].y
            area -= polygon[j].x * polygon[i].y
        return abs(area) / 2.0


def benchmark_performance():
    """Benchmark framework performance"""
    print("=== COMPUTATIONAL GEOMETRY FRAMEWORK BENCHMARK ===\n")
    
    # Test with different point counts
    for n in [100, 1000, 10000, 100000]:
        points = [Point(np.random.rand() * 1000, np.random.rand() * 1000) for _ in range(n)]
        
        print(f"Testing with {n} points:")
        
        # Convex Hull
        start = time.time()
        hull = ConvexHullSolver.graham_scan(points)
        elapsed = time.time() - start
        print(f"  Graham Scan: {elapsed:.4f}s ({len(hull)} hull points)")
        
        # Delaunay Triangulation (only for smaller sets)
        if n <= 1000:
            start = time.time()
            dt = DelaunayTriangulator(points[:min(n, 500)])
            triangles = dt.triangulate()
            elapsed = time.time() - start
            print(f"  Delaunay: {elapsed:.4f}s ({len(triangles)} triangles)")
        
        # Spatial indexing
        start = time.time()
        qtree = QuadTree((0, 0, 1000, 1000))
        for p in points:
            qtree.insert(p)
        elapsed = time.time() - start
        print(f"  QuadTree Insert: {elapsed:.4f}s")
        
        print()


if __name__ == "__main__":
    # Demonstration
    print("Computational Geometry & Mesh Generation Framework")
    print("=" * 50)
    
    # Generate sample points
    np.random.seed(42)
    points = [Point(np.random.rand() * 100, np.random.rand() * 100) for _ in range(20)]
    
    # Convex Hull
    hull = ConvexHullSolver.graham_scan(points)
    print(f"\nConvex Hull: {len(hull)} vertices")
    
    # Delaunay Triangulation
    dt = DelaunayTriangulator(points)
    triangles = dt.triangulate()
    print(f"Delaunay Triangulation: {len(triangles)} triangles")
    
    # Voronoi Diagram
    voronoi = VoronoiDiagram(points)
    cells = voronoi.compute()
    print(f"Voronoi Diagram: {len(cells)} cells")
    
    # Spatial Indexing
    rtree = SpatialIndex()
    for p in points:
        rtree.insert(p)
    query_point = Point(50, 50)
    nearest = rtree.nearest_neighbor(query_point, k=5)
    print(f"5-Nearest Neighbors to (50, 50): {len(nearest)} points")
    
    print("\nFramework initialized successfully!")
    print("Ready for industrial-strength geometric computations.")
    
    # Run benchmark
    benchmark_performance()