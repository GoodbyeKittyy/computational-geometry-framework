/**
 * Computational Geometry & Mesh Generation Framework
 * Functional programming approach with immutable data structures
 * Scala implementation for JVM-based distributed processing
 */

package com.geometry.framework

import scala.math._
import scala.collection.mutable
import scala.annotation.tailrec

// ============================================================================
// CORE DATA STRUCTURES
// ============================================================================

case class Point(x: Double, y: Double) {
  def distanceTo(other: Point): Double = {
    val dx = x - other.x
    val dy = y - other.y
    sqrt(dx * dx + dy * dy)
  }
  
  def angleTo(other: Point): Double = {
    atan2(other.y - y, other.x - x)
  }
  
  def +(other: Point): Point = Point(x + other.x, y + other.y)
  def -(other: Point): Point = Point(x - other.x, y - other.y)
  def *(scalar: Double): Point = Point(x * scalar, y * scalar)
}

case class Triangle(p1: Point, p2: Point, p3: Point) {
  def circumcenter: Point = {
    val ax = p1.x
    val ay = p1.y
    val bx = p2.x
    val by = p2.y
    val cx = p3.x
    val cy = p3.y
    
    val d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    
    if (abs(d) < 1e-9) {
      Point((ax + bx + cx) / 3.0, (ay + by + cy) / 3.0)
    } else {
      val ux = ((ax * ax + ay * ay) * (by - cy) +
                (bx * bx + by * by) * (cy - ay) +
                (cx * cx + cy * cy) * (ay - by)) / d
      val uy = ((ax * ax + ay * ay) * (cx - bx) +
                (bx * bx + by * by) * (ax - cx) +
                (cx * cx + cy * cy) * (bx - ax)) / d
      Point(ux, uy)
    }
  }
  
  def circumradius: Double = {
    val center = circumcenter
    center.distanceTo(p1)
  }
  
  def circumcircleContains(p: Point): Boolean = {
    val center = circumcenter
    val radius = circumradius
    center.distanceTo(p) < radius + 1e-9
  }
  
  def containsVertex(p: Point): Boolean = {
    p == p1 || p == p2 || p == p3
  }
  
  def area: Double = {
    val ax = p2.x - p1.x
    val ay = p2.y - p1.y
    val bx = p3.x - p1.x
    val by = p3.y - p1.y
    abs(ax * by - ay * bx) / 2.0
  }
  
  def vertices: List[Point] = List(p1, p2, p3)
}

case class Edge(p1: Point, p2: Point) {
  override def equals(obj: Any): Boolean = obj match {
    case Edge(a, b) => (p1 == a && p2 == b) || (p1 == b && p2 == a)
    case _ => false
  }
  
  override def hashCode(): Int = {
    p1.hashCode() + p2.hashCode()
  }
}

// ============================================================================
// CONVEX HULL ALGORITHMS
// ============================================================================

object ConvexHullSolver {
  
  def crossProduct(o: Point, a: Point, b: Point): Double = {
    (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x)
  }
  
  /**
   * Graham Scan algorithm for convex hull
   * Time Complexity: O(n log n)
   */
  def grahamScan(points: List[Point]): List[Point] = {
    if (points.length < 3) return points
    
    // Find lowest point (leftmost if tie)
    val pivot = points.minBy(p => (p.y, p.x))
    
    // Sort by polar angle
    val sorted = points.filter(_ != pivot)
      .sortBy(p => pivot.angleTo(p))
    
    // Build hull
    @tailrec
    def buildHull(remaining: List[Point], hull: List[Point]): List[Point] = {
      remaining match {
        case Nil => hull
        case point :: rest =>
          @tailrec
          def popInvalidPoints(currentHull: List[Point]): List[Point] = {
            currentHull match {
              case a :: b :: tail if crossProduct(a, b, point) <= 0 =>
                popInvalidPoints(a :: tail)
              case _ => currentHull
            }
          }
          
          val validHull = popInvalidPoints(hull)
          buildHull(rest, point :: validHull)
      }
    }
    
    val hull = buildHull(sorted, List(pivot))
    hull.reverse
  }
  
  /**
   * Quickhull algorithm for convex hull
   * Time Complexity: O(n log n) average case
   */
  def quickhull(points: List[Point]): List[Point] = {
    if (points.length < 3) return points
    
    val left = points.minBy(_.x)
    val right = points.maxBy(_.x)
    
    val upper = points.filter(p => p != left && p != right && crossProduct(left, right, p) > 0)
    val lower = points.filter(p => p != left && p != right && crossProduct(left, right, p) < 0)
    
    def findHull(p1: Point, p2: Point, points: List[Point]): List[Point] = {
      if (points.isEmpty) return List()
      
      val furthest = points.maxBy(p => abs(crossProduct(p1, p2, p)))
      
      val left = points.filter(p => p != furthest && crossProduct(p1, furthest, p) > 0)
      val right = points.filter(p => p != furthest && crossProduct(furthest, p2, p) > 0)
      
      findHull(p1, furthest, left) ++ List(furthest) ++ findHull(furthest, p2, right)
    }
    
    List(left) ++ findHull(left, right, upper) ++ List(right) ++ findHull(right, left, lower)
  }
}

// ============================================================================
// DELAUNAY TRIANGULATION
// ============================================================================

class DelaunayTriangulator(val points: List[Point]) {
  private var triangles: List[Triangle] = List()
  
  /**
   * Bowyer-Watson algorithm for Delaunay triangulation
   * Time Complexity: O(n log n) expected
   */
  def triangulate(): List[Triangle] = {
    if (points.length < 3) return List()
    
    // Find bounding box
    val minX = points.minBy(_.x).x
    val maxX = points.maxBy(_.x).x
    val minY = points.minBy(_.y).y
    val maxY = points.maxBy(_.y).y
    
    val dx = maxX - minX
    val dy = maxY - minY
    val deltaMax = max(dx, dy)
    val midX = (minX + maxX) / 2.0
    val midY = (minY + maxY) / 2.0
    
    // Create super-triangle
    val p1 = Point(midX - 20 * deltaMax, midY - deltaMax)
    val p2 = Point(midX, midY + 20 * deltaMax)
    val p3 = Point(midX + 20 * deltaMax, midY - deltaMax)
    
    triangles = List(Triangle(p1, p2, p3))
    
    // Add points incrementally
    points.foreach { point =>
      val badTriangles = triangles.filter(_.circumcircleContains(point))
      
      val polygon = badTriangles.flatMap { tri =>
        List(
          Edge(tri.p1, tri.p2),
          Edge(tri.p2, tri.p3),
          Edge(tri.p3, tri.p1)
        )
      }.groupBy(identity)
        .filter(_._2.length == 1)
        .keys
        .toList
      
      triangles = triangles.filterNot(badTriangles.contains)
      
      polygon.foreach { edge =>
        triangles = Triangle(edge.p1, edge.p2, point) :: triangles
      }
    }
    
    // Remove triangles with super-triangle vertices
    triangles = triangles.filterNot { tri =>
      tri.containsVertex(p1) || tri.containsVertex(p2) || tri.containsVertex(p3)
    }
    
    triangles
  }
  
  def getTriangles: List[Triangle] = triangles
}

// ============================================================================
// VORONOI DIAGRAM
// ============================================================================

class VoronoiDiagram(val points: List[Point]) {
  private var cells: Map[Point, List[Point]] = Map()
  
  def compute(): Map[Point, List[Point]] = {
    val triangulator = new DelaunayTriangulator(points)
    val triangles = triangulator.triangulate()
    
    // Calculate circumcenters
    val cellMap = mutable.Map[Point, mutable.ListBuffer[Point]]()
    
    triangles.foreach { tri =>
      val circumcenter = tri.circumcenter
      tri.vertices.foreach { vertex =>
        cellMap.getOrElseUpdate(vertex, mutable.ListBuffer()) += circumcenter
      }
    }
    
    // Sort vertices counter-clockwise
    cells = cellMap.map { case (point, vertices) =>
      val centerX = vertices.map(_.x).sum / vertices.length
      val centerY = vertices.map(_.y).sum / vertices.length
      val center = Point(centerX, centerY)
      
      val sorted = vertices.toList.sortBy(v => center.angleTo(v))
      (point, sorted)
    }.toMap
    
    cells
  }
  
  def getCells: Map[Point, List[Point]] = cells
}

// ============================================================================
// SPATIAL INDEXING - QUADTREE
// ============================================================================

case class BoundingBox(minX: Double, minY: Double, maxX: Double, maxY: Double) {
  def contains(p: Point): Boolean = {
    p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY
  }
  
  def intersects(other: BoundingBox): Boolean = {
    !(maxX < other.minX || other.maxX < minX || 
      maxY < other.minY || other.maxY < minY)
  }
}

class QuadTree(val boundary: BoundingBox, val capacity: Int = 4) {
  private var points: List[Point] = List()
  private var divided: Boolean = false
  private var children: Option[(QuadTree, QuadTree, QuadTree, QuadTree)] = None
  
  private def subdivide(): Unit = {
    val midX = (boundary.minX + boundary.maxX) / 2.0
    val midY = (boundary.minY + boundary.maxY) / 2.0
    
    val nw = new QuadTree(BoundingBox(boundary.minX, midY, midX, boundary.maxY), capacity)
    val ne = new QuadTree(BoundingBox(midX, midY, boundary.maxX, boundary.maxY), capacity)
    val sw = new QuadTree(BoundingBox(boundary.minX, boundary.minY, midX, midY), capacity)
    val se = new QuadTree(BoundingBox(midX, boundary.minY, boundary.maxX, midY), capacity)
    
    children = Some((nw, ne, sw, se))
    divided = true
    
    points.foreach { p =>
      nw.insert(p) || ne.insert(p) || sw.insert(p) || se.insert(p)
    }
    points = List()
  }
  
  def insert(p: Point): Boolean = {
    if (!boundary.contains(p)) return false
    
    if (points.length < capacity && !divided) {
      points = p :: points
      return true
    }
    
    if (!divided) subdivide()
    
    children match {
      case Some((nw, ne, sw, se)) =>
        nw.insert(p) || ne.insert(p) || sw.insert(p) || se.insert(p)
      case None => false
    }
  }
  
  def query(range: BoundingBox): List[Point] = {
    if (!boundary.intersects(range)) return List()
    
    val found = points.filter(range.contains)
    
    children match {
      case Some((nw, ne, sw, se)) =>
        found ++ nw.query(range) ++ ne.query(range) ++ sw.query(range) ++ se.query(range)
      case None => found
    }
  }
}

// ============================================================================
// GEOMETRIC UTILITIES
// ============================================================================

object GeometryUtils {
  
  /**
   * Point-in-polygon test using ray casting
   * Time Complexity: O(n)
   */
  def pointInPolygon(point: Point, polygon: List[Point]): Boolean = {
    @tailrec
    def raycast(i: Int, j: Int, inside: Boolean): Boolean = {
      if (i >= polygon.length) return inside
      
      val xi = polygon(i).x
      val yi = polygon(i).y
      val xj = polygon(j).x
      val yj = polygon(j).y
      
      val intersect = ((yi > point.y) != (yj > point.y)) &&
                     (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi)
      
      raycast(i + 1, i, if (intersect) !inside else inside)
    }
    
    if (polygon.isEmpty) false
    else raycast(0, polygon.length - 1, inside = false)
  }
  
  /**
   * Line segment intersection
   * Time Complexity: O(1)
   */
  def lineIntersection(p1: Point, p2: Point, p3: Point, p4: Point): Option[Point] = {
    val x1 = p1.x; val y1 = p1.y
    val x2 = p2.x; val y2 = p2.y
    val x3 = p3.x; val y3 = p3.y
    val x4 = p4.x; val y4 = p4.y
    
    val denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    
    if (abs(denom) < 1e-9) return None
    
    val t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
    val u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom
    
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
      Some(Point(x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
    } else {
      None
    }
  }
  
  /**
   * Polygon area using shoelace formula
   * Time Complexity: O(n)
   */
  def polygonArea(polygon: List[Point]): Double = {
    if (polygon.length < 3) return 0.0
    
    val area = polygon.indices.map { i =>
      val j = (i + 1) % polygon.length
      polygon(i).x * polygon(j).y - polygon(j).x * polygon(i).y
    }.sum
    
    abs(area) / 2.0
  }
  
  /**
   * K-nearest neighbors search
   * Time Complexity: O(n log k)
   */
  def kNearestNeighbors(queryPoint: Point, points: List[Point], k: Int): List[Point] = {
    points
      .map(p => (p, queryPoint.distanceTo(p)))
      .sortBy(_._2)
      .take(k)
      .map(_._1)
  }
  
  /**
   * Pairwise distances
   * Time Complexity: O(n^2)
   */
  def pairwiseDistances(points: List[Point]): Map[(Point, Point), Double] = {
    (for {
      i <- points.indices
      j <- i + 1 until points.length
    } yield {
      val p1 = points(i)
      val p2 = points(j)
      ((p1, p2), p1.distanceTo(p2))
    }).toMap
  }
}

// ============================================================================
// BENCHMARK & DEMO
// ============================================================================

object GeometryFramework {
  
  def benchmark(): Unit = {
    println("=== COMPUTATIONAL GEOMETRY FRAMEWORK BENCHMARK ===\n")
    
    val sizes = List(100, 1000, 5000, 10000)
    val random = new scala.util.Random(42)
    
    sizes.foreach { n =>
      val points = List.fill(n)(Point(random.nextDouble() * 1000, random.nextDouble() * 1000))
      
      println(s"Testing with $n points:")
      
      // Convex Hull
      val t1 = System.nanoTime()
      val hull = ConvexHullSolver.grahamScan(points)
      val t2 = System.nanoTime()
      println(f"  Graham Scan: ${(t2 - t1) / 1e6}%.4f ms (${hull.length} hull points)")
      
      // Delaunay Triangulation (only for smaller sets)
      if (n <= 1000) {
        val t3 = System.nanoTime()
        val triangulator = new DelaunayTriangulator(points.take(500))
        val triangles = triangulator.triangulate()
        val t4 = System.nanoTime()
        println(f"  Delaunay: ${(t4 - t3) / 1e6}%.4f ms (${triangles.length} triangles)")
      }
      
      // K-Nearest Neighbors
      val queryPoint = Point(500, 500)
      val t5 = System.nanoTime()
      val neighbors = GeometryUtils.kNearestNeighbors(queryPoint, points, 5)
      val t6 = System.nanoTime()
      println(f"  5-NN: ${(t6 - t5) / 1e6}%.4f ms")
      
      println()
    }
  }
  
  def demo(): Unit = {
    println("Computational Geometry & Mesh Generation Framework")
    println("===================================================\n")
    
    // Generate sample points
    val random = new scala.util.Random(42)
    val points = List.fill(20)(Point(random.nextDouble() * 100, random.nextDouble() * 100))
    
    // Convex Hull
    val hull = ConvexHullSolver.grahamScan(points)
    println(s"Convex Hull: ${hull.length} vertices")
    
    // Delaunay Triangulation
    val triangulator = new DelaunayTriangulator(points)
    val triangles = triangulator.triangulate()
    println(s"Delaunay Triangulation: ${triangles.length} triangles")
    
    // Voronoi Diagram
    val voronoi = new VoronoiDiagram(points)
    val cells = voronoi.compute()
    println(s"Voronoi Diagram: ${cells.size} cells")
    
    // Point-in-polygon test
    val testPoint = Point(50, 50)
    val inside = GeometryUtils.pointInPolygon(testPoint, hull)
    println(s"Point (50,50) inside convex hull: $inside")
    
    // Spatial indexing
    val quadtree = new QuadTree(BoundingBox(0, 0, 100, 100))
    points.foreach(quadtree.insert)
    val rangePoints = quadtree.query(BoundingBox(25, 25, 75, 75))
    println(s"Points in range [25,25] to [75,75]: ${rangePoints.length}")
    
    println("\nFramework initialized successfully!")
    println("Ready for functional, immutable geometric computations.")
  }
  
  def main(args: Array[String]): Unit = {
    demo()
    println()
    benchmark()
  }
}