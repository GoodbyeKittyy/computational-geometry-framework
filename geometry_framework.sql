-- ============================================================================
-- Computational Geometry & Mesh Generation Framework - PostgreSQL Extension
-- Spatial database schema and geometric functions
-- Optimized for PostGIS integration and spatial queries
-- ============================================================================

-- Enable PostGIS extension
CREATE EXTENSION IF NOT EXISTS postgis;
CREATE EXTENSION IF NOT EXISTS postgis_topology;

-- ============================================================================
-- CORE TABLES
-- ============================================================================

-- Points table with spatial indexing
CREATE TABLE IF NOT EXISTS geometry_points (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255),
    geom GEOMETRY(Point, 4326),
    attributes JSONB,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_points_geom ON geometry_points USING GIST(geom);
CREATE INDEX idx_points_attrs ON geometry_points USING GIN(attributes);

-- Triangles table for Delaunay triangulation
CREATE TABLE IF NOT EXISTS geometry_triangles (
    id SERIAL PRIMARY KEY,
    geom GEOMETRY(Polygon, 4326),
    p1_id INTEGER REFERENCES geometry_points(id),
    p2_id INTEGER REFERENCES geometry_points(id),
    p3_id INTEGER REFERENCES geometry_points(id),
    area DOUBLE PRECISION,
    circumcenter GEOMETRY(Point, 4326),
    circumradius DOUBLE PRECISION,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_triangles_geom ON geometry_triangles USING GIST(geom);
CREATE INDEX idx_triangles_circumcenter ON geometry_triangles USING GIST(circumcenter);

-- Voronoi cells table
CREATE TABLE IF NOT EXISTS geometry_voronoi_cells (
    id SERIAL PRIMARY KEY,
    point_id INTEGER REFERENCES geometry_points(id),
    geom GEOMETRY(Polygon, 4326),
    area DOUBLE PRECISION,
    perimeter DOUBLE PRECISION,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_voronoi_geom ON geometry_voronoi_cells USING GIST(geom);
CREATE INDEX idx_voronoi_point ON geometry_voronoi_cells(point_id);

-- Convex hull results table
CREATE TABLE IF NOT EXISTS geometry_convex_hulls (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255),
    geom GEOMETRY(Polygon, 4326),
    point_count INTEGER,
    hull_point_count INTEGER,
    area DOUBLE PRECISION,
    perimeter DOUBLE PRECISION,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_hulls_geom ON geometry_convex_hulls USING GIST(geom);

-- Spatial index table (R-tree simulation)
CREATE TABLE IF NOT EXISTS geometry_spatial_index (
    id SERIAL PRIMARY KEY,
    bbox GEOMETRY(Polygon, 4326),
    level INTEGER,
    parent_id INTEGER REFERENCES geometry_spatial_index(id),
    point_ids INTEGER[],
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE INDEX idx_spatial_bbox ON geometry_spatial_index USING GIST(bbox);
CREATE INDEX idx_spatial_level ON geometry_spatial_index(level);

-- ============================================================================
-- GEOMETRIC COMPUTATION FUNCTIONS
-- ============================================================================

-- Function: Calculate cross product of vectors
CREATE OR REPLACE FUNCTION cross_product(
    o_x DOUBLE PRECISION, o_y DOUBLE PRECISION,
    a_x DOUBLE PRECISION, a_y DOUBLE PRECISION,
    b_x DOUBLE PRECISION, b_y DOUBLE PRECISION
) RETURNS DOUBLE PRECISION AS $$
BEGIN
    RETURN (a_x - o_x) * (b_y - o_y) - (a_y - o_y) * (b_x - o_x);
END;
$$ LANGUAGE plpgsql IMMUTABLE;

-- Function: Compute convex hull using Graham Scan
CREATE OR REPLACE FUNCTION compute_convex_hull(point_set_name VARCHAR)
RETURNS TABLE(id INTEGER, x DOUBLE PRECISION, y DOUBLE PRECISION) AS $$
DECLARE
    hull_geom GEOMETRY;
BEGIN
    -- Use PostGIS ST_ConvexHull for efficient computation
    SELECT ST_ConvexHull(ST_Collect(geom)) INTO hull_geom
    FROM geometry_points
    WHERE name = point_set_name OR point_set_name IS NULL;
    
    -- Extract hull points
    RETURN QUERY
    SELECT 
        ROW_NUMBER() OVER ()::INTEGER as id,
        ST_X(geom)::DOUBLE PRECISION as x,
        ST_Y(geom)::DOUBLE PRECISION as y
    FROM (
        SELECT (ST_DumpPoints(hull_geom)).geom
    ) AS points;
END;
$$ LANGUAGE plpgsql;

-- Function: Compute Delaunay triangulation
CREATE OR REPLACE FUNCTION compute_delaunay_triangulation(point_set_name VARCHAR)
RETURNS VOID AS $$
DECLARE
    triangulation GEOMETRY;
    tri GEOMETRY;
    tri_points GEOMETRY[];
BEGIN
    -- Delete existing triangulations for this point set
    DELETE FROM geometry_triangles 
    WHERE p1_id IN (SELECT id FROM geometry_points WHERE name = point_set_name);
    
    -- Compute Delaunay triangulation using PostGIS
    SELECT ST_DelaunayTriangles(ST_Collect(geom), 0.0, 0) INTO triangulation
    FROM geometry_points
    WHERE name = point_set_name OR point_set_name IS NULL;
    
    -- Insert triangles
    FOR tri IN SELECT (ST_Dump(triangulation)).geom
    LOOP
        tri_points := ARRAY[
            ST_PointN(ST_ExteriorRing(tri), 1),
            ST_PointN(ST_ExteriorRing(tri), 2),
            ST_PointN(ST_ExteriorRing(tri), 3)
        ];
        
        INSERT INTO geometry_triangles (geom, area, circumcenter, circumradius)
        VALUES (
            tri,
            ST_Area(tri),
            ST_Centroid(tri),
            ST_Distance(ST_Centroid(tri), tri_points[1])
        );
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Function: Compute Voronoi diagram
CREATE OR REPLACE FUNCTION compute_voronoi_diagram(point_set_name VARCHAR)
RETURNS VOID AS $$
DECLARE
    voronoi GEOMETRY;
    cell GEOMETRY;
    pt_id INTEGER;
BEGIN
    -- Delete existing Voronoi cells
    DELETE FROM geometry_voronoi_cells 
    WHERE point_id IN (SELECT id FROM geometry_points WHERE name = point_set_name);
    
    -- Compute Voronoi diagram using PostGIS
    SELECT ST_VoronoiPolygons(ST_Collect(geom), 0.0) INTO voronoi
    FROM geometry_points
    WHERE name = point_set_name OR point_set_name IS NULL;
    
    -- Insert cells
    FOR cell IN SELECT (ST_Dump(voronoi)).geom
    LOOP
        -- Find the point this cell belongs to
        SELECT id INTO pt_id
        FROM geometry_points
        WHERE ST_Contains(cell, geom)
        LIMIT 1;
        
        IF pt_id IS NOT NULL THEN
            INSERT INTO geometry_voronoi_cells (point_id, geom, area, perimeter)
            VALUES (
                pt_id,
                cell,
                ST_Area(cell),
                ST_Perimeter(cell)
            );
        END IF;
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- Function: Point-in-polygon test
CREATE OR REPLACE FUNCTION point_in_polygon_test(
    test_point GEOMETRY,
    polygon_id INTEGER
) RETURNS BOOLEAN AS $$
DECLARE
    poly GEOMETRY;
    result BOOLEAN;
BEGIN
    SELECT geom INTO poly FROM geometry_convex_hulls WHERE id = polygon_id;
    RETURN ST_Contains(poly, test_point);
END;
$$ LANGUAGE plpgsql;

-- Function: Line intersection
CREATE OR REPLACE FUNCTION line_intersection(
    p1_x DOUBLE PRECISION, p1_y DOUBLE PRECISION,
    p2_x DOUBLE PRECISION, p2_y DOUBLE PRECISION,
    p3_x DOUBLE PRECISION, p3_y DOUBLE PRECISION,
    p4_x DOUBLE PRECISION, p4_y DOUBLE PRECISION
) RETURNS GEOMETRY AS $$
DECLARE
    line1 GEOMETRY;
    line2 GEOMETRY;
BEGIN
    line1 := ST_MakeLine(ST_MakePoint(p1_x, p1_y), ST_MakePoint(p2_x, p2_y));
    line2 := ST_MakeLine(ST_MakePoint(p3_x, p3_y), ST_MakePoint(p4_x, p4_y));
    RETURN ST_Intersection(line1, line2);
END;
$$ LANGUAGE plpgsql;

-- Function: K-nearest neighbors
CREATE OR REPLACE FUNCTION k_nearest_neighbors(
    query_point GEOMETRY,
    k INTEGER
) RETURNS TABLE(
    point_id INTEGER,
    distance DOUBLE PRECISION,
    x DOUBLE PRECISION,
    y DOUBLE PRECISION
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        gp.id,
        ST_Distance(gp.geom, query_point) as dist,
        ST_X(gp.geom) as x,
        ST_Y(gp.geom) as y
    FROM geometry_points gp
    ORDER BY gp.geom <-> query_point
    LIMIT k;
END;
$$ LANGUAGE plpgsql;

-- Function: Range query (spatial window)
CREATE OR REPLACE FUNCTION range_query(
    min_x DOUBLE PRECISION, min_y DOUBLE PRECISION,
    max_x DOUBLE PRECISION, max_y DOUBLE PRECISION
) RETURNS TABLE(
    point_id INTEGER,
    x DOUBLE PRECISION,
    y DOUBLE PRECISION,
    attributes JSONB
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        id,
        ST_X(geom) as x,
        ST_Y(geom) as y,
        attributes
    FROM geometry_points
    WHERE ST_Intersects(
        geom,
        ST_MakeEnvelope(min_x, min_y, max_x, max_y, 4326)
    );
END;
$$ LANGUAGE plpgsql;

-- Function: Compute polygon area
CREATE OR REPLACE FUNCTION compute_polygon_area(polygon_id INTEGER)
RETURNS DOUBLE PRECISION AS $$
DECLARE
    poly GEOMETRY;
BEGIN
    SELECT geom INTO poly FROM geometry_convex_hulls WHERE id = polygon_id;
    RETURN ST_Area(poly);
END;
$$ LANGUAGE plpgsql;

-- Function: Build spatial index (R-tree)
CREATE OR REPLACE FUNCTION build_spatial_index()
RETURNS VOID AS $$
DECLARE
    bbox GEOMETRY;
    point_array INTEGER[];
BEGIN
    -- Clear existing index
    DELETE FROM geometry_spatial_index;
    
    -- Build root node
    SELECT ST_Extent(geom), ARRAY_AGG(id)
    INTO bbox, point_array
    FROM geometry_points;
    
    INSERT INTO geometry_spatial_index (bbox, level, point_ids)
    VALUES (bbox, 0, point_array);
    
    -- Recursively subdivide (simplified for demonstration)
    -- In production, implement full R-tree algorithm
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- ANALYSIS AND REPORTING FUNCTIONS
-- ============================================================================

-- Function: Generate mesh quality report
CREATE OR REPLACE FUNCTION generate_mesh_quality_report()
RETURNS TABLE(
    metric VARCHAR,
    value DOUBLE PRECISION,
    unit VARCHAR
) AS $$
BEGIN
    RETURN QUERY
    SELECT 'Total Triangles'::VARCHAR, COUNT(*)::DOUBLE PRECISION, 'count'::VARCHAR
    FROM geometry_triangles
    UNION ALL
    SELECT 'Average Triangle Area'::VARCHAR, AVG(area), 'sq units'::VARCHAR
    FROM geometry_triangles
    UNION ALL
    SELECT 'Min Triangle Area'::VARCHAR, MIN(area), 'sq units'::VARCHAR
    FROM geometry_triangles
    UNION ALL
    SELECT 'Max Triangle Area'::VARCHAR, MAX(area), 'sq units'::VARCHAR
    FROM geometry_triangles
    UNION ALL
    SELECT 'Total Points'::VARCHAR, COUNT(*)::DOUBLE PRECISION, 'count'::VARCHAR
    FROM geometry_points
    UNION ALL
    SELECT 'Average Voronoi Cell Area'::VARCHAR, AVG(area), 'sq units'::VARCHAR
    FROM geometry_voronoi_cells;
END;
$$ LANGUAGE plpgsql;

-- Function: Find closest triangle to point
CREATE OR REPLACE FUNCTION find_closest_triangle(query_point GEOMETRY)
RETURNS TABLE(
    triangle_id INTEGER,
    distance DOUBLE PRECISION,
    area DOUBLE PRECISION
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        id,
        ST_Distance(geom, query_point) as dist,
        area
    FROM geometry_triangles
    ORDER BY geom <-> query_point
    LIMIT 1;
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- PERFORMANCE OPTIMIZATION FUNCTIONS
-- ============================================================================

-- Function: Vacuum and analyze spatial tables
CREATE OR REPLACE FUNCTION optimize_spatial_tables()
RETURNS VOID AS $$
BEGIN
    VACUUM ANALYZE geometry_points;
    VACUUM ANALYZE geometry_triangles;
    VACUUM ANALYZE geometry_voronoi_cells;
    VACUUM ANALYZE geometry_convex_hulls;
    REINDEX TABLE geometry_points;
    REINDEX TABLE geometry_triangles;
END;
$$ LANGUAGE plpgsql;

-- Function: Get spatial index statistics
CREATE OR REPLACE FUNCTION get_spatial_statistics()
RETURNS TABLE(
    table_name VARCHAR,
    row_count BIGINT,
    index_size TEXT,
    table_size TEXT
) AS $$
BEGIN
    RETURN QUERY
    SELECT 
        'geometry_points'::VARCHAR,
        COUNT(*)::BIGINT,
        pg_size_pretty(pg_total_relation_size('idx_points_geom')),
        pg_size_pretty(pg_total_relation_size('geometry_points'))
    FROM geometry_points
    UNION ALL
    SELECT 
        'geometry_triangles'::VARCHAR,
        COUNT(*)::BIGINT,
        pg_size_pretty(pg_total_relation_size('idx_triangles_geom')),
        pg_size_pretty(pg_total_relation_size('geometry_triangles'))
    FROM geometry_triangles
    UNION ALL
    SELECT 
        'geometry_voronoi_cells'::VARCHAR,
        COUNT(*)::BIGINT,
        pg_size_pretty(pg_total_relation_size('idx_voronoi_geom')),
        pg_size_pretty(pg_total_relation_size('geometry_voronoi_cells'))
    FROM geometry_voronoi_cells;
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- SAMPLE DATA GENERATION
-- ============================================================================

-- Function: Generate random points for testing
CREATE OR REPLACE FUNCTION generate_random_points(
    num_points INTEGER,
    min_x DOUBLE PRECISION DEFAULT 0,
    max_x DOUBLE PRECISION DEFAULT 1000,
    min_y DOUBLE PRECISION DEFAULT 0,
    max_y DOUBLE PRECISION DEFAULT 1000,
    point_set_name VARCHAR DEFAULT 'random_set'
) RETURNS VOID AS $$
DECLARE
    i INTEGER;
    x DOUBLE PRECISION;
    y DOUBLE PRECISION;
BEGIN
    FOR i IN 1..num_points LOOP
        x := min_x + (random() * (max_x - min_x));
        y := min_y + (random() * (max_y - min_y));
        
        INSERT INTO geometry_points (name, geom, attributes)
        VALUES (
            point_set_name,
            ST_SetSRID(ST_MakePoint(x, y), 4326),
            jsonb_build_object('index', i, 'generated', true)
        );
    END LOOP;
END;
$$ LANGUAGE plpgsql;

-- ============================================================================
-- MATERIALIZED VIEWS FOR PERFORMANCE
-- ============================================================================

-- Materialized view: Convex hull summary
CREATE MATERIALIZED VIEW IF NOT EXISTS mv_convex_hull_summary AS
SELECT 
    ch.id,
    ch.name,
    ch.point_count,
    ch.hull_point_count,
    ch.area,
    ch.perimeter,
    ST_AsText(ch.geom) as geom_wkt
FROM geometry_convex_hulls ch
ORDER BY ch.area DESC;

CREATE INDEX idx_mv_hull_area ON mv_convex_hull_summary(area);

-- Materialized view: Triangle quality metrics
CREATE MATERIALIZED VIEW IF NOT EXISTS mv_triangle_quality AS
SELECT 
    id,
    area,
    circumradius,
    area / (circumradius * circumradius) as quality_ratio,
    ST_AsText(geom) as geom_wkt
FROM geometry_triangles
WHERE area > 0;

CREATE INDEX idx_mv_triangle_quality ON mv_triangle_quality(quality_ratio);

-- ============================================================================
-- TRIGGERS FOR AUTOMATIC UPDATES
-- ============================================================================

-- Trigger: Update timestamp on point modification
CREATE OR REPLACE FUNCTION update_timestamp()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trg_update_points_timestamp
BEFORE UPDATE ON geometry_points
FOR EACH ROW
EXECUTE FUNCTION update_timestamp();

-- ============================================================================
-- USAGE EXAMPLES (COMMENTED SQL)
-- ============================================================================

/*
-- Example 1: Generate random points
SELECT generate_random_points(100, 0, 1000, 0, 1000, 'test_set');

-- Example 2: Compute convex hull
SELECT * FROM compute_convex_hull('test_set');

-- Example 3: Compute Delaunay triangulation
SELECT compute_delaunay_triangulation('test_set');

-- Example 4: Compute Voronoi diagram
SELECT compute_voronoi_diagram('test_set');

-- Example 5: Find k-nearest neighbors
SELECT * FROM k_nearest_neighbors(ST_MakePoint(500, 500), 5);

-- Example 6: Range query
SELECT * FROM range_query(100, 100, 500, 500);

-- Example 7: Generate quality report
SELECT * FROM generate_mesh_quality_report();

-- Example 8: Get spatial statistics
SELECT * FROM get_spatial_statistics();

-- Example 9: Optimize tables
SELECT optimize_spatial_tables();

-- Example 10: Refresh materialized views
REFRESH MATERIALIZED VIEW mv_convex_hull_summary;
REFRESH MATERIALIZED VIEW mv_triangle_quality;
*/