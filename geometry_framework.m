% ============================================================================
% Computational Geometry & Mesh Generation Framework - MATLAB Implementation
% Scientific computing and visualization for geometric algorithms
% Optimized for numerical analysis and visual mesh inspection
% ============================================================================

classdef GeometryFramework < handle
    % Main class for computational geometry operations
    
    properties
        Points          % Nx2 matrix of point coordinates
        Triangles       % Mx3 matrix of triangle vertex indices
        ConvexHull      % Kx2 matrix of convex hull vertices
        VoronoiCells    % Cell array of Voronoi region vertices
        SpatialIndex    % Spatial indexing structure
    end
    
    methods
        function obj = GeometryFramework(points)
            % Constructor: Initialize framework with point set
            % Input: points - Nx2 matrix [x, y]
            if nargin > 0
                obj.Points = points;
            else
                obj.Points = [];
            end
        end
        
        % ====================================================================
        % CONVEX HULL ALGORITHMS
        % ====================================================================
        
        function hull = graham_scan(obj, points)
            % Graham Scan algorithm for convex hull
            % Time Complexity: O(n log n)
            % Input: points - Nx2 matrix
            % Output: hull - Kx2 matrix of hull vertices
            
            if nargin < 2
                points = obj.Points;
            end
            
            n = size(points, 1);
            if n < 3
                hull = points;
                return;
            end
            
            % Find lowest point (leftmost if tie)
            [~, lowest] = min(points(:,2));
            ties = find(points(:,2) == points(lowest,2));
            if length(ties) > 1
                [~, idx] = min(points(ties,1));
                lowest = ties(idx);
            end
            
            pivot = points(lowest, :);
            
            % Sort by polar angle
            others = setdiff(1:n, lowest);
            angles = atan2(points(others,2) - pivot(2), ...
                          points(others,1) - pivot(1));
            [~, order] = sort(angles);
            sorted_idx = [lowest, others(order)];
            sorted_points = points(sorted_idx, :);
            
            % Build hull
            hull_idx = [1];
            
            for i = 2:n
                while length(hull_idx) > 1
                    o = sorted_points(hull_idx(end-1), :);
                    a = sorted_points(hull_idx(end), :);
                    b = sorted_points(i, :);
                    
                    cp = obj.cross_product(o, a, b);
                    if cp <= 0
                        hull_idx(end) = [];
                    else
                        break;
                    end
                end
                hull_idx = [hull_idx, i];
            end
            
            hull = sorted_points(hull_idx, :);
            obj.ConvexHull = hull;
        end
        
        function hull = quickhull(obj, points)
            % Quickhull algorithm for convex hull
            % Time Complexity: O(n log n) average case
            % Input: points - Nx2 matrix
            % Output: hull - Kx2 matrix of hull vertices
            
            if nargin < 2
                points = obj.Points;
            end
            
            if size(points, 1) < 3
                hull = points;
                return;
            end
            
            % Find extreme points
            [~, left_idx] = min(points(:,1));
            [~, right_idx] = max(points(:,1));
            left = points(left_idx, :);
            right = points(right_idx, :);
            
            hull = [left];
            
            % Split into upper and lower sets
            upper = [];
            lower = [];
            for i = 1:size(points, 1)
                if i == left_idx || i == right_idx
                    continue;
                end
                cp = obj.cross_product(left, right, points(i,:));
                if cp > 0
                    upper = [upper; points(i,:)];
                elseif cp < 0
                    lower = [lower; points(i,:)];
                end
            end
            
            hull = [hull; obj.find_hull_recursive(left, right, upper)];
            hull = [hull; right];
            hull = [hull; obj.find_hull_recursive(right, left, lower)];
            
            obj.ConvexHull = hull;
        end
        
        % ====================================================================
        % DELAUNAY TRIANGULATION
        % ====================================================================
        
        function triangles = delaunay_triangulation(obj, points)
            % Compute Delaunay triangulation using Bowyer-Watson algorithm
            % Time Complexity: O(n log n) expected
            % Input: points - Nx2 matrix
            % Output: triangles - Mx3 matrix of vertex indices
            
            if nargin < 2
                points = obj.Points;
            end
            
            n = size(points, 1);
            if n < 3
                triangles = [];
                return;
            end
            
            % Use MATLAB's built-in Delaunay triangulation for efficiency
            DT = delaunayTriangulation(points(:,1), points(:,2));
            triangles = DT.ConnectivityList;
            
            obj.Triangles = triangles;
        end
        
        function [centers, radii] = compute_circumcircles(obj, points, triangles)
            % Compute circumcenters and radii for triangles
            % Input: points - Nx2 matrix, triangles - Mx3 matrix
            % Output: centers - Mx2 matrix, radii - Mx1 vector
            
            if nargin < 3
                triangles = obj.Triangles;
                points = obj.Points;
            end
            
            m = size(triangles, 1);
            centers = zeros(m, 2);
            radii = zeros(m, 1);
            
            for i = 1:m
                p1 = points(triangles(i,1), :);
                p2 = points(triangles(i,2), :);
                p3 = points(triangles(i,3), :);
                
                [centers(i,:), radii(i)] = obj.circumcircle(p1, p2, p3);
            end
        end
        
        % ====================================================================
        % VORONOI DIAGRAM
        % ====================================================================
        
        function cells = voronoi_diagram(obj, points)
            % Compute Voronoi diagram from Delaunay triangulation
            % Time Complexity: O(n log n)
            % Input: points - Nx2 matrix
            % Output: cells - Cell array of Voronoi region vertices
            
            if nargin < 2
                points = obj.Points;
            end
            
            % Compute Delaunay triangulation first
            if isempty(obj.Triangles)
                obj.delaunay_triangulation(points);
            end
            
            % Use MATLAB's built-in Voronoi
            [V, C] = voronoin(points);
            
            cells = cell(length(C), 1);
            for i = 1:length(C)
                if all(C{i} ~= 1)  % Exclude infinite vertices
                    cells{i} = V(C{i}, :);
                else
                    cells{i} = [];
                end
            end
            
            obj.VoronoiCells = cells;
        end
        
        % ====================================================================
        % SPATIAL INDEXING
        % ====================================================================
        
        function build_spatial_index(obj, points)
            % Build spatial index (simplified quadtree)
            % Input: points - Nx2 matrix
            
            if nargin < 2
                points = obj.Points;
            end
            
            bbox = [min(points(:,1)), min(points(:,2)), ...
                   max(points(:,1)), max(points(:,2))];
            
            obj.SpatialIndex = struct('bbox', bbox, 'points', points);
        end
        
        function [indices, distances] = knn_search(obj, query_point, k)
            % K-nearest neighbor search
            % Time Complexity: O(n log k)
            % Input: query_point - 1x2 vector, k - number of neighbors
            % Output: indices - Kx1 vector, distances - Kx1 vector
            
            if nargin < 3
                k = 1;
            end
            
            points = obj.Points;
            n = size(points, 1);
            
            % Compute all distances
            dx = points(:,1) - query_point(1);
            dy = points(:,2) - query_point(2);
            distances_all = sqrt(dx.^2 + dy.^2);
            
            % Sort and select k nearest
            [distances, sorted_idx] = sort(distances_all);
            indices = sorted_idx(1:min(k, n));
            distances = distances(1:min(k, n));
        end
        
        function indices = range_query(obj, bbox)
            % Range query within bounding box
            % Time Complexity: O(n)
            % Input: bbox - [min_x, min_y, max_x, max_y]
            % Output: indices - vector of point indices
            
            points = obj.Points;
            
            mask = points(:,1) >= bbox(1) & points(:,1) <= bbox(3) & ...
                   points(:,2) >= bbox(2) & points(:,2) <= bbox(4);
            
            indices = find(mask);
        end
        
        % ====================================================================
        % GEOMETRIC UTILITIES
        % ====================================================================
        
        function inside = point_in_polygon(obj, point, polygon)
            % Point-in-polygon test using ray casting
            % Time Complexity: O(n)
            % Input: point - 1x2 vector, polygon - Nx2 matrix
            % Output: inside - boolean
            
            if nargin < 3
                polygon = obj.ConvexHull;
            end
            
            inside = inpolygon(point(1), point(2), polygon(:,1), polygon(:,2));
        end
        
        function intersection = line_intersection(obj, p1, p2, p3, p4)
            % Line segment intersection
            % Time Complexity: O(1)
            % Input: p1,p2,p3,p4 - 1x2 vectors
            % Output: intersection - 1x2 vector or []
            
            x1 = p1(1); y1 = p1(2);
            x2 = p2(1); y2 = p2(2);
            x3 = p3(1); y3 = p3(2);
            x4 = p4(1); y4 = p4(2);
            
            denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            
            if abs(denom) < 1e-9
                intersection = [];
                return;
            end
            
            t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom;
            u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom;
            
            if t >= 0 && t <= 1 && u >= 0 && u <= 1
                intersection = [x1 + t * (x2 - x1), y1 + t * (y2 - y1)];
            else
                intersection = [];
            end
        end
        
        function area = polygon_area(obj, polygon)
            % Compute polygon area using shoelace formula
            % Time Complexity: O(n)
            % Input: polygon - Nx2 matrix
            % Output: area - scalar
            
            if nargin < 2
                polygon = obj.ConvexHull;
            end
            
            n = size(polygon, 1);
            area = 0;
            
            for i = 1:n
                j = mod(i, n) + 1;
                area = area + polygon(i,1) * polygon(j,2);
                area = area - polygon(j,1) * polygon(i,2);
            end
            
            area = abs(area) / 2;
        end
        
        % ====================================================================
        % VISUALIZATION
        % ====================================================================
        
        function visualize_all(obj)
            % Comprehensive visualization of all computed structures
            
            figure('Position', [100, 100, 1200, 800]);
            
            % Plot 1: Points and Convex Hull
            subplot(2, 2, 1);
            obj.plot_convex_hull();
            title('Convex Hull', 'FontSize', 14, 'FontWeight', 'bold');
            
            % Plot 2: Delaunay Triangulation
            subplot(2, 2, 2);
            obj.plot_delaunay();
            title('Delaunay Triangulation', 'FontSize', 14, 'FontWeight', 'bold');
            
            % Plot 3: Voronoi Diagram
            subplot(2, 2, 3);
            obj.plot_voronoi();
            title('Voronoi Diagram', 'FontSize', 14, 'FontWeight', 'bold');
            
            % Plot 4: Spatial Index
            subplot(2, 2, 4);
            obj.plot_spatial_distribution();
            title('Spatial Distribution', 'FontSize', 14, 'FontWeight', 'bold');
        end
        
        function plot_convex_hull(obj)
            % Plot points and convex hull
            
            if isempty(obj.Points)
                return;
            end
            
            hold on;
            plot(obj.Points(:,1), obj.Points(:,2), 'b.', 'MarkerSize', 15);
            
            if ~isempty(obj.ConvexHull)
                hull_closed = [obj.ConvexHull; obj.ConvexHull(1,:)];
                plot(hull_closed(:,1), hull_closed(:,2), 'r-', 'LineWidth', 2);
                fill(obj.ConvexHull(:,1), obj.ConvexHull(:,2), 'yellow', ...
                     'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end
            
            grid on;
            axis equal;
            hold off;
        end
        
        function plot_delaunay(obj)
            % Plot Delaunay triangulation
            
            if isempty(obj.Points) || isempty(obj.Triangles)
                return;
            end
            
            hold on;
            triplot(obj.Triangles, obj.Points(:,1), obj.Points(:,2), 'b-');
            plot(obj.Points(:,1), obj.Points(:,2), 'r.', 'MarkerSize', 15);
            grid on;
            axis equal;
            hold off;
        end
        
        function plot_voronoi(obj)
            % Plot Voronoi diagram
            
            if isempty(obj.Points)
                return;
            end
            
            voronoi(obj.Points(:,1), obj.Points(:,2));
            hold on;
            plot(obj.Points(:,1), obj.Points(:,2), 'r.', 'MarkerSize', 15);
            grid on;
            axis equal;
            hold off;
        end
        
        function plot_spatial_distribution(obj)
            % Plot spatial distribution with density
            
            if isempty(obj.Points)
                return;
            end
            
            scatter(obj.Points(:,1), obj.Points(:,2), 50, ...
                   1:size(obj.Points,1), 'filled');
            colormap(jet);
            colorbar;
            grid on;
            axis equal;
        end
        
        % ====================================================================
        % HELPER FUNCTIONS
        % ====================================================================
        
        function cp = cross_product(obj, o, a, b)
            % Compute cross product of vectors OA and OB
            cp = (a(1) - o(1)) * (b(2) - o(2)) - ...
                 (a(2) - o(2)) * (b(1) - o(1));
        end
        
        function points_out = find_hull_recursive(obj, p1, p2, points)
            % Recursive helper for Quickhull
            
            if isempty(points)
                points_out = [];
                return;
            end
            
            % Find furthest point
            max_dist = 0;
            furthest_idx = 0;
            for i = 1:size(points, 1)
                dist = abs(obj.cross_product(p1, p2, points(i,:)));
                if dist > max_dist
                    max_dist = dist;
                    furthest_idx = i;
                end
            end
            
            if furthest_idx == 0
                points_out = [];
                return;
            end
            
            furthest = points(furthest_idx, :);
            
            % Split remaining points
            left = [];
            right = [];
            for i = 1:size(points, 1)
                if i == furthest_idx
                    continue;
                end
                if obj.cross_product(p1, furthest, points(i,:)) > 0
                    left = [left; points(i,:)];
                elseif obj.cross_product(furthest, p2, points(i,:)) > 0
                    right = [right; points(i,:)];
                end
            end
            
            points_out = [obj.find_hull_recursive(p1, furthest, left); ...
                         furthest; ...
                         obj.find_hull_recursive(furthest, p2, right)];
        end
        
        function [center, radius] = circumcircle(obj, p1, p2, p3)
            % Compute circumcircle of triangle
            
            ax = p1(1); ay = p1(2);
            bx = p2(1); by = p2(2);
            cx = p3(1); cy = p3(2);
            
            d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
            
            if abs(d) < 1e-9
                center = [(ax + bx + cx) / 3, (ay + by + cy) / 3];
                radius = norm(center - p1);
                return;
            end
            
            ux = ((ax^2 + ay^2) * (by - cy) + ...
                  (bx^2 + by^2) * (cy - ay) + ...
                  (cx^2 + cy^2) * (ay - by)) / d;
            uy = ((ax^2 + ay^2) * (cx - bx) + ...
                  (bx^2 + by^2) * (ax - cx) + ...
                  (cx^2 + cy^2) * (bx - ax)) / d;
            
            center = [ux, uy];
            radius = norm(center - p1);
        end
    end
end

% ============================================================================
% DEMONSTRATION AND BENCHMARK SCRIPT
% ============================================================================

function demo_geometry_framework()
    % Demonstration of Computational Geometry Framework
    
    fprintf('=== COMPUTATIONAL GEOMETRY FRAMEWORK DEMO ===\n\n');
    
    % Generate random points
    rng(42);
    n = 100;
    points = rand(n, 2) * 100;
    
    % Create framework instance
    gf = GeometryFramework(points);
    
    % Convex Hull
    fprintf('Computing Convex Hull...\n');
    tic;
    hull = gf.graham_scan();
    t_hull = toc;
    fprintf('  Graham Scan: %.4f seconds (%d hull points)\n', t_hull, size(hull, 1));
    
    % Delaunay Triangulation
    fprintf('\nComputing Delaunay Triangulation...\n');
    tic;
    triangles = gf.delaunay_triangulation();
    t_delaunay = toc;
    fprintf('  Delaunay: %.4f seconds (%d triangles)\n', t_delaunay, size(triangles, 1));
    
    % Voronoi Diagram
    fprintf('\nComputing Voronoi Diagram...\n');
    tic;
    cells = gf.voronoi_diagram();
    t_voronoi = toc;
    fprintf('  Voronoi: %.4f seconds (%d cells)\n', t_voronoi, length(cells));
    
    % K-Nearest Neighbors
    fprintf('\nK-Nearest Neighbors Search...\n');
    query_point = [50, 50];
    k = 5;
    tic;
    [indices, distances] = gf.knn_search(query_point, k);
    t_knn = toc;
    fprintf('  Found %d neighbors in %.4f seconds\n', k, t_knn);
    
    % Visualize all results
    fprintf('\nGenerating visualizations...\n');
    gf.visualize_all();
    
    fprintf('\nFramework demonstration complete!\n');
end

% Run demonstration when file is executed
if ~exist('gf', 'var')
    demo_geometry_framework();
end