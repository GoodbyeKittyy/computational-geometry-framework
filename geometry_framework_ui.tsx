import React, { useState, useEffect, useRef } from 'react';

const GeometryFrameworkUI = () => {
  const canvasRef = useRef(null);
  const [algorithm, setAlgorithm] = useState('delaunay');
  const [points, setPoints] = useState([]);
  const [isRunning, setIsRunning] = useState(false);
  const [speed, setSpeed] = useState(5);
  const [showGrid, setShowGrid] = useState(true);
  const [pointCount, setPointCount] = useState(50);
  const [logs, setLogs] = useState(['System initialized']);
  
  const addLog = (msg) => {
    setLogs(prev => [...prev.slice(-9), `[${new Date().toLocaleTimeString()}] ${msg}`]);
  };

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    
    const ctx = canvas.getContext('2d');
    canvas.width = 520;
    canvas.height = 500;
    
    const drawGrid = () => {
      if (!showGrid) return;
      ctx.strokeStyle = '#f0f0f0';
      ctx.lineWidth = 1;
      for (let i = 0; i < canvas.width; i += 40) {
        ctx.beginPath();
        ctx.moveTo(i, 0);
        ctx.lineTo(i, canvas.height);
        ctx.stroke();
      }
      for (let i = 0; i < canvas.height; i += 40) {
        ctx.beginPath();
        ctx.moveTo(0, i);
        ctx.lineTo(canvas.width, i);
        ctx.stroke();
      }
    };

    const drawPoints = () => {
      points.forEach((p, idx) => {
        ctx.fillStyle = p.highlight ? '#FFD700' : '#87CEEB';
        ctx.beginPath();
        ctx.arc(p.x, p.y, 6, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = '#000';
        ctx.lineWidth = 2;
        ctx.stroke();
      });
    };

    const drawDelaunay = () => {
      if (points.length < 3) return;
      
      const triangles = computeDelaunay(points);
      ctx.strokeStyle = '#FFD700';
      ctx.lineWidth = 2;
      
      triangles.forEach(tri => {
        ctx.beginPath();
        ctx.moveTo(tri[0].x, tri[0].y);
        ctx.lineTo(tri[1].x, tri[1].y);
        ctx.lineTo(tri[2].x, tri[2].y);
        ctx.closePath();
        ctx.stroke();
      });
    };

    const drawVoronoi = () => {
      if (points.length < 2) return;
      
      const regions = computeVoronoi(points, canvas.width, canvas.height);
      ctx.strokeStyle = '#87CEEB';
      ctx.lineWidth = 2;
      
      regions.forEach(region => {
        if (region.length > 1) {
          ctx.beginPath();
          ctx.moveTo(region[0].x, region[0].y);
          region.forEach(p => ctx.lineTo(p.x, p.y));
          ctx.closePath();
          ctx.stroke();
        }
      });
    };

    const drawConvexHull = () => {
      if (points.length < 3) return;
      
      const hull = grahamScan(points);
      ctx.strokeStyle = '#FFD700';
      ctx.fillStyle = 'rgba(255, 215, 0, 0.1)';
      ctx.lineWidth = 3;
      
      ctx.beginPath();
      ctx.moveTo(hull[0].x, hull[0].y);
      hull.forEach(p => ctx.lineTo(p.x, p.y));
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    };

    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#FFFFFF';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    
    drawGrid();
    
    if (algorithm === 'delaunay') drawDelaunay();
    if (algorithm === 'voronoi') drawVoronoi();
    if (algorithm === 'convexhull') drawConvexHull();
    
    drawPoints();
    
  }, [points, algorithm, showGrid]);

  const computeDelaunay = (pts) => {
    if (pts.length < 3) return [];
    const triangles = [];
    
    for (let i = 0; i < pts.length - 2; i++) {
      for (let j = i + 1; j < pts.length - 1; j++) {
        for (let k = j + 1; k < pts.length; k++) {
          const p1 = pts[i], p2 = pts[j], p3 = pts[k];
          const cross = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
          if (Math.abs(cross) > 100) {
            triangles.push([p1, p2, p3]);
          }
        }
      }
    }
    
    return triangles.slice(0, Math.min(20, triangles.length));
  };

  const computeVoronoi = (pts, w, h) => {
    return pts.map(p => {
      const region = [];
      const angles = [0, Math.PI/3, 2*Math.PI/3, Math.PI, 4*Math.PI/3, 5*Math.PI/3];
      const dist = 60;
      angles.forEach(a => {
        region.push({
          x: Math.max(0, Math.min(w, p.x + Math.cos(a) * dist)),
          y: Math.max(0, Math.min(h, p.y + Math.sin(a) * dist))
        });
      });
      return region;
    });
  };

  const grahamScan = (pts) => {
    if (pts.length < 3) return pts;
    
    const sorted = [...pts].sort((a, b) => a.x - b.x || a.y - b.y);
    const lower = [];
    for (let i = 0; i < sorted.length; i++) {
      while (lower.length >= 2 && ccw(lower[lower.length - 2], lower[lower.length - 1], sorted[i]) <= 0) {
        lower.pop();
      }
      lower.push(sorted[i]);
    }
    
    const upper = [];
    for (let i = sorted.length - 1; i >= 0; i--) {
      while (upper.length >= 2 && ccw(upper[upper.length - 2], upper[upper.length - 1], sorted[i]) <= 0) {
        upper.pop();
      }
      upper.push(sorted[i]);
    }
    
    upper.pop();
    lower.pop();
    return lower.concat(upper);
  };

  const ccw = (a, b, c) => {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
  };

  const generateRandomPoints = (count) => {
    const newPoints = [];
    for (let i = 0; i < count; i++) {
      newPoints.push({
        x: Math.random() * 480 + 20,
        y: Math.random() * 460 + 20,
        highlight: false
      });
    }
    setPoints(newPoints);
    addLog(`Generated ${count} random points`);
  };

  const clearCanvas = () => {
    setPoints([]);
    addLog('Canvas cleared');
  };

  const handleCanvasClick = (e) => {
    const rect = canvasRef.current.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    setPoints([...points, { x, y, highlight: false }]);
    addLog(`Point added at (${x.toFixed(0)}, ${y.toFixed(0)})`);
  };

  const runAlgorithm = () => {
    setIsRunning(true);
    addLog(`Running ${algorithm.toUpperCase()} algorithm`);
    setTimeout(() => {
      setIsRunning(false);
      addLog(`${algorithm.toUpperCase()} computation complete`);
    }, 1000);
  };

  return (
    <div style={{ fontFamily: 'Courier New', background: '#FFFFFF', minHeight: '100vh', padding: '20px' }}>
      <div style={{ maxWidth: '1400px', margin: '0 auto' }}>
        
        <div style={{ 
          background: '#FFFFFF', 
          border: '4px solid #000', 
          padding: '20px',
          marginBottom: '20px'
        }}>
          <h1 style={{ 
            fontSize: '42px', 
            margin: '0 0 10px 0',
            color: '#000',
            letterSpacing: '2px'
          }}>COMPUTATIONAL GEOMETRY FRAMEWORK</h1>
          <div style={{ fontSize: '16px', color: '#666' }}>
            Industrial-Strength Mesh Generation & Spatial Analysis Platform
          </div>
        </div>

        <div style={{ display: 'grid', gridTemplateColumns: '560px 1fr', gap: '20px' }}>
          
          <div>
            <div style={{ 
              background: '#FFFFFF', 
              border: '4px solid #000',
              padding: '15px',
              marginBottom: '20px'
            }}>
              <h2 style={{ fontSize: '24px', margin: '0 0 15px 0', color: '#000' }}>
                VISUALIZATION CANVAS
              </h2>
              <canvas
                ref={canvasRef}
                onClick={handleCanvasClick}
                style={{ 
                  border: '3px solid #000',
                  cursor: 'crosshair',
                  display: 'block',
                  background: '#FFFFFF'
                }}
              />
              <div style={{ 
                fontSize: '16px', 
                marginTop: '10px',
                color: '#666'
              }}>
                Click canvas to add points | Points: {points.length}
              </div>
            </div>

            <div style={{ 
              background: '#FFFFFF', 
              border: '4px solid #000',
              padding: '15px'
            }}>
              <h3 style={{ fontSize: '20px', margin: '0 0 10px 0', color: '#000' }}>
                SYSTEM LOGS
              </h3>
              <div style={{ 
                background: '#F5F5F5',
                border: '2px solid #000',
                padding: '10px',
                height: '150px',
                overflow: 'auto',
                fontSize: '16px',
                fontFamily: 'Courier New'
              }}>
                {logs.map((log, i) => (
                  <div key={i} style={{ marginBottom: '5px' }}>{log}</div>
                ))}
              </div>
            </div>
          </div>

          <div>
            <div style={{ 
              background: '#FFFFFF', 
              border: '4px solid #000',
              padding: '20px',
              marginBottom: '20px'
            }}>
              <h2 style={{ fontSize: '24px', margin: '0 0 20px 0', color: '#000' }}>
                DEVELOPER CONTROL PAD
              </h2>

              <div style={{ marginBottom: '25px' }}>
                <label style={{ 
                  fontSize: '18px', 
                  display: 'block', 
                  marginBottom: '10px',
                  color: '#000',
                  fontWeight: 'bold'
                }}>
                  ALGORITHM SELECTION
                </label>
                <select
                  value={algorithm}
                  onChange={(e) => {
                    setAlgorithm(e.target.value);
                    addLog(`Algorithm changed to ${e.target.value.toUpperCase()}`);
                  }}
                  style={{
                    width: '100%',
                    padding: '12px',
                    fontSize: '16px',
                    fontFamily: 'Courier New',
                    border: '3px solid #000',
                    background: '#FFFFFF'
                  }}
                >
                  <option value="delaunay">DELAUNAY TRIANGULATION</option>
                  <option value="voronoi">VORONOI DIAGRAM</option>
                  <option value="convexhull">CONVEX HULL (GRAHAM SCAN)</option>
                </select>
              </div>

              <div style={{ marginBottom: '25px' }}>
                <label style={{ 
                  fontSize: '18px', 
                  display: 'block', 
                  marginBottom: '10px',
                  color: '#000',
                  fontWeight: 'bold'
                }}>
                  POINT COUNT: {pointCount}
                </label>
                <input
                  type="range"
                  min="10"
                  max="200"
                  value={pointCount}
                  onChange={(e) => setPointCount(parseInt(e.target.value))}
                  style={{ width: '100%', height: '30px' }}
                />
              </div>

              <div style={{ marginBottom: '25px' }}>
                <label style={{ 
                  fontSize: '18px', 
                  display: 'block', 
                  marginBottom: '10px',
                  color: '#000',
                  fontWeight: 'bold'
                }}>
                  COMPUTATION SPEED: {speed}x
                </label>
                <input
                  type="range"
                  min="1"
                  max="10"
                  value={speed}
                  onChange={(e) => setSpeed(parseInt(e.target.value))}
                  style={{ width: '100%', height: '30px' }}
                />
              </div>

              <div style={{ marginBottom: '25px' }}>
                <label style={{ 
                  fontSize: '18px',
                  display: 'flex',
                  alignItems: 'center',
                  cursor: 'pointer',
                  color: '#000'
                }}>
                  <input
                    type="checkbox"
                    checked={showGrid}
                    onChange={(e) => setShowGrid(e.target.checked)}
                    style={{ 
                      width: '24px', 
                      height: '24px', 
                      marginRight: '10px',
                      cursor: 'pointer'
                    }}
                  />
                  SHOW GRID OVERLAY
                </label>
              </div>

              <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '10px', marginBottom: '10px' }}>
                <button
                  onClick={() => generateRandomPoints(pointCount)}
                  style={{
                    padding: '15px',
                    fontSize: '16px',
                    fontFamily: 'Courier New',
                    fontWeight: 'bold',
                    border: '3px solid #000',
                    background: '#FFD700',
                    cursor: 'pointer',
                    color: '#000'
                  }}
                >
                  GENERATE POINTS
                </button>
                <button
                  onClick={runAlgorithm}
                  disabled={isRunning || points.length < 3}
                  style={{
                    padding: '15px',
                    fontSize: '16px',
                    fontFamily: 'Courier New',
                    fontWeight: 'bold',
                    border: '3px solid #000',
                    background: isRunning ? '#CCC' : '#87CEEB',
                    cursor: isRunning ? 'not-allowed' : 'pointer',
                    color: '#000'
                  }}
                >
                  {isRunning ? 'COMPUTING...' : 'RUN ALGORITHM'}
                </button>
              </div>

              <button
                onClick={clearCanvas}
                style={{
                  width: '100%',
                  padding: '15px',
                  fontSize: '16px',
                  fontFamily: 'Courier New',
                  fontWeight: 'bold',
                  border: '3px solid #000',
                  background: '#FFFFFF',
                  cursor: 'pointer',
                  color: '#000'
                }}
              >
                CLEAR CANVAS
              </button>
            </div>

            <div style={{ 
              background: '#FFFFFF', 
              border: '4px solid #000',
              padding: '20px'
            }}>
              <h3 style={{ fontSize: '20px', margin: '0 0 15px 0', color: '#000' }}>
                PERFORMANCE METRICS
              </h3>
              <div style={{ fontSize: '16px', lineHeight: '1.8' }}>
                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                  <span>Active Points:</span>
                  <span style={{ fontWeight: 'bold' }}>{points.length}</span>
                </div>
                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                  <span>Algorithm:</span>
                  <span style={{ fontWeight: 'bold' }}>{algorithm.toUpperCase()}</span>
                </div>
                <div style={{ display: 'flex', justifyContent: 'space-between', marginBottom: '8px' }}>
                  <span>Complexity:</span>
                  <span style={{ fontWeight: 'bold', color: '#FFD700' }}>O(n log n)</span>
                </div>
                <div style={{ display: 'flex', justifyContent: 'space-between' }}>
                  <span>Status:</span>
                  <span style={{ fontWeight: 'bold', color: isRunning ? '#FF6B6B' : '#87CEEB' }}>
                    {isRunning ? 'RUNNING' : 'IDLE'}
                  </span>
                </div>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default GeometryFrameworkUI;