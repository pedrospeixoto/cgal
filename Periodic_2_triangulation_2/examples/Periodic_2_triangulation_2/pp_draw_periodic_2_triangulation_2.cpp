#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_periodic_2_triangulation_2.h>

#include <fstream>
#include <map>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       PDT;

typedef PDT::Point                                          Point;
typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Offset                                         Offset;

// Convert CGAL triangulation to points and cells
void triangulation_to_mesh(const PDT& T,
                           std::vector<std::array<double, 3>>& points,
                           std::vector<std::array<int, 3>>& cells)
{
  // Clear output vectors
  points.clear();
  cells.clear();
  
  // Use map to track unique vertex+offset combinations
  std::map<std::pair<Vertex_handle, Offset>, int> vertex_map;
  
  // Get domain size
  auto domain = T.domain();
  double domain_width = domain.xmax() - domain.xmin();
  double domain_height = domain.ymax() - domain.ymin();
  
  int vertex_id = 0;
  
  // Iterate through all faces in the triangulation (single sheet only)
  for(auto fit = T.finite_faces_begin(); fit != T.finite_faces_end(); ++fit)
  {
    std::array<int, 3> triangle_vertices;
    
    // Process each of the 3 vertices of the triangle
    for(int i = 0; i < 3; i++)
    {
      auto vertex_handle = fit->vertex(i);
      Offset offset = T.get_offset(fit, i);  // Get actual offset from triangulation
      
      auto key = std::make_pair(vertex_handle, offset);
      
      // Check if this vertex+offset combination already exists
      auto it = vertex_map.find(key);
      if(it == vertex_map.end())
      {
        // New vertex - add it to the map and points list
        vertex_map[key] = vertex_id;
        
        // Get the base point
        Point base_point = vertex_handle->point();
        
        // No periodic offset applied (single sheet)
        double x = base_point.x();
        double y = base_point.y();
        double z = 0.0; // 2D grid, so z=0
        
        // if offset is not (0,0), apply periodic shift
        x += offset.x() * domain_width;
        y += offset.y() * domain_height;
        
        // print point x,y
        std::cout << "Point ID " << vertex_id << ": (" << x << ", " << y << ")\n";
        points.push_back({x, y, z});
        triangle_vertices[i] = vertex_id;
        vertex_id++;
      }
      else
      {
        // Existing vertex - reuse its ID
        triangle_vertices[i] = it->second;
      }
    }
    
    cells.push_back(triangle_vertices);
  }
}

// Write mesh data to VTU XML format
void write_vtu_xml(const std::vector<std::array<double, 3>>& points,
                   const std::vector<std::array<int, 3>>& cells,
                   const std::string& filename)
{
  std::ofstream out(filename);
  
  // Write VTU file header
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <UnstructuredGrid>\n";
  out << "    <Piece NumberOfPoints=\"" << points.size() 
      << "\" NumberOfCells=\"" << cells.size() << "\">\n";
  
  // Write points
  out << "      <Points>\n";
  out << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(const auto& p : points)
  {
    out << "          " << p[0] << " " << p[1] << " " << p[2] << "\n";
  }
  out << "        </DataArray>\n";
  out << "      </Points>\n";
  
  // Write cells
  out << "      <Cells>\n";
  
  // Connectivity
  out << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
  for(const auto& cell : cells)
  {
    out << "          " << cell[0] << " " << cell[1] << " " << cell[2] << "\n";
  }
  out << "        </DataArray>\n";
  
  // Offsets
  out << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  for(size_t i = 0; i < cells.size(); i++)
  {
    out << "          " << (i + 1) * 3 << "\n";
  }
  out << "        </DataArray>\n";
  
  // Types (5 = VTK_TRIANGLE)
  out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for(size_t i = 0; i < cells.size(); i++)
  {
    out << "          5\n";
  }
  out << "        </DataArray>\n";
  
  out << "      </Cells>\n";
  out << "    </Piece>\n";
  out << "  </UnstructuredGrid>\n";
  out << "</VTKFile>\n";
  
  out.close();
  std::cout << "Wrote " << points.size() << " points and " 
            << cells.size() << " triangles to " << filename << std::endl;
}

// Combined function for convenience
void write_vtu(const PDT& T, const std::string& filename)
{
  std::vector<std::array<double, 3>> points;
  std::vector<std::array<int, 3>> cells;
  
  triangulation_to_mesh(T, points, cells);
  write_vtu_xml(points, cells, filename);
}

int main(int argc, char* argv[])
{
  // Declare periodic triangulation 2D
  PDT T;

  // Read points and insert in T
  Point p;
  std::ifstream ifs((argc > 1) ? argv[1] : "data/data1.dt.cin");
  if (ifs)
  {
    while (ifs >> p)
    { T.insert(p);
      std::cout << "Inserted point: " << p << std::endl;
    }

    if( T.is_triangulation_in_1_sheet())
    { T.convert_to_9_sheeted_covering(); }

    // Write to VTU file
    std::string vtu_filename = "periodic_triangulation.vtu";
    write_vtu(T, vtu_filename);

    // Draw the periodic triangulation
    CGAL::draw(T);
  }

  return EXIT_SUCCESS;
}
