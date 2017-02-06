#ifndef _SCENE_H_
#define _SCENE_H_

// std
#include <fstream>
#include <algorithm>
#include <list>

// Qt
#include <QtOpenGL>

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "dt.h"

// Delaunay triangulation
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Triangulation_vertex_base_2<Kernel> Vertex_base;
typedef My_vertex_base<Kernel> Vb;
typedef CGAL::Triangulation_face_base_2<Kernel> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef DT<Kernel, TDS> Triangulation;

#undef min
#undef max


class Scene
{
public:
  typedef Kernel::FT FT;
  typedef Kernel::Point_2 Point;


private:
  // Delaunay triangulation
  Triangulation m_dt;
  bool m_view_delaunay_edges;
  bool m_view_voronoi_edges;
  bool m_view_convex_hull;
  bool m_view_face_under_mouse;
  bool m_view_line_walk;
  bool m_view_nearest_vertex;
  bool m_view_draw_centroid;
  double lloyd_speed = 100;
  Point m_mouse_pos;
  bool write_monte_carlo;
  bool write_analytic_energy;
  bool write_max_vec;
  bool control_oscilation;

public:    
  int frame_count;
  Scene()
  {
    m_view_delaunay_edges = true;
    m_view_voronoi_edges = true;
    m_view_convex_hull = true;
    m_view_face_under_mouse = true;
    m_view_line_walk = true;
    m_view_nearest_vertex = true;
    m_view_draw_centroid = true;
    frame_count=0;
  }

  ~Scene()
  {
  }

  // accessors
  Point& mouse_pos() {return m_mouse_pos; }

  void clear()
  {
    m_dt.clear();
  }

  void add_point(const Point& point)
  {
    m_dt.insert_with_check(point);
  }



  void render()
  {
    m_dt.draw_border();


    if(m_view_line_walk){
      m_dt.draw_line_walk(0,255,0);
    }

    // render generators
    m_dt.gl_draw_generators(4.0f, 0, 0, 255); // blue dots

    // render Delaunay edges
    if(m_view_delaunay_edges)
      m_dt.gl_draw_delaunay_edges(1.0f, 0, 0, 128); // dark blue line segments

    // render Voronoi edges
    if(m_view_voronoi_edges)
      m_dt.gl_draw_voronoi_edges(1.0f, 0, 128, 0); // dark green line segments

    if(m_view_convex_hull)
      m_dt.draw_convex_hull(2.0f, 255, 0, 0); // red thick line segments

    if(m_view_face_under_mouse)
      m_dt.draw_face_under_mouse(m_mouse_pos,255,0,0);

    if(m_view_nearest_vertex){
      m_dt.draw_nearest_vertex(m_mouse_pos);
    }

    if(m_view_draw_centroid){
      m_dt.draw_centroid();
    }
  }

  void set_view_Delaunay_edges(bool value) { m_view_delaunay_edges = value; }
  void set_view_Voronoi_edges(bool value) { m_view_voronoi_edges = value; }
  void set_view_Convex_hull(bool value) { m_view_convex_hull = value; }
  void set_view_Face_under_mouse(bool value){ m_view_face_under_mouse = value; }
  void set_view_Line_walk(bool value){ m_view_line_walk = value; }
  void set_view_Nearest_vertex(bool value){ m_view_nearest_vertex = value; }
  void set_view_Cell_centroid(bool value){ m_view_draw_centroid = value; }
  void set_monte_carlo(bool value){ write_monte_carlo = value; }
  void set_analytic_energy(bool value){ write_analytic_energy = value; }
  void set_max_vec(bool value){write_max_vec=value;}
  void set_control_oscillation(bool value){control_oscilation=value;}
  void set_border_shape(Triangulation::BorderType type){
    m_dt.clear();
    m_dt.change_border_shape(type);
  }

  void set_step_param(StepParam sp){
    m_dt.set_step_param(sp);
  }

  void change_lloyd_speed(double speed){this->lloyd_speed=speed;}
  void initial_unit_circle(){
    m_dt.initial_unit_circle();
  }

  void initial_uniform(){
    m_dt.initial_uniform();
  }

  void lloyd()
  {
    // to complete
    // m_dt.lloyd();
    m_dt.lloyd(this->lloyd_speed,write_monte_carlo,write_analytic_energy,write_max_vec,control_oscilation);
    frame_count++;
  }

};

#endif // _SCENE_H_
