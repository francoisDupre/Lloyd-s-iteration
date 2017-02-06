#ifndef TRIANGULATION_H
#define TRIANGULATION_H
#include <algorithm>
#include <vector>
#include <CGAL/basic.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/number_utils.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h> //CGAL::bounded_side_2()
#include <stdexcept>
#include <list>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


// a base vertex with attributes
template< class Kernel, class Vbb = CGAL::Triangulation_vertex_base_2<Kernel> >
class My_vertex_base : public Vbb
{
public:
  typedef Vbb Vertex_base;
  typedef typename Vertex_base::Point Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vbb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef My_vertex_base<Kernel,Vb2>                           Other;
  };

private:
  // additional member data
  bool m_input;

public:
  bool& input() { return m_input; }

public:
  double last_x,last_y;
  double speed;
  My_vertex_base()
    : Vertex_base()
  {
    speed=150.0;
  }
  My_vertex_base(const Point & p, void* f)
    : Vertex_base(p,f)
  {
    speed=150.0;
  }
  My_vertex_base(const Point & p)
    : Vertex_base(p)
  {
    speed=150.0;
  }
};

struct StepParam{
  double default_speed;
  double gain_step;
  double moment;
  double max_speed;
};

template <class Kernel, class TDS>
class DT : public CGAL::Delaunay_triangulation_2<Kernel,TDS>
{
public:

  // typedefs for basic primitives
  typedef typename Kernel::FT FT;
  typedef typename Kernel::Ray_2      Ray;
  typedef typename Kernel::Line_2     Line;
  typedef typename Kernel::Point_2    Point;
  typedef typename Kernel::Vector_2   Vector;
  typedef typename Kernel::Segment_2  Segment;
  typedef typename Kernel::Triangle_2 Triangle;

  typedef typename CGAL::Delaunay_triangulation_2<Kernel,TDS> Dt2;

  // handles
  typedef typename Dt2::Face_handle         Face_handle;
  typedef typename Dt2::Vertex_handle       Vertex_handle;

  // iterators
  typedef typename Dt2::Face_iterator       Face_iterator;
  typedef typename Dt2::Edge_iterator       Edge_iterator;
  typedef typename Dt2::Finite_vertices_iterator Vertex_iterator; //

  // circulators
  typedef typename Dt2::Edge_circulator      Edge_circulator;
  typedef typename Dt2::Vertex_circulator    Vertex_circulator;
  typedef typename Dt2::Face_circulator      Face_circulator;
  typedef typename Dt2::Line_face_circulator Line_face_circulator;

  // exact kernel and polygon
  typedef typename CGAL::Exact_predicates_exact_constructions_kernel Kernel_exact;
  typedef typename Kernel_exact::Point_2 Point_exact;
  typedef typename CGAL::Polygon_2<Kernel_exact> Polygon_exact;
  typedef typename CGAL::Polygon_with_holes_2<Kernel_exact> Polygon_with_holes;
  typedef typename std::list<Polygon_with_holes> Pwh_list;
private:
  Polygon_exact border;
  FT border_surface;
  StepParam step_parameter;
public:
  enum BorderType{
    HOURGLASS,HOURGLASS_SMALL,HOURGLASS_LARGE,HEXAGON,SNAKE,SOFTSNAKE
  };
  std::vector<Point_exact> hourglass,hourglass_small,hourglass_large,hexagon,snake_soft,snake;

  void change_border_shape(BorderType type){
    std::vector<Point_exact>* p;
    switch(type){
    case DT::HOURGLASS:
      p=&hourglass;
      break;
    case DT::HOURGLASS_SMALL:
      p=&hourglass_small;
      break;
    case DT::HOURGLASS_LARGE:
      p=&hourglass_large;
      break;
    case DT::HEXAGON:
      p=&hexagon;
      break;
    case DT::SNAKE:
      p=&snake;
      break;
    case DT::SOFTSNAKE:
      p=&snake_soft;
      break;

    }
    this->border = Polygon_exact(p->begin(),p->end());
    if(!border.is_simple()){
      throw "border not simple";
    }
    std::vector<Point> points;
    for(auto i=border.vertices_begin();i!=border.vertices_end();i++){
      points.push_back(convert_point_inexact(*i));
    }
    this->border_surface = this->surface_polygon(points);

    std::cout<<"The surface of this shape is : "<<border_surface<<std::endl;
  }

  DT()
  {
    step_parameter.default_speed = 100;
    step_parameter.gain_step = 10;
    step_parameter.moment = 0.6;
    step_parameter.max_speed = 300;

    std::vector<Point_exact> _hourglass = {
      Point_exact(  8,  4),
      Point_exact(- 8,  4),
      Point_exact(- 1,- 4),
      Point_exact(- 8,-12),
      Point_exact(  8,-12),
      Point_exact(  1,- 4)
    };
    std::vector<Point_exact> _hourglass_small = {
      Point_exact(  8,  4),
      Point_exact(- 8,  4),
      Point_exact(- 0.5,- 4),
      Point_exact(- 8,-12),
      Point_exact(  8,-12),
      Point_exact(  0.5,- 4)
    };
    std::vector<Point_exact> _hourglass_large = {
      Point_exact(  8,  4),
      Point_exact(- 8,  4),
      Point_exact(- 1.5,- 4),
      Point_exact(- 8,-12),
      Point_exact(  8,-12),
      Point_exact(  1.5,- 4)
    };

    std::vector<Point_exact> _hexagon = {
      Point_exact(  4,  4),
      Point_exact( -4,  4),
      Point_exact( -8,- 4),
      Point_exact( -4,-12),
      Point_exact(  4,-12),
      Point_exact(  8,- 4)
    };
    std::vector<Point_exact> _snake_soft = {
      Point_exact(  2,  1),
      Point_exact(  1,  2),
      Point_exact(- 1,  2),
      Point_exact(- 2,  1),
      Point_exact(- 2,- 9),
      Point_exact(- 1,-10),
      Point_exact(  9,-10),
      Point_exact( 10,- 9),
      Point_exact( 10,- 3),
      Point_exact( 11,- 2),
      Point_exact( 13,- 2),
      Point_exact( 14,- 3),
      Point_exact( 14,- 9),
      Point_exact( 15,-10),
      Point_exact( 17,-10),
      Point_exact( 18,- 9),
      Point_exact( 18,  1),
      Point_exact( 17,  2),
      Point_exact(  7,  2),
      Point_exact(  6,  1),
      Point_exact(  6,- 5),
      Point_exact(  5,- 6),
      Point_exact(  3,- 6),
      Point_exact(  2,- 5),
    };
    std::vector<Point_exact> _snake = {
      Point_exact(  2,  2),
      Point_exact(- 2,  2),
      Point_exact(- 2,-10),
      Point_exact( 10,-10),
      Point_exact( 10,- 2),
      Point_exact( 14,- 2),
      Point_exact( 14,-10),
      Point_exact( 18,-10),
      Point_exact( 18,  2),
      Point_exact(  6,  2),
      Point_exact(  6,- 6),
      Point_exact(  2,- 6),
    };
    hourglass=_hourglass;
    hourglass_small=_hourglass_small;
    hourglass_large=_hourglass_large;
    hexagon=_hexagon;
    snake=_snake;
    snake_soft=_snake_soft;
    change_border_shape(HOURGLASS);

  }

  void set_step_param(StepParam sp){
    step_parameter=sp;
  }

  Point convert_point_inexact(Point_exact p){
    FT x = CGAL::to_double(p.hx()/p.hw());
    FT y = CGAL::to_double(p.hy()/p.hw());
    return Point(x,y);
  }

  virtual ~DT()
  {
    Dt2::clear();
  }

public:

  // random (uniform)
  void generators_random(unsigned int nb_generators)
  {
    Dt2::clear();
    for(unsigned int i = 0; i < nb_generators; i++)
      insert(Point(r(), r()));
  }

  // random number between zero and max
  FT r(FT max = 1.0) { return max * (FT)rand() / (FT)RAND_MAX; }

  // OPENGL DRAWINGS FUNCTIONS

  // draw generators
  void gl_draw_generators(const float point_size,
			  const unsigned char red,
			  const unsigned char green,
			  const unsigned char blue)
  {
    ::glColor3ub(red, green, blue);
    ::glPointSize(point_size);

    ::glBegin(GL_POINTS);
    typename Dt2::Point_iterator it;
    for(it = Dt2::points_begin();
	it != Dt2::points_end();
	it++)
    {
      const Point& p = *it;
      ::glVertex2f(p.x(), p.y());
    }
    ::glEnd();
  }

  // draw delaunay edges
  void gl_draw_delaunay_edges(const float line_width,
			      const unsigned char red,
			      const unsigned char green,
			      const unsigned char blue)
  {
    ::glColor3ub(red, green, blue);
    ::glLineWidth(line_width);
    ::glBegin(GL_LINES);
    Edge_iterator hEdge;
    for(hEdge  = Dt2::edges_begin();
	hEdge != Dt2::edges_end();
	hEdge++)
    {
      const Point& p1 = (*hEdge).first->vertex(Dt2::ccw((*hEdge).second))->point();
      const Point& p2 = (*hEdge).first->vertex(Dt2::cw((*hEdge).second))->point();
      ::glVertex2d(p1.x(), p1.y());
      ::glVertex2d(p2.x(), p2.y());
    }
    ::glEnd();
  }

  void gl_draw_voronoi_edges(const float line_width,
			     const unsigned char red,
			     const unsigned char green,
			     const unsigned char blue)
  {
    ::glColor3ub(red, green, blue);
    ::glLineWidth(line_width);
    Edge_iterator hEdge;

    ::glBegin(GL_LINES);
    for(hEdge = this->Dt2::edges_begin();hEdge!=this->Dt2::edges_end();hEdge++){
      CGAL::Object object = Dt2::dual(hEdge);
      Segment segment;
      Ray ray;
      Point source, target;
      if(CGAL::assign(segment,object))
      {
	source = segment.source();
	target = segment.target();
      }
      else if(CGAL::assign(ray,object))
      {
	source = ray.source();

	Vector v = ray.to_vector();
	FT min_coord = std::min(std::abs(v.x()),std::abs(v.y()))+0.1;
	Vector v_long(v.x()/min_coord*100.0,v.y()/min_coord*100.0);
	target = ray.point(0) + v_long;
      }
      ::glVertex2f(source.x(),source.y());
      ::glVertex2f(target.x(),target.y());
    }
    ::glEnd();
    // TO COMPLETE
  }
  void draw_convex_hull(const float line_width,
			const unsigned char red,
			const unsigned char green,
			const unsigned char blue){
    if( this->Dt2::dimension() != 2) return;
    ::glColor3ub(red, green, blue);
    ::glLineWidth(line_width);
    ::glBegin(GL_LINES);
    Vertex_handle inf_vertex = this->Dt2::infinite_vertex();
    Face_circulator fc =this->Dt2::incident_faces(inf_vertex);
    Face_circulator start(fc);
    do{
      Vertex_handle v1 = fc->vertex(this->tds().ccw(fc->index(inf_vertex)));
      Vertex_handle v2 = fc->vertex(this->tds().cw(fc->index(inf_vertex)));
      ::glVertex2f(v1->point().x(),v1->point().y());
      ::glVertex2f(v2->point().x(),v2->point().y());
    }while(++fc != start);
    ::glEnd();
  }
  void gl_draw_face(Face_handle f,const unsigned char red,
		    const unsigned char green,
		    const unsigned char blue){
    if(f==NULL) return;
    if(f->vertex(0) == this->Dt2::infinite_vertex()) return;
    if(f->vertex(1) == this->Dt2::infinite_vertex()) return;
    if(f->vertex(2) == this->Dt2::infinite_vertex()) return;
    ::glColor3ub(red, green, blue);
    ::glLineWidth(1.0);
    ::glBegin(GL_TRIANGLES);
    Vertex_handle v0 = f->vertex(0);
    Vertex_handle v1 = f->vertex(1);
    Vertex_handle v2 = f->vertex(2);
    ::glVertex2f(v0->point().x(),v0->point().y());
    ::glVertex2f(v1->point().x(),v1->point().y());
    ::glVertex2f(v2->point().x(),v2->point().y());
    ::glEnd();
  }

  void draw_face_under_mouse(const Point& mouse_pos,const unsigned char red,
			     const unsigned char green,
			     const unsigned char blue){
    if (this->Dt2::dimension() != 2) return;
    Face_handle f = Dt2::locate(mouse_pos);
    this->gl_draw_face(f,red,green,blue);
  }
  void draw_line_walk(const unsigned char red,
		      const unsigned char green,
		      const unsigned char blue){
    if (this->Dt2::dimension() != 2) return;
    Line_face_circulator circ = this->Dt2::line_walk(Point(0.0,0.0),Point(1.0,1.0));
    Line_face_circulator end = circ;
    ::glColor3ub(red, green, blue);
    ::glLineWidth(1.0);
    ::glBegin(GL_LINES);
    CGAL_For_all(circ,end){
      this->gl_draw_face(circ,red,green,blue);
    }

    ::glEnd();
  }
  void draw_nearest_vertex(const Point& mouse_pos){
    Vertex_handle v = this->Dt2::nearest_vertex(mouse_pos);
    if(v==NULL) return;
    ::glColor3ub(128, 0, 128);
    ::glPointSize( 9.0 );
    ::glBegin(GL_POINTS);
    ::glVertex2f(v->point().x(),v->point().y());
    ::glEnd();
  }

  FT inertia_polygon(std::vector<Vector> points){
    std::vector<std::pair<Vector,Vector> > edges;
    for(int i=0;i<points.size()-1;i++){
      edges.push_back(std::make_pair(points[i],points[i+1]));
    }
    edges.push_back(std::make_pair(points.back(),points[0]));
    FT I=0;
    for(auto e : edges){
      FT x1 = e.first.x(), x2 = e.second.x(), y1 = e.first.y(), y2 = e.second.y();
      FT Ii = 1/12.*(y1*y1+y1*y2+y2*y2);
      Ii += 1/12.*(x1*x1+x1*x2+x2*x2);
      Ii *= x1*y2-x2*y1;
      I+=Ii;
    }
    return I;
  }

  FT surface_polygon(std::vector<Point> points){
    std::vector<std::pair<Point,Point> > edges;
    for(int i=0;i<points.size()-1;i++){
      edges.push_back(std::make_pair(points[i],points[i+1]));
    }
    edges.push_back(std::make_pair(points.back(),points[0]));
    FT A=0;
    for(auto e : edges){
      FT area = 0.5*(e.first.x()*e.second.y()-e.second.x()*e.first.y());
      A += area;
    }
    return A;
  }

  Point centroid_polygon(std::vector<Point> points){
    std::vector<std::pair<Point,Point> > edges;
    for(int i=0;i<points.size()-1;i++){
      edges.push_back(std::make_pair(points[i],points[i+1]));
    }
    edges.push_back(std::make_pair(points.back(),points[0]));
    FT Cx=0,Cy=0,A=0;
    for(auto e : edges){
      FT area = 0.5*(e.first.x()*e.second.y()-e.second.x()*e.first.y());
      A += area;
      Cx += 1.0/3*(e.first.x()+e.second.x())*area;
      Cy += 1.0/3*(e.first.y()+e.second.y())*area;
    }
    Cx/=A;
    Cy/=A;
    return Point(Cx,Cy);
  }

  std::pair<Point,FT> centroid_voronoi_cell_in_domain(Vertex_handle v,
						  double lloyd_speed=100,
						  bool calculate_inertia=false,
						  bool control_oscillation=false,
						  FT* new_speed = NULL){
    Face_circulator fc =this->Dt2::incident_faces(v);
    Face_circulator start(fc);
    std::vector<Point> points;
    do{
      if(this->Dt2::is_infinite(fc)){
	typename TDS::Edge e = std::make_pair(fc,fc->index(this->Dt2::infinite_vertex()));
	CGAL::Object obj = this->Dt2::dual(e);
	Ray ray;
	if(CGAL::assign(ray,obj))
	{
	  Vector v = ray.to_vector();
	  FT min_coord = std::min(std::abs(v.x()),std::abs(v.y()))+0.1;
	  Vector v_long(v.x()/min_coord*1000.0,v.y()/min_coord*1000.0);
	  Point far_away = ray.point(0) + v_long;
	  points.push_back(far_away);
	}else{
	  throw std::logic_error("Should be a ray!");
	}
      }else{
	Point p = this->Dt2::dual(fc);
	points.push_back(p);
      }
    }while(++fc != start);

    Polygon_exact Q;
    for(auto p : points){
      Q.push_back(Point_exact(p.x(),p.y()));
    }

    if(!Q.is_simple()){
      for(int i=0;i<Q.size();i++){
	std::cerr<<Q.vertex(i).x()<<' '<<Q.vertex(i).y()<<std::endl;
      }
      throw std::logic_error("Q should be a simple polygon!");
    }
    if (!(CGAL::do_intersect(this->border, Q))){
      for(int i=0;i<Q.size();i++){
	std::cerr<<Q.vertex(i).x()<<' '<<Q.vertex(i).y()<<std::endl;
      }
      throw std::logic_error("They should intersect!");
    }

    Pwh_list intR;
    CGAL::intersection(this->border, Q, std::back_inserter(intR));

    auto selected_poly = intR.end();
    if(intR.size()>1){
      //throw std::logic_error("Shouldn't have more than one polygon");
      double maxSqDist=1e10;
      for(auto poly=intR.begin();poly!=intR.end();poly++){

	Polygon_exact boundary = poly->outer_boundary();
	double cx=0,cy=0;
	for(std::size_t i=0;i<boundary.size();i++){
	  Point p = convert_point_inexact(boundary.vertex(i));
	  cx+=p.x();cy+=p.y();
	}
	cx/=boundary.size();
	cy/=boundary.size();
	FT x = v->point().x();
	FT y = v->point().y();
	FT sqDist = (x - cx)*(x - cx)+(y - cy)*(y - cy);
	if(sqDist<maxSqDist){
	  selected_poly=poly;
	  maxSqDist=sqDist;
	}
      }
    }else{
      selected_poly=intR.begin();
    }
    assert(selected_poly!=intR.end());


    Polygon_exact Q_cut = selected_poly->outer_boundary();
    std::vector<Point> cut_points;
    for(std::size_t i=0;i<Q_cut.size();i++){
      Point p = convert_point_inexact(Q_cut.vertex(i));
      cut_points.push_back(p);
    }
    Point target = centroid_polygon(cut_points);

    FT x0 = v->last_x;
    FT y0 = v->last_y;
    FT x = v->point().x();
    FT y = v->point().y();
    FT cx = target.x();
    FT cy = target.y();
    Vector lastmove(x-x0,y-y0);
    Vector move(cx-x,cy-y);
    if(control_oscillation){
      if(lastmove*move > 0.0){
	lloyd_speed=std::min(step_parameter.max_speed,v->speed+step_parameter.gain_step);
      }else{
	lloyd_speed=step_parameter.default_speed;
      }
      *new_speed = lloyd_speed;
      target = v->point() + move * lloyd_speed / 100.;
    }else{
      target = v->point() + move * lloyd_speed / 100.;
    }
    target = target + lloyd_speed / 100. * step_parameter.moment * lastmove;


    // if centroid is outside the domain
    int tries=0;
    while(CGAL::ON_BOUNDED_SIDE !=
	  CGAL::bounded_side_2(cut_points.begin(),cut_points.end(),target,Kernel())){
      target=Point( (target.x()+v->point().x())/2. , (target.y()+v->point().y())/2. ) ;
      if(tries++ > 10) break;
    }

    FT inertia = 0;
    if(calculate_inertia){
      std::vector<Vector> vector_relative;
      for(auto i = cut_points.begin(); i!=cut_points.end() ; i++){
	vector_relative.push_back(*i - v->point());
      }
      inertia += inertia_polygon(vector_relative);
    }
    return std::make_pair(target,inertia);
  }

  void initial_unit_circle(){
    this->Dt2::clear();
    int N1=30;
    for(int i=0;i<N1;i++){
      FT a= r()*2*3.1415926;
      this->Dt2::insert(Point(cos(a),sin(a)));
    }
    int i=0;
    int N2=150;
    while(i<N2){
      FT x=r()*2-1;
      FT y=r()*2-1;
      if(x*x+y*y>1)continue;
      i++;
      this->Dt2::insert(Point(x,y));
    }
  }

  void initial_uniform(){
    this->Dt2::clear();
    double xmin=1e10,xmax=-1e10,ymin=1e10,ymax=-1e10;
    for(int i=0;i<this->border.size();i++){
      Point p = convert_point_inexact(border.vertex(i));
      FT x = p.x();
      FT y = p.y();
      xmin=std::min(xmin,x);
      xmax=std::max(xmax,x);
      ymin=std::min(ymin,y);
      ymax=std::max(ymax,y);
    }

    int N=180;
    int i=0;
    while(i<N){
      double x=r()*(xmax-xmin)+xmin;
      double y=r()*(ymax-ymin)+ymin;
      if(this->border.bounded_side(Point_exact(x,y))==CGAL::ON_BOUNDED_SIDE){
	this->Dt2::insert(Point(x,y));
	i++;
      }
    }
  }

  void draw_centroid(){
    if (this->Dt2::dimension() != 2) return;
    ::glColor3ub(255, 0, 0);
    ::glLineWidth(1.0);
    ::glBegin(GL_LINES);
    for(auto v = this->Dt2::finite_vertices_begin();v!=this->Dt2::finite_vertices_end();v++){
      Point center = this->centroid_voronoi_cell_in_domain(v).first;
      //::glVertex2f(center.x(),center.y());

      ::glVertex2f(center.x()-0.12,center.y());
      ::glVertex2f(center.x()+0.12,center.y());
      ::glVertex2f(center.x(),center.y()-0.12);
      ::glVertex2f(center.x(),center.y()+0.12);
    }
    ::glEnd();
  }

  void energy_monte_carlo(FT analytic_answer=0,bool monte_carlo = false, FT* maxLengthToWrite = NULL){
    static std::string file = "./energy.txt";
    static bool first_call = true;
    static double xmin=1e10,xmax=-1e10,ymin=1e10,ymax=-1e10;
    if(first_call){
      first_call=false;
      for(int i=0;i<this->border.size();i++){
	Point p = convert_point_inexact(border.vertex(i));
	FT x = p.x();
	FT y = p.y();
	xmin=std::min(xmin,x);
	xmax=std::max(xmax,x);
	ymin=std::min(ymin,y);
	ymax=std::max(ymax,y);
      }
    }

    std::fstream fs;
    fs.open(file, std::fstream::out | std::fstream::app);
    if(monte_carlo){
      for(int count=0;count<1;count++){
	int N=10000;
	double energy=0;
	int i=0;
	while(i<N){
	  double x=r()*(xmax-xmin)+xmin;
	  double y=r()*(ymax-ymin)+ymin;
	  if(this->border.bounded_side(Point_exact(x,y))==CGAL::ON_BOUNDED_SIDE){
	    i++;
	    double minSqDist = 1e10;
	    for(auto v = this->Dt2::finite_vertices_begin();v!=this->Dt2::finite_vertices_end();v++){
	      double xp = v->point().x();
	      double yp = v->point().y();
	      double sqDist = (x-xp)*(x-xp)+(y-yp)*(y-yp);
	      if(sqDist<minSqDist){
		minSqDist=sqDist;
	      }
	    }
	    energy+=minSqDist;
	  }
	}
	energy*=this->border_surface/N;

	fs<<energy<<"\t";
      }
    }
    fs<<analytic_answer<<"\t";
    if(maxLengthToWrite){
      fs<<*maxLengthToWrite;
    }
    fs<<std::endl;
    fs.close();
  }

  void lloyd(double lloyd_speed,
	     bool monte_carlo = false,
	     bool analytic_energy = true,
	     bool max_vec = true,
	     bool control_oscillation = false){

    if (this->Dt2::dimension() != 2) return;

    FT inertia=0;
    FT new_speed;
    std::vector<Point> newPoints;
    std::vector<FT> new_speeds;
    std::vector<Point> oldPoints;
    for(auto v = this->Dt2::finite_vertices_begin();v!=this->Dt2::finite_vertices_end();v++){
      std::pair<Point,FT> p = this->centroid_voronoi_cell_in_domain(v,lloyd_speed,true,control_oscillation,&new_speed);
      Point center = p.first;
      inertia += p.second;
      newPoints.push_back(center);
      new_speeds.push_back(new_speed);
      oldPoints.push_back(v->point());
    }

    FT maxLengthSq=0;
    FT maxLength;
    FT* pLength=NULL;
    if(max_vec){
      int i=0;
      for(auto v = this->Dt2::finite_vertices_begin();v!=this->Dt2::finite_vertices_end();v++){
	FT dx = newPoints[i].x()-v->point().x();
	FT dy = newPoints[i].y()-v->point().y();
	maxLengthSq=std::max(maxLengthSq, dx*dx+dy*dy);
	i++;
      }
      maxLength = std::sqrt(maxLengthSq);
      pLength=&maxLength;
    }

    this->energy_monte_carlo(inertia,monte_carlo,pLength);

    this->Dt2::clear();
    int i=0;
    for(auto p : newPoints){
      auto v = this->Dt2::insert(p);
      v->last_x = oldPoints[i].x();
      v->last_y = oldPoints[i].y();
      v->speed = new_speeds[i];
      i++;
    }


  }

  void draw_border(){
    ::glColor3ub(50, 200, 50);
    ::glLineWidth(2.0);
    ::glBegin(GL_LINE_LOOP);
    for(std::size_t i=0;i<border.size();i++){
      Point p = convert_point_inexact(border.vertex(i));
      ::glVertex2f(p.x(),p.y());
    }
    ::glEnd();
  }

  void insert_with_check(const Point& p){
    if(border.bounded_side(Point_exact(p.x(),p.y()))==CGAL::ON_BOUNDED_SIDE){
      this->insert(p);
    }
    else
      return;
  }
};

#endif // TRIANGULATION_H
