#include <QtGui>
#include "glviewer.h"
#include <QTime>

GlViewer::GlViewer(QWidget *pParent)
	: QGLWidget(QGLFormat(QGL::SampleBuffers), pParent)
{
	m_scene = NULL;
	m_center_x = m_center_y = 0.0;
	m_scale = 0.3;
        is_ready = true;
        //Will tick as soon as all the events are finished processing
        startTimer(0);
	setAutoFillBackground(false);
	setMouseTracking(true);
}
//Qt5 stacks the events, so you need this to avoid being flooded by the mouseEvents, which makes you lagg a lot.
void GlViewer::timerEvent(QTimerEvent* /*event*/)
{
  is_ready = true;
}
void GlViewer::resizeGL(int width, int height) {
	glViewport(0, 0, width, height);
	double aspect_ratio = double(height) / double(width);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.0, 1.0, -aspect_ratio, aspect_ratio, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void GlViewer::initializeGL() {
	glClearColor(1., 1., 1., 0.);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_SMOOTH);
}

void GlViewer::paintGL() {
	glClear(GL_COLOR_BUFFER_BIT);
	if (!m_scene) return;

	glPushMatrix();
	glScaled(m_scale, m_scale, m_scale);
	glTranslated(-m_center_x, -m_center_y, 0.0);
	m_scene->render();
	glPopMatrix();
}

void GlViewer::wheelEvent(QWheelEvent *event) {
	if (!m_scene) return;
	double angle = event->delta() / 120;
	m_scale*=1+angle*0.1;
	//m_scale += 0.05 * (event->delta() / 120);
	if (m_scale <= 0.0) m_scale = 0.0;
	updateGL();
}

void GlViewer::mousePressEvent(QMouseEvent *event) {
  if(is_ready)
  {
      is_ready = false;
	if (!m_scene) return;
	m_mouse_click = event->pos();

	if (event->button() == Qt::LeftButton) 
	{
		setCursor(QCursor(Qt::PointingHandCursor));
		sample_mouse_path(m_mouse_click);
	} 
	else 
	{
		setCursor(QCursor(Qt::ClosedHandCursor));
	}
  }
}

void GlViewer::mouseMoveEvent(QMouseEvent *event)
{
    if(is_ready)
    {
        is_ready = false;
	if(!m_scene) return;    
	m_mouse_move = event->pos();

	if (event->buttons() == Qt::LeftButton)
	{
		if (m_mouse_move != m_mouse_click)
			sample_mouse_path(m_mouse_move);
	}
	if (event->buttons() == Qt::RightButton)
	{        
		move_camera(m_mouse_click, m_mouse_move);
	}

	// store mouse pos in scene
	double x, y;
	convert_to_world_space(m_mouse_move, x, y);
	m_scene->mouse_pos() = Scene::Point(x, y);
	std::cout << "Mousemove event (" << x << "; " << y << ")" << std::endl;

	m_mouse_click = m_mouse_move;
	updateGL();
    }
}

void GlViewer::mouseReleaseEvent(QMouseEvent *event) {
	if (!m_scene) return;
	m_mouse_move = event->pos();

	if (event->button() == Qt::LeftButton) 
	{
		if (m_mouse_move != m_mouse_click)
			sample_mouse_path(m_mouse_move);
	} 
	else 
	{
		move_camera(m_mouse_click, m_mouse_move);
	}    



	m_mouse_click = m_mouse_move;

	setCursor(QCursor(Qt::ArrowCursor));
	updateGL();
}

void GlViewer::sample_mouse_path(const QPoint& point)
{
	double x, y;
	convert_to_world_space(point, x, y);

	m_scene->add_point(Scene::Point(x, y));
}

void GlViewer::move_camera(const QPoint& p0, const QPoint& p1)
{
	m_center_x -= double(p1.x() - p0.x()) / double(width())/m_scale;
	m_center_y += double(p1.y() - p0.y()) / double(height())/m_scale;
}

void GlViewer::convert_to_world_space(const QPoint& point, double &x, double &y)
{
	double aspect_ratio = double(height()) / double(width());

	x = double(point.x()) / double(width());
	x = (2.0*x - 1.0) / m_scale;
	x += m_center_x;

	y = 1.0 - double(point.y()) / double(height());
	y = (2.0*y - 1.0) * aspect_ratio / m_scale;
	y += m_center_y;
}
