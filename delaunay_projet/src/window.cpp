// Qt
#include <QtGui>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QClipboard>
#include <QPixmap>
#include <QString>

#include "window.h"

#include <fstream>
#include <iostream>
MainWindow::MainWindow() : 
  QMainWindow(), Ui_MainWindow(),
  maxNumRecentFiles(15), recentFileActs(15)
{
  setupUi(this);

  // init scene
  m_scene = new Scene;
  viewer->set_scene(m_scene);
  button_state_update();

  // Handling actions
  connect(actionQuit, SIGNAL(triggered()), this, SLOT(close()));

  lloyd_count=1;
}

MainWindow::~MainWindow()
{
  delete m_scene;
}

void MainWindow::update()
{
  viewer->repaint();
}

void MainWindow::button_state_update()
{
  m_scene->set_view_Delaunay_edges(actionView_Delaunay_edges->isChecked()) ;
  m_scene->set_view_Voronoi_edges(actionView_Voronoi_edges->isChecked());
  m_scene->set_view_Convex_hull(actionView_Convex_hull->isChecked());
  m_scene->set_view_Face_under_mouse(actionView_Face_under_mouse->isChecked());
  m_scene->set_view_Line_walk(actionView_Line_walk->isChecked());
  m_scene->set_view_Nearest_vertex(actionView_Nearest_vertex->isChecked());
  m_scene->set_view_Cell_centroid(actionView_Cell_centroid->isChecked());
  m_scene->change_lloyd_speed(lloyd_speed->value());
  write_image = actionWrite_Image->isChecked();
  m_scene->set_monte_carlo(actionMonteCarlo_Energy->isChecked());
  m_scene->set_analytic_energy(actionAnalytic_Energy->isChecked());
  m_scene->set_max_vec(actionWrite_Max_Movement->isChecked());
  m_scene->set_control_oscillation(actionControl_Oscillation->isChecked());

  bool control_oscillation = actionControl_Oscillation->isChecked();
  lloyd_speed->setEnabled(!control_oscillation);
  default_speed->setEnabled(control_oscillation);
  max_speed->setEnabled(control_oscillation);
  speed_step->setEnabled(control_oscillation);

  StepParam sp;
  sp.default_speed = default_speed->value();
  sp.max_speed = max_speed->value();
  sp.gain_step = speed_step->value();
  sp.moment = speed_moment->value();
  m_scene->set_step_param(sp);

}


void MainWindow::closeEvent(QCloseEvent *event)
{
  event->accept();
}



void MainWindow::on_actionClear_triggered()
{
  m_scene->clear();
  update();
}

void MainWindow::on_actionView_Delaunay_edges_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Voronoi_edges_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Convex_hull_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Face_under_mouse_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Line_walk_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Nearest_vertex_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionView_Cell_centroid_toggled()
{
  button_state_update();
  update();
}

void MainWindow::on_actionWrite_Image_toggled()
{
  button_state_update();
}

void MainWindow::on_actionMonteCarlo_Energy_toggled()
{
  button_state_update();
}

void MainWindow::on_actionAnalytic_Energy_toggled()
{
  button_state_update();
}

void MainWindow::on_actionWrite_Max_Movement_toggled()
{
  button_state_update();
}

void MainWindow::on_actionControl_Oscillation_toggled()
{
  button_state_update();
}

void MainWindow::on_lloyd_speed_valueChanged(double value)
{
  button_state_update();
}

void MainWindow::on_default_speed_valueChanged(double value)
{
  button_state_update();
}

void MainWindow::on_max_speed_valueChanged(double value)
{
  button_state_update();
}

void MainWindow::on_speed_step_valueChanged(double value)
{
  button_state_update();
}

void MainWindow::on_speed_moment_valueChanged(double value)
{
  button_state_update();
}

void MainWindow::on_actionLloyd_triggered()
{

  m_scene->lloyd();
  update();
  if(write_image){
    QPixmap pixmap = viewer->renderPixmap();
    QString filename = QString("renderings/image_%1.jpg").arg(this->lloyd_count, 4, 10, QChar('0'));
    if(!pixmap.save(filename,"JPG")){
      std::cout<<"Failed to save image to : "<<filename.toStdString()<<std::endl;
    }
  }
  this->lloyd_count++;
}

void MainWindow::on_actionUnitCircle_triggered()
{
  m_scene->initial_unit_circle();
  update();
}

void MainWindow::on_actionUniform_points_triggered()
{
  m_scene->initial_uniform();
  update();
}

void MainWindow::on_actionHexagon_triggered()
{
  m_scene->set_border_shape(Triangulation::HEXAGON);
  update();
}

void MainWindow::on_actionHourglass_large_triggered()
{
  m_scene->set_border_shape(Triangulation::HOURGLASS_LARGE);
  update();
}

void MainWindow::on_actionHourglass_small_triggered()
{
  m_scene->set_border_shape(Triangulation::HOURGLASS_SMALL);
  update();
}

void MainWindow::on_actionHourglass_triggered()
{
  m_scene->set_border_shape(Triangulation::HOURGLASS);
  update();
}

void MainWindow::on_actionSnake_triggered()
{
  m_scene->set_border_shape(Triangulation::SNAKE);
  update();
}

void MainWindow::on_actionSnake_soft_triggered()
{
  m_scene->set_border_shape(Triangulation::SOFTSNAKE);
  update();
}

void MainWindow::on_Lloyd1_clicked()
{
  on_actionLloyd_triggered();
}

void MainWindow::on_Lloyd10_clicked()
{
  for(int i=0;i<10;i++)on_actionLloyd_triggered();
}

void MainWindow::on_Lloyd100_clicked()
{
  for(int i=0;i<100;i++)on_actionLloyd_triggered();
}

void MainWindow::on_LloydTo1000_clicked()
{
  while(this->lloyd_count<=1000){
    on_actionLloyd_triggered();
  }
}

