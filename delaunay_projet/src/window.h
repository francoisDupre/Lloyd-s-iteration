#ifndef _WINDOW_
#define _WINDOW_

// std
#include <list>
#include <QWidget>

// Qt
#include <QString>
#include "scene.h"
#include "ui_delaunay.h"

class MainWindow : public QMainWindow, public Ui_MainWindow
{
  Q_OBJECT

private:
  Scene* m_scene;

  unsigned int maxNumRecentFiles;
  QAction* recentFilesSeparator;
  QVector<QAction*> recentFileActs;
  bool write_image = true;
  int lloyd_count;
public:
  MainWindow();
  ~MainWindow();

  void update();
  void button_state_update();
protected:

protected slots:

  // drag & drop
  void closeEvent(QCloseEvent *event);

public slots:

  // Data
  void on_actionClear_triggered();

  // View
  void on_actionView_Delaunay_edges_toggled();
  void on_actionView_Voronoi_edges_toggled();
  void on_actionView_Convex_hull_toggled();
  void on_actionView_Face_under_mouse_toggled();
  void on_actionView_Line_walk_toggled();
  void on_actionView_Nearest_vertex_toggled();
  void on_actionView_Cell_centroid_toggled();
  void on_actionWrite_Image_toggled();
  void on_actionMonteCarlo_Energy_toggled();
  void on_actionAnalytic_Energy_toggled();
  void on_actionWrite_Max_Movement_toggled();

  void on_actionControl_Oscillation_toggled();

  void on_lloyd_speed_valueChanged(double value);
  void on_default_speed_valueChanged(double value);
  void on_max_speed_valueChanged(double value);
  void on_speed_step_valueChanged(double value);
  void on_speed_moment_valueChanged(double value);

  // algorithms
  void on_actionLloyd_triggered();
  void on_actionUnitCircle_triggered();
  void on_actionUniform_points_triggered();

  //border
  void on_actionHexagon_triggered();
  void on_actionHourglass_large_triggered();
  void on_actionHourglass_small_triggered();
  void on_actionHourglass_triggered();
  void on_actionSnake_triggered();
  void on_actionSnake_soft_triggered();

  //buttoms
  void on_Lloyd1_clicked();
  void on_Lloyd10_clicked();
  void on_Lloyd100_clicked();
  void on_LloydTo1000_clicked();
};

#endif // _WINDOW_
