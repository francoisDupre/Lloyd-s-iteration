#include <QtGui>
#include "window.h"

int main(int argv, char **args)
{	
	QApplication app(argv, args);
	app.setApplicationName("2D Delaunay triangulation");
	MainWindow window;
	window.show();
	return app.exec();
}
