/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


#ifndef _GLWIDGET_H
#define _GLWIDGET_H

/*
#include <QtGui>
#include <QtGui/QMouseEvent>
#include <QtGui/QImage>
#include <QtGui/QLabel>
#include <QtOpenGL/QGLWidget>

#include <QImage>

#include <stdio.h>

#include <vector>

using namespace std;

extern class dr_group_image_buffer;


class GLWidget : public QGLWidget {

	Q_OBJECT // must include this if you use Qt signals/slots


	public:
		QImage data, gldata;
		//GLWidget(vector <QImage> *images, QWidget *parent = NULL);
		GLWidget(dr_group_image_buffer *img_buffer, QWidget *parent = NULL);
	private:
		//vector <QImage> *images;
		dr_group_image_buffer *img_buffer;

	protected:
		void initializeGL();
		void resizeGL(int w, int h);
		void paintGL();
		void mousePressEvent(QMouseEvent *event);
		void mouseMoveEvent(QMouseEvent *event);
		void keyPressEvent(QKeyEvent *event);

signals:
		void mouseMoved(QString str);

		public slots:
};






class video_analysis:public QWidget {
	private:
		GLWidget *video;
		QLabel *label;
	public:
		//video_analysis(vector <QImage> *images) 
		video_analysis(dr_group_image_buffer *img_buffer) 
		{
			label = new QLabel("Coordinates");
			label->resize(200,100);
			video = new GLWidget(img_buffer, NULL);

			QObject::connect(video, SIGNAL(mouseMoved(QString)), label, SLOT(setText(QString)));
			label->show();
			video->show();
		}
};



*/


#endif  /* _GLWIDGET_H */
