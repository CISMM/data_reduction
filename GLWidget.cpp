/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


#include <QtGui/QMouseEvent>
#include <GL/glu.h>

#include <stdio.h>
#include <iostream>
#include <dr_data.h>

//extern QLabel* label;

using namespace std;

//GLWidget::GLWidget(vector <QImage> *images, QWidget *parent) : QGLWidget(parent) {
GLWidget::GLWidget(dr_group_image_buffer *img_buffer, QWidget *parent) : QGLWidget(parent) {
    cout<<"GLWidget constructor..."<<endl;

    setMouseTracking(true);
    //data.load(filename);

    //this->images = images; //install the data set
    //data = this->images->at(0);
    this->img_buffer = img_buffer;
    this->img_buffer->get_frame(0, &data);


    gldata = QGLWidget::convertToGLFormat(data);
    resize(data.size());
}

void GLWidget::initializeGL() {
    cout<<"initializeGL..."<<endl;
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(0, 0, 0, 0);
}

void GLWidget::resizeGL(int w, int h) {
    cout<<"resizeGL..."<<endl;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void GLWidget::paintGL() {
    cout<<"paintGL..."<<endl;
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1,0,0);
    /*
    glBegin(GL_POLYGON);
    glVertex2f(0,0);
    glVertex2f(100,500);
    glVertex2f(500,100);
    glEnd();
    */
    glDrawPixels(data.width(), data.height(), GL_RGBA, GL_UNSIGNED_BYTE, gldata.bits());
    glFlush();
}

void GLWidget::mousePressEvent(QMouseEvent *event) {
    cout<<"mousePressEvent..."<<endl;
    emit mousePressed();
}

void GLWidget::mouseMoveEvent(QMouseEvent *event) {
    //printf("%d, %d\n", event->x(), event->y());

    char buf[1024];
    sprintf(buf, "%d, %d\n", event->x(), event->y());
    QString str(buf);
    emit mouseMoved(str);
}


void GLWidget::keyPressEvent(QKeyEvent* event) {
    switch(event->key()) {
    case Qt::Key_Escape:
        close();
        break;
    default:
        event->ignore();
        break;
    }
}



