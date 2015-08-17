/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


#include "dr_data.h"
#include <iostream>
#include <string>
#include <QMap>
#include <QString>


void main_helper::read_config_files(QMap<QString, QString> * key_val, QString * file_name) {
    //  this function reads config file and fill in the QMap that suppose to contain the key value pairs.
    QFile file(*file_name);
    if (file.open(QIODevice::ReadOnly)) {
        //  QMessageBox::information(0, "error in reading config file: ", file.errorString());

        QTextStream in(&file);
      //  cout << "in is at end:" << in.atEnd() << endl;
        while(!in.atEnd()) {
            QString line = in.readLine();
            // analysis the line,
            //   if it starts with a "#", ignore this line
            //   if it starts with a ">", have not decided what to do; right now just ignore it
            //   otherwise, parse
            if (line.at(0) == '#' || line.at(0) == '>' ) {
              //  cout << "ignore" << endl;
            }
            else {
                QStringList fields = line.split(" = ");
        //        cout << "inseted xxx:" << " " << fields.length() << endl;
                key_val->insert(fields.at(0), fields.at(1));}
        }
    }
    file.close();
}

