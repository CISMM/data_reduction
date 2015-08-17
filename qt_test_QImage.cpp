/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


#include <iostream>
#include <fstream>
#include <stdint.h>
#include <stdio.h>

#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "dr_data.h"
#include "file_list.h"
#include <QApplication>

#include <QLabel>
#include <QtGui/QApplication>

using namespace std;

int NUMBER_OF_SEGS = 1;
double KS_ALPHA = 0.05;
int DILATE_DIM = 7;
int BUFFER_SIZE = 100;
bool ONLINE_CALC_MEAN_VAR = true;
int THRESHOLD_FOR_CORRELATION = 95;

bool DO_SLIDING_WINDOWING = false;
int WINDOW_SIZE = 1;
//int WINDOW_SIZE = 1;
int HARD_CODED_THRESHOLD = 40;

int WINDOW_METHOD_VAR = 1;

// video with identical frames, theoretical limit
// 0           1                 2   3                              4                 5                           6                       7                      8               9
// test_QImage "save file name"  -i "image name"                    "dilate diameter" "threshold for correlation" "do sliding windowing?" "hard_coded_threshold" "save_file_dir" "debugging_image_dir"
// added a directory containing video frames
// test_QImage "save file name"  -d "image directory with datalist" "dilate diameter" "threshold for correlation" "do sliding windowing?" "hard_coded_threshold" "save_file_dir" "debugging_image_dir"

// notice that "threshold_for_correlation" is not used when 'USE_HARD_CODED_THRESHOLD' is true
int main(int argc, char *argv[])
{
    int ret;
    QCoreApplication app(argc, argv);

    // *****************************
    //    parsing input arguments
    // *****************************
    if (argc<4) {
        cerr<<"Usage: " << argv[0] << " [save_file] -i [image_file] [dilate_diameter] [threshold_for_correlation_in_percentage] [sliding_window_size] [hard_coded_threshold] [save_file_dir] [debugging_image_dir]"<<endl;
        cerr<<"or: ./test [save_file] -d [video directory]  [dilate_diameter] [threshold_for_correlation_in_percentage] [sliding_window_size] [hard_coded_threshold] [save_file_dir] [debugging_image_dir]"<<endl;
        return -1;
    }

    int buffer_size;

    if (argc >= 5) {
        cout << "dilation diameter was defined in command line opt \n";
        DILATE_DIM = atoi(argv[4]);
        cout << DILATE_DIM << "\n";
    }
    if (argc >= 6) {
        int input_correlation_threshold = atoi(argv[5]);
        if (input_correlation_threshold >= 0 && input_correlation_threshold <= 100) {
            cout << "set correlation threshold: " << static_cast<double>(input_correlation_threshold) / 100.0 << "% percentile" << endl;
            THRESHOLD_FOR_CORRELATION = input_correlation_threshold;
        }
        else {
            cout << "the input correlation threshold is not valid. use the default setting: " << THRESHOLD_FOR_CORRELATION << endl;
        }
    }

    if (argc >= 7) {
        int is_doing_sliding_windowing = atoi(argv[6]);
        cout << "window size:" << is_doing_sliding_windowing << endl;
        if (is_doing_sliding_windowing >= 2) { // do a sliding windowing
            if (ONLINE_CALC_MEAN_VAR) { // we cannot do sliding windowing and online at the same time, so this should be a conflict
                cout << "previously online methods was chosen, since the sliding window size is larger than 1, switch into using the sliding window method.";
            }
            else {
                cout << "using sliding windowing: true\n";

            }
            DO_SLIDING_WINDOWING = true;
            WINDOW_SIZE = is_doing_sliding_windowing;
            ONLINE_CALC_MEAN_VAR = false;
        }
        else {
            DO_SLIDING_WINDOWING = false;
            ONLINE_CALC_MEAN_VAR = true;
        }
    }

    if (argc >= 8) {
        HARD_CODED_THRESHOLD = atoi(argv[7]);
        cout << "hard coded threshold: " << HARD_CODED_THRESHOLD << endl;
    }
    QString * save_file_dir;
    if (argc >= 9) {
        save_file_dir = new QString(argv[8]);
        cout << "save file dir: " << save_file_dir->toStdString() << endl;
    }
    QString * debug_image_name;
    if (argc >= 10) {
        debug_image_name = new QString(argv[9]);
        cout << "debug image name: " << debug_image_name->toStdString() << endl;
    }

    if ( strcmp(argv[2],"-i") == 0 ) {
        if (argc >= 5) {
            buffer_size = atoi(argv[4]);
        } else if (argc == 4) {
            buffer_size = 100;
        } else {
            cerr<<"Usage: ./test [save_file] -i [image_file] <number>"<<endl;
            return -1;
        }
    } else if ( strcmp(argv[2],"-d") == 0 ) {
        if (argc <= 5) {
            cerr<<"Usage: ./test [save_file] -d [video directory] [number]"<<endl;
            return -1;
        }

        buffer_size = atoi(argv[4]);
    } else {
        cerr<<"Usage: ./test [save_file] -i [image_file] <number>"<<endl;
        cerr<<"or: ./test [save_file] -d [video directory] [number]"<<endl;
        return -1;
    }

    main_helper helpp ;
    QMap<QString, QString> * key_val_map = new QMap<QString, QString>();
    QString * config_file_name = new QString("sample_config.txt");
    helpp.read_config_files(key_val_map, config_file_name);

    dr_group_image_buffer buf(buffer_size, 0);
    QImage * mean_img = NULL;
    QVector<QVector<double> *> * mean_img_double = NULL;
    QImage * variance_img = NULL;
    QImage * new_img = NULL;
    //identical frame test
    QVector<string> * datalist = new QVector<string>();
    QVector<QVector<double> *> * var_img_double = new QVector<QVector<double> *>();

    // check if input image sets are OK
    if ( strcmp(argv[2],"-i") == 0 ) {
        string fileroot(argv[3]);

        if(!file_list(fileroot, datalist)) {
            cerr<<"Attempt to get list of files was unsuccessful"<<endl;
            return -1;
        }

        QImage *image=NULL;

        buf.set_data_list(datalist);

        // *********************
        // without sliding window
        // *********************
        if (ONLINE_CALC_MEAN_VAR) {
            // calculate mean and variance online
            // for n times, n we know,
            // call the function in the inner loop.
            // finally calculate the variance image
            // here we read one image (first one) from the data list
            // and initialize all the stats (mean, variance) based on that image
            // so here we basically don't fill up the image array in dr_group_buffer, instead we calculate the
            // mean and variance online .
            if( strcmp(datalist->at(0).c_str(),"")!=0 ) {
                if (image != NULL)
                    delete(image);
                image = new QImage(datalist->at(0).c_str());
                if( image->isNull() ) {
                    cerr<<datalist->at(0).c_str()<<" doesn't exist."<<endl;
                    return -1;
                }

                int n = 0;
                mean_img = new QImage(*image);
                mean_img_double = new QVector<QVector<double> *> ();
                variance_img = new QImage(*image);
                // M2 is used for updating mean and variance at the same time.
                QVector<QVector<double> *>* M2_img = new QVector<QVector<double> *>();
                int w = image->width();
                int h = image->height();
                int x, y;

                for (y = 0; y < h; y++) {
                    mean_img_double->append(new QVector<double>());
                    M2_img->append(new QVector<double>());
                    for (x = 0; x < w; x++) {
                        mean_img_double->at(y)->append(0.0);
                        variance_img->setPixel(x,y,0);
                        M2_img->at(y)->append(0.0);
                    }
                }

                int i;
                for ( i = 0; i < datalist->size(); i++) {
                    new_img = new QImage(datalist->at(i).c_str());
                    // repeat n time to get the result
                    n = buf.online_mean_and_variance(n, mean_img_double, M2_img, new_img);
                }

                int curr_mean_pixel = 0;
                double curr_var_pixel_double = 0.0;
                for (y = 0; y < h; y++) {
                    var_img_double->append(new QVector<double>());
                    for (x = 0; x < w; x++) {
                        curr_mean_pixel = static_cast<int>(mean_img_double->at(y)->at(x));
                        mean_img->setPixel(x,y,curr_mean_pixel);
                        curr_var_pixel_double = M2_img->at(y)->at(x) / static_cast<double>(n-1);
                        var_img_double->at(y)->append(curr_var_pixel_double);
                    }
                }
                mean_img->save(*(debug_image_name) + "mean_img.pgm", 0, 100);
                QImage * scaled_variance_img = new QImage(*(variance_img));
                buf.scale_image(var_img_double, scaled_variance_img);

                scaled_variance_img->save(*(debug_image_name) + "variance_img.pgm", 0, 100);

                // set the mean_img and variance_img into the buffer object
                buf.need_update_state = false;
                buf.set_mean_img(mean_img);

                buf.set_variance_img(scaled_variance_img);
            }
        }

        else {
            if (DO_SLIDING_WINDOWING) {
                // we now have a buffer
                //   -> involves rewriting dr_group_buffer class,
                //   -> in a for loop, call (update for / background) method
                //      -> by doing this we can determine the time point where movement starts
                // for the compression. we now have the following choices:
                //   -> we can do a second pass to make the fore / background mapping for static time perioids and don't compress the moving period
                //   -> or we can just compute the correlation in the window, determine the fore/background, and save the thrown image based on the result
                //   -> or in the begging we save the first half window frames as indicated, and when go one step, we save the old middle image

                // pseudo code:
                //   -> setting the N for for loop
                //   -> reading first m/2 images,
                //   -> start poping out image until reach m
                //   -> in the for loop read and pop until at the end
                //   -> in the last loop  (m / 2) pop
                //   -> DONE

                buf.set_do_sliding_windowing(true);
                buf.set_sliding_window_size(WINDOW_SIZE);
            }
        }
    }
    else {
        cout << "should never reach here! check input parameter settings!" << endl;
    }

    //getting statistics
    //FROM ME: Obtain mean, min, max, and difference images.
    cout<<"Computing statistics of the group..."<<endl;
    dr_group_stat group_stat;
    if (!ONLINE_CALC_MEAN_VAR ) {
        if (!DO_SLIDING_WINDOWING) {
            // DO NOTHING
            //     buf.need_update_state = true;
            //      buf.get_group_stat(&group_stat);
            cout << "do nothing, not calculating stats.." << endl;
        }
        else { // do it with a sliding windowing
            buf.need_update_state = false;
            /// THERE ARE TWO VERSIONS OF IT:
            /// THE FIRST PART WILL BE: (WINDOW_METHOD_VAR = 0)
            ///     for the boundaries, pop every old frame once reading another new frame when window center goes after the 1/2 window size in the begining / end
            ///
            /// THE SECOND VARIATION IS: (WINDOW_METHOD_VAR = 1) <- currently we should use this one
            ///     when window arrives the first 1/2 window size, pop all the frames using the statistsics in the window
            ///     when right side of the window arrives the end, pop all the frames

            // we now have a buffer
            //   -> involves rewriting dr_group_buffer class,
            //   -> in a for loop, call (update for / background) method
            //      -> by doing this we can determine the time point where movement starts
            // for the compression. we now have the following choices:
            //   -> we can do a second pass to make the fore / background mapping for static time perioids and don't compress the moving period
            //   -> or we can just compute the correlation in the window, determine the fore/background, and save the thrown iamge basesd on the result
            //   -> or in the begging we save the first half window frames as indicated, and when go one step, we save the old middle image

            // pseudo code:
            //   -> first online compute mean and variance ( does not differ from the online verion )
            //   -> setting the N for for loop
            //   -> reading first m/2 images,
            //   -> start poping out image until reach m
            //   -> in the for loop read and pop until at the end
            //   -> in the last loop  ( m / 2) pop
            //   -> DONE

            // at this time the buffer contains no image

            //     int number_of_iterations = datalist->size() - WINDOW_SIZE;

            // recommend: WINDOW_SIZE to be an even number
            // reading first m/2 images

            if (WINDOW_METHOD_VAR == 0) {

                int i;
                for (i = 0; i < WINDOW_SIZE / 2; i++) {
                    buf.fill_window_buffer(new QImage(datalist->at(i).c_str()));
                }
                /// TODO: here print things like:
                ///     taking image no .xxx
                ///     poping out image no .xxx
                ///
                ///
                /// ALSO save the bin_map for every window movement
                buf.init_window_stats();
                // then we can use the statistics to compute the correlation, save the image and read new image
                // I/O:
                //     at where ( which frame is readed )
                //     what has been saved ( which frame has been saved )
                //     at this time, the foreground / background ratio
                //

                // after you fill in half of the buffers, you can then use the image to do statistics
                // before read N/2 + 1 image, you can save 0th imabe by the N/2-statistics
                // the 1th image is based on the N/2 +1 statistics
                // the N/2'th image is based on N/2 image statistics <- some subtle things here, consider later, ( now we assume the frame number >> window size )
                QImage * out_img = new QImage(datalist->at(i).c_str());
                QString * out_img_file_name = new QString(*(save_file_dir)+"compressed_");
                QString * in_img_file_name;
                int j = 0;
                for (i = WINDOW_SIZE / 2; i < WINDOW_SIZE; i++) {

                    in_img_file_name = new QString(datalist->at(i).c_str());
                    out_img_file_name->append(QString::number(j));
                    out_img_file_name->append(".pgm");
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, j, debug_image_name);
                    out_img->save(*out_img_file_name,0,100);
                    out_img_file_name->remove(".pgm");
                    // out_img_file_name->remove(QString::number(j));
                    out_img_file_name->remove(out_img_file_name->size() - QString::number(j).size(), QString::number(j).size() );

                    buf.fill_window_buffer(new QImage(*in_img_file_name));
                    // no deque, has enque
                    buf.update_window_stats(false, true);
                    j++;
                }
                for (i = WINDOW_SIZE; i < datalist->size(); i++) {
                    in_img_file_name= new QString(datalist->at(i).c_str());
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, j, debug_image_name);
                    out_img_file_name->append(QString::number(j));
                    out_img_file_name->append(".pgm");
                    out_img->save(*out_img_file_name,0,100);
                    out_img_file_name->remove(".pgm");
                    //out_img_file_name->remove(QString::number(j));
                    out_img_file_name->remove(out_img_file_name->size() - QString::number(j).size(), QString::number(j).size() );
                    buf.fill_window_buffer(new QImage(*in_img_file_name));
                    // both deque and enque
                    buf.update_window_stats(true, true);
                    buf.deque_buffer();
                    j++;
                }
                // here; no more images, we don't fil lthe buffer
                for (i = 0; i < WINDOW_SIZE / 2; i++) {
                    //        in_img_file_name = new QString(datalist->at(i).c_str());
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, j,debug_image_name);
                    out_img_file_name->append(QString::number(j));
                    out_img_file_name->append(".pgm");
                    out_img->save(*out_img_file_name,0,100);
                    out_img_file_name->remove(".pgm");
                    //out_img_file_name->remove(QString::number(j));
                    out_img_file_name->remove(out_img_file_name->size() - QString::number(j).size(), QString::number(j).size() );
                    // only deque
                    buf.deque_buffer();
                    buf.update_window_stats(true, false);
                    j++;
                }
                // here; we don't deque anymore since then there will be not enough images for doing stats
                for (i = 0; i < WINDOW_SIZE / 2; i++) {
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, j, debug_image_name);
                    out_img_file_name->append(QString::number(j));
                    out_img_file_name->append(".pgm");
                    out_img->save(*out_img_file_name,0,100);
                    out_img_file_name->remove(".pgm");
                    //out_img_file_name->remove(QString::number(j));
                    out_img_file_name->remove(out_img_file_name->size() - QString::number(j).size(), QString::number(j).size() );
                    j++;
                }
                // DONE.
            }
            else if (WINDOW_METHOD_VAR == 1) {
                int i;
                cout << "---------------------------------- \n" << "datalist size : " << datalist->size() << "\n--------------------------\n";
                // first fill the buffer
                for (i = 0; i < WINDOW_SIZE; i++) {
                    buf.fill_window_buffer(new QImage(datalist->at(i).c_str()));
                    cout << "reading image no: " << i << endl;
                }

                buf.init_window_stats2(debug_image_name);
                // then we can use the statistics to compute the correlation, save the image and read new image

                // ** Avoid memory glitch
                //QImage * out_img = new QImage(datalist->at(i).c_str());
                QImage * out_img;

                QString * out_img_file_name = new QString(*(save_file_dir)+"compressed_");
                QString * in_img_file_name;

                // save all the half frames, not pop them
                QString * foo = new QString("foo");
                for (i = 0; i < WINDOW_SIZE / 2; i++) {
                    out_img = new QImage(datalist->at(i).c_str());
                    QString s = foo->sprintf("%04d",i);
                    out_img_file_name->append(s);
                    out_img_file_name->append(".pgm");

                    // ** Apply a binary threshold on the image and return a modified compressed version. 
                    buf.save_compressed_image_in_queue(out_img,3,DILATE_DIM, i, debug_image_name);
                    out_img->save(*out_img_file_name,0,100);

                    cout << "saving image: " << out_img_file_name->toStdString() << endl;
                    cout << "buffer size: " << buf.get_window_buffer_size() << endl;

                    out_img_file_name->remove(".pgm");
                    // out_img_file_name->remove(s);
                    out_img_file_name->remove(out_img_file_name->size()-s.size(), s.size() );
                    delete out_img;
                }
                cout << "buffer size: " << buf.get_window_buffer_size() << endl;
                // read one, pop one
                for (i = WINDOW_SIZE / 2; i < datalist->size() - WINDOW_SIZE / 2; i++) {
              //      cout << "WINDOW_SIZE / 2: "  << WINDOW_SIZE /2 << endl;
                    out_img = new QImage(datalist->at(i).c_str());
                    in_img_file_name= new QString(datalist->at(i+WINDOW_SIZE/2).c_str());
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, i, debug_image_name);
                    QString s = foo->sprintf("%04d",i);
                    out_img_file_name->append(s);
                    out_img_file_name->append(".pgm");
                    out_img->save(*out_img_file_name,0,100);
                    // TODO: also save the floating point version of the compressed image.


                    cout << "saving image: " << out_img_file_name->toStdString() << endl;
                    cout << "buffer size: " << buf.get_window_buffer_size() << endl;

                    out_img_file_name->remove(".pgm");
                    //out_img_file_name->remove(s);
                    out_img_file_name->remove(out_img_file_name->size()-s.size(), s.size() );
                    buf.fill_window_buffer(new QImage(*in_img_file_name));
                    cout << "reading image: " << in_img_file_name->toStdString() << endl;
                    cout << "buffer size: " << buf.get_window_buffer_size() << endl;

                    // both deque and enque
                    buf.update_window_stats2(i, debug_image_name);
                    buf.deque_buffer();
                    cout << "buffer size: " << buf.get_window_buffer_size() << endl;
                    delete out_img;
                }

                // save everything left
                for (i = datalist->size()-WINDOW_SIZE/2; i < datalist->size(); i++) {
                    out_img = new QImage(datalist->at(i).c_str());
                    buf.save_compressed_image_in_queue(out_img, 3, DILATE_DIM, i, debug_image_name);
                    QString s = foo->sprintf("%04d",i);
                    out_img_file_name->append(s);
                    out_img_file_name->append(".pgm");
                    out_img->save(*out_img_file_name,0,100);
                    // TODO: also save the floating point version of the compressed image.

                    cout << "saving image: " << out_img_file_name->toStdString() << endl;
                    cout << "buffer size: " << buf.get_window_buffer_size() << endl;

                    out_img_file_name->remove(".pgm");
                    //out_img_file_name->remove(s);
                    out_img_file_name->remove(out_img_file_name->size()-s.size(), s.size() );
                    delete out_img;
                }
                // DONE.
            }
        }
    }
    else {  // online calc mean variance
        group_stat.mean_images = new QVector<QImage *>();
        group_stat.mean_images->append(mean_img);
        group_stat.var_images = new QVector<QImage *>();
        group_stat.var_images->append(variance_img);
        // compute correlation image without calling get_group_stat()
        int i;
        QImage * new_img = NULL;
        QVector<QVector<QVector<double> *> *> * old_corr_values = new QVector<QVector<QVector<double> *> *>();
        int h = mean_img->height();
        int w = mean_img->width();
        for ( int y = 0; y < h; y++) {
            old_corr_values->append(new QVector<QVector <double> *>());
            for ( int x =0; x < w; x++) {
                old_corr_values->at(y)->append(new QVector<double>());
                for (int z = 0; z < 9; z++) {
                    old_corr_values->at(y)->at(x)->append(0.0);
                }
            }
        }
        for ( i = 0; i < datalist->size(); i++) {
            new_img = new QImage(datalist->at(i).c_str());
            //  cout << "image: " << datalist->at(i).c_str() << endl;
            buf.compute_corr_image_online(new_img, old_corr_values, mean_img_double);
        }
        // finally we have a vector which contains 8 vectors for 1 pixel
        // we process it to get the correlation image
        QImage * corr_img = new QImage(*new_img);
        for(int y =0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double maximum_corr_rvalue = 0.0;
                double curr_var = var_img_double->at(y)->at(x);
                //        int curr_var = qGray(variance_img->pixel(x,y));
                double corr1, corr2, corr3, corr4, corr5, corr6, corr7, corr8;
                if (curr_var <= 0.001) {
                    curr_var = 1.0;
                }
                // 1 2 3
                // 8 C 4
                // 7 6 5
                corr1 = 0;
                corr2 = 0;
                corr3 = 0;
                corr4 = 0;
                corr5 = 0;
                corr6 = 0;
                corr7 = 0;
                corr8 = 0;
                double other_var = 0.0;
                if (x != 0) {
                    other_var = var_img_double->at(y)->at(x-1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr8 = fabs(old_corr_values->at(y)->at(x)->at(8) / sqrt( curr_var * other_var ) );
                }
                if (y != 0) {
                    other_var = var_img_double->at(y-1)->at(x);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr2 = fabs(old_corr_values->at(y)->at(x)->at(2) / sqrt( curr_var * other_var ) );
                }
                if (x != w - 1) {
                    other_var = var_img_double->at(y)->at(x+1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr4 = fabs(old_corr_values->at(y)->at(x)->at(4) / sqrt(curr_var * other_var ) );
                }
                if (y != h - 1) {
                    other_var = var_img_double->at(y+1)->at(x);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr6 = fabs( old_corr_values->at(y)->at(x)->at(6) / sqrt( curr_var * other_var )) ;
                }
                if (x != 0 && y != 0) {
                    other_var = var_img_double->at(y-1)->at(x-1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr1 = fabs(old_corr_values->at(y)->at(x)->at(1) / sqrt( curr_var * other_var ) );
                }
                if (x != 0 && y != h - 1) {
                    other_var = var_img_double->at(y+1)->at(x-1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr3 = fabs( old_corr_values->at(y)->at(x)->at(3) / sqrt( curr_var * other_var ) ) ;
                }
                if (x != w-1 && y != 0) {
                    other_var = var_img_double->at(y-1)->at(x+1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr7 = fabs( old_corr_values->at(y)->at(x)->at(7) / sqrt( curr_var * other_var ) );
                }
                if (x != w-1 && y != h-1) {
                    other_var = var_img_double->at(y+1)->at(x+1);
                    if (other_var <= 0.001) {
                        other_var = 1.0;
                    }
                    corr5 = fabs( old_corr_values->at(y)->at(x)->at(5) / sqrt( curr_var * other_var ) );
                }
                QVector<double> * corrs = new QVector<double>();
                double d_n = static_cast<double>(datalist->size());
                corrs->append(corr1/ d_n);
                corrs->append(corr2/ d_n);
                corrs->append(corr3/ d_n);
                corrs->append(corr4/ d_n);
                corrs->append(corr5/ d_n);
                corrs->append(corr6/ d_n);
                corrs->append(corr7/ d_n);
                corrs->append(corr8/ d_n);
                if ( x == 0 && y ==2) {
                    cout << "eight neighbors: " << corrs->at(0) << " " << corrs->at(1) << " " << corrs->at(2)<< " " << corrs->at(3) << " " << corrs->at(4) << " " << corrs->at(5) << " " << corrs->at(6) << " " << corrs->at(7) << endl;
                }
                for( int z = 0; z < 8; z ++) {
                    if (corrs->at(z) > maximum_corr_rvalue) {
                        maximum_corr_rvalue = corrs->at(z);
                    }
                }

                maximum_corr_rvalue = maximum_corr_rvalue;

                if (maximum_corr_rvalue > 1.0 or maximum_corr_rvalue < -1.0){
                    cout << "shut! \n";
                    cout << "maximum corr rvalue: " << maximum_corr_rvalue << endl;
                    maximum_corr_rvalue = 1.0;
                }
                //                if (x==0 && y == 0) {
                //                    cout << "first pixel: " << maximum_corr_rvalue * 255.0 << endl;
                //                }
                corr_img->setPixel(x,y,static_cast<int>(maximum_corr_rvalue * 255.0));
            }
        }
        group_stat.corr_image = corr_img;
        buf.set_corr_img(corr_img);
        //       QImage * scaled_corr_image = new QImage(*corr_img);
        //       buf.scale_image2(corr_img, scaled_corr_image);
        corr_img->save(*(debug_image_name)+"corr_img.pgm", 0, 100);
    }

    // because in the sliding window version, everthing has been saved before this point, so we dont' execute any code
    //  below this point in the sliding window version.
    if (!DO_SLIDING_WINDOWING) {  // is not doing sliding windowing, means doing it online

        if (group_stat.mean_images->at(0) != NULL && group_stat.max_image != NULL && group_stat.min_image != NULL && group_stat.var_images->at(0) != NULL) {
            if (!ONLINE_CALC_MEAN_VAR ) { // not doing online, means do nothing
                cout << "neigther online nor windowing, do nothing..." << endl;
            }     
        }

        // run it here to compute slope and offset
        // extreme value test is no longer a plan
        //if (DO_EVT_LEARNING || DO_EVT_LEARNING2 || DO_KS_LEARNING || DO_SK_LEARNING || DO_CORRELATION || DO_CORRELATIONS) {
        if (DO_SK_LEARNING || DO_CORRELATION || DO_CORRELATIONS) {
            // new: do not show GUI, still call this method but chagne the things inside create_video_display();
            /// TODO: change the function name later.
            buf.create_video_display();
        }
        // get group stat again to retrieve the slope and offset...
        /// TODO: feel nervous here, this call to get_group_stat should go somewhere else
        buf.get_group_stat(&group_stat);
        /*
    if (DO_EVT_LEARNING || DO_EVT_LEARNING2 ) {
        printf("outer group stat slope: %f, offset: %f\n", group_stat.slopes->last(), group_stat.offsets->last());
    }
    */
        QVector<QImage*> * bin_diff_imgs = new QVector<QImage*>();
        if (DO_SK_LEARNING) {
            buf.two_levels_an_image_other(&group_stat, bin_diff_imgs, 1);
        }
        else if (DO_KS_LEARNING) {
            buf.two_levels_an_image_other(&group_stat, bin_diff_imgs, 2);
            buf.open_bin_image(bin_diff_imgs,3,DILATE_DIM, debug_image_name);
            buf.set_bin_diff_images_into_video_display(bin_diff_imgs);
        }
        else if (DO_CORRELATION) {
            // here: use the online version of the two_levels_an_image(); etc
            if (!DO_SLIDING_WINDOWING) {
                if (ONLINE_CALC_MEAN_VAR) {
                    buf.two_levels_an_image_other_online(bin_diff_imgs, 2, THRESHOLD_FOR_CORRELATION);
                    bin_diff_imgs->at(0)->save(*(debug_image_name)+"old_bin_diff_img.pgm", 0,100);
                    buf.open_bin_image(bin_diff_imgs,3,DILATE_DIM, debug_image_name);
                    bin_diff_imgs->at(0)->save(*(debug_image_name)+"bin_diff_img.pgm",0,100);
                    buf.set_bin_diff_images_into_video_display_online(bin_diff_imgs);
                }
                else {
                    buf.two_levels_an_image_other(&group_stat, bin_diff_imgs, 3);
                    buf.open_bin_image(bin_diff_imgs,3,DILATE_DIM, debug_image_name);  // mathematical morphology open operation
                    buf.set_bin_diff_images_into_video_display(bin_diff_imgs); // set message to video_analysis object
                }
            }
        }
        else if (DO_CORRELATIONS) {
            buf.two_levels_an_image_windowed_correlations(&group_stat, bin_diff_imgs);
            cout << "generated binary images based on correlaton map" << endl;

            buf.open_bin_image(bin_diff_imgs,3, DILATE_DIM, debug_image_name);  // mathematical morphology open operation
            buf.set_bin_diff_images_into_video_display(bin_diff_imgs); // set message to video_analysis object
        }

        if ( DO_LEARNING) {

            if (DO_KS_LEARNING || DO_SK_LEARNING) {

                cout << "EVT or Normality test is used in determing foreground/background" << endl;
                buf.output_pixels_intensity_over_time(true);

                buf.learning_from_evt_background(bin_diff_imgs); // set some useful variables in image buffer object
                buf.clean_multiple_background(bin_diff_imgs); // replace background with mean, this might be OK here for evt

                // save the (compressed) images
                char savename[1024];
                sprintf(savename, "Data_out/out_");
                strcat(savename, argv[1]);
                // set the save name

                ret = buf.save_buffer_images(savename); // save the compressed images
                //save the buffer in disk
                ret = buf.save_buffer(argv[1]);

                sprintf(savename, "Data_out_backempty/out_");
                strcat(savename, argv[1]);
                ret = buf.save_empty_back_buffer_images(savename);

                char saved_by_frame_filename[1024];
                saved_by_frame_filename[0] = '\0';
                strcat(saved_by_frame_filename, argv[1]);
                strcat(saved_by_frame_filename, ".frames");

                ret = buf.save_buffer_by_frame(saved_by_frame_filename);

                if(ret != 0)
                    cerr<<"Error: Saving buffer fail"<<endl;

                //compress using bzip2
                cout<<"Calling bzip2 to compress the video stored by pixel time series..."<<endl;
                string compress_command("bzip2 -k -9 ");
                compress_command += "\"";
                compress_command += argv[1];
                compress_command += "\"";
                ret = system(compress_command.c_str());

                //calculate the file size compression rate
                cout<<"Computing compression ratio..."<<endl;
                struct stat filestat_before;
                if( stat(argv[1], &filestat_before)<0 )
                    cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                else {
                    string compressed_filename(argv[1]);
                    compressed_filename += ".bz2";

                    //test if file exits, if so, delete it

                    struct stat filestat_after;
                    if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                        cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                    else {
                        cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                        cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                        cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                    }
                }

                //compress using bzip2
                cout<<"Calling bzip2 to compress the video stored by frames..."<<endl;
                string compress_command2("bzip2 -k -9 ");
                compress_command2 += "\"";
                compress_command2 += saved_by_frame_filename;
                compress_command2 += "\"";
                ret = system(compress_command2.c_str());

                //calculate the file size compression rate
                cout<<"Computing compression ratio..."<<endl;
                if( stat(saved_by_frame_filename, &filestat_before)<0 )
                    cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                else {
                    string compressed_filename(saved_by_frame_filename);
                    compressed_filename += ".bz2";

                    //test if file exits, if so, delete it

                    struct stat filestat_after;
                    if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                        cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                    else {
                        cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                        cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                        cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                    }
                }
            }
            else if (DO_CORRELATION) {
                cout << "Correlation is used in determing foreground/background" << endl;
                // set this to false for now.

                buf.output_pixels_intensity_over_time(false);

                //    buf.learning_from_corr_background(bin_diff_imgs); // set some useful variables in image buffer object
                // provide an option to do this online: write image files on the way
                // and  then disable the bzip codes.

                if (ONLINE_CALC_MEAN_VAR) {
                    // save images on the way
                    int number_of_frames = datalist->size();
                    cout << "number of frames: " << number_of_frames << endl;
                    buf.clean_multiple_background_online(bin_diff_imgs, number_of_frames, save_file_dir );
                }
                else {
                    buf.clean_multiple_background(bin_diff_imgs); // replace background with mean, this might be OK here for evt
                    // save the (compressed) images
                    char savename[1024];
                    sprintf(savename, "Data_out/out_");
                    strcat(savename, argv[1]);
                    // set the save name

                    ret = buf.save_buffer_images(savename); // save the compressed images
                    //save the buffer in disk
                    ret = buf.save_buffer(argv[1]);

                    // this does not happen in online version
                    sprintf(savename, "Data_out_backempty/out_");
                    strcat(savename, argv[1]);
                    ret = buf.save_empty_back_buffer_images(savename);

                    char saved_by_frame_filename[1024];
                    saved_by_frame_filename[0] = '\0';
                    strcat(saved_by_frame_filename, argv[1]);
                    strcat(saved_by_frame_filename, ".frames");

                    ret = buf.save_buffer_by_frame(saved_by_frame_filename);

                    if(ret != 0)
                        cerr<<"Error: Saving buffer fail"<<endl;

                    //compress using bzip2
                    cout<<"Calling bzip2 to compress the video stored by pixel time series..."<<endl;
                    string compress_command("bzip2 -k -9 ");
                    compress_command += "\"";
                    compress_command += argv[1];
                    compress_command += "\"";
                    ret = system(compress_command.c_str());

                    //calculate the file size compression rate
                    cout<<"Computing compression ratio..."<<endl;
                    struct stat filestat_before;
                    if( stat(argv[1], &filestat_before)<0 )
                        cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                    else {
                        string compressed_filename(argv[1]);
                        compressed_filename += ".bz2";

                        //test if file exits, if so, delete it

                        struct stat filestat_after;
                        if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                            cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                        else {
                            cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                            cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                            cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                        }
                    }

                    //compress using bzip2
                    cout<<"Calling bzip2 to compress the video stored by frames..."<<endl;
                    string compress_command2("bzip2 -k -9 ");
                    compress_command2 += "\"";
                    compress_command2 += saved_by_frame_filename;
                    compress_command2 += "\"";
                    ret = system(compress_command2.c_str());

                    //calculate the file size compression rate
                    cout<<"Computing compression ratio..."<<endl;
                    if( stat(saved_by_frame_filename, &filestat_before)<0 )
                        cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                    else {
                        string compressed_filename(saved_by_frame_filename);
                        compressed_filename += ".bz2";

                        //test if file exits, if so, delete it

                        struct stat filestat_after;
                        if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                            cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                        else {
                            cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                            cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                            cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                        }
                    }
                }
            }
            else if (DO_CORRELATIONS) {
                cout << "Windowed correlation is used in determing foreground/background" << endl;
                buf.output_pixels_intensity_over_time(true);
                cout << "before replacing images with mean" << endl;

                buf.clean_multiple_background(bin_diff_imgs); // replace background with mean, this might be OK here for evt

                // save the (compressed) images
                char savename[1024];
                sprintf(savename, "Data_out/out_");
                strcat(savename, argv[1]);
                // set the save name

                ret = buf.save_buffer_images(savename); // save the compressed images
                //save the buffer in disk
                ret = buf.save_buffer(argv[1]);

                sprintf(savename, "Data_out_backempty/out_");
                strcat(savename, argv[1]);
                ret = buf.save_empty_back_buffer_images(savename);

                char saved_by_frame_filename[1024];
                saved_by_frame_filename[0] = '\0';
                strcat(saved_by_frame_filename, argv[1]);
                strcat(saved_by_frame_filename, ".frames");

                ret = buf.save_buffer_by_frame(saved_by_frame_filename);

                if(ret != 0)
                    cerr<<"Error: Saving buffer fail"<<endl;

                //compress using bzip2
                cout<<"Calling bzip2 to compress the video stored by pixel time series..."<<endl;
                string compress_command("bzip2 -k -9");
                compress_command += "\"";
                compress_command += argv[1];
                compress_command += "\"";
                ret = system(compress_command.c_str());

                //calculate the file size compression rate
                cout<<"Computing compression ratio..."<<endl;
                struct stat filestat_before;
                if( stat(argv[1], &filestat_before)<0 )
                    cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                else {
                    string compressed_filename(argv[1]);
                    compressed_filename += ".bz2";

                    //test if file exists, if so, delete it

                    struct stat filestat_after;
                    if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                        cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                    else {
                        cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                        cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                        cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                    }
                }

                //compress using bzip2
                cout<<"Calling bzip2 to compress the video stored by frames..."<<endl;
                string compress_command2("bzip2 -k -9 ");
                compress_command2 += "\"";
                compress_command2 += saved_by_frame_filename;
                compress_command2 += "\"";
                ret = system(compress_command2.c_str());

                //calculate the file size compression rate
                cout<<"Computing compression ratio..."<<endl;
                if( stat(saved_by_frame_filename, &filestat_before)<0 )
                    cerr<<"Error: can't get statistics of original file \""<<argv[1]<<"\"."<<endl;
                else {
                    string compressed_filename(saved_by_frame_filename);
                    compressed_filename += ".bz2";

                    //test if file exits, if so, delete it

                    struct stat filestat_after;
                    if( stat(compressed_filename.c_str(), &filestat_after)<0 )
                        cerr<<"Error: can't get statistics of compressed file \""<<argv[1]<<"\"."<<endl;
                    else {
                        cout<<"Original data size: "<<filestat_before.st_size<<" bytes"<<endl;
                        cout<<"Compressed data size: "<<filestat_after.st_size<<" bytes"<<endl;
                        cout<<"compressed/original="<<(double)filestat_after.st_size/(double)filestat_before.st_size*100<<"%"<<endl;
                    }
                }
            }

#ifdef DISPLAY_SAMPLE_TIME_SERIES
            //show trace of a pixel, explicitly plotting
            int time_steps = buf.getSize();
            double pixel_time_series[time_steps];

            //show trace of a pixel using method
            cout<<"Plotting one pixel's time series..."<<endl;
            ret = buf.trace_one_pixel_with_plot(20, 20, pixel_time_series);
            if ( strcmp(argv[2],"-d") == 0 ) {
                ret = buf.trace_one_pixel_with_distribution_plot(20, 20, pixel_time_series,100);
            }

            cout<<"Plotting one pixel's time series..."<<endl;
            ret = buf.trace_one_pixel_with_plot(10, 30, pixel_time_series);
            if ( strcmp(argv[2],"-d") == 0 ) {
                ret = buf.trace_one_pixel_with_distribution_plot(10, 30, pixel_time_series,100);
            }

            cout<<"Plotting one pixel's time series..."<<endl;
            ret = buf.trace_one_pixel_with_plot(80, 80, pixel_time_series);
            if ( strcmp(argv[2],"-d") == 0 ) {
                ret = buf.trace_one_pixel_with_distribution_plot(80, 80, pixel_time_series,100);
            }
#endif
            cout<<"Done."<<endl;
            app.exit();
            return 0;
        }
    }
    else {
        cout<<"Done."<<endl;
        app.exit(2);
        return 2;
    }
}

