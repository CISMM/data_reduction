/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


/* Data reduction Library
 *
 */

/* Implementation plan
 * 1. Finish adding image to the buffer
 * 2. Finish functions that updates statistic information of the buffer
 * 	Updating is need-based
 * 3. Test adding image and statistics
 */


#ifndef DR_DATA_H
#define DR_DATA_H


//#include <QtGui/QImage>
//#include <QtGui/QLabel

//#include <QImage>

//tnt
#include <tnt.h>
#include <tnt_array1d.h>
#include <tnt_array2d.h>

//jama
#include <jama_lu.h>

//qt
#include <QtCore/QCoreApplication>
#include <QImage>
#include <QtOpenGL/QGLWidget>
#include <QMap>
#include <QTime>
#include <QtCore/qmath.h>

//qwt
//#include <qwt_series_data.h>

//gsl
#include <stdio.h>
#include <stdlib.h>

//cluster library
#include <cluster.h>

//c++
#include <vector>

//c
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <GLWidget.h>


#include "tnt.h"
#include <jama_qr.h>
#define MAX_PRESET_SIZE 36028797018963968ll

#include <cmath>
#include <QApplication>
#include <QSize>
#include <math.h>
#include <QFile>
#include <QCoreApplication>
#include <QScopedPointer>
#include <QTextStream>

// #include <random>

#include <qglobal.h>

//error code
#define BUFFER_FULL -1
#define RESOLUTION_NOT_MATCH -2

#define CANNOT_OPEN_FILE -3
#define EMPTY_BUFFER -4
#define MALLOC_RETURN_NULL -5

#define DO_LEARNING true

//#define DO_CLUSTER_LEARNING false
//#define DO_EVT_LEARNING false
//#define DO_EVT_LEARNING2 false
#define DO_KS_LEARNING false
#define DO_SK_LEARNING false
#define DO_CORRELATION true
#define DO_CORRELATIONS false // set this to be true will enable the correlation computatoin with windowing.

#define BACKGROUND_THRESHOLD 50
//  four methods for fitting a line to the mean over variance scatter plot
// if we don't use that line, we don't need them ( if all three  are false we don't use the line and
//   don't display the scatter plots )
#define USE_LEAST_SQUARES false
#define USE_RANSAC false
#define USE_MAXMIN_AND_MINMAX false
#define USE_THRESHOLD_BASED_TRIAL_NUMBER false

#define USE_HARD_CODED_THRESHOLD true

#define SHOW_GUI false

#define WRITE_XANDY_INTO_TXT false
// #define TRIAL_THRESHOLD 90
// ransac metric
extern int RANSAC_METRIC;
extern int NUMBER_OF_SEGS;
extern int HARD_CODED_THRESHOLD;
extern int WINDOW_SIZE;
using namespace std;



class dr_group_image_buffer;
/*
typedef enum {
    DR_FLOAT,
    DR_UINT32,
    DR_UINT8
    //DR_UINT16;
} dr_pixel_depth;
*/

/*
typedef struct {
    unsigned resolution_x;
    unsigned resolution_y;
    unsigned n_pixels;
    dr_pixel_depth pixel_depth;
    void* dataPtr;
} dr_image;
*/

typedef struct {
    int64_t gid;
    QImage *mean_image;
    QImage *max_image;
    QImage *min_image;
    QImage *diff_image;
    uint64_t preset_size;

    int resolution_x; //so it is consistent with Qt
    int resolution_y;
    bool is_full;

    QImage *var_image;
    double *var_double;
    QVector<QVector<double>*> * var_doubles;
    // here we may store a array of QImages representing the segments variance
    //  and in the code calculates the group stats, there is an if-else condition
    //  which for some case, calculate the variance over segments.
    QVector<QImage*> * var_images;
    QVector<QImage*> * mean_images;

    int number_of_seg;
    double slope;
    double offset;

    QVector<double> * slopes;
    QVector<double> * offsets;

    QImage * corr_image;
    QVector<QImage*> * window_corr_images;
} dr_group_stat;


//This is like a page table
//Keeping the statistics of which buffer (or called group) is in the main memory
//DB is like the execution engine, keeping which buffer to go to disk and which to load to main memory
//I am not entirely clear about this two yet
//I will incrementally test the buffer and build up this system
/*
class cache_table_entry {
    private:
        int pid;
        int score;
    public:
        cache_table_entry(int id, int init_score):pid(id), score(init_score) {}
        //bool
};

class dr_db {
    //accounting
    QString dir;
    QString dir_full_data;
    QString dir_reduced_data;
    int cache_size;
    int group_loaded;
    vector<int> all_pid_list;

    //keep stat of all groups
    vector<dr_group_stat> group_stats;

    //only keep cache of some groups' data
};
*/

/*
typedef struct {
    QString dir;
} dr_image_database;
*/

//int init_dr_image(dr_image *img)
//


class GLWidget : public QGLWidget {

    Q_OBJECT // must include this if you use Qt signals/slots


public:
    QImage data, gldata;
    //GLWidget(vector <QImage> *images, QWidget *parent = NULL);
    GLWidget(dr_group_image_buffer *img_buffer, QWidget *parent = NULL);
    //vector <QImage> *images;
    dr_group_image_buffer *img_buffer;
private:

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);

signals:
    void mouseMoved(QString str);
    void mousePressed();
public slots:
};


class main_helper {
public:
    void read_config_files(QMap<QString, QString> * key_val, QString * file_name);

private:

protected:


};

class video_analysis:public QWidget {

    Q_OBJECT // must include this if you use Qt signals/slots

private:
    GLWidget *video;

    dr_group_stat groupStat;

    QVector<QImage*> * bin_diff_imgs;
    unsigned h_bins;

public:
    bool show_time_series_mean;
    bool show_time_series_stdev;

    bool show_mean_images;
    bool show_bin_diff_images;
    bool show_one_mean_image;
    bool show_variance_image;
    bool show_correlation_image;
    bool show_frame;
    bool show_time_series;
    bool show_pixel_hist;
    bool show_least_square_fitted_line;
    bool show_least_square_fitted_lines;
    bool show_windowed_correlation_image;
    bool show_pixel_cumulative_hist;
    bool show_normalized_pixel_cumulative_hist;
    bool show_variance_vs_mean_intensity_plot;
    bool show_least_square_fitted_line_on_variance_vs_mean_intensity_plot;

    int current_x,current_y;

    int one_variance_index;

    //QLabel * corr_map_label;

    //video_analysis(vector <QImage> *images)
    video_analysis(dr_group_image_buffer *img_buffer) ;
    void plot_time_series();
    void set_show_time_series_mean(bool b);
    void set_show_time_series_stdev(bool b);

    void plot_histogram();
    void set_hist_bins(unsigned n);

    void plot_pixel_cumulative_distribution_histogram();
    void plot_normalized_pixel_cumulative_distribution_histogram();

    void plot_pixel_variance_vs_mean_intensity_scatter_plot(double & slope, double & offset, QVector<double> * ys_out, QVector<double> * xs_out, int plot_index);

    void plot_line_fit_plot(double slope, double offset, const QVector<double> * ys, const QVector<double> * xs, int plot_index);

    void ransac(QVector<double> * lineParams, QVector<double> * x, QVector<double> * y, QVector<double> * weights);
    void ransac_with_dynamic_threshold(QVector<double> * lineParams, QVector<double> * x, QVector<double> * y, QVector<double> * weights);

    void maxmin_and_minmax(double & slope, double & offset, TNT::Array2D<double> * x, TNT::Array2D<double> * y);

    int rand_int(int low, int high);
    int calculate_current_number_of_trails(double threshold, double inlierRatio, int m);
    void output_x_and_y_into_txt(QVector<double> * xs, QVector<double> * ys);

    void set_bin_diff_imgs(QVector<QImage*> * bin_diff_imgs);

    void plot_correlation_r_value_image();

    QRgb gray_scale_to_blackbody_rgb(int gray_scale_intensity);

    void  plot_windowed_corr_timeseries_plot();

public slots:
    void pixel_timeseries_update(QString str);
    void pixel_histogram_update(QString str);
    // cumulative update method
    void pixel_cumulative_histogram_update(QString str);
    void pixel_normalized_histogram_update(QString str);
    // corr image update method
    void windowed_corr_pixel_update(QString str);

};


typedef struct {
    double max_slope;
    double min_slope;
    double max_sd; //max standard deviation
    double max_of__max_min_diff;
    int n_backgrounds;
} background_features;


// function names in dr_group_image_buffer
class dr_group_image_buffer {

private:
    int64_t gid;		//buffer (group) ID

    vector <QImage> images;

    vector <QImage> empty_bgd_images;

    //statistics
    QImage *one_mean_image;
    QVector<QImage*> * mean_images;
    QImage *max_image;             //max image of the video
    QImage *min_image;             //min image of the video
    QImage *diff_image;             //min image of the video
    uint64_t preset_size;   //Upper limit of number of images in the buffer

    QVector<QImage*> * var_images;
    QVector<QImage*> * windowed_corr_images;
    double *var_double; // old stuff, try remove
    QVector<QVector<double>*> * var_doubles;
    // double **var_doubles;
    int number_of_seg;

    int resolution_x;       //row resolution
    int resolution_y;       //column resolution_x
    bool is_full;		//The buffer is full or not

    video_analysis *v;
    //void setupVideoData();

    background_features bg_features;

    double slope;
    double offset;

    QVector<double> * slopes;
    QVector<double> * offsets;

    QVector<QImage*> * bin_diff_images;

    QImage * corr_image;

    QVector<string> * datalist;

    QVector<QImage*> * buffer_for_sliding_windowing;

    // -- window stats --
    QVector<QVector<double> *> * curr_intensity_sum_in_window;
    QVector<QVector<double> *> * curr_mean_double_in_window;
    QVector<QVector<double> *> * curr_variance_double_in_window;
    QVector<QVector<double> *> * curr_sum_for_variance_double_in_window;
    QVector<QVector<QVector<double> *> *> * curr_sum_for_correlation_double_in_window;
    QVector<QVector<double> *> * curr_correlation_double_in_window;

    bool do_sliding_windowing;

    int sliding_window_buffer_size;

public:
    bool need_update_state;    //whenever a image is added or removed, it needs update

    dr_group_image_buffer(int preset_size, int id); //constructor, given size limit and ID
    int get_size();				//how many images are currently in the buffer
    int get_width();
    int get_height();
    int get_preset_size();			//limit of the buffer
    int add_image(QImage img);			//add an image to the buffer
    void get_group_stat(dr_group_stat* stat);	//get group statistic, update if needed
    int save_buffer(string filename);
    int save_buffer_by_frame(string filename);
    int save_buffer_images(string filename);

    int save_empty_back_buffer_images(string filename);

  //  int trace_one_pixel(int x, int y, double* data_out); //return array, the users need to free data_out after using
 //   int trace_one_pixel_with_plot(int x, int y, double* data_out); //return array, the users need to free data_out after using
  //  int trace_one_pixel_with_distribution_plot(int x, int y, double* data_out, int bins);
    int create_video_display();
    int get_frame(int k, QImage *img);
    //int playBuffer(int frames, QLabel* video);
    int two_levels_an_image(QImage img_in, QImage *img_out);

    int learning_from_kmeans_background(QImage *bin_diff_img);
    int clean_background(QImage *bin_diff_img);

    int clean_multiple_background(QVector<QImage*> * bin_diff_imgs);

    bool is_background(double slope, double sd, double max_min_diff);
    int plot_intensity_by_frame();
    int dilate_bin_image(QVector<QImage*> *imgs, int diameter);
    int erode_bin_image(QVector<QImage *> *imgs, int diameter);
    int open_bin_image(QVector<QImage *> *imgs, int erode_diameter, int dilate_diameter, QString * debug_image_name);

    void output_pixels_intensity_over_time(bool use_full_image);
    int plot_diff_curve();

    QVector<double> * learning_from_background_line(double slope, double offset, double mean);
    void divide_video_into_segments(int n);
    bool is_background_evt(QVector<double> * background_distribution_param, int max_pixel_value, int min_pixel_value, int number_of_samples);
    bool is_background_evt2(QVector<double> * background_distribution_param, int max_pixel_value, int number_of_samples);
    bool is_background_skew_kurto(QVector<double> * intensity_over_time, double skewness_threshold, double kurtosis_threshold);
    bool is_background_ks_test(QVector<double> * intensity_over_time, QVector<double> * standard_gaussian_cumulative_table, double threshold);
    double compute_skewness(QVector<double> * intensity_over_time, double * m2_out);
    double compute_kurtosis(QVector<double> * intensity_over_time, double * m2_out);

    double phi(double x1, double x2);

    int two_levels_an_image_evt(dr_group_stat * group_stat, QImage *img_out);
    int two_levels_an_image_evt2(dr_group_stat * group_stat, QImage *img_out);
    int two_levels_an_image_other(dr_group_stat * group_stat, QVector<QImage*> *imgs_out, int method);
    int two_levels_an_image_windowed_correlations(dr_group_stat * group_stat, QVector<QImage*> *imgs_out);

    int learning_from_evt_background(QVector<QImage*> * bin_diff_imgs);
    int learning_from_corr_background(QVector<QImage*> * bin_diff_imgs);
    void set_slope_and_offset(double theSlope, double theOffset);

    void compute_one_mean_image(QImage *one_mean_image);

    void set_bin_diff_images_into_video_display(QVector<QImage*> * bin_diff_imgs);


    void compute_corr_image(QImage * img);

    void compute_corr_image_with_windowing(QVector<QImage*> * out_corr_images, int window_size);

    void compute_corr_image_memory_friendly(QString * image_file_name_before_index, int number_of_frames, bool start_with_zero, QImage * image);

    void compute_corr_image_online(QImage * new_img, QVector<QVector<QVector<double> *> *> * old_corr_values, QVector<QVector<double> *> * mean_img_double);

    void compute_corr_image_with_windowing_memory_friendly();
    double compute_correlation_between_two_vectors(QVector<int> * vector1, QVector<int> * vector2);

    double compute_correlation_between_two_pixel_time_series(QVector<int> * vector1, QVector<int> * vector2, int x1, int y1, int x2, int y2);


    double anderson_darling_test(QVector<double> * x, double alpha);

    int clean_windowed_correlation_background(QVector<QImage*> * bin_diff_imgs);

    int online_mean_and_variance(int n, QVector<QVector<double> *> * mean_img_double,  QVector<QVector<double> *> *  M2_img, QImage * new_img);
    void set_mean_img(QImage * mean_img);
    void set_variance_img(QImage * variance_img);

    int two_levels_an_image_other_online(QVector<QImage*> * imgs_out, int method, int threshold);

    void set_bin_diff_images_into_video_display_online(QVector<QImage*> * bin_diff_imgs);

    int clean_multiple_background_online(QVector<QImage *> * bin_diff_imgs, int number_of_frames, QString * save_dir);

    void set_data_list(QVector<string> * datalist);

    void set_corr_img(QImage * corr_img);

    void scale_image(QVector<QVector<double> *> * img, QImage * new_img);

    void scale_image2(QImage * img, QImage * new_img);

    void save_compressed_image_in_queue(QImage * out_img, int erosion_size, int dilation_size, int index, QString * debugging_image_dir);

    void fill_window_buffer(QImage * in_img);

    void deque_buffer();

    void update_window_stats(bool has_deque, bool has_enque);

    void init_window_stats();

    void set_do_sliding_windowing(bool doit);

    void set_sliding_window_size(int size);

    int get_window_buffer_size();

    void init_window_stats2( QString * debugging_image_dir);

    void update_window_stats2(int index, QString * debugging_image_dir);

};

#endif
