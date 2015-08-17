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
#include <stdio.h>
#include <stdlib.h>


// ======================== stats functions ============================
double asl_stat_mean(double *data, int n) {
    int i;
    double s=0;
    for (i=0; i<n; i++) {
        s += data[i];
    }
    return s/(double)n;
}


double asl_stat_sd(double *data, double mean, int n) {
    int i;
    double s=0;
    for (i=0; i<n; i++) {
        s += pow( (data[i]-mean), 2);
    }
    return sqrt(s/(double)n);
}


void asl_stat_max(double* data, int n, double *max_value, int *max_index) {
    int i;
    *max_value = 1.0e-19;
    for (i=0; i<n; i++) {
        if (data[i] > *max_value) {
            *max_value = data[i];
            *max_index = i;
        }
    }
}


void asl_stat_max_min(double* data, int n, double *max_value, int *max_index, double *min_value, int *min_index) {
    int i;
    *max_value = data[0];
    *min_value = data[0];
    for (i=0; i<n; i++) {
        if (data[i] > *max_value) {
            *max_value = data[i];
            *max_index = i;
        }
        if (data[i] < *min_value) {
            *min_value = data[i];
            *min_index = i;
        }
    }
}


// ============================ matrix functions ===========================
template<class T>
TNT::Array2D<T> transpose(const TNT::Array2D<T> &M)
{
    TNT::Array2D<T> tran(M.dim2(), M.dim1() );
    for(int r=0; r<M.dim1(); ++r)
        for(int c=0; c<M.dim2(); ++c)
            tran[c][r] = M[r][c];
    return tran;
}


template<class T>
TNT::Array2D<T> invert(const TNT::Array2D<T> &M)
{
    assert(M.dim1() == M.dim2()); // square matrices only please

    // solve for inverse with LU decomposition
    JAMA::LU<T> lu(M);

    // create identity matrix
    TNT::Array2D<T> id(M.dim1(), M.dim2(), (T)0);
    for (int i = 0; i < M.dim1(); i++) id[i][i] = 1;

    // solves A * A_inv = Identity
    return lu.solve(id);
}


template <class T>
Array2D <T> aslpp_ml_linear_fit_normal_equation(Array2D <T> x, Array2D <T> y)
{
    return matmult ( matmult ( invert( matmult(transpose(x), x) ), transpose(x) ), y );
}


//only works for ncolumns == 1
//this function will return true if one cluster is merged with other cluster
bool merge_clusters(int *clusterID, int data_size, int n_clusters, double *means, double *sds) {
    int i,j;
    int id_of_larger, id_of_smaller;
    double scale=1.5;
    //double mean1, mean2, sd1, sd2;

    int clusterID_backup[data_size];
    //backup original assignment
    for (i=0; i<data_size; i++) {
        clusterID_backup[i] = clusterID[i];
    }

    //cout<<"hi 16"<<endl;
    for (i=0; i<n_clusters; i++) {
        for (j=i+1; j<n_clusters; j++) {
            id_of_larger = (means[i]>means[j])?i:j;
            id_of_smaller = (means[i]>means[j])?j:i;
            //overlap
            if (means[id_of_larger] - scale*sds[id_of_larger] <  means[id_of_smaller] + scale*sds[id_of_smaller]) {
                return true;
                //change last cluster to j
            }
        }
    }
    return false;
}


//  =========================== dr_group_image_buffer class methods ===========================
dr_group_image_buffer::dr_group_image_buffer(int preset_size, int id) {
    gid = id;
    images.clear();
    this->preset_size = preset_size;
    this->need_update_state = true;
    resolution_x = 0; //so it is consistent with Qt
    resolution_y = 0;
    is_full = false;
    one_mean_image = NULL;
    max_image = NULL;
    min_image = NULL;
    diff_image = NULL;

    var_double = NULL;
    var_doubles = NULL;

    slope = 0.0;
    offset = 0.0;

    this->var_images = new QVector<QImage*>();
    this->mean_images = new QVector<QImage*>();
    this->windowed_corr_images = new QVector<QImage*>();

    slopes = new QVector<double>();
    offsets = new QVector<double>();

    this->bin_diff_images = new QVector<QImage*>(); // new!
    v = NULL;

    corr_image = NULL;

    // buffer for sliding windowing
    do_sliding_windowing = false;
}


void dr_group_image_buffer::set_slope_and_offset(double theSlope, double theOffset) {
    slopes->push_back(theSlope);
    offsets->push_back(theOffset);
    return;
}


int dr_group_image_buffer::create_video_display() {
    if(images.size() == 0)
        return -1;


    v = new video_analysis(this);
    return 0;
}


void dr_group_image_buffer::set_bin_diff_images_into_video_display(QVector<QImage*> * binDiffImgs)
{
    if (v == NULL) {
        printf("error: video display not initialized.\n");
    }
    else {
        v->set_bin_diff_imgs(binDiffImgs);
        //  generate frames background replaced by 0, and the last image is background image with foreground replaced by 0.
        int i,h,w;
        for (i = 0; i < images.size(); i++) {
            QImage * emt_bgd_img = new QImage(images[0]);
            for (h= 0; h< emt_bgd_img->height(); h++) {
                for (w = 0; w < emt_bgd_img->width(); w++) {
                    // if is foreground, use that pixel;
                    // if is background, use 0;
                    if ( qGray((binDiffImgs->at(0))->pixel(w,h)) > 0 ) {  // foreground
                        emt_bgd_img->setPixel(w,h,qGray(((QImage)images[i]).pixel(w,h)));
                    }
                    else {
                        emt_bgd_img->setPixel(w,h,qRgb(0,0,0));
                    }
                }
            }
            this->empty_bgd_images.push_back(*emt_bgd_img);
        }
        QImage * bgd_img = new QImage(images[0]);
        for (h = 0; h < ((QImage)images[0]).height(); h++) {
            for (w = 0; w < ((QImage)images[0]).width(); w++) {
                if ( qGray((binDiffImgs->at(0))->pixel(w,h)) > 0 ) {  // foreground
                    bgd_img->setPixel(w,h,qRgb(0,0,0));
                }
                else {
                    bgd_img->setPixel(w,h,((QImage)images[0]).pixel(w,h));
                }
            }
        }
        this->empty_bgd_images.push_back(*bgd_img);
    }
    return;

}


// =========================== getters ===========================
int dr_group_image_buffer::get_size() {
    return images.size(); //images is the vector that holds images
}


int dr_group_image_buffer::get_preset_size() {
    return preset_size;
}


int dr_group_image_buffer::get_width() {
    return resolution_x;
}


int dr_group_image_buffer::get_height() {
    return resolution_y;
}


// ==========================  group stats ==========================
void dr_group_image_buffer::get_group_stat(dr_group_stat* stat) {
    if(need_update_state == false) {  // does not need to recalculate
        stat->gid = gid;

        stat->max_image = max_image;
        stat->min_image = min_image;
        stat->diff_image = diff_image;
        stat->preset_size = preset_size;
        stat->resolution_x = resolution_x;
        stat->resolution_y = resolution_y;
        stat->is_full = is_full;

        stat->var_doubles = var_doubles;
        stat->mean_images = mean_images;
        stat->var_images = var_images;

        stat->slopes = slopes;
        stat->offsets = offsets;
        stat->corr_image = corr_image;
        stat->window_corr_images = windowed_corr_images;
        return;
    }
    else { 	//need recalculation
        // need to cut video into pieces
        uint i, j, k,n, last_n, h, w;
        int x;

        // n is number of frames for each segment, assume it can be divided without remainder
        //   also assume each segment have the same length
        //   there are problems in calculating the variance over two time periods
        //   right now these two variances are the same..
        //   thinking about changeing the data structure implementation? (QVector)
        n = (images.size() / NUMBER_OF_SEGS); // number of segments minus 1
        last_n = images.size() - n * (NUMBER_OF_SEGS-1);  // handle the remainder, the last segment's number of
        //  frames will always be greater or equal to the previous ones
        h = images[0].height();
        w = images[0].width();
        var_doubles = new QVector<QVector<double> *>();
        //      new double*[NUMBER_OF_SEGS];

        // for each segment, compute the stats
        if (corr_image != NULL) {
            free(corr_image);
            corr_image = NULL;
        }
        corr_image = new QImage(images[0]);


        for (x = 0; x < NUMBER_OF_SEGS - 1; x++) {
            int64_t *sum = new int64_t[w*h];
            int64_t *max = new int64_t[w*h];
            int64_t *min = new int64_t[w*h];
            int64_t *diff = new int64_t[w*h];

            int diff_from_mean;
            int64_t max_var_sum;
            int64_t *mean = new int64_t[w*h];
            int64_t *var_sum = new int64_t[w*h];

            QVector<double> * curr_var_double = new QVector<double>();

            //initialize memory
            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {
                    if (SHOW_GUI) {
                        max[i*w + j] = 0ll;
                        min[i*w + j] = 255ll; //max 256
                        diff[i*w + j] = 0ll;  //max 256

                        var_sum[i*w + j] = 0ll;
                    }
                    sum[i*w + j] = 0ll;
                    mean[i*w + j] = 0ll;
                }
            }

            //FROM ME: For each pixel in the images, find and store
            //the sum of intensity values across all images, as well
            //as the most and least intense pixel among all images.
            //Pixel sums are used to calculate the average intensity (mean)
            //per pixel below.
            for (k=x*n; k<x*n+n; k++) {
                for (i=0; i<h; i++) {  //row is same as y
                    // parallellize the chosen segment of code
#pragma omp parallel for private(j)
                    for (j=0; j<w; j++) {  //column is same as x
                        //printf("image %d pixel(%d,%d) = %d\n",k, j, i, qGray(images[k].pixel(j, i)));
                        sum[i*w + j] += qGray(images[k].pixel(j, i));
                        if (SHOW_GUI) {
                            if ( max[i*w + j] <  qGray( images[k].pixel(j, i) ) )
                                max[i*w + j] =  qGray( images[k].pixel(j, i) );
                            if ( min[i*w + j] > qGray( images[k].pixel(j, i) ) )
                                min[i*w + j] = qGray( images[k].pixel(j, i) );
                        }
                    }
                }
            }

            // calculate the mean here in order to compute the variance
            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {
                    mean[i*w + j] = sum[i*w + j]/n;
                }
            }
            max_var_sum = 0;

            // calculate variance
            if (SHOW_GUI) {
                for (i=0; i<h; i++) {
                    for (j=0; j<w; j++) {
                        for (k=x*n; k<x*n+n; k++) {
                            diff_from_mean = qGray(images[k].pixel(j,i))- mean[i*w + j];
                            var_sum[i*w + j] += diff_from_mean * diff_from_mean;
                        }
                        if (var_sum[i*w+j] > max_var_sum) {
                            max_var_sum = var_sum[i*w+j];  // record the maximum variance sum over image
                        }
                    }
                }

                printf("======maximum variance sum: %ld======\n",max_var_sum);
                //FROM ME: For each pixel in the images, find and store the
                //difference of the max and min values found above.
                for (i=0; i<h; i++) {  //row is same as y
                    for (j=0; j<w; j++) {  //column is same as x
                        diff[i*w + j] = max[i*w + j]  - min[i*w + j];
                    }
                }
            }

            int current_meanpx, current_maxpx, current_minpx;
            double current_varpx;
            //allocate memory for mean_image, max_image, min_image

            if (max_image != NULL) {
                free(max_image);
                max_image = NULL;
            }
            if (min_image != NULL) {
                free(min_image);
                min_image = NULL;
            }
            if (diff_image != NULL) {
                free(diff_image);
                diff_image = NULL;
            }

            if (one_mean_image != NULL) {
                free(one_mean_image);
                one_mean_image = NULL;
            }



            max_image = new QImage(images[0]);
            min_image = new QImage(images[0]);
            diff_image = new QImage(images[0]);

            mean_images->push_back(new QImage(images[x]));
            var_images->push_back(new QImage(images[x]));

            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {
                    current_meanpx = mean[i*w + j];
                    mean_images->at(x)->setPixel(j,i, qRgb(current_meanpx,current_meanpx,current_meanpx));
                    if (SHOW_GUI) {
                        current_maxpx = max[i*w + j];
                        max_image->setPixel(j,i, qRgb(current_maxpx,current_maxpx,current_maxpx));
                        current_minpx = min[i*w + j];
                        min_image->setPixel(j,i,qRgb(current_minpx,current_minpx,current_minpx));
                        diff_image->setPixel(j,i,qRgb(diff[i*w+j], diff[i*w+j], diff[i*w+j]));

                        current_varpx = static_cast<double>(var_sum[i*w + j]) * 255.0 / static_cast<double>(max_var_sum);
                        //   var_double[i*w+j] = current_varpx;
                        curr_var_double->append(current_varpx);
                        var_images->at(x)->setPixel(j,i,qRgb(static_cast<int>(current_varpx),static_cast<int>(current_varpx),static_cast<int>(current_varpx)));
                    }
                }
            }
            if (SHOW_GUI) {
                var_doubles->append(curr_var_double);
            }
            // does this related to segmentation fault?
            delete [] sum;
            delete [] max;
            delete [] min;
            delete [] diff;
            delete [] var_sum;
            delete [] mean;
        }

        int64_t *sum = new int64_t[w*h];
        int64_t *max = new int64_t[w*h];
        int64_t *min = new int64_t[w*h];
        int64_t *diff = new int64_t[w*h];

        int diff_from_mean;
        int64_t max_var_sum;
        int64_t *mean = new int64_t[w*h];
        int64_t *var_sum = new int64_t[w*h];

        QVector<double> * curr_var_double = new QVector<double>();

        //initialize memory
        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                sum[i*w + j] = 0ll;
                mean[i*w + j] = 0ll;

                if (SHOW_GUI) {
                    max[i*w + j] = 0ll;
                    min[i*w + j] = 255ll; //max 256
                    diff[i*w + j] = 0ll;  //max 256

                    var_sum[i*w + j] = 0ll;
                }
            }
        }

        for (k=(NUMBER_OF_SEGS-1) * n; k<images.size(); k++) {
            for (i=0; i<h; i++) {  //row is same as y
                // parallellize the chosen segment of code
                for (j=0; j<w; j++) {  //column is same as x
                    sum[i*w + j] += qGray( images[k].pixel(j, i) );
                    if (SHOW_GUI) {
                        if (max[i*w + j] <  qGray( images[k].pixel(j, i) ))
                            max[i*w + j] =  qGray( images[k].pixel(j, i) );
                        if (min[i*w + j] > qGray( images[k].pixel(j, i) ))
                            min[i*w + j] = qGray( images[k].pixel(j, i) );
                    }
                }
            }
        }

        // calculate the mean here in order to compute the variance
        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                mean[i*w + j] = sum[i*w + j]/last_n;
            }
        }
        max_var_sum = 0;

        // calculate variance
        if (SHOW_GUI) {
            for (i=0; i<h; i++) {
                for (j=0; j<w; j++) {
                    for (k=(NUMBER_OF_SEGS-1) * n; k<images.size(); k++) {
                        diff_from_mean = qGray(images[k].pixel(j,i))- mean[i*w + j];
                        var_sum[i*w + j] += diff_from_mean * diff_from_mean;
                    }
                    if (var_sum[i*w+j] > max_var_sum) {
                        max_var_sum = var_sum[i*w+j];  // record the maximum variance sum over image
                    }
                }
            }
            printf("======maximum variance sum: %ld======\n",max_var_sum);
            //FROM ME: For each pixel in the images, find and store the
            //difference of the max and min values found above.
            for (i=0; i<h; i++) {  //row is same as y
                for (j=0; j<w; j++) {  //column is same as x
                    diff[i*w + j] = max[i*w + j]  - min[i*w + j];
                }
            }
        }

        int current_meanpx, current_maxpx, current_minpx;
        double current_varpx;
        //allocate memory for mean_image, max_image, min_image
        if (max_image != NULL) {
            free(max_image);
            max_image = NULL;
        }
        if (min_image != NULL) {
            free(min_image);
            min_image = NULL;
        }
        if (diff_image != NULL) {
            free(diff_image);
            diff_image = NULL;
        }

        max_image = new QImage(images[0]);
        min_image = new QImage(images[0]);
        diff_image = new QImage(images[0]);

        mean_images->push_back(new QImage(images[0]));
        var_images->push_back(new QImage(images[0]));

        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                current_meanpx = mean[i*w + j];
                mean_images->last()->setPixel(j,i, qRgb(current_meanpx,current_meanpx,current_meanpx));
                if (SHOW_GUI) {
                    current_maxpx = max[i*w + j];
                    max_image->setPixel(j,i, qRgb(current_maxpx,current_maxpx,current_maxpx));
                    current_minpx = min[i*w + j];
                    min_image->setPixel(j,i,qRgb(current_minpx,current_minpx,current_minpx));
                    diff_image->setPixel(j,i,qRgb(diff[i*w+j], diff[i*w+j], diff[i*w+j]));

                    // set one pixel of variance image
                    current_varpx = static_cast<double>(var_sum[i*w + j]) * 255.0 / static_cast<double>(max_var_sum);
                    //   var_double[i*w+j] = current_varpx;
                    curr_var_double->append(current_varpx);
                    // new! commented the line above
                    //      printf("curr variance px: %f\n", current_varpx);
                    var_images->last()->setPixel(j,i,qRgb(static_cast<int>(current_varpx),static_cast<int>(current_varpx),static_cast<int>(current_varpx)));
                }
            }
        }
        if (SHOW_GUI) {
            var_doubles->append(curr_var_double);
        }
        // does is related to segmentation fault?
        delete [] sum;
        delete [] max;
        delete [] min;
        delete [] diff;
        delete [] var_sum;
        delete [] mean;
    }
    // correlation is at last so we can use the pre-computed mean image(s).
    if (DO_CORRELATION) {
        compute_corr_image(corr_image);
    }
    if (DO_CORRELATIONS) {
        compute_corr_image_with_windowing(windowed_corr_images,20);
    }
    need_update_state = false;
    stat->gid = gid;
    stat->max_image = max_image;
    stat->min_image = min_image;
    stat->diff_image = diff_image;
    stat->preset_size = preset_size;
    stat->resolution_x = resolution_x;
    stat->resolution_y = resolution_y;
    stat->is_full = is_full;
    stat->mean_images = mean_images;
    stat->var_images = var_images;

    stat->var_doubles = var_doubles;

    stat->slopes = slopes;
    stat->offsets = offsets;
    stat->corr_image = corr_image;
    stat->window_corr_images = windowed_corr_images;
    //  printf("===========================\nvar images size: %d\n==============\n", var_images->size());
    //  printf("===========================\ngroup stat var images size: %d\n==============\n", stat->var_images->size());

    return ;
}


int dr_group_image_buffer::get_frame(int k, QImage *img) {
    unsigned kk = (unsigned)k;
    if(k<0 && kk>images.size()-1) {
        cerr<<"Error: dr_group_image_buffer::get_frame: index "<<k<<" out of range"<<endl;
        return -1;
    }
    *img = images[k];
    return 0;
}


// =========================== buffer operation methods ===========================
int dr_group_image_buffer::add_image(QImage img) {

    if(is_full)
        return BUFFER_FULL;

    if (images.size() == 0) {
        resolution_x = img.width();
        resolution_y = img.height();
    } else { //in a group, images must be of the same size
        if (resolution_x != img.width() || resolution_y != img.height() )
            return RESOLUTION_NOT_MATCH;
    }
    images.push_back(img);
    need_update_state = true;
    if( preset_size == (uint64_t)(images.size()) )
        is_full = true;
    return 0;
}


// =========================== I/O methods ===========================
//column store
int dr_group_image_buffer::save_buffer(string filename) {

    int n = images.size();
    if( n == 0 ) {
        cerr<<"Error: Trying to save an empty buffer."<<endl;
        return EMPTY_BUFFER;
    }

    FILE* f = fopen(filename.c_str(), "w");
    if (f == NULL) {
        cerr<<"Error: Can(not) open file \""<<filename<<"\"."<<endl;
        return CANNOT_OPEN_FILE;
    }

    char byte;
    int column, row, k;

    //FROM ME: Store header info (frame count, frame width, frame height)
    fprintf(f, "%d %d %d\n", static_cast<int>(images.size()), resolution_x, resolution_y);

    //for each coordinate, store all images' pixel value together, they are closed in values, good for compression
    //image index
    for (row=0; row<resolution_y; row++) {  //row corresponds to y
        for (column=0; column<resolution_x; column++) {  //column corresponds to x
            for (k=0; k<n; k++) {
                byte = (char)qGray( images[k].pixel(column, row));
                fputc(byte, f);
            }
        }
    }
    fclose(f);
    return 0;
}


int  dr_group_image_buffer::save_buffer_by_frame(string filename) {
    int n = images.size();
    if( n == 0 ) {
        cerr<<"Error: Trying to save an empty buffer."<<endl;
        return EMPTY_BUFFER;
    }

    FILE* f = fopen(filename.c_str(), "w");
    if (f == NULL) {
        cerr<<"Error: Can(not) open file \""<<filename<<"\"."<<endl; // should this "can" be replaced by "cannot" ??? -- Chong
        return CANNOT_OPEN_FILE;
    }

    char byte;
    int column, row, k;

    //FROM ME: Store header info (frame count, frame width, frame height)
    fprintf(f, "%d %d %d\n", static_cast<int>(images.size()), resolution_x, resolution_y);

    //for each coordinate, store all images' pixel value together, they are closed in values, good for compression
    //image index
    for (k=0; k<n; k++) {
        for (row=0; row<resolution_y; row++) {  //row - y
            for (column=0; column<resolution_x; column++) {  //column - x
                byte = (char)qGray( images[k].pixel(column, row));
                fputc(byte, f);
            }
        }
    }
    fclose(f);
    return 0;
}


int  dr_group_image_buffer::save_buffer_images(string filename) {
    int n = get_size();
    int i;

    char name[1024];
    char number[1024];
    for (i=0; i<n; i++) {
        sprintf(number,"%04d",i);
        strcpy(name, filename.c_str());
        strcat(name, number);
        strcat(name, ".pgm");
        if (images[i].save(name, 0, 100) == false) {
            printf("Error: Can't save %s!\n", name);
        }
    }

    return 0;
}


int dr_group_image_buffer::save_empty_back_buffer_images(string filename) {
    int n = get_size();
    int i;

    char name[1024];
    char number[1024];
    for (i=0; i<n+1; i++) {
        sprintf(number,"%04d",i);
        strcpy(name, filename.c_str());
        strcat(name, number);
        strcat(name, ".pgm");
        if (empty_bgd_images[i].save(name, 0, 100) == false) {
            printf("Error: Can't save %s!\n", name);
        }
    }

    return 0;
}


// ======================= background decision methods ==========================
int dr_group_image_buffer::two_levels_an_image_evt(dr_group_stat *group_stat, QImage *img_out) {
    printf("Performing evt background computation...\n");
    int x,y,z;
    //convert img_in to a column vector
    int h,w,k;
    k = images.size();
    QImage * sample_frame = &images[0];
    h = sample_frame->height();
    w = sample_frame->width();

    if (h != img_out->height() || w != img_out->width()) {
        cout << "There are problems with image size! may be the output image as parameter does not have the same size as the video frame!\n";
        return 0;
    }
    double slope = (*group_stat).slope;
    double offset = (*group_stat).offset;
    // for each pixel location, calculate its probability of being foreground
    int white, black; // count the number of pixels being foreground / background
    white = 0;
    black = 0;
    for (x = 0; x < h; x++) {
        for (y = 0; y < w; y++) {
            // calculate the parameter needed by is_foreground()
            double sum = 0.0;
            // calculate the maximum / minimum pixel value
            double max_pixel_value = static_cast<double>(qGray(sample_frame->pixel(y,x)));
            double min_pixel_value = 255.0;
            for (z = 0; z < k; z++) {
                QImage *curr_frame  = &images[z];
                double curr_pixel_double = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                sum = sum + curr_pixel_double;
                if (curr_pixel_double > max_pixel_value) {
                    max_pixel_value = curr_pixel_double;
                }
                if (curr_pixel_double < min_pixel_value) {
                    min_pixel_value = curr_pixel_double;
                }
            }
            double mean = sum / static_cast<double>(k); // the mean of that pixel position over time
            QVector<double> * background_distribution_param;
            // printf("max, min, k: %f, %f, %f\n", max_pixel_value, min_pixel_value, mean);
            // currently get the background distribution param from minimum values
            background_distribution_param = learning_from_background_line(slope, offset, mean); // WORK_ING

            //   printf("background distribution param: %f, %f\n", *background_distribution_param[0][0], *background_distribution_param[0][1]);
            if ( is_background_evt(background_distribution_param, max_pixel_value, min_pixel_value, k) ) {
                (*img_out).setPixel(y,x,qRgb(0,0,0));
                black += 1;
            }
            else {
                (*img_out).setPixel(y,x,qRgb(255,255,255));
                white += 1;
            }
        }
    }
    //remember to free space using delete later
    printf("number of foreground pixels: %d, number of background pixels: %d\n", white, black);
    return 0;
}


int dr_group_image_buffer::two_levels_an_image_evt2(dr_group_stat *group_stat, QImage *img_out) {
    // this function only check the maximum value rather than both max and min like in evt1.
    printf("Performing evt2 background computation...\n");
    int x,y,z;
    //convert img_in to a column vector
    int h,w,k;
    k = images.size();
    QImage * sample_frame = &images[0];
    h = sample_frame->height();
    w = sample_frame->width();

    if (h != img_out->height() || w != img_out->width()) {
        cout << "There are problems with image size! may be the output image as parameter does not have the same size as the video frame!\n";
        return 0;
    }
    double slope = (*group_stat).slope;
    double offset = (*group_stat).offset;
    // for each pixel location, calculate its probability of being foreground
    int white, black; // count the number of pixels being foreground / background
    white = 0;
    black = 0;
    for (x = 0; x < h; x++) {
        for (y = 0; y < w; y++) {
            // calculate the parameter needed by is_foreground()
            double sum = 0.0;
            // calculate the maximum / minimum pixel value
            double max_pixel_value = static_cast<double>(qGray(sample_frame->pixel(y,x)));
            for (z = 0; z < k; z++) {
                QImage *curr_frame  = &images[z];
                double curr_pixel_double = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                sum = sum + curr_pixel_double;
                if (curr_pixel_double > max_pixel_value) {
                    max_pixel_value = curr_pixel_double;
                }
            }
            double mean = sum / static_cast<double>(k); // the mean of that pixel position over time
            QVector<double> * background_distribution_param;

            // get the background mean and variance from the slope and offset of the line
            background_distribution_param = learning_from_background_line(slope, offset, mean);

            //   printf("background distribution param: %f, %f\n", *background_distribution_param[0][0], *background_distribution_param[0][1]);
            if ( is_background_evt2(background_distribution_param, max_pixel_value, k) ) {
                (*img_out).setPixel(y,x,qRgb(0,0,0));
                black += 1;
            }
            else {
                (*img_out).setPixel(y,x,qRgb(255,255,255));
                white += 1;
            }
        }
    }
    //remember to free space using delete later
    printf("number of foreground pixels: %d, number of background pixels: %d\n", white, black);
    return 0;
}


int dr_group_image_buffer::two_levels_an_image_other(dr_group_stat * group_stat, QVector<QImage*> *imgs_out, int method ) {
    // method 1: skewness and kurtosis
    // method 2: KS test
    // method 3: correlation

    printf("Performing m = %d background computation... m = 1 for skewness and kurtosis, 2 for KS test, 3 for correlation\n", method);
    int x,y,z;
    //convert img_in to a column vector
    int h,w,k;
    k = images.size();
    QImage * sample_frame = &images[0];
    QImage * img_out = new QImage(images[0]);
    h = sample_frame->height();
    w = sample_frame->width();

    if (h != img_out->height() || w != img_out->width()) {
        cout << "There are problems with image size! may be the output image as parameter does not have the same size as the video frame!\n";
        return 0;
    }
    int white = 0, black = 0;
    int n = k / NUMBER_OF_SEGS;
    double slope = 0.0;
    double offset = 0.0;
    // for correlation, temporal information has already been computed somewhere else. so we just use a single correlation map and set the threshold on it
    // we precompute the threshold, for example the 10% smallest protion of the correlation map value.
    QImage * corr_map = (*group_stat).corr_image;
    int corr_threshold = 0.0;
    if (method == 3) {
        QVector<int> * corr_values = new QVector<int>();
        for (int i = 0; i < corr_map->height(); i++) {
            for (int j = 0; j < corr_map->width(); j++) {
                int pixel = qGray(corr_map->pixel(j,i));
                if (corr_values->size() == 0) {
                    corr_values->append(pixel);
                }
                else {
                    /// TODO: here uses insertion sort... change it later.
                    for (int k = 0; k < corr_values->size(); k++) {
                        if (pixel > corr_values->at(k)) {
                            corr_values->insert(k,pixel);
                            break;
                        }
                    }
                }
            }
        }
        // take the 25%th correlation value as the threshold
        int threshold_index = corr_values->size() / 4;
        corr_threshold = corr_values->at(static_cast<int>(threshold_index));
    }
    printf("corr threshold: %d \n", corr_threshold);

    for (int i = 0; i < NUMBER_OF_SEGS - 1; i++) {
        if (method == 1 || method == 2) { // normality test need slope and offset
            slope = (*group_stat).slopes->at(i);
            offset = (*group_stat).offsets->at(i);
        }
        white = 0;
        black = 0;
        img_out = new QImage(images[0]);
        // find the threshold based on alpha value

        for (x = 0; x < h; x++) {
            for (y = 0; y < w; y++) {
                // calculate the parameter needed by is_foreground()
                double sum = 0.0;
                // calculate the maximum / minimum pixel value
                double max_pixel_value = static_cast<double>(qGray(sample_frame->pixel(y,x)));
                for (z = n * i ; z < n * i + n; z++) {
                    QImage *curr_frame  = &images[z];
                    double curr_pixel_double = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                    sum = sum + curr_pixel_double;
                    if (curr_pixel_double > max_pixel_value) {
                        max_pixel_value = curr_pixel_double;
                    }
                }
                double mean = sum / static_cast<double>(n); // the mean of that pixel position over time
                QVector<double> * background_distribution_param;

                // get the background mean and variance from the slope and offset of the line
                background_distribution_param = learning_from_background_line(slope, offset, mean);

                // we need to provide a QVector of the normalized value of a pixel intensities ( substract the mean and divide by the standard derivation )
                // compute the intensity over time there

                QVector<double> *intensity_over_time = new QVector<double>(n,0.0);

                double sum_for_standard_deviation = 0.0;

                for (z = n * i ; z < n * i + n; z++) {
                    QImage *curr_frame  = &images[z];
                    double curr_diff = static_cast<double>(qGray(curr_frame->pixel(y,x))) - mean;
                    sum_for_standard_deviation += curr_diff * curr_diff;
                }

                double st_dev = sqrt(sum_for_standard_deviation / n);

                //  double st_dev = sqrt(background_distribution_param->at(1));
                int index = 0;
                for (z = n * i ; z < n * i + n; z++) {
                    QImage *curr_frame  = &images[z];
                    double normalized_pixel = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                    normalized_pixel = normalized_pixel - mean;
                    normalized_pixel = normalized_pixel / st_dev;
                    (*intensity_over_time)[index] = normalized_pixel;
                    index++;
                }

                if (method == 1) {
                    //   printf("background distribution param: %f, %f\n", *background_distribution_param[0][0], *background_distribution_param[0][1]);
                    if ( is_background_skew_kurto(intensity_over_time, 1, 3) ) {
                        (*img_out).setPixel(y,x,qRgb(0,0,0));
                        black += 1;
                    }
                    else {
                        (*img_out).setPixel(y,x,qRgb(255,255,255));
                        white += 1;
                    }
                }
                else if (method == 3) {
                    // since there are only one  map, we simply generate the same binary image for each segment, a variation (multiple segments) is
                    //   not corporated with the code here. we use a window size. It is implemented in a separate function.
                    if (qGray(corr_map->pixel(y,x)) < corr_threshold) {
                        (*img_out).setPixel(y,x, qRgb(0,0,0));
                        black += 1;
                    }
                    else {
                        (*img_out).setPixel(y,x,qRgb(255,255,255));
                        white+= 1;
                    }
                }
            }
        }
        imgs_out->push_back(img_out); // push the bin diff image
        this->bin_diff_images->push_back(img_out);
    }


    img_out = new QImage(images[0]);
    // handle the last segment
    int last_n = k - n * (NUMBER_OF_SEGS - 1);
    for (x = 0; x < h; x++) {
        for (y = 0; y < w; y++) {
            double sum = 0.0;
            // calculate the maximum / minimum pixel value
            double max_pixel_value = static_cast<double>(qGray(sample_frame->pixel(y,x)));
            for (z = n * (NUMBER_OF_SEGS - 1); z < images.size(); z++) {
                QImage *curr_frame  = &images[z];
                double curr_pixel_double = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                sum = sum + curr_pixel_double;
                if (curr_pixel_double > max_pixel_value) {
                    max_pixel_value = curr_pixel_double;
                }
            }
            double mean = sum / static_cast<double>(last_n); // the mean of that pixel position over time
            QVector<double> * background_distribution_param;

            // get the background mean and variance from the slope and offset of the line
            background_distribution_param = learning_from_background_line(slope, offset, mean);
            QVector<double> *intensity_over_time = new QVector<double>(last_n,0.0);


            //     double st_dev = sqrt(background_distribution_param->at(1));
            double sum_for_standard_deviation = 0.0;

            for (z = n * (NUMBER_OF_SEGS - 1); z < images.size(); z++) {
                QImage *curr_frame  = &images[z];
                double curr_diff = static_cast<double>(qGray(curr_frame->pixel(y,x))) - mean;
                sum_for_standard_deviation += curr_diff * curr_diff;
            }

            double st_dev = sqrt(sum_for_standard_deviation / last_n);

            int index = 0;
            for (z = n * (NUMBER_OF_SEGS - 1); z < images.size(); z++) {
                QImage *curr_frame  = &images[z];
                double normalized_pixel = static_cast<double>(qGray(curr_frame->pixel(y,x)));
                normalized_pixel = normalized_pixel - mean;
                normalized_pixel = normalized_pixel / st_dev;
                (*intensity_over_time)[index] = normalized_pixel;
                index++;
            }
            if (method == 1) {
                //   printf("background distribution param: %f, %f\n", *background_distribution_param[0][0], *background_distribution_param[0][1]);
                if ( is_background_skew_kurto(intensity_over_time, 1, 3) ) {
                    (*img_out).setPixel(y,x,qRgb(0,0,0));
                    black += 1;
                }
                else {
                    (*img_out).setPixel(y,x,qRgb(255,255,255));
                    white += 1;
                }
            }
            else if (method == 3) {
                if (qGray(corr_map->pixel(y,x)) < corr_threshold) {
                    (*img_out).setPixel(y,x, qRgb(0,0,0));
                    black += 1;
                }
                else {
                    //   cout << "found one foreground!" << endl;
                    (*img_out).setPixel(y,x,qRgb(255,255,255));
                    white+= 1;
                }
            }
        }
    }

    imgs_out->push_back(img_out); // push the bin diff image
    this->bin_diff_images->push_back(img_out);
    //remember to free space using delete later
    printf("number of foreground pixels: %d, number of background pixels: %d\n", white, black);
    return 0;
}


int dr_group_image_buffer::two_levels_an_image_windowed_correlations(dr_group_stat * group_stat, QVector<QImage*> *imgs_out) {
 //   printf("Performing m = %d background computation using windowed correlation maps");
    int x,y;
    //convert img_in to a column vector
    int h,w;
    QImage * sample_frame = &images[0];
    QImage * img_out = new QImage(images[0]);
    h = sample_frame->height();
    w = sample_frame->width();

    if (h != img_out->height() || w != img_out->width()) {
        cout << "There are problems with image size! may be the output image as parameter does not have the same size as the video frame!\n";
        return 0;
    }
    int white = 0, black = 0;
    // for correlation, temporal information has already been computed somewhere else. so we just use a single correlation map and set the threshold on it
    // we precompute the threshold, for example the 10% smallest portion of the correlation map value.
    QVector<QImage*> * corr_maps = (*group_stat).window_corr_images;
    int n = corr_maps->size();
    int corr_threshold = 37;
    // currently set the magic number threshold to be 37

    for (int i = 0; i < n; i++) {
        white = 0;
        black = 0;
        img_out = new QImage(images[0]);
        // find the threshold based on alpha value

        for (x = 0; x < h; x++) {
            for (y = 0; y < w; y++) {
                // since there are only one correlation map, we simply generate the same binary image for each segment, a variation (multiple segments) is
                //   not corporated with the code here. we use a window size. It is implemented in a separate function.
                if (qGray(corr_maps->at(i)->pixel(y,x)) < corr_threshold) {
                    (*img_out).setPixel(y,x, qRgb(0,0,0));
                    black += 1;
                }
                else {
                    (*img_out).setPixel(y,x,qRgb(255,255,255));
                    white+= 1;
                }
            }
        }
        imgs_out->push_back(img_out); // push the bin diff image
        this->bin_diff_images->push_back(img_out);
    }
    //remember to free space using delete later
    printf("number of foreground pixels: %d, number of background pixels: %d\n", white, black);
    return 0;
}


void plotResults(double* xData, double* yData, int dataSize, char style) {  
    FILE *gnuplotPipe,*tempDataFile;
    const char *tempDataFileName;
    double x,y;
    int i;
    tempDataFileName = "tempData";
    //gnuplotPipe = popen("c:\\gnuplot\\bin\\pgnuplot -persist","w");
    gnuplotPipe = popen("gnuplot -persist","w");
    if (gnuplotPipe) {
        if(style == '-')
            fprintf(gnuplotPipe,"plot \"%s\" with lines\n",tempDataFileName);
        else if  (style == 'p') {
            fprintf(gnuplotPipe,"plot \"%s\" with linespoints\n",tempDataFileName);
        }
        else if  (style == 'x')
            fprintf(gnuplotPipe,"plot \"%s\" with points\n",tempDataFileName);
        else if  (style == 'b')
            fprintf(gnuplotPipe,"plot \"%s\" with boxes\n",tempDataFileName);

        fflush(gnuplotPipe);
        tempDataFile = fopen(tempDataFileName,"w");
        for (i=0; i <= dataSize; i++) {
            x = xData[i];
            y = yData[i];
            fprintf(tempDataFile,"%lf %lf\n",x,y);
        }
        fclose(tempDataFile);
        //printf("press enter to continue...");
        //getchar();
        remove(tempDataFileName);
        fprintf(gnuplotPipe,"exit \n");
    } else {
        printf("gnuplot not found...");
    }
}   


// ======================= background methods ==========================
// get 2 figures, mean and standard deviation
int dr_group_image_buffer::learning_from_kmeans_background(QImage *bin_diff_img) {
    int n = images.size();
    double time_series[n];

    //count how many "backgrounds" in diff_image
    int i=0,j=0,k=0;
    int w,h;
    w = images[0].width();
    h = images[0].height();
    int n_backgrounds = 0;
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            //printf("%d \n",qGray(diff_image->pixel(j,i)));
            if ( qGray(bin_diff_img->pixel(j,i)) == 0 ) {
                n_backgrounds ++;
            }
        }
    }
    printf("n_backgrounds=%d\n",n_backgrounds);
    printf("n_signal=%d\n",w*h-n_backgrounds);
    printf("background percentage=%f %%\n", (double)n_backgrounds/(w*h)*100);


    //double mean[n_backgrounds];
    double mean_a[n_backgrounds];
    //double sd[n_backgrounds];
    double sd_a[n_backgrounds];
    double max_of__max_min_diff = 0;

    //FROM ME: So among pixels that are considered background, we
    //are getting the average intensity value, standard deviation for
    //that pixel among all the images and we are also finding the
    //largest max-min difference seen among all background pixels across
    //all the images.
    k = 0;
    // looks like this for loop is a little bit slow, but it may be impossible to speed up since it calls trace_one_pixel() -- Chong
    double mean;
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            //printf("%d \n",qGray(diff_image->pixel(j,i)));
            if ( qGray(bin_diff_img->pixel(j,i)) == 0 ) {
          //      trace_one_pixel(j,i,time_series); // trace that pixel in time series. fill in time_series array
                mean_a[k] = asl_stat_mean(time_series,n); // compute mean from that array of that traced pixel
                mean = mean_a[k];
                //mean[k] = gsl_stats_mean(time_series,1,n);
                sd_a[k++] = asl_stat_sd(time_series, mean, n); // compute std from that array of that traced pixel
                //sd[k++] = gsl_stats_sd(time_series,1,n);
                //max min difference
                double temp_max, temp_min;
                int temp_max_index, temp_min_index;
                asl_stat_max_min(time_series, n, &temp_max, &temp_max_index, &temp_min, &temp_min_index); // compute max/min from that time series of that pixel
                if(k==0)
                    max_of__max_min_diff = temp_max - temp_min;
                else {
                    double temp_max_min_diff = temp_max-temp_min;
                    if (temp_max_min_diff > max_of__max_min_diff)
                        max_of__max_min_diff = temp_max_min_diff;
                }
            }
        }
    }
    printf("k=%d\n",k);

    //for (i=0; i<n_backgrounds; i++)
    // looks like the reason for 3 is that k is too large -- Chong
    for (i=0; i<3; i++) {
        printf("mean_a=%lf, std_a=%lf\n", mean_a[i], sd_a[i]);
    }

    /*
    FILE* f = fopen("bg_pixel_trace","r");
    if (f == NULL) {
        fprintf(stderr,"Can't open file to save background pixel trace\n");
    }
    bool done=false;
    */

    //fitting straight lines to each pixel's time series
    //frame index
    Array2D< double > X(n, 2);
    for (int i=0; i<n; i++)
    {
        X[i][0] = 1;
        X[i][1] = i;
    }

    Array2D< double > y(n, 1); // <- looks like then this is just 1D array -- Chong

    k=0;
    //double theta0[n_backgrounds];
    double theta1[n_backgrounds];
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            //printf(" (%d,%d) \n", j,i);
            if ( qGray(bin_diff_img->pixel(j,i)) == 0 ) {



                Array2D< double > y(n, 1);
                for (int l=0; l<n; l++) {
                    y[l][0] = time_series[l];
                }
                Array2D <double> Theta = aslpp_ml_linear_fit_normal_equation(X, y);
                //printf("k=%d\n",k);
                //theta0[k] = Theta[0][0];
                theta1[k] = Theta[1][0];
                //printf("theta0[%d]=%lf\n",k,theta0[k]);
                //printf("theta1[%d]=%lf\n",k,theta1[k]);
                k++;
            }
        }
    }
    cout<<"k="<<k<<endl;

    //plot sda vs. slop
    //plotResults(theta1, sd_a, n_backgrounds, 'x');

    //find max standard deviation
    double max_sd;
    int max_sd_i;
    asl_stat_max(sd_a, n, &max_sd, &max_sd_i);
    printf("max time series standard deviation: %lf\n", max_sd);

    double max_slope;
    double min_slope;
    int max_slope_i;
    int min_slope_i;
    asl_stat_max_min(theta1, n, &max_slope, &max_slope_i, &min_slope, &min_slope_i);
    printf("max time series background slope: %lf\n", max_slope);
    printf("min time series background slope: %lf\n", min_slope);

    //FROM ME: Setting these in object attribute
    bg_features.max_slope = max_slope;
    bg_features.min_slope = min_slope;
    bg_features.max_sd = max_sd;
    bg_features.n_backgrounds = n_backgrounds;
    bg_features.max_of__max_min_diff = max_of__max_min_diff;
    return 0;
}


int dr_group_image_buffer::clean_background(QImage *bin_diff_img) {
    int n = images.size();
    int i=0,j=0,k=0,l=0,m=0;
    int w,h;
    w = images[0].width();
    h = images[0].height();
    double time_series[n];

    //for each background pixel, replace with its mean
    double mean;
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            //printf("%d \n",qGray(diff_image->pixel(j,i)));
            if ( qGray(bin_diff_img->pixel(j,i)) == 0 ) {
                /*
                trace_one_pixel(j,i,time_series);
                mean = asl_stat_mean(time_series,n);
                */
                mean = qGray(mean_images->at(0)->pixel(j,i));

                for (k=0; k<n; k++) {
                    images[k].setPixel(j,i,qRgb(mean,mean,mean));
                }
            }
        }
    }

    //k-means with 3 levels on each pixel that contains non-background
    int col = 1;
    int row = n;
    //const int nclusters = 3;
    const int transpose = 0;
    const char dist = 'e';
    const char method = 'a'; //arithmetic mean, 'm' for median
    int npass = 3;
    int ifound = 0;
    double error;

    int **mask;
    //allocate rows
    mask = new int*[row];
    //allocate columns
    //cout<<"hi -5"<<endl;
    for (i=0; i<row; i++)
        mask[i] = new int;

    //cout<<"hi -4"<<endl;
    for (i=0; i<row; i++) {
        for (j=0; j<w; j++) {
            mask[i][0] = 1;
        }
    }

    //cout<<"hi -3"<<endl;
    double **data;
    //allocate rows
    data = new double*[row];
    //allocate columns
    for (i=0; i<row; i++)
        data[i] = new double;

    //cout<<"hi -2"<<endl;
    double* weight = (double*)malloc(col*sizeof(double));
    for (i = 0; i < col; i++)
        weight[i] = 1.0;

    //cout<<"hi -1"<<endl;
    int* clusterid = (int*)malloc(row*sizeof(int));

    int max_nclusters = 3;
    int final_n_clusters;
    double final_means[max_nclusters];
    double final_sds[max_nclusters];
    int elements[max_nclusters];
    int X_size;

    //FROM ME: Visits every pixel in image looking for level 1 pixels.
    //Clusters the trace of the pixel into 3 groups and checks to see
    //if the line fit for any of the clusters falls within the bounds
    //of the background that was determined earlier.
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            //cout<<"hi 0"<<endl;
            //printf("%d \n",qGray(diff_image->pixel(j,i)));
            if ( qGray(bin_diff_img->pixel(j,i)) != 0 ) {
                //trace_one_pixel(j,i,time_series);
                //cout<<"hi 1"<<endl;
                for (k=0; k<n; k++) {
                    data[k][0] = time_series[k] ;
                }
#ifdef DB_DISPLAY
                printf("Pixel(%d,%d) : %d levels\n",j,i,final_n_clusters);
#endif
                for (k=0; k<final_n_clusters; k++) {
#ifdef DB_DISPLAY
                    printf("mean[%d] = %lf  sd[%d] = %lf element[%d]=%d,",k, final_means[k], k, final_sds[k], k, elements[k]);
#endif
                    X_size=0;

                    Array2D <double> X(elements[k],2);
                    Array2D <double> y(elements[k],1);
                    double cluster_time_series[elements[k]];

                    for(l=0;l<n;l++) { //check the whole time series
                        if(clusterid[l] == k) {
                            X_size++;
                        }
                    }
                    if (X_size != elements[k]) {
                        fprintf(stderr,"Error: number of elements doesn't match!!!!\n");
                        exit(-1);
                    }

                    m = 0;
                    for(l=0;l<n;l++) { //check the whole time series
                        if(clusterid[l] == k) {
                            X[m][0] = 1;
                            X[m][1] = m;
                            y[m][0] = time_series[l];
                            cluster_time_series[m] = time_series[l];
                            m++;
                        }
                    }
                    if (m != elements[k]) {
                        fprintf(stderr,"Error: number of elements doesn't match!\n");
                        exit(4);
                    }
                    Array2D <double> theta =  aslpp_ml_linear_fit_normal_equation(X, y);

                    if(theta.dim1() == 2 && theta.dim2() == 1) {
#ifdef DB_DISPLAY
                        printf("slope=%lf\n", theta[1][0]);
#endif

                        double cluster_time_series_max, cluster_time_series_min;
                        int cluster_time_series_max_index, cluster_time_series_min_index;
                        asl_stat_max_min(cluster_time_series, n, &cluster_time_series_max, &cluster_time_series_max_index, &cluster_time_series_min, &cluster_time_series_min_index);


                        if (is_background(theta[1][0], final_sds[k], (cluster_time_series_max-cluster_time_series_min)))
                        {
                            printf("Found partial background\n");
                            //replace background by its mean value
                            for (int p=0; p<n; p++) {
                                if (clusterid[p] == k) {
                                    images[p].setPixel(j,i, qRgb(final_means[k], final_means[k], final_means[k]));
                                }
                            }

                        }
                    }
                }
            }
        }
    }
    //merge similar levels
    //find which level matches features of "background" and clean it.
    return 0;
}


int dr_group_image_buffer::clean_multiple_background(QVector<QImage*> * bin_diff_imgs){
    // change this to support multiple segments
    int n = images.size() / NUMBER_OF_SEGS;
    //  int last_n = images.size() - n * (NUMBER_OF_SEGS - 1);

    int i=0,j=0,k=0;
    int w,h;
    w = images[0].width();
    h = images[0].height();
    // double time_series[n];
    double mean;
    QImage * mean_img = this->one_mean_image;
    QImage * curr_bin_diff_img;
    //for each background pixel, replace with its mean
    // do it NUMBER_OF_SEGS - 1 times
    for (int x = 0; x < NUMBER_OF_SEGS - 1; x++) {
        curr_bin_diff_img = bin_diff_imgs->at(x);
        for (i = 0; i < h; i++) {
            for (j = 0; j < w; j++) {
                if( qGray(curr_bin_diff_img->pixel(j,i)) == 0 ) {
                    mean = qGray(mean_img->pixel(j,i));
                    for (k = x * n; k < x * n + n; k++) {
                        images[k].setPixel(j,i,qRgb(mean,mean,mean));
                    }
                }
            }
        }
    }
    // handle last segment
    curr_bin_diff_img = bin_diff_imgs->last();
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            if (qGray(curr_bin_diff_img->pixel(j,i)) == 0) {
                mean = qGray(mean_img->pixel(j,i));
                for (k = n * (NUMBER_OF_SEGS-1); k < images.size(); k++) {
                    images[k].setPixel(j,i,qRgb(mean,mean,mean));
                }
            }
        }
    }
    return 0;
}


bool dr_group_image_buffer::is_background(double slope, double sd, double max_min_diff) {
    if (slope < 0 && slope < bg_features.min_slope) {
#ifdef DB_DISPLAY
        printf("slope=%lf, background min slope=%lf, too small, not bg\n", slope, bg_features.min_slope);
#endif
        return false;
    }
    if (slope > bg_features.max_slope) {
#ifdef DB_DISPLAY
        printf("slope=%lf, background max slope=%lf, too big, not bg\n", slope, bg_features.max_slope);
#endif
        return false;
    }
    if (sd > bg_features.max_sd) {
#ifdef DB_DISPLAY
        printf("sd=%lf, background max sd=%lf, too big, not bg\n", sd, bg_features.max_sd);
#endif
        return false;
    }
    if (max_min_diff > bg_features.max_of__max_min_diff) {
#ifdef DB_DISPLAY
        printf("max_min_diff=%lf, background max max_min_diff=%lf, too big, not bg\n", max_min_diff , bg_features.max_of__max_min_diff);
#endif
        return false;
    }

    return true;
}


int dr_group_image_buffer::plot_intensity_by_frame() {
    int n = images.size();
    double intensity[n];
    double frames[n];
    int i,j,k;
    int h = resolution_y, w = resolution_x;
    for (k=0; k<n; k++) {
        intensity[k] = 0;
        frames[k] = k;
        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                intensity[k]+=qGray(images[k].pixel(j,i));
            }
        }
        intensity[k]/=(h*w);
        printf("%lf, %lf\n", frames[k], intensity[k]);
    }

    plotResults(frames, intensity, n, 'x');
    return 0;
}


int dr_group_image_buffer::dilate_bin_image(QVector<QImage*> *imgs, int diameter) {
    for (int k = 0; k < imgs->size(); k++) {
        QImage *img2 = new QImage(*imgs->at(k));

        int radius = diameter/2;
        int h = imgs->at(0)->height();
        int w = imgs->at(0)->width();
        int i;
        int j;
        for (i=0; i<h; i++) {
            int j;
            for (j=0; j<w; j++) {
                if (  qGray(imgs->at(k)->pixel(j,i)) != 0 ) {
                    //scan from height
                    int x,y, x1, x2;
                    for (y = i-radius; y<= i+radius; y++) {
                        //calculate x range
                        x1 = (int)( j - sqrt( pow(radius,2) - pow((y-i),2))  );
                        x2 = (int)( j + sqrt( pow(radius,2) - pow((y-i),2))  );
                        for (x = x1; x<= x2; x++) {
                            if(x>=0 && x<w && y>=0 && y<h) {
                                img2->setPixel(x,y, 255);
                            }
                        }
                    }
                }
            }
        }
        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                imgs->at(k)->setPixel(j,i, qGray(img2->pixel(j,i)));
            }
        }
        delete img2;
    }
    return 0;
}


int dr_group_image_buffer::erode_bin_image(QVector<QImage *> *imgs, int diameter) {
    for (int k = 0; k <imgs->size(); k++) {
        QImage *img2 = new QImage(*imgs->at(k));
        int radius = diameter/2;
        int h = imgs->at(0)->height();
        int w = imgs->at(0)->width();
        int i;
        int j;
        bool flag = false;
        for (i=0; i<h; i++) {
            int j;
            for (j=0; j<w; j++) {
                flag = false;
                //scan from height
                int x,y, x1, x2;
                for (y = i-radius; y<= i+radius; y++) {
                    if (flag) {
                        break;
                    }
                    //calculate x range
                    x1 = (int)( j - sqrt( pow(radius,2) - pow((y-i),2))  );
                    x2 = (int)( j + sqrt( pow(radius,2) - pow((y-i),2))  );
                    for (x = x1; x<= x2; x++) {
                        if(x<0 || x>=w || y<0 || y>=h || qGray(imgs->at(k)->pixel(x,y)) == 0) {
                            img2->setPixel(j,i,0);
                            flag = true;
                            break;
                        }
                    }
                }
                if (!flag) {
                    img2->setPixel(j,i,255);
                }
            }
        }
        for (i=0; i<h; i++) {
            for (j=0; j<w; j++) {
                imgs->at(k)->setPixel(j,i, qGray(img2->pixel(j,i)));
            }
        }
        delete img2;
    }
    return 0;
}


int dr_group_image_buffer::open_bin_image(QVector<QImage *> *imgs, int erode_diameter, int dilate_diameter, QString * debug_image_name) {
    erode_bin_image(imgs, erode_diameter);
    if (*(debug_image_name) != "no") {
    imgs->at(0)->save(*(debug_image_name)+"erode_img.pgm",0,100); }
    dilate_bin_image(imgs, dilate_diameter);
    return 0;
}


int dr_group_image_buffer::plot_diff_curve() {
    int h = resolution_y, w = resolution_x;
    int n = h*w;
    int i,j,k;

    int data[n];
    k = 0;
    for (i=0; i<h; i++) {
        for (j=0; j<w; j++) {
            data[k++]=qGray(diff_image->pixel(j,i));
        }
    }
    if (n!=k)
        return 1;
    return 0;
}


QVector<double> * dr_group_image_buffer::learning_from_background_line(double slope, double offset, double mean) {
    // this is where we use stats in determing which pixel is foreground and which is background.
    // the offset of the fitted line is the variance when itensity is zero, which is the variance of the gaussian noise
    // the mean is the mean intensity of the pixel you want to compute on, which is the mean of the background distribution (possion + gaussian noise)..
    // on Russ' notes, sensitivity is defined as photons/count, thus the slope is 1 over sensitivity here
    QVector<double> * noiseGaussianParam = new QVector<double>(2,0.0);
    double g_mean = mean;
    double g_variance = offset + slope * mean;
    //   printf("offset, slope, mean: %f %f %f\n", offset, slope, mean);
    //  printf("in learing from background, g_variance value: %f\n", g_variance);
    (*noiseGaussianParam)[0] = g_mean;
    (*noiseGaussianParam)[1] = g_variance;
    return noiseGaussianParam;
}


bool dr_group_image_buffer::is_background_evt(QVector<double> * background_distribution_param, int max_pixel_value, int min_pixel_value, int number_of_samples) {
    // first center the values by shifting with the mean
    double max_pixel_value_double = static_cast<double>(max_pixel_value);
    double min_pixel_value_double = static_cast<double>(min_pixel_value);
    // printf("background distribution param dimension : 1: %d, 2: %d \n", background_distribution_param->dim1(), background_distribution_param->dim2());
    double mean = (*background_distribution_param)[0];
    double centered_max_pixel_value = max_pixel_value_double - mean;
    double centered_min_pixel_value = min_pixel_value_double - mean;
    // then normalize it by dividing it over the variance
    double variance = (*background_distribution_param)[1];
    double centered_normalized_max_pixel_value = centered_max_pixel_value / sqrt(variance);
    double centered_normalized_min_pixel_value = centered_min_pixel_value / sqrt(variance);

    double probability = phi(centered_normalized_min_pixel_value, centered_normalized_max_pixel_value);
    double extreme_pr = 1 - pow(probability, number_of_samples);

    if (extreme_pr >= static_cast<double>(BACKGROUND_THRESHOLD)/100.0) {
        return true;
    }
    return false;
}


bool dr_group_image_buffer::is_background_evt2(QVector<double> * background_distribution_param, int max_pixel_value, int number_of_samples) {
    // first center the values by shifting with the mean
    double max_pixel_value_double = static_cast<double>(max_pixel_value);
    double mean = (*background_distribution_param)[0];
    double centered_max_pixel_value = max_pixel_value_double - mean;
    // then normalize it by dividing it over the variance
    double variance = (*background_distribution_param)[1];
    double centered_normalized_max_pixel_value = centered_max_pixel_value / sqrt(variance);
    // there might be other problems..

    double probability = 0.5 + 0.5 * erf(centered_normalized_max_pixel_value / sqrt(2.0));
    double extreme_pr = 1- pow((probability), number_of_samples);

    if (extreme_pr >= static_cast<double>(BACKGROUND_THRESHOLD)/100.0) {
        return true;
    }
    return false;
}


bool dr_group_image_buffer::is_background_skew_kurto(QVector<double> * intensity_over_time, double skewness_threshold, double kurtosis_threshold) {
    // assume the intensity over time vector is already normalized by subtracting the mean and dividing by the variance
    double m2 = 0;
    for(int i = 0; i < intensity_over_time->size(); i++) {
    }
    double skewness = compute_skewness(intensity_over_time, &m2);
    double kurtosis = compute_kurtosis(intensity_over_time, &m2);
    if (abs(skewness) > skewness_threshold || abs(abs(kurtosis) - 3) > kurtosis_threshold) {
        return false;
    }
    return true;
}


// the standard_gaussian_cumulative_table should also have 601 items
bool dr_group_image_buffer::is_background_ks_test(QVector<double> * intensity_over_time, QVector<double> * standard_gaussian_cumulative_table, double threshold) {
    // assume  the intensity over time vector is already normalized by subtracting the mean and dividing by the standard deviation
    double max_diff = 0.0;
    // assume this histogram is from -6 to 6, 601 bins
    QVector<double> pixel_distribution(601,0.0);
    // fill in this distribution ....
    int n = intensity_over_time->size();
    for (int i = 0; i < n; i++) {
        int bin = static_cast<int>(100.0 * (*intensity_over_time)[i]);
        if (bin < -300) {
            pixel_distribution[0] += 1;
        }
        else if (bin > 300) {
            pixel_distribution[600] += 1;
        }
        else {
            pixel_distribution[bin+300] += 1;
        }
    }
    // normalize the pixel distribution so that it integrates up to 1
    for (int i = 0; i < 601; i++) {
        pixel_distribution[i] = pixel_distribution[i] / static_cast<double>(n);
    }

    double current_cumulative_val = 0.0;
    QVector<double> pixel_cumulative_distribution(601,0.0);
    // calculate the cumulative distribution
    for (int i = 0; i < 601; i++) {
        current_cumulative_val += pixel_distribution[i];
        pixel_cumulative_distribution[i] = current_cumulative_val;
        // get current diff, update the max ciff
        double curr_diff = abs(pixel_cumulative_distribution[i] - (*standard_gaussian_cumulative_table)[i]);
        if (curr_diff > max_diff) {
            max_diff = curr_diff;
        }
    }

    if (max_diff < threshold) {
        return true;
    }
    return false;
}


double dr_group_image_buffer::phi(double x1, double x2) {
    return (erf(x2/std::sqrt(2.0)) - erf(x1/std::sqrt(2.0))) / 2.0 ;
}


void dr_group_image_buffer::divide_video_into_segments(int n) {
    n = 1;
    return;
}


double dr_group_image_buffer::compute_skewness(QVector<double> * intensity_over_time, double * m2_out) {
    // the mean of the data should be zero
    int n = intensity_over_time->size();
    double sum_for_m3 = 0.0;
    double sum_for_m2 = 0.0;
    for (int i = 0; i < n; i++) {
        sum_for_m3 += pow( (*intensity_over_time)[i], 3.0 );
        sum_for_m2 += pow( (*intensity_over_time)[i], 2.0 );
    }
    double m3 = sum_for_m3 / static_cast<double>(n);
    double m2 = sum_for_m2 / static_cast<double>(n);
    double g1 = m3 / pow(m2, 1.5);
    double G1 = g1 * sqrt(static_cast<double>(n * (n - 1))) / static_cast<double>(n-2);
    *(m2_out) = m2;
    return G1;
}


double dr_group_image_buffer::compute_kurtosis(QVector<double> * intensity_over_time, double * m2_out) {
    int n = intensity_over_time->size();
    double sum_for_m4 = 0.0;
    for (int i = 0; i < n; i++) {
        sum_for_m4 += pow( (*intensity_over_time)[i], 4.0 );
    }
    double m4 = sum_for_m4 / static_cast<double>(n);
    double a4 = m4 / pow(*(m2_out), 2.0);
    double g2 = a4 - 3;
    double G2 = (n-1)*( (n-1)*g2 + 6 ) / ( (n-2) * (n-3) );
    return G2;
}


void dr_group_image_buffer::compute_one_mean_image(QImage *one_mean_image) {
    if (this->bin_diff_images->size() != NUMBER_OF_SEGS) {
        cout << "the one mean image should not be computed until all the binary mask images have been generated!" << endl;
        return;
    }
    cout << "compute the one mean image!" << endl;

    int h = images[0].height();
    int w = images[0].width();
    int64_t *sum = new int64_t[w*h];
    int64_t *frame_count = new int64_t[w*h];
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            sum[i*w + j] = 0ll;
            frame_count[i*w + j] = 0ll;
        }
    }
    printf("in one mean image, h and w: %d %d\n", h, w);
    int n = images.size() / NUMBER_OF_SEGS;
    for (int k = 0; k < NUMBER_OF_SEGS-1; k++) {
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                if (qGray(bin_diff_images->at(k)->pixel(j,i)) == 0) { // background
                    for (int x = k * n; x < k*n + n; x++) {
                        sum[i*w + j] += qGray(images[x].pixel(j,i));
                        frame_count[i*w+j] += 1;
                    }
                }
            }
        }
    }

    // handle the last segment
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            if (qGray(bin_diff_images->last()->pixel(j,i)) == 0) { // background
                for (int x = n * (NUMBER_OF_SEGS-1); x < images.size(); x++) {
                    sum[i*w + j] += qGray(images[x].pixel(j,i));

                    frame_count[i*w + j] += 1;
                }
            }
        }
    }

    // calculate the one mean image
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            int curr_one_mean_image_pixel;
            if (frame_count[i*w+j] > 0) {
                curr_one_mean_image_pixel = sum[i*w + j] / frame_count[i*w+j];// just try using integer division
            }
            else {
                curr_one_mean_image_pixel = 0;
            }
            //           printf("sum and frame count i w + j: %d %d\n", sum[i*w+j], frame_count[i*w+j]);
            one_mean_image->setPixel(j,i, qRgb(curr_one_mean_image_pixel,curr_one_mean_image_pixel,curr_one_mean_image_pixel));
        }
    }
    this->one_mean_image = one_mean_image;
    return;
}


void dr_group_image_buffer::output_pixels_intensity_over_time(bool use_full_image) {
    QFile file("pixelOverTime.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    int dim1; // y
    int dim2; // x
    if (!use_full_image) {
        dim1 = 200;
        dim2 = 200;
    }
    else {
        QImage * first_image = &(images[0]);
        dim1 = first_image->height();
        dim2 = first_image->width();
    }
    for (int i = 0; i < dim1; i++) {
        for (int j = 0; j < dim2; j++) {
            for (int k = 0; k < images.size(); k++) {
                QImage * currImage = &(images[k]);
                out << qGray(currImage->pixel(j,i)) << " ";
            }
            out << "\n";
        }
    }
    file.close();
    return;
}


void dr_group_image_buffer::compute_corr_image(QImage * image) {
    int i ,j, k;
    int video_length = images.size();
    int width = images[0].width(); // x
    int height = images[0].height(); // y

    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            QVector<int> * curr_pixel_values = new QVector<int>();
            // 8 neighbors
            QVector<int> * neighbor_pixel1_values = new QVector<int>(); // upper left
            QVector<int> * neighbor_pixel2_values = new QVector<int>(); // upper
            QVector<int> * neighbor_pixel3_values = new QVector<int>(); // upper right
            QVector<int> * neighbor_pixel4_values = new QVector<int>(); // right
            QVector<int> * neighbor_pixel5_values = new QVector<int>(); // lower right
            QVector<int> * neighbor_pixel6_values = new QVector<int>(); // lower
            QVector<int> * neighbor_pixel7_values = new QVector<int>(); // lower left
            QVector<int> * neighbor_pixel8_values = new QVector<int>(); // left

            // fill in corr_pixel_values
            for (k = 0; k < video_length; k++) {
                // we can change this line to support more images...
                QImage * curr_image = &(images[k]);
                curr_pixel_values->append(qGray(curr_image->pixel(i, j)));
                if (i != 0)  { // current pixel not on left most column
                    neighbor_pixel8_values->append(qGray(curr_image->pixel(i-1,j)));
                }
                if (j != 0) { // current pixel not on first row
                    neighbor_pixel2_values->append(qGray(curr_image->pixel(i, j-1)));
                }
                if (i != width - 1) { // current pixel not on right most column
                    neighbor_pixel4_values->append(qGray(curr_image->pixel(i+1, j)));
                }
                if (j != height - 1) { // current pixel not on last row
                    neighbor_pixel6_values->append(qGray(curr_image->pixel(i, j+1)));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel1_values->append(qGray(curr_image->pixel(i-1, j-1)));
                }
                if (i != 0 && j != height -1) {
                    neighbor_pixel3_values->append(qGray(curr_image->pixel(i-1, j+1)));
                }
                if (i != width -1 && j != 0) {
                    neighbor_pixel7_values->append(qGray(curr_image->pixel(i+1, j-1)));
                }
                if (i != width -1 && j != height -1 ) {
                    neighbor_pixel5_values->append(qGray(curr_image->pixel(i+1, j+1)));
                }
            }
            QVector<double> corr_values;
            // commpute correlation between center pixel and each neighbor
            if (i != 0 && j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
            }
            if (j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
            }
            if (i != 0 && j != height -1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel3_values)));
            }
            if (i != width - 1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel4_values)));
            }
            if (i != width -1 && j != height - 1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel5_values)));
            }
            if (j != height -1 ) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel6_values)));
            }
            if (i != width -1 && j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel7_values)));
            }
            if (i != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel8_values)));
            }
            double max_corr = 0.0;
            bool changed = false;
            for (int l = 0; l < corr_values.size(); l++) {
                if (fabs(corr_values.at(l)) > max_corr) {
                    max_corr = corr_values.at(l);
                    changed = true;
                }
            }
            /**
            if (max_corr == 0.0) {
                printf("corr_values: ");
                for (int l = 0; l < corr_values.size(); l++) {
                    printf("%f ", corr_values.at(l));
                }
                printf("\n");
             }
            */
            int corr_image_pixel_intensity = static_cast<int>(max_corr*255);
            image->setPixel(i, j, qRgb(corr_image_pixel_intensity, corr_image_pixel_intensity, corr_image_pixel_intensity));
        }
    }

}


void dr_group_image_buffer::compute_corr_image_memory_friendly(QString * image_file_name_before_index, int number_of_frames, bool start_with_zero, QImage * image) {
    // ths
    int i ,j, k;
    int video_length = 100;
    int width = images[0].width(); // x
    int height = images[0].height(); // y

    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            QVector<int> * curr_pixel_values = new QVector<int>();
            // 8 neighbors
            QVector<int> * neighbor_pixel1_values = new QVector<int>(); // upper left
            QVector<int> * neighbor_pixel2_values = new QVector<int>(); // upper
            QVector<int> * neighbor_pixel3_values = new QVector<int>(); // upper right
            QVector<int> * neighbor_pixel4_values = new QVector<int>(); // right
            QVector<int> * neighbor_pixel5_values = new QVector<int>(); // lower right
            QVector<int> * neighbor_pixel6_values = new QVector<int>(); // lower
            QVector<int> * neighbor_pixel7_values = new QVector<int>(); // lower left
            QVector<int> * neighbor_pixel8_values = new QVector<int>(); // left

            QImage * curr_image = NULL;

            for (k = 0; k < video_length; k++) {
                // we can change this line to support more images...
                curr_image = &(images[k]);
                curr_pixel_values->append(qGray(curr_image->pixel(i, j)));
                if (i != 0)  { // current pixel not on left most column
                    neighbor_pixel8_values->append(qGray(curr_image->pixel(i-1,j)));
                }
                if (j != 0) { // current pixel not on first row
                    neighbor_pixel2_values->append(qGray(curr_image->pixel(i, j-1)));
                }
                if (i != width - 1) { // current pixel not on right most column
                    neighbor_pixel4_values->append(qGray(curr_image->pixel(i+1, j)));
                }
                if (j != height - 1) { // current pixel not on last row
                    neighbor_pixel6_values->append(qGray(curr_image->pixel(i, j+1)));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel1_values->append(qGray(curr_image->pixel(i-1, j-1)));
                }
                if (i != 0 && j != height -1) {
                    neighbor_pixel3_values->append(qGray(curr_image->pixel(i-1, j+1)));
                }
                if (i != width -1 && j != 0) {
                    neighbor_pixel7_values->append(qGray(curr_image->pixel(i+1, j-1)));
                }
                if (i != width -1 && j != height -1 ) {
                    neighbor_pixel5_values->append(qGray(curr_image->pixel(i+1, j+1)));
                }
            }

            // handle the images out of the first 100 images
            int number_of_frames_remaining = number_of_frames - 100;

            for (k = 0; k < number_of_frames_remaining; k++) {
                QString * digit_str = new QString();

                if (start_with_zero) {
                    digit_str->append(digit_str->sprintf("%4d", k+100));
                }
                else { // start with 1
                    digit_str->append(digit_str->sprintf("%4d", k+101));
                }

                QString * image_file_name = &(image_file_name_before_index->append(digit_str));

                curr_image = new QImage(*image_file_name);

                curr_pixel_values->append(qGray(curr_image->pixel(i, j)));

                if (i != 0)  { // current pixel not on left most column
                    neighbor_pixel8_values->append(qGray(curr_image->pixel(i-1,j)));
                }
                if (j != 0) { // current pixel not on first row
                    neighbor_pixel2_values->append(qGray(curr_image->pixel(i, j-1)));
                }
                if (i != width - 1) { // current pixel not on right most column
                    neighbor_pixel4_values->append(qGray(curr_image->pixel(i+1, j)));
                }
                if (j != height - 1) { // current pixel not on last row
                    neighbor_pixel6_values->append(qGray(curr_image->pixel(i, j+1)));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel1_values->append(qGray(curr_image->pixel(i-1, j-1)));
                }
                if (i != 0 && j != height -1) {
                    neighbor_pixel3_values->append(qGray(curr_image->pixel(i-1, j+1)));
                }
                if (i != width -1 && j != 0) {
                    neighbor_pixel7_values->append(qGray(curr_image->pixel(i+1, j-1)));
                }
                if (i != width -1 && j != height -1 ) {
                    neighbor_pixel5_values->append(qGray(curr_image->pixel(i+1, j+1)));
                }
            }
            QVector<double> corr_values;
            // compute correlation between center pixel and each neighbor
            if (i != 0 && j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
            }
            if (j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
            }
            if (i != 0 && j != height -1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel3_values)));
            }
            if (i != width - 1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel4_values)));
            }
            if (i != width -1 && j != height - 1) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel5_values)));
            }
            if (j != height -1 ) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel6_values)));
            }
            if (i != width -1 && j != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel7_values)));
            }
            if (i != 0) {
                corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel8_values)));
            }
            double max_corr = 0.0;
            bool changed = false;
            for (int l = 0; l < corr_values.size(); l++) {
                if (fabs(corr_values.at(l)) > max_corr) {
                    max_corr = corr_values.at(l);
                    changed = true;
                }
            }
            int corr_image_pixel_intensity = static_cast<int>(max_corr*255);
            image->setPixel(i, j, qRgb(corr_image_pixel_intensity, corr_image_pixel_intensity, corr_image_pixel_intensity));
        }
    }

}


void dr_group_image_buffer::compute_corr_image_online(QImage * new_img, QVector<QVector<QVector<double> *> *> * old_corr_values, QVector<QVector<double> *> * mean_img_double) {
    // in this function, we basically make use of the previous computed 2 images: mean and variance.
    //  everytime we see an image, we load it and update the correlation image.
    // therefore we don't actually need the buffer: memory for just 3 images is fine.
    // passes: 2 passes
    //
    // we may change it to editable with and height numbers,
    // for now we just do 3x3

    // update the neighbors

    int width = new_img->width();
    int height = new_img->height();
    int i,j;
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            // 1 2 3
            // 8 C 4
            // 7 6 5
            double inc1, inc2, inc3, inc4, inc5, inc6, inc7, inc8;
            int curr_pixel_val = qGray(new_img->pixel(i, j));
            double x_minus_miuX = static_cast<double>(curr_pixel_val) - mean_img_double->at(j)->at(i);
            if (i != 0)  { // current pixel not on left most column
                inc8 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i-1,j))) - mean_img_double->at(j)->at(i-1));
                old_corr_values->at(j)->at(i)->replace(8, old_corr_values->at(j)->at(i)->at(8) + inc8);
            }
            if (j != 0) { // current pixel not on first row
                inc2 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i,j-1))) - mean_img_double->at(j-1)->at(i));
                old_corr_values->at(j)->at(i)->replace(2, old_corr_values->at(j)->at(i)->at(2) + inc2);
            }
            if (i != width - 1) { // current pixel not on right most column
                inc4 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i+1,j))) - mean_img_double->at(j)->at(i+1));
                old_corr_values->at(j)->at(i)->replace(4, old_corr_values->at(j)->at(i)->at(4) + inc4);
            }
            if (j != height - 1) { // current pixel not on last row
                inc6 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i,j+1))) - mean_img_double->at(j+1)->at(i));
                old_corr_values->at(j)->at(i)->replace(6, old_corr_values->at(j)->at(i)->at(6) + inc6);
            }
            if (i != 0 && j != 0) {
                inc1 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i-1,j-1))) - mean_img_double->at(j-1)->at(i-1) );
                old_corr_values->at(j)->at(i)->replace(1, old_corr_values->at(j)->at(i)->at(1) + inc1);
            }
            if (i != 0 && j != height -1) {
                inc3 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i-1,j+1))) - mean_img_double->at(j+1)->at(i-1));
                old_corr_values->at(j)->at(i)->replace(3, old_corr_values->at(j)->at(i)->at(3) + inc3);
            }
            if (i != width -1 && j != 0) {
                inc7 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i+1,j-1))) - mean_img_double->at(j-1)->at(i+1) );
                old_corr_values->at(j)->at(i)->replace(7, old_corr_values->at(j)->at(i)->at(7) + inc7);
            }
            if (i != width -1 && j != height -1 ) {
                inc5 = x_minus_miuX * (static_cast<double>(qGray(new_img->pixel(i+1,j+1))) - mean_img_double->at(j+1)->at(i+1) );
                old_corr_values->at(j)->at(i)->replace(5, old_corr_values->at(j)->at(i)->at(5) + inc5);
            }
            old_corr_values->at(j)->at(i)->replace(0, 0.0);
        }
    }
}


void dr_group_image_buffer::dr_group_image_buffer::compute_corr_image_with_windowing_memory_friendly() {
    // in this function. we do the windowing

}


void dr_group_image_buffer::compute_corr_image_with_windowing(QVector<QImage*> *out_corr_images, int window_size) {
    int i, j, k, h;  // h is for iteration inside the window
    int video_length = images.size();
    int width = images[0].width(); // x
    int height = images[0].height(); // y

    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {

            for (k = 0; k < video_length-window_size; k += window_size / 2) {
                QVector<int> * curr_pixel_values = new QVector<int>();
                // 8 neighbors
                QVector<int> * neighbor_pixel1_values = new QVector<int>(); // upper left
                QVector<int> * neighbor_pixel2_values = new QVector<int>(); // upper
                QVector<int> * neighbor_pixel3_values = new QVector<int>(); // upper right
                QVector<int> * neighbor_pixel4_values = new QVector<int>(); // right
                QVector<int> * neighbor_pixel5_values = new QVector<int>(); // lower right
                QVector<int> * neighbor_pixel6_values = new QVector<int>(); // lower
                QVector<int> * neighbor_pixel7_values = new QVector<int>(); // lower left
                QVector<int> * neighbor_pixel8_values = new QVector<int>(); // left

                QImage * curr_corr_image = &(images[0]);
                out_corr_images->append(curr_corr_image);
                // for every k (starting point of the window), compute a correlation value
                // fill in corr_pixel_values
                for (h = 0; h < window_size; h++) {
                    if (k+h >= video_length) {
                        break;
                    }
                    QImage * curr_image = &(images[k+h]);
                    curr_pixel_values->append(qGray(curr_image->pixel(i, j)));
                    if (i != 0)  { // current pixel not on left most column
                        neighbor_pixel8_values->append(qGray(curr_image->pixel(i-1,j)));
                    }
                    if (j != 0) { // current pixel not on first row
                        neighbor_pixel2_values->append(qGray(curr_image->pixel(i, j-1)));
                    }
                    if (i != width - 1) { // current pixel not on right most column
                        neighbor_pixel4_values->append(qGray(curr_image->pixel(i+1, j)));
                    }
                    if (j != height - 1) { // current pixel not on last row
                        neighbor_pixel6_values->append(qGray(curr_image->pixel(i, j+1)));
                    }
                    if (i != 0 && j != 0) {
                        neighbor_pixel1_values->append(qGray(curr_image->pixel(i-1, j-1)));
                    }
                    if (i != 0 && j != height -1) {
                        neighbor_pixel3_values->append(qGray(curr_image->pixel(i-1, j+1)));
                    }
                    if (i != width -1 && j != 0) {
                        neighbor_pixel7_values->append(qGray(curr_image->pixel(i+1, j-1)));
                    }
                    if (i != width -1 && j != height -1 ) {
                        neighbor_pixel5_values->append(qGray(curr_image->pixel(i+1, j+1)));
                    }
                }
                QVector<double> corr_values;
                // commpute correlation between center pixel and each neighbor
                if (i != 0 && j != 0) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
                }
                if (j != 0) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel2_values)));
                }
                if (i != 0 && j != height -1) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel3_values)));
                }
                if (i != width - 1) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel4_values)));
                }
                if (i != width -1 && j != height - 1) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel5_values)));
                }
                if (j != height -1 ) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel6_values)));
                }
                if (i != width -1 && j != 0) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel7_values)));
                }
                if (i != 0) {
                    corr_values.append(fabs(compute_correlation_between_two_vectors(curr_pixel_values, neighbor_pixel8_values)));
                }
                double max_corr = 0.0;
                for (int l = 0; l < corr_values.size(); l++) {
                    if (fabs(corr_values.at(l)) > max_corr) {
                        max_corr = corr_values.at(l);
                    }
                }
                int corr_image_pixel_intensity = static_cast<int>(max_corr*255);
                curr_corr_image->setPixel(i, j, corr_image_pixel_intensity);
            }
        }
    }
}


double dr_group_image_buffer::compute_correlation_between_two_vectors(QVector<int> * vector1, QVector<int> * vector2) {
    // compute mean of v1 and of v2
    // corr(v1, v2) = E[(v1 - mean(v1))(v2 - mean(v2)] / (std(v1) * std(v2))
    long long int sum_v1 = 0;
    long long int sum_v2 = 0;
    int i;
    for (i=0; i<vector1->size(); i++) {
        sum_v1 += vector1->at(i);
        sum_v2 += vector2->at(i);
    }
    double mean_v1 = static_cast<double>(sum_v1) / static_cast<double>(vector1->size());
    double mean_v2 = static_cast<double>(sum_v2) / static_cast<double>(vector2->size());

    // compute the standard dev and convariance
    long int sum_of_diff_from_mean_v1 = 0;
    long int sum_of_diff_from_mean_v2 = 0;
    long int sum_of_diff_v1_times_v2 = 0;
    for (i = 0; i<vector1->size(); i++) {
        sum_of_diff_from_mean_v1 += (static_cast<double>(vector1->at(i))-mean_v1)*(static_cast<double>(vector1->at(i))-mean_v1);
        sum_of_diff_from_mean_v2 += (static_cast<double>(vector2->at(i))-mean_v2)*(static_cast<double>(vector2->at(i))-mean_v2);
        sum_of_diff_v1_times_v2 += (static_cast<double>(vector1->at(i))-mean_v1) * (static_cast<double>(vector2->at(i))-mean_v2);
    }
    double std_dev_v1 = sqrt(static_cast<double>(sum_of_diff_from_mean_v1) / static_cast<double>(vector1->size()));
    double std_dev_v2 = sqrt(static_cast<double>(sum_of_diff_from_mean_v2) / static_cast<double>(vector2->size()));
    if (std_dev_v1 <= 0.001) {
        std_dev_v1 = 1.0;
    }
    if (std_dev_v2 <= 0.001) {
        std_dev_v2 = 1.0;
    }
    double cov = static_cast<double>(sum_of_diff_v1_times_v2) / static_cast<double>(vector1->size());
    return cov / (std_dev_v1 * std_dev_v2);
}


double dr_group_image_buffer::compute_correlation_between_two_pixel_time_series(QVector<int> * vector1, QVector<int> * vector2, int x1, int y1, int x2, int y2) {
    // compute mean of v1 and of v2
    // corr(v1, v2) = E[(v1 - mean(v1))(v2 - mean(v2)] / (std(v1) * std(v2))

    double mean_v1 = static_cast<double>(qGray(this->mean_images->at(0)->pixel(x1,y1)));
    double mean_v2 = static_cast<double>(qGray(this->mean_images->at(0)->pixel(x2,y2)));

    // compute the standard dev and convariance
    long int sum_of_diff_from_mean_v1 = 0;
    long int sum_of_diff_from_mean_v2 = 0;
    long int sum_of_diff_v1_times_v2 = 0;
    int i;
    for (i = 0; i<vector1->size(); i++) {
        sum_of_diff_from_mean_v1 += (static_cast<double>(vector1->at(i))-mean_v1)*(static_cast<double>(vector1->at(i))-mean_v1);
        sum_of_diff_from_mean_v2 += (static_cast<double>(vector2->at(i))-mean_v2)*(static_cast<double>(vector2->at(i))-mean_v2);
        sum_of_diff_v1_times_v2 += (static_cast<double>(vector1->at(i))-mean_v1) * (static_cast<double>(vector2->at(i))-mean_v2);
    }
    double std_dev_v1 = sqrt(static_cast<double>(sum_of_diff_from_mean_v1) / static_cast<double>(vector1->size()));
    double std_dev_v2 = sqrt(static_cast<double>(sum_of_diff_from_mean_v2) / static_cast<double>(vector2->size()));
    if (std_dev_v1 <= 0.001) {
        std_dev_v1 = 1.0;
    }
    if (std_dev_v2 <= 0.001) {
        std_dev_v2 = 1.0;
    }
    double cov = static_cast<double>(sum_of_diff_v1_times_v2) / static_cast<double>(vector1->size());
    return cov / (std_dev_v1 * std_dev_v2);
}


int compare_two_doubles(double a, double b) {
    return a - b;
}


double dr_group_image_buffer::anderson_darling_test(QVector<double> * x, double alpha) {
    if (x->size() < 7) {
        cout << "sample size must be greater than 7";
        return -1.0;
    }
    // convert x into a normal cumulative distribution function
    // first compute the mean and stdev of x
    double mean_x = 0.0;
    double sum_x = 0.0;
    for (int i = 0; i < x->size(); i++) {
        sum_x += x->at(i);
    }
    mean_x = sum_x / static_cast<double>(x->size());
    double sum_dev_x = 0.0;
    for (int i = 0; i < x->size(); i++) {
        sum_dev_x += pow((x->at(i) - mean_x) ,2);
    }
    double stdev_x = sum_dev_x / x->size();
    // standardize x
    QVector<double> * standardized_x = new QVector<double>();
    for (int i = 0; i < x->size(); i++) {
        standardized_x->append((x->at(i) - mean_x)/stdev_x);
    }
    // sort standardized x
    qSort(standardized_x->begin(), standardized_x->end(), compare_two_doubles);
    // construct the fx vector
    QVector<double> *fx = new QVector<double>(x->size(),0.0);
    for (int i = 0; i < x->size(); i++) {
        fx->replace(i,0.5 + 0.5 * erf(standardized_x->at(i)/sqrt(2.0)));
    }
    // compute S
    double S = 0.0;
    for (int i = 1; i < x->size()+1; i++) {
        S += ((2.0 * static_cast<double>(i)) - 1.0 ) *(log(fx->at(static_cast<double>(i)-1.0)) + log(fx->at(x->size()-static_cast<double>(i))) )/ static_cast<double>(x->size());
    }
    double AD2 = -1.0 * (static_cast<double>(x->size()) + S);
    // correction factor for small smaples size, case normal
    double AD2a = AD2 * (1.0 + 0.75 / (static_cast<double>(x->size())) + 2.25 / (static_cast<double>(x->size()) * static_cast<double>(x->size())));
    double P = 0.0;
    if (AD2a >= 0.0 && AD2a < 0.2) {
        P = 1 - exp(-13.436 + 101.14*AD2a + 223.73*AD2a*AD2a);
    }
    else {
        if (AD2a >= 0.2 && AD2a < 0.34) {
            P = 1 - exp(-8.318 + 42.796*AD2a - 59.938*AD2a*AD2a);
        }
        else {
            if (AD2a >= 0.340 && AD2a < 0.600) {
                P = exp(0.9177 - 4.279*AD2a - 1.38*AD2a*AD2a);
            }
            else {
                P = exp(1.2937 - 5.709*AD2a + 0.0186*AD2a*AD2a);
            }
        }
    }
    alpha = 0.0;
    return P;
}


int dr_group_image_buffer::online_mean_and_variance(int n, QVector<QVector<double> *> * mean_img_double, QVector<QVector<double> *> *  M2_img, QImage * new_img){
    //  n = 0;
    //  mean = 0;
    //  M2 = 0;
    //  for x in data:
    //      n = n+1;
    //      delta = x - mean
    //      mean = mean + delta/n
    //      M2 = M2 + delta*(x-mean)

    //  variance = M2 / (n-1)
    //  return variance

    // this function returns updated n as return value..
    int out_n = n+1;
    int h = new_img->height();
    int w = new_img->width();

    int x,y;
    double delta;
    double new_mean_pixel, new_M2_pixel;
    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            delta = static_cast<double>(qGray(new_img->pixel(x,y))) - mean_img_double->at(y)->at(x);
            new_mean_pixel = mean_img_double->at(y)->at(x) + ( delta / static_cast<double>(out_n) );
            mean_img_double->at(y)->replace(x, new_mean_pixel);
            new_M2_pixel = M2_img->at(y)->at(x) + delta * (static_cast<double>(qGray(new_img->pixel(x,y)))  - new_mean_pixel);

            M2_img->at(y)->replace(x, new_M2_pixel);
        }
    }
    return out_n;
}


void dr_group_image_buffer::set_mean_img(QImage * mean_img) {
    this->mean_images->replace(0,mean_img);
    this->one_mean_image = mean_img;
}


void dr_group_image_buffer::set_variance_img(QImage * variance_img) {
    this->var_images->replace(0, variance_img);
}


int dr_group_image_buffer::two_levels_an_image_other_online(QVector<QImage*> *imgs_out, int method, int threshold) {
    // method 1: KS test
    // method 2: correlation
    /// TODO: DONE on this
    if (method == 1) {
        printf("Performing K-S test background computation...\n");
        printf("nothing implemented for this right now, quit....\n");
        return -1;
    }
    if (method == 2) {
        printf("Performing correlation background computation...\n");
    }

    int x,y;
    //convert img_in to a column vector
    int h, w, k;
    QImage * corr_map = this->corr_image;
    h = corr_map->height();
    w = corr_map->width();
    int corr_threshold;
    if (method == 2) {
        if (!USE_HARD_CODED_THRESHOLD)  { // then we compute a threshold
            QVector<int> * corr_values = new QVector<int>();
            for (int i = 0; i < h; i++) {
                for (int j = 0; j < w; j++) {
                    int pixel = qGray(corr_map->pixel(j,i));
                    if (corr_values->size() == 0) {
                        corr_values->append(pixel);
                    }
                    else { /// TODO: inserting sort here..
                        // this is to sort the values decreasing order
                        for (k = 0; k < corr_values->size(); k++) {
                            if (pixel > corr_values->at(k)) {
                                corr_values->insert(k,pixel);
                                break;
                            }
                        }
                    }
                }
            }
            // threshold the correlation value based on in the input parameter
            // default: 75 % percent quantile
            int threshold_index = (corr_values->size()-1) * (100-threshold) / 100;
            cout << "correlation threshold index: " << threshold_index << endl;
            corr_threshold = corr_values->at(static_cast<int>(threshold_index));
        }
        else {
            corr_threshold = HARD_CODED_THRESHOLD;
            cout << "use hard coded threshold: ";
        }
    }
    printf("corr threshold: %d \n", corr_threshold);

    // for now we do not do windowing
    int white = 0;
    int black = 0;
    QImage * img_out = new QImage(*(this->corr_image));
    // find the threshold based on alpha value
    for (x = 0; x < h; x++) {
        for (y = 0; y < w; y++) {
            if (method == 2) {
                // since there are only one  map, we simply generate the same binary image for each segment, a variation (multiple segments) is
                //   not corporated with the code here. we use a window size. It is implemented in a separate function.
                if (qGray(corr_map->pixel(y,x)) < corr_threshold) {
                    (*img_out).setPixel(y,x,0);
                    black += 1;
                }
                else {
                    (*img_out).setPixel(y,x,255);
                    white+= 1;
                }
            }
        }
    }
    imgs_out->push_back(img_out); // push the bin diff image
    this->bin_diff_images->push_back(img_out);
    //  }

    //remember to free space using delete later
    printf("number of foreground pixels: %d, number of background pixels: %d\n", white, black);
    return 0;
}


int dr_group_image_buffer::clean_multiple_background_online(QVector<QImage *> * bin_diff_imgs, int number_of_frames, QString * save_dir) {
    // now we are dealing with no windowing version of this since windowing
    //  requires... something rather than 3 passes and a memory space for a few images .
    //  int n = images.size() / NUMBER_OF_SEGS;
    //  int last_n = images.size() - n * (NUMBER_OF_SEGS - 1);
    QImage * bin_diff_image = bin_diff_imgs->at(0);

    int i=0,x=0,y=0;
    int w,h;
    w = bin_diff_image->width();
    h = bin_diff_image->height();
    // double time_series[n];
    int mean;
    QImage * mean_img = this->mean_images->at(0);
    QImage * new_img;
    //for each background pixel, replace with its mean
    // do it NUMBER_OF_SEGS - 1 times
    for(i = 0; i < number_of_frames; i++) { // number of frames
        // read one image
        new_img = new QImage(*(new QString(datalist->at(i).c_str())));
        for (y = 0; y < h; y++) {
            for (x = 0; x < w; x++) {
                if( qGray(bin_diff_image->pixel(x,y)) == 0 ) {
                    mean = qGray(mean_img->pixel(x,y));
                    new_img->setPixel(x,y,mean);
                }
            }
        }
        // write one image
        QString * new_file_name = new QString(*(save_dir) + "compressed_");
        QString * foo = new QString("foo");
        QString s = foo->sprintf("%04d",i);
        new_file_name->append(s);
        new_file_name->append(".pgm");
        new_img->save(*new_file_name, 0 ,100);
    }
    return 0;
}


void dr_group_image_buffer::set_bin_diff_images_into_video_display_online(QVector<QImage*> * bin_diff_imgs) {
    // this function does two things:
    //   1. set the bin_diff_imgs into video_display class object
    //   2. set a list of frames with background pixels being replaced by 0 to a new list, for further storage
    // Point 2 is not memory friendly, so we comment it out.
    if (v == NULL) {
        printf("error: video display not initialized.\n");
    }
    else {
        v->set_bin_diff_imgs(bin_diff_imgs);
    }
    return;
}


void dr_group_image_buffer::set_data_list(QVector<string> * datalist) {
    this->datalist = datalist;
}


void dr_group_image_buffer::set_corr_img(QImage * corr_img) {
    this->corr_image = corr_img;
}


void dr_group_image_buffer::scale_image(QVector<QVector<double> *> * img, QImage * new_img) {
    int height = img->size();
    int width = img->at(0)->size();
    double maximum_intensity = 0.0;
    int x,y;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            if (img->at(y)->at(x) > maximum_intensity) {
                maximum_intensity = img->at(y)->at(x);
            }
        }
    }
   // cout << "maximum intensity: " << maximum_intensity << endl;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            int new_pixel_value = static_cast<int>(img->at(y)->at(x) * 255.0 / maximum_intensity);
            new_img->setPixel(x,y,new_pixel_value);
        }
    }
}


void dr_group_image_buffer::scale_image2(QImage * img, QImage * new_img) {
    int height = img->height();
    int width = img->width();
    int maximum_intensity = 0;
    int x,y;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            if (qGray(img->pixel(x,y)) > maximum_intensity) {
                maximum_intensity = qGray(img->pixel(x,y));
            }
        }
    }
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            int new_pixel_value = qGray(img->pixel(x,y)) * 255 / maximum_intensity;
            new_img->setPixel(x,y,new_pixel_value);
        }
    }
}


void dr_group_image_buffer::save_compressed_image_in_queue(QImage * out_img, int erosion_size, int dilation_size, int index, QString * debugging_image_dir) {
    // two_levels an image based oh the correlation threshold
    /// TODO: two choices:
    /// (1)   set a hard coded threshold: given in parameter
    /// (2)   or we can do Otsu's thresholding: ...

    /// currently we choose the first one.

    // first make a bin_diff_image
    QImage * bin_diff_image = new QImage(*(out_img));
    QImage * corr_image = new QImage(*(out_img));
    int i,j,h,w;
    h = bin_diff_image->height();
    w = bin_diff_image->width();
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            // Notice: the correlation value in double is before dividing by number of frames and times by 255.
            // So we do it now here
            if (curr_correlation_double_in_window->at(j)->at(i) * 255.0  >= HARD_CODED_THRESHOLD) {
                bin_diff_image->setPixel(i,j,255);
            }
            else {
                bin_diff_image->setPixel(i,j,0);
            }
            corr_image->setPixel(i,j,static_cast<int>(curr_correlation_double_in_window->at(j)->at(i) * 255.0));
        }
    }
    int n_of_background = 0;

    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            if (qGray(bin_diff_image->pixel(i,j)) == 0) {
                n_of_background++;
            }

        }
    }
    cout << "number of foreground pixels: " << h * w - n_of_background << endl;
    cout << "number of background pixels: " << n_of_background << endl;

    // then we "open" the binary image
    // this funtion require the image in qvector,
    // XXX memory glitch
    QVector<QImage*> * binary_image_vec = new QVector<QImage*>();
    binary_image_vec->append(bin_diff_image);

    // XXX memory glitch
    // ** Apply erode and dilate to the bin image
    this->open_bin_image(binary_image_vec, erosion_size, dilation_size, new QString("no"));
    // then (hopefully) the QImage which has its pointer stored in the vector will be passed into the function and the
    // image will eventually be "opened"
    QString * foo = new QString("foo");
    QString s = foo->sprintf("%04d",index);

    // also save the floating point version of the compressed image.
    QString * frame_f_file_name = new QString(*debugging_image_dir);
    frame_f_file_name->append("/frame_f_");
    frame_f_file_name->append(s);
    frame_f_file_name->append(".dat");
    QFile file(*frame_f_file_name);
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);

    n_of_background = 0;
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            if (qGray(binary_image_vec->at(0)->pixel(i,j)) == 0) {
                int mean_background_pixel = static_cast<int>(curr_mean_double_in_window->at(j)->at(i));
                out << curr_mean_double_in_window->at(j)->at(i);
                out_img->setPixel(i,j,mean_background_pixel);
                n_of_background++;
            }
            else{
                out << qGray( out_img->pixel(i,j) ); // for the unchanged foreground pixel, just store the integer value.
            }
        }
    }
    cout << "number of foreground pixels: " << h * w - n_of_background << endl;
    cout << "number of background pixels: " << n_of_background << endl;
    QString * bin_map_file_name = new QString(*debugging_image_dir);
    bin_map_file_name->append("/bin_map_");
    bin_map_file_name->append(s);
    bin_map_file_name->append(".pgm");
    binary_image_vec->at(0)->save(*bin_map_file_name,0,100);

    QString * corr_image_file_name = new QString(*debugging_image_dir);
    corr_image_file_name->append("/corr_img_");
    corr_image_file_name->append(s);
    corr_image_file_name->append(".pgm");
    corr_image->save(*corr_image_file_name,0,100);
}


void dr_group_image_buffer::fill_window_buffer(QImage * in_img){
    this->buffer_for_sliding_windowing->append(in_img);
}


void dr_group_image_buffer:: deque_buffer(){
    // remove the oldest image in buffer
    QImage * removed = this->buffer_for_sliding_windowing->at(0);
     this->buffer_for_sliding_windowing->remove(0);
    delete removed;
}


void dr_group_image_buffer::update_window_stats2(int index, QString * debugging_image_dir){
    // init the mean, variance in the buffer
    QImage * sample_image;
    sample_image = buffer_for_sliding_windowing->at(0);
    int h, w;
    h = sample_image->height();
    w = sample_image->width();

    int i, j;
    for (j = 0; j < h; j++) {

        for (i = 0; i < w; i++) {
             curr_intensity_sum_in_window->at(j)->replace(i,0.0);

             curr_sum_for_variance_double_in_window->at(j)->replace(i,0.0);

             curr_correlation_double_in_window->at(j)->replace(i,0.0);
        }
    }
    // sum all intensity values in buffer pixelwise
    int k;
    for (k = 0; k < this->sliding_window_buffer_size; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                curr_intensity_sum_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) + static_cast<double>( qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)) ) );
            }
        }
    }
    // compute the mean...
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            curr_mean_double_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) / static_cast<double>(sliding_window_buffer_size));
        }
    }
    // compute the sum for variance...
    double diff_from_mean = 0.0;
    for(k = 0; k < this->sliding_window_buffer_size; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                diff_from_mean = static_cast<double>( qGray( this->buffer_for_sliding_windowing->at(k)->pixel(i,j) ) - curr_mean_double_in_window->at(j)->at(i));
                curr_sum_for_variance_double_in_window->at(j)->replace(i, curr_sum_for_variance_double_in_window->at(j)->at(i) + diff_from_mean * diff_from_mean);
            }
        }
    }
    // compute the variance...
    QString * foo = new QString("foo");
    QString s = foo->sprintf("%04d",index);

    QString * variance_f_file_name = new QString(*debugging_image_dir);
    variance_f_file_name->append("/variance_f_");
    variance_f_file_name->append(s);
    variance_f_file_name->append(".dat");
    QFile file(*variance_f_file_name);
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);
    double window_size = static_cast<double>(sliding_window_buffer_size);
    for(j = 0; j < h; j++) {
        for(i=0; i<w;i++) {
            double curr_variance_value = curr_sum_for_variance_double_in_window->at(j)->at(i) / window_size;
            curr_variance_double_in_window->at(j)->replace(i, curr_variance_value);
            out << curr_variance_value;
        }
    }

    // compute the correlation
    double curr_pixel;
    double neighbor_pixel;

    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            double corr_sum_1 = 0.0;
            double corr_sum_2 = 0.0;
            double corr_sum_3 = 0.0;
            double corr_sum_4 = 0.0;
            double corr_sum_5 = 0.0;
            double corr_sum_6 = 0.0;
            double corr_sum_7 = 0.0;
            double corr_sum_8 = 0.0;
            for(k = 0; k < this->sliding_window_buffer_size; k++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                curr_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)));
                if (i != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j)));
                    corr_sum_8 = corr_sum_8 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i-1));
                }
                if (j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j-1)));
                    corr_sum_2 = corr_sum_2 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i));
                }
                if (i != w-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j)));
                    corr_sum_4 = corr_sum_4 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i+1));
                }
                if (j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j+1)));
                    corr_sum_6 = corr_sum_6 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j-1)));
                    corr_sum_1 = corr_sum_1 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i-1));
                }
                if (i != 0 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j+1)));
                    corr_sum_3 = corr_sum_3 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i-1));
                }
                if (i != w-1 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j-1)));
                    corr_sum_7 = corr_sum_7 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i+1));
                }
                if (i != w-1 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j+1)));
                    corr_sum_5 = corr_sum_5 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i+1));
                }

            }
           // cout << "curr sums: " << corr_sum_1 << " " <<  corr_sum_2 << " " <<  corr_sum_3 << " " <<  corr_sum_4 << " " <<  corr_sum_5 << " " <<  corr_sum_6 << " " <<  corr_sum_7 << " " <<   corr_sum_8 << endl;
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(0,corr_sum_1);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(1,corr_sum_2);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(2,corr_sum_3);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(3,corr_sum_4);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(4,corr_sum_5);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(5,corr_sum_6);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(6,corr_sum_7);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(7,corr_sum_8);
        }
    }
    double corr1, corr2, corr3, corr4, corr5, corr6, corr7, corr8;
    corr1 = 0.0;
    corr2 = 0.0;
    corr3 = 0.0;
    corr4 = 0.0;
    corr5 = 0.0;
    corr6 = 0.0;
    corr7 = 0.0;
    corr8 = 0.0;
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            // 1 2 3
            // 8 C 4
            // 7 6 5
            double max_corr = 0.0;
            double this_var = curr_variance_double_in_window->at(j)->at(i);
            double other_var = 0.0;
            if (i != 0) {
                other_var = curr_variance_double_in_window->at(j)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr8 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
              //  cout << "bar: " << curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) << endl;
              //  cout << "corr8: " << corr8 << endl;
                if (corr8 > max_corr) {
                    max_corr = corr8;
                }
            }
            if (j != 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr2 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr2 > max_corr) {
                    max_corr = corr2;
                }
            }
            if (i != w-1) {
                other_var = curr_variance_double_in_window->at(j)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr4 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (WINDOW_SIZE* qSqrt(other_var * this_var )));
                if (corr4 > max_corr) {
                    max_corr = corr4;
                }
            }
            if (j != h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr6 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
               if (corr6 > max_corr) {
                        max_corr = corr6;
                    }
            }
            if (i != 0 && j != 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr1 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr1 > max_corr) {
                    max_corr = corr1;
                }
            }
            if (i != 0 && j!= h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr3 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr3 > max_corr) {
                    max_corr = corr3;
                }
            }
            if ( i != w-1 && j!= 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr7 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
              //  cout << "WINDOW SIZE: " << WINDOW_SIZE << endl;
                if (corr7 > max_corr) {
                    max_corr = corr7;
                }
            }
            if (i!= w-1 && j!= h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr5 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr5 > max_corr) {
                    max_corr = corr5;
                }
            }
            this->curr_correlation_double_in_window->at(j)->replace(i, max_corr);
        }
    }
}


void dr_group_image_buffer::update_window_stats(bool has_deque, bool has_enque){

    // if has_deque = true,
    // then currently in the buffer we have buffer_size + 1 images,
    //  update based on the diff between the two imgs,
    //
    // if has_deque = false;
    //  then we update the stats based on the new coming image
    if (has_deque && has_enque) {
        QImage * oldest_image = buffer_for_sliding_windowing->first();
        QImage * newest_image = buffer_for_sliding_windowing->at(buffer_for_sliding_windowing->size()-1);
        int i, j;
        int h, w;
        h = oldest_image->height();
        w = oldest_image->width();

        double new_mean, old_mean, old_var, var_diff, xn, x0, sum_1_to_nM1_xi,n;
 //       double old_var1, old_var2, old_var3, old_var4, old_var5, old_var6, old_var7, old_var8;
        double x1_mean,x2_mean,x3_mean,x4_mean,x5_mean,x6_mean,x7_mean,x8_mean;
 //       double x1_sum,x2_sum,x3_sum,x4_sum,x5_sum,x6_sum,x7_sum,x8_sum;
        double x1,x2,x3,x4,x5,x6,x7,x8;
        double old_sum, new_sum;
        // new mean = old_mean + (new_x - old_x) / n
        // new_variance = old_variance  + diff
        // diff = [ (xn^2 - x0^2) - (2 * sum_(1 to n-1)[xi] * ( new_mean - old_mean ) )  - 2 * ( new_mean * xn - old_mean * x0 ) ] / n  + ( new_mean^2 - old_mean^2 )
        // new_corr = max( new_neighbors )
        //  new_neighbor = new_corr_sum / (new_neighbor_var * new_var)
        //   new_corr_sum = old_corr_sum - old_stuff + new_stuff
        //   old_stuff: calculate before updating mean/var
        //   new_sutff: calculate after updating mean/var

        QVector<QVector<QVector<double> *> *> * diff_for_corr_sum = new QVector<QVector<QVector<double> *> *>();
        // initialize diff_for_corr
        int x;
        for (j = 0; j < h; j++) {
            diff_for_corr_sum->append(new QVector<QVector<double>*>());
            for (i = 0; i < w; i++){
                diff_for_corr_sum->at(j)->append(new QVector<double>());
                for (x = 0; x < 8; x++) {
                    diff_for_corr_sum->at(j)->at(i)->append(0.0);
                }
            }
        }
        // minus the old stuff from the diff for corr
        n = static_cast<double>(WINDOW_SIZE);
        double old_corr_sum = curr_sum_for_correlation_double_in_window->at(200)->at(200)->at(7);
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                old_sum = curr_intensity_sum_in_window->at(j)->at(i);
                x0 = static_cast<double>( qGray(oldest_image->pixel(i,j)) );
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    x8 = static_cast<double>( qGray(oldest_image->pixel(i-1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(7, 0.0 + x8_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j)->at(i-1) - n * old_mean * x8_mean - x8*x0);
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    x2 = static_cast<double>( qGray(oldest_image->pixel(i,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(1, 0.0 + x2_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j-1)->at(i) - n * old_mean * x2_mean - x2*x0);
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    x4 = static_cast<double>( qGray(oldest_image->pixel(i+1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(3, 0.0 + x4_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j)->at(i+1) - n * old_mean * x4_mean - x4*x0);
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    x6 = static_cast<double>( qGray(oldest_image->pixel(i,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(5, 0.0 + x6_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j+1)->at(i) - n * old_mean * x6_mean - x6*x0);
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    x1 = static_cast<double>( qGray(oldest_image->pixel(i-1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(0, 0.0 + x1_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) - n * old_mean * x1_mean - x1*x0);
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    x7 = static_cast<double>( qGray(oldest_image->pixel(i-1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(6, 0.0 + x7_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) - n * old_mean * x7_mean - x7*x0);
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    x3 = static_cast<double>( qGray(oldest_image->pixel(i+1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(2, 0.0 + x3_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) - n * old_mean * x3_mean - x3*x0);
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    x5 = static_cast<double>( qGray(oldest_image->pixel(i+1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(4, 0.0 + x5_mean * old_sum + old_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) - n * old_mean * x5_mean - x5*x0);
                }
            }
        }

        // repalce the old mean and variance by the new ones
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                xn = static_cast<double>( qGray(newest_image->pixel(i,j)) );
                x0 = static_cast<double>( qGray(oldest_image->pixel(i,j)) );
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                double n = static_cast<double>(buffer_for_sliding_windowing->size()-1);
                new_mean = curr_mean_double_in_window->at(j)->at(i) +  (xn-x0) / (buffer_for_sliding_windowing->size()-1);
                curr_mean_double_in_window->at(j)->replace(i, new_mean);
                //cout << "oldmean minus new mean: " << old_mean - new_mean << endl;
                sum_1_to_nM1_xi = curr_intensity_sum_in_window->at(j)->at(i) - x0;
                var_diff = ( (qPow(xn,2) - qPow(x0,2)) - (2.0 * sum_1_to_nM1_xi * (new_mean - old_mean)) - 2.0 * (new_mean * xn - old_mean * x0) ) / n + (qPow(new_mean,2) - qPow(old_mean, 2));
                old_var = curr_variance_double_in_window->at(j)->at(i);
                // here: keep neighbor's old_var
                curr_intensity_sum_in_window->at(j)->replace(i, sum_1_to_nM1_xi+xn);
                curr_variance_double_in_window->at(j)->replace(i, old_var + var_diff);
            }
        }
        double curr_sum_for_correlation;
        // plus the new stuff from the diff for corr
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                double n = static_cast<double>(10);
                new_mean = curr_mean_double_in_window->at(j)->at(i);
                new_sum = curr_intensity_sum_in_window->at(j)->at(i);
                x0 = static_cast<double>( qGray(newest_image->pixel(i,j)) );
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    x8 = static_cast<double>( qGray(newest_image->pixel(i-1,j)) );

                    diff_for_corr_sum->at(j)->at(i)->replace(7, diff_for_corr_sum->at(j)->at(i)->at(7) - x8_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j)->at(i-1) + n * new_mean * x8_mean + x8*x0);
             //       cout << "diff for corr sum: " << diff_for_corr_sum->at(j)->at(i)->at(7) / 100.0 << endl;
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7);
             //       cout << "corr sum: " << curr_sum_for_correlation << endl;
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(7, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(7) / 255.0 );
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    x2 = static_cast<double>( qGray(newest_image->pixel(i,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(1, diff_for_corr_sum->at(j)->at(i)->at(1) - x2_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j-1)->at(i) + n * new_mean * x2_mean + x2*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(1, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(1)  /  255.0);
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    x4 = static_cast<double>( qGray(newest_image->pixel(i+1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(3, diff_for_corr_sum->at(j)->at(i)->at(3) - x4_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j)->at(i+1) + n * new_mean * x4_mean + x4*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(3, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(3)  /  255.0);
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    x6 = static_cast<double>( qGray(newest_image->pixel(i,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(5, diff_for_corr_sum->at(j)->at(i)->at(5) - x6_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j+1)->at(i) + n * new_mean * x6_mean + x6*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(5, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(5)  /  255.0);
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    x1 = static_cast<double>( qGray(newest_image->pixel(i-1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(0, diff_for_corr_sum->at(j)->at(i)->at(0) - x1_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) + n * new_mean * x1_mean + x1*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(0, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(0)  /  255.0);
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    x7 = static_cast<double>( qGray(newest_image->pixel(i-1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(6, diff_for_corr_sum->at(j)->at(i)->at(6) - x7_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) + n * new_mean * x7_mean + x7*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(6, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(6) /  255.0 );
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    x3 = static_cast<double>( qGray(newest_image->pixel(i+1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(2, diff_for_corr_sum->at(j)->at(i)->at(2) - x3_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) + n * new_mean * x3_mean + x3*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(2, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(2)  /  255.0);
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    x5 = static_cast<double>( qGray(newest_image->pixel(i+1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(4, diff_for_corr_sum->at(j)->at(i)->at(4) - x5_mean * new_sum - new_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) + n * new_mean * x5_mean + x5*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(4, curr_sum_for_correlation + diff_for_corr_sum->at(j)->at(i)->at(4) / 255.0 );
                }
            }
        }
        double new_corr_sum = curr_sum_for_correlation_double_in_window->at(200)->at(200)->at(7);
        cout << "old and new: " << old_corr_sum << " " << new_corr_sum << endl;
        double neigh_var = 0.0;
        double this_var = 0.0;
        double x1_corr, x2_corr,x3_corr,x4_corr,x5_corr,x6_corr,x7_corr,x8_corr;
        // update correlation values
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                double max_corr = 0.0;
                // 1 2 3
                // 8 C 4
                // 7 6 5
                //
                // 8 2 4 6 1 7 3 5
                this_var = curr_variance_double_in_window->at(j)->at(i);
                if (this_var < 0.001) {
                    this_var = 1.0;
                }

                if (i != 0) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x8_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (10* sqrt(neigh_var * this_var)));
                    if (x8_corr > max_corr) {
                        max_corr = x8_corr;
                    }
                }
                if (j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x2_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (10* sqrt(neigh_var * this_var)));
                    if (x2_corr > max_corr) {
                        max_corr = x2_corr;
                    }
                }
                if (i != w-1) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x4_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (10* sqrt(neigh_var * this_var)));
                    if (x4_corr > max_corr) {
                        max_corr = x4_corr;
                    }
                }
                if (j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x6_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (10* sqrt(neigh_var * this_var)));
                    if (x6_corr > max_corr) {
                        max_corr = x6_corr;
                    }
                }
                if (i != 0 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x1_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (10* sqrt(neigh_var * this_var)));
                    if (x1_corr > max_corr) {
                        max_corr = x1_corr;
                    }
                }
                if (i != 0 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x7_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (10* sqrt(neigh_var * this_var)));
                    if (x7_corr > max_corr) {
                        max_corr = x7_corr;
                    }
                }
                if (i != w-1 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x3_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (10* sqrt(neigh_var * this_var)));
                    if (x3_corr > max_corr) {
                        max_corr = x3_corr;
                    }
                }
                if (i != w-1 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x5_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (10* sqrt(neigh_var * this_var)));
                    if (x5_corr > max_corr) {
                        max_corr = x5_corr;
                    }
                }

          //      cout << "old corr: " << curr_correlation_double_in_window->at(j)->at(i) << endl;
                curr_correlation_double_in_window->at(j)->replace(i,max_corr);
          //      cout << "new corr: " << curr_correlation_double_in_window->at(j)->at(i) << endl;
            }
        }
        for (int foo = 0; foo < diff_for_corr_sum->size(); foo++) {
            QVector<QVector<double> *> * currxx = diff_for_corr_sum->at(foo);
            for (int bar = 0; bar < currxx->size(); bar ++) {
                QVector<double> * curroo = currxx->at(bar);
                delete curroo;
            }
            delete currxx;
        }
        delete diff_for_corr_sum;
    }
    else if (has_enque) {
        // if we don't deque:
        // then we don't  need the oldest image information
        QImage * newest_image = buffer_for_sliding_windowing->last();
        int i, j;
        int h, w;
        h = newest_image->height();
        w = newest_image->width();

        double new_mean, old_mean, old_var, xn, x0, sum_0_to_nM1_xi,n;
 //       double old_var1, old_var2, old_var3, old_var4, old_var5, old_var6, old_var7, old_var8;
        double x1_mean,x2_mean,x3_mean,x4_mean,x5_mean,x6_mean,x7_mean,x8_mean;
 //       double x1_sum,x2_sum,x3_sum,x4_sum,x5_sum,x6_sum,x7_sum,x8_sum;
        double x1,x2,x3,x4,x5,x6,x7,x8;
        double old_sum, new_sum;
        // new mean = old_mean + (new_x - old_x) / n

        // new_variance = old_variance  + diff
        // diff = [ (xn^2 - x0^2) - (2 * sum_(1 to n-1)[xi] * ( new_mean - old_mean ) )  - 2 * ( new_mean * xn - old_mean * x0 ) ] / n  + ( new_mean^2 - old_mean^2 )

        // new_corr = max( new_neighbors )
        //  new_neighbor = new_corr_sum / (new_neighbor_var * new_var)
        //   new_corr_sum = old_corr_sum - old_stuff + new_stuff
        //   old_stuff: calculate before updating mean/var
        //   new_sutff: calculate after updating mean/var

        QVector<QVector<QVector<double> *> *> * diff_for_corr_sum = new QVector<QVector<QVector<double> *> *>();
        // initialize diff_for_corr
        int x;
        for (j = 0; j < h; j++) {
            diff_for_corr_sum->append(new QVector<QVector<double>*>());
            for (i = 0; i < w; i++){
                diff_for_corr_sum->at(j)->append(new QVector<double>());
                for (x = 0; x < 8; x++) {
                    diff_for_corr_sum->at(j)->at(i)->append(0.0);
                }
            }
        }
        // minus the old stuff from the diff for corr
        n = static_cast<double>(buffer_for_sliding_windowing->size()-1);

        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                old_sum = curr_intensity_sum_in_window->at(j)->at(i);
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(7, 0.0 - x8_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j)->at(i-1) - n * old_mean * x8_mean);
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    diff_for_corr_sum->at(j)->at(i)->replace(1, 0.0 - x2_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i) - n * old_mean * x2_mean);
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(3, 0.0 - x4_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j)->at(i+1) - n * old_mean * x4_mean);
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    diff_for_corr_sum->at(j)->at(i)->replace(5, 0.0 - x6_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i) - n * old_mean * x6_mean);
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(0, 0.0 - x1_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) - n * old_mean * x1_mean);
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(6, 0.0 - x7_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) - n * old_mean * x7_mean);
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(2, 0.0 - x3_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) - n * old_mean * x3_mean);
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(4, 0.0 - x5_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) - n * old_mean * x5_mean);
                }
            }
        }
        n = static_cast<double>(buffer_for_sliding_windowing->size()-1);
        double new_var;
        double var_diff_on_numerator;
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
             //   int hh = newest_image->height();
             //   cout << "================\ntesting: \n this image in the window buffer size: \n===========\n" << hh << endl;
                // newest image is OK
                xn = static_cast<double>( qGray(newest_image->pixel(i,j)) );
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                new_mean = (old_mean * n + xn) / (n+1);
                curr_mean_double_in_window->at(j)->replace(i, new_mean);

                sum_0_to_nM1_xi = curr_intensity_sum_in_window->at(j)->at(i);
                var_diff_on_numerator = qPow(xn,2) - (2.0 * sum_0_to_nM1_xi * (new_mean - old_mean)) - 2.0 * new_mean * xn + (n+1.0) * qPow(new_mean,2) - n * qPow(old_mean, 2);
                old_var = curr_variance_double_in_window->at(j)->at(i);

                new_var = ( (old_var * n) + var_diff_on_numerator ) / (n+1.0);
                curr_variance_double_in_window->at(j)->replace(i, new_var);
            }
        }

        double curr_sum_for_correlation;
        // plus the new stuff from the diff for corr
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                new_mean = curr_mean_double_in_window->at(j)->at(i);
                new_sum = curr_intensity_sum_in_window->at(j)->at(i);
                x0 = static_cast<double>( qGray(newest_image->pixel(i,j)) );
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    x8 = static_cast<double>( qGray(newest_image->pixel(i-1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(7, diff_for_corr_sum->at(j)->at(i)->at(7) + x8_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j)->at(i-1) + new_mean * x8_mean + x8*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(7, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(7) );
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    x2 = static_cast<double>( qGray(newest_image->pixel(i,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(1, diff_for_corr_sum->at(j)->at(i)->at(1) + x2_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i) + n * new_mean * x2_mean + x2*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(1, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(1) );
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    x4 = static_cast<double>( qGray(newest_image->pixel(i+1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(3, diff_for_corr_sum->at(j)->at(i)->at(3) + x4_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j)->at(i+1) + n * new_mean * x4_mean + x4*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(3, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(3) );
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    x6 = static_cast<double>( qGray(newest_image->pixel(i,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(5, diff_for_corr_sum->at(j)->at(i)->at(5) + x6_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i) + n * new_mean * x6_mean + x6*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(5, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(5) );
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    x1 = static_cast<double>( qGray(newest_image->pixel(i-1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(0, diff_for_corr_sum->at(j)->at(i)->at(0) + x1_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) + n * new_mean * x1_mean + x1*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(0, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(0) );
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    x7 = static_cast<double>( qGray(newest_image->pixel(i-1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(6, diff_for_corr_sum->at(j)->at(i)->at(6) + x7_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) + n * new_mean * x7_mean + x7*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(6, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(6) );
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    x3 = static_cast<double>( qGray(newest_image->pixel(i+1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(2, diff_for_corr_sum->at(j)->at(i)->at(2) + x3_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) + n * new_mean * x3_mean + x3*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(2, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(2) );
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    x5 = static_cast<double>( qGray(newest_image->pixel(i+1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(4, diff_for_corr_sum->at(j)->at(i)->at(4) + x5_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) + n * new_mean * x5_mean + x5*x0);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(4, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(4) );
                }
            }
        }

        double max_corr = 0.0;
        double neigh_var = 0.0;
        double this_var = 0.0;
        double x1_corr, x2_corr,x3_corr,x4_corr,x5_corr,x6_corr,x7_corr,x8_corr;

        // update correlation values
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                //
                // 8 2 4 6 1 7 3 5
                this_var = curr_variance_double_in_window->at(j)->at(i);
                if (this_var < 0.001) {
                    this_var = 1.0;
                }

                if (i != 0) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x8_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (neigh_var * this_var));
                    if (x8_corr > max_corr) {
                        max_corr = x8_corr;
                    }
                }
                if (j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x2_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (neigh_var * this_var));
                    if (x2_corr > max_corr) {
                        max_corr = x2_corr;
                    }
                }
                if (i != w-1) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x4_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (neigh_var * this_var));
                    if (x4_corr > max_corr) {
                        max_corr = x4_corr;
                    }
                }
                if (j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x6_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (neigh_var * this_var));
                    if (x6_corr > max_corr) {
                        max_corr = x6_corr;
                    }
                }
                if (i != 0 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x1_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (neigh_var * this_var));
                    if (x1_corr > max_corr) {
                        max_corr = x1_corr;
                    }
                }
                if (i != 0 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x7_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (neigh_var * this_var));
                    if (x7_corr > max_corr) {
                        max_corr = x7_corr;
                    }
                }
                if (i != w-1 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x3_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (neigh_var * this_var));
                    if (x3_corr > max_corr) {
                        max_corr = x3_corr;
                    }
                }
                if (i != w-1 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x5_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (neigh_var * this_var));
                    if (x5_corr > max_corr) {
                        max_corr = x5_corr;
                    }
                }
                cout << "old corr: " << curr_sum_for_correlation_double_in_window->at(j)->at(i);
                curr_correlation_double_in_window->at(j)->replace(i,max_corr);
                cout << "new corr: " << curr_sum_for_correlation_double_in_window->at(j)->at(i);

            }
        }
    }
    else if (has_deque) { // if we only deque
        QImage * oldest_image = buffer_for_sliding_windowing->first();
        int i, j;
        int h, w;
        h = oldest_image->height();
        w = oldest_image->width();

        double new_mean, old_mean, var_diff, x0, sum_1_to_nM1_xi,n;
        double x1_mean,x2_mean,x3_mean,x4_mean,x5_mean,x6_mean,x7_mean,x8_mean;
        double x1,x2,x3,x4,x5,x6,x7,x8;
        double old_sum, new_sum;
        // before doing this loop to update mean and variance, firstly make use all the old_mean, old_sum vals to compute what should be subtracted in making the corr diff
        QVector<QVector<QVector<double> *> *> * diff_for_corr_sum = new QVector<QVector<QVector<double> *> *>();
        // initialize diff_for_corr
        int x;
        for (j = 0; j < h; j++) {
            diff_for_corr_sum->append(new QVector<QVector<double>*>());
            for (i = 0; i < w; i++){
                diff_for_corr_sum->at(j)->append(new QVector<double>());
                for (x = 0; x < 8; x++) {
                    diff_for_corr_sum->at(j)->at(i)->append(0.0);
                }
            }
        }
        // minus the old stuff from the diff for corr
        n = static_cast<double>(buffer_for_sliding_windowing->size()-1);

        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                old_sum = curr_intensity_sum_in_window->at(j)->at(i);
                x0 = static_cast<double>( qGray(oldest_image->pixel(i,j)) );
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    x8 = static_cast<double>( qGray(oldest_image->pixel(i-1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(7, 0.0 - x8_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j)->at(i-1) - n * old_mean * x8_mean - x8*x0);
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    x2 = static_cast<double>( qGray(oldest_image->pixel(i,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(1, 0.0 - x2_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i) - n * old_mean * x2_mean - x2*x0);
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    x4 = static_cast<double>( qGray(oldest_image->pixel(i+1,j)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(3, 0.0 - x4_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j)->at(i+1) - n * old_mean * x4_mean - x4*x0);
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    x6 = static_cast<double>( qGray(oldest_image->pixel(i,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(5, 0.0 - x6_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i) - n * old_mean * x6_mean - x6*x0);
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    x1 = static_cast<double>( qGray(oldest_image->pixel(i-1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(0, 0.0 - x1_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) - n * old_mean * x1_mean - x1*x0);
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    x7 = static_cast<double>( qGray(oldest_image->pixel(i-1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(6, 0.0 - x7_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) - n * old_mean * x7_mean - x7*x0);
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    x3 = static_cast<double>( qGray(oldest_image->pixel(i+1,j-1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(2, 0.0 - x3_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) - n * old_mean * x3_mean - x3*x0);
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    x5 = static_cast<double>( qGray(oldest_image->pixel(i+1,j+1)) );
                    diff_for_corr_sum->at(j)->at(i)->replace(4, 0.0 - x5_mean * old_sum - old_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) - n * old_mean * x5_mean - x5*x0);
                }
            }
        }
        n = static_cast<double>(buffer_for_sliding_windowing->size()-1);
        /// TODO: update variance sum here later so we can make use of it (? may not be necessary)
        double old_var_numerator = 0.0;
        double new_var = 0.0;
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                x0 = static_cast<double>( qGray(oldest_image->pixel(i,j)) );
                old_mean = curr_mean_double_in_window->at(j)->at(i);
                new_mean = (old_mean * n - x0) / (n-1.0);
                curr_mean_double_in_window->at(j)->replace(i, new_mean);

                sum_1_to_nM1_xi = curr_intensity_sum_in_window->at(j)->at(i) - x0;
                var_diff = ( 0.0 - qPow(x0,2) - (2.0 * sum_1_to_nM1_xi * (new_mean - old_mean)) + ( 2.0 * x0 * old_mean) + (n-1) * qPow(new_mean,2) - n * qPow(old_mean, 2) );
                old_var_numerator = (curr_variance_double_in_window->at(j)->at(i)) * n;
                // here: keep neighbor's old_var
                new_var = (old_var_numerator + var_diff) / (n-1);
                curr_variance_double_in_window->at(j)->replace(i, new_var);
            }
        }
        double curr_sum_for_correlation;
        // plus the new stuff from the diff for corr
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                new_mean = curr_mean_double_in_window->at(j)->at(i);
                new_sum = curr_intensity_sum_in_window->at(j)->at(i);
                if (i != 0) {
                    x8_mean = curr_mean_double_in_window->at(j)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(7, diff_for_corr_sum->at(j)->at(i)->at(7) + x8_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j)->at(i-1) + (n-1) * new_mean * x8_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(7, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(7) );
                }
                if (j != 0) {
                    x2_mean = curr_mean_double_in_window->at(j-1)->at(i);
                    diff_for_corr_sum->at(j)->at(i)->replace(1, diff_for_corr_sum->at(j)->at(i)->at(1) + x2_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i) + n * new_mean * x2_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(1, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(1) );
                }
                if (i != w-1) {
                    x4_mean = curr_mean_double_in_window->at(j)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(3, diff_for_corr_sum->at(j)->at(i)->at(3) + x4_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j)->at(i+1) + n * new_mean * x4_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(3, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(3) );
                }
                if (j != h-1) {
                    x6_mean = curr_mean_double_in_window->at(j+1)->at(i);
                    diff_for_corr_sum->at(j)->at(i)->replace(5, diff_for_corr_sum->at(j)->at(i)->at(5) + x6_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i) + n * new_mean * x6_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(5, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(5) );
                }
                if (i != 0 && j != 0) {
                    x1_mean = curr_mean_double_in_window->at(j-1)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(0, diff_for_corr_sum->at(j)->at(i)->at(0) + x1_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i-1) + n * new_mean * x1_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(0, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(0) );
                }
                if (i != 0 && j != h-1) {
                    x7_mean = curr_mean_double_in_window->at(j+1)->at(i-1);
                    diff_for_corr_sum->at(j)->at(i)->replace(6, diff_for_corr_sum->at(j)->at(i)->at(6) + x7_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i-1) + n * new_mean * x7_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(6, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(6) );
                }
                if (i != w-1 && j != 0) {
                    x3_mean = curr_mean_double_in_window->at(j-1)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(2, diff_for_corr_sum->at(j)->at(i)->at(2) + x3_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j-1)->at(i+1) + n * new_mean * x3_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(2, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(2) );
                }
                if (i != w-1 && j != h-1) {
                    x5_mean = curr_mean_double_in_window->at(j+1)->at(i+1);
                    diff_for_corr_sum->at(j)->at(i)->replace(4, diff_for_corr_sum->at(j)->at(i)->at(4) + x5_mean * new_sum + new_mean * curr_intensity_sum_in_window->at(j+1)->at(i+1) + n * new_mean * x5_mean);
                    curr_sum_for_correlation = curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4);
                    curr_sum_for_correlation_double_in_window->at(j)->at(i)->replace(4, curr_sum_for_correlation - diff_for_corr_sum->at(j)->at(i)->at(4) );
                }
            }
        }

        double max_corr = 0.0;
        double neigh_var = 0.0;
        double this_var = 0.0;
        double x1_corr, x2_corr,x3_corr,x4_corr,x5_corr,x6_corr,x7_corr,x8_corr;
        // update correlation values
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                //
                // 8 2 4 6 1 7 3 5
                this_var = curr_variance_double_in_window->at(j)->at(i);
                if (this_var < 0.001) {
                    this_var = 1.0;
                }

                if (i != 0) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x8_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (neigh_var * this_var));
                    if (x8_corr > max_corr) {
                        max_corr = x8_corr;
                    }
                }
                if (j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x2_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (neigh_var * this_var));
                    if (x2_corr > max_corr) {
                        max_corr = x2_corr;
                    }
                }
                if (i != w-1) {
                    neigh_var = curr_variance_double_in_window->at(j)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x4_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (neigh_var * this_var));
                    if (x4_corr > max_corr) {
                        max_corr = x4_corr;
                    }
                }
                if (j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x6_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (neigh_var * this_var));
                    if (x6_corr > max_corr) {
                        max_corr = x6_corr;
                    }
                }
                if (i != 0 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x1_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (neigh_var * this_var));
                    if (x1_corr > max_corr) {
                        max_corr = x1_corr;
                    }
                }
                if (i != 0 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i-1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x7_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (neigh_var * this_var));
                    if (x7_corr > max_corr) {
                        max_corr = x7_corr;
                    }
                }
                if (i != w-1 && j != 0) {
                    neigh_var = curr_variance_double_in_window->at(j-1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x3_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (neigh_var * this_var));
                    if (x3_corr > max_corr) {
                        max_corr = x3_corr;
                    }
                }
                if (i != w-1 && j != h-1) {
                    neigh_var = curr_variance_double_in_window->at(j+1)->at(i+1);
                    if (neigh_var < 0.001) {
                        neigh_var = 1.0;
                    }
                    x5_corr = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (neigh_var * this_var));
                    if (x5_corr > max_corr) {
                        max_corr = x5_corr;
                    }
                }
            }
        }
    }
}


void dr_group_image_buffer::init_window_stats() {
    // init the mean, variance in the buffer
    QImage * sample_image;
    sample_image = buffer_for_sliding_windowing->at(0);
    curr_intensity_sum_in_window = new QVector<QVector<double> *>();
    curr_mean_double_in_window = new QVector<QVector<double> *>();
    curr_variance_double_in_window = new QVector<QVector<double> *>();
    curr_sum_for_variance_double_in_window = new QVector<QVector<double> *>();

    curr_sum_for_correlation_double_in_window = new QVector<QVector<QVector<double> *> *>();
    curr_correlation_double_in_window = new QVector<QVector<double> *>();

    int h, w;
    h = sample_image->height();
    w = sample_image->width();

    int i, j;
    for (j = 0; j < h; j++) {
        curr_intensity_sum_in_window->append( new QVector<double>() );
        curr_mean_double_in_window->append( new QVector<double>() );
        curr_variance_double_in_window->append(new QVector<double>());
        curr_sum_for_variance_double_in_window->append(new QVector<double>());

        curr_sum_for_correlation_double_in_window->append(new QVector<QVector<double> *>());
        curr_correlation_double_in_window->append(new QVector<double>());

        for (i = 0; i < w; i++) {
             curr_intensity_sum_in_window->at(j)->append(0.0);
             curr_mean_double_in_window->at(j)->append(0.0);
             curr_variance_double_in_window->at(j)->append(0.0);

             curr_sum_for_variance_double_in_window->at(j)->append(0.0);

             curr_sum_for_correlation_double_in_window->at(j)->append(new QVector<double>());
             curr_correlation_double_in_window->at(j)->append(0.0);
        }
    }
    // sum all intensity values in buffer pixelwise
    // ** Sum up the pixels from the first half of input frames for each pixel.
    // ** Save the result to curr_intensity_sum_in_window
    int k;
    for (k = 0; k < this->sliding_window_buffer_size / 2; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                curr_intensity_sum_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) + static_cast<double>( qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)) ) );
            }
        }
    }
    /// TODO: this can also be replaced by the online version, even though I am not sure it is necessary to do so.
    // compute the mean...
    // ** Devide each pixel sum in 'curr_mean_double_in_window' by half window size to get the mean
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            curr_mean_double_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) / static_cast<double>(sliding_window_buffer_size / 2));
        }
    }
    // compute the sum for variance...
    double diff_from_mean = 0.0;
    for(k = 0; k < this->sliding_window_buffer_size / 2; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                diff_from_mean = static_cast<double>( qGray( this->buffer_for_sliding_windowing->at(k)->pixel(i,j) ) - curr_mean_double_in_window->at(j)->at(i));
                curr_sum_for_variance_double_in_window->at(j)->replace(i, curr_sum_for_variance_double_in_window->at(j)->at(i) + diff_from_mean * diff_from_mean);
            }
        }
    }
    // compute the variance...

    for(j = 0; j < h; j++) {
        for(i=0; i<w;i++) {
            curr_variance_double_in_window->at(j)->replace(i, curr_sum_for_variance_double_in_window->at(j)->at(i) / static_cast<double>(sliding_window_buffer_size / 2));
        }
    }
    // compute the correlation
    double curr_pixel;
    double neighbor_pixel;

    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            double corr_sum_1 = 0.0;
            double corr_sum_2 = 0.0;
            double corr_sum_3 = 0.0;
            double corr_sum_4 = 0.0;
            double corr_sum_5 = 0.0;
            double corr_sum_6 = 0.0;
            double corr_sum_7 = 0.0;
            double corr_sum_8 = 0.0;
            for(k = 0; k < this->sliding_window_buffer_size / 2; k++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                curr_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)));
                if (i != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j)));
                    corr_sum_8 = corr_sum_8 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i-1));
                }
                if (j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j-1)));
                    corr_sum_2 = corr_sum_2 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i));
                }
                if (i != w-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j)));
                    corr_sum_4 = corr_sum_4 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i+1));
                }
                if (j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j+1)));
                    corr_sum_6 = corr_sum_6 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j-1)));
                    corr_sum_1 = corr_sum_1 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i-1));
                }
                if (i != 0 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j+1)));
                    corr_sum_3 = corr_sum_3 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i-1));
                }
                if (i != w-1 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j-1)));
                    corr_sum_7 = corr_sum_7 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i+1));
                }
                if (i != w-1 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j+1)));
                    corr_sum_5 = corr_sum_5 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i+1));
                }

            }
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_1);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_2);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_3);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_4);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_5);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_6);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_7);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_8);
        }
    }
    double corr1, corr2, corr3, corr4, corr5, corr6, corr7, corr8;
    corr1 = 0.0;
    corr2 = 0.0;
    corr3 = 0.0;
    corr4 = 0.0;
    corr5 = 0.0;
    corr6 = 0.0;
    corr7 = 0.0;
    corr8 = 0.0;
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            // 1 2 3
            // 8 C 4
            // 7 6 5
            double max_corr = 0.0;
            double this_var = curr_variance_double_in_window->at(j)->at(i);
            double other_var = 0.0;
            if (i != 0) {
                other_var = curr_variance_double_in_window->at(j)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr8 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (10* qSqrt(other_var * this_var)));
                if (corr8 > max_corr) {
                    max_corr = corr8;
                }
            }
            if (j != 0) {
                other_var = curr_sum_for_variance_double_in_window->at(j-1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr2 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (10* qSqrt(other_var * this_var)));
                if (corr2 > max_corr) {
                    max_corr = corr2;
                }
            }
            if (i != w-1) {
                other_var = curr_sum_for_variance_double_in_window->at(j)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr4 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (10* qSqrt(other_var * this_var )));
                if (corr4 > max_corr) {
                    max_corr = corr4;
                }
            }
            if (j != h-1) {
                other_var = curr_sum_for_variance_double_in_window->at(j+1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr6 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (10* qSqrt(other_var * this_var)));
               if (corr6 > max_corr) {
                        max_corr = corr6;
                    }
            }
            if (i != 0 && j != 0) {
                other_var = curr_sum_for_variance_double_in_window->at(j-1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr1 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (10* qSqrt(other_var * this_var)));
                if (corr1 > max_corr) {
                    max_corr = corr1;
                }
            }
            if (i != 0 && j!= h-1) {
                other_var = curr_sum_for_variance_double_in_window->at(j+1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr3 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (10* qSqrt(other_var * this_var)));
                if (corr3 > max_corr) {
                    max_corr = corr3;
                }
            }
            if ( i != w-1 && j!= 0) {
                other_var = curr_sum_for_variance_double_in_window->at(j-1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr7 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (10* qSqrt(other_var * this_var)));
                if (corr7 > max_corr) {
                    max_corr = corr7;
                }
            }
            if (i!= w-1 && j!= h-1) {
                other_var = curr_sum_for_variance_double_in_window->at(j+1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr5 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (10* qSqrt(other_var * this_var)));
                if (corr5 > max_corr) {
                    max_corr = corr5;
                }
            }
            this->curr_correlation_double_in_window->at(j)->replace(i, max_corr);
        }
    }
}


// ** Compute a correlation image from the images from the first half of a window
void dr_group_image_buffer::init_window_stats2( QString * debugging_image_dir) {
    // init the mean, variance in the buffer
    QImage * sample_image;
    sample_image = buffer_for_sliding_windowing->at(0);
    curr_intensity_sum_in_window = new QVector<QVector<double> *>();
    curr_mean_double_in_window = new QVector<QVector<double> *>();
    curr_variance_double_in_window = new QVector<QVector<double> *>();
    curr_sum_for_variance_double_in_window = new QVector<QVector<double> *>();

    curr_sum_for_correlation_double_in_window = new QVector<QVector<QVector<double> *> *>();
    curr_correlation_double_in_window = new QVector<QVector<double> *>();

    int h, w;
    h = sample_image->height();
    w = sample_image->width();

    int i, j;
    for (j = 0; j < h; j++) {
        curr_intensity_sum_in_window->append( new QVector<double>() );
        curr_mean_double_in_window->append( new QVector<double>() );
        curr_variance_double_in_window->append(new QVector<double>());
        curr_sum_for_variance_double_in_window->append(new QVector<double>());

        curr_sum_for_correlation_double_in_window->append(new QVector<QVector<double> *>());
        curr_correlation_double_in_window->append(new QVector<double>());

        for (i = 0; i < w; i++) {
             curr_intensity_sum_in_window->at(j)->append(0.0);
             curr_mean_double_in_window->at(j)->append(0.0);
             curr_variance_double_in_window->at(j)->append(0.0);

             curr_sum_for_variance_double_in_window->at(j)->append(0.0);

             curr_sum_for_correlation_double_in_window->at(j)->append(new QVector<double>());

             curr_correlation_double_in_window->at(j)->append(0.0);
        }
    }
    // sum all intensity values in buffer pixelwise
    int k;
    for (k = 0; k < this->sliding_window_buffer_size; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                curr_intensity_sum_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) + static_cast<double>( qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)) ) );
            }
        }
    }
    // compute the mean...
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            curr_mean_double_in_window->at(j)->replace(i, curr_intensity_sum_in_window->at(j)->at(i) / static_cast<double>(sliding_window_buffer_size));
        }
    }
    // compute the sum for variance...
    double diff_from_mean = 0.0;
    for(k = 0; k < this->sliding_window_buffer_size; k++) {
        for (j = 0; j < h; j++) {
            for (i = 0; i < w; i++) {
                diff_from_mean = static_cast<double>( qGray( this->buffer_for_sliding_windowing->at(k)->pixel(i,j) ) - curr_mean_double_in_window->at(j)->at(i));
                curr_sum_for_variance_double_in_window->at(j)->replace(i, curr_sum_for_variance_double_in_window->at(j)->at(i) + diff_from_mean * diff_from_mean);
            }
        }
    }
    // compute the variance...
    QString * foo = new QString("foo");
    QString s = foo->sprintf("%04d",-1);
    QString * variance_f_file_name = new QString(*debugging_image_dir);
    variance_f_file_name->append("/variance_f_");
    variance_f_file_name->append(s);
    variance_f_file_name->append(".dat");
    QFile file(*variance_f_file_name);
    file.open(QIODevice::WriteOnly);
    QDataStream out(&file);
    double window_size_in_double = static_cast<double>(sliding_window_buffer_size);
    for(j = 0; j < h; j++) {
        for(i = 0; i < w; i++) {
            double curr_variance_value = curr_sum_for_variance_double_in_window->at(j)->at(i) / window_size_in_double;
            curr_variance_double_in_window->at(j)->replace(i, curr_variance_value);
            out << curr_variance_value;
        }
    }
    // compute the correlation
    double curr_pixel;
    double neighbor_pixel;

    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            double corr_sum_1 = 0.0;
            double corr_sum_2 = 0.0;
            double corr_sum_3 = 0.0;
            double corr_sum_4 = 0.0;
            double corr_sum_5 = 0.0;
            double corr_sum_6 = 0.0;
            double corr_sum_7 = 0.0;
            double corr_sum_8 = 0.0;
            for(k = 0; k < this->sliding_window_buffer_size; k++) {
                // 1 2 3
                // 8 C 4
                // 7 6 5
                curr_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j)));
                if (i != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j)));
                    corr_sum_8 = corr_sum_8 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i-1));
                }
                if (j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j-1)));
                    corr_sum_2 = corr_sum_2 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i));
                }
                if (i != w-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j)));
                    corr_sum_4 = corr_sum_4 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j)->at(i+1));
                }
                if (j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i,j+1)));
                    corr_sum_6 = corr_sum_6 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i));
                }
                if (i != 0 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j-1)));
                    corr_sum_1 = corr_sum_1 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i-1));
                }
                if (i != 0 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i-1,j+1)));
                    corr_sum_3 = corr_sum_3 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i-1));
                }
                if (i != w-1 && j != 0) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j-1)));
                    corr_sum_7 = corr_sum_7 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j-1)->at(i+1));
                }
                if (i != w-1 && j != h-1) {
                    neighbor_pixel = static_cast<double>(qGray(this->buffer_for_sliding_windowing->at(k)->pixel(i+1,j+1)));
                    corr_sum_5 = corr_sum_5 + (curr_pixel - curr_mean_double_in_window->at(j)->at(i)) * (neighbor_pixel - curr_mean_double_in_window->at(j+1)->at(i+1));
                }

            }
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_1);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_2);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_3);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_4);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_5);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_6);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_7);
            curr_sum_for_correlation_double_in_window->at(j)->at(i)->append(corr_sum_8);
        }
    }

    // XXX could be rewritten in CUDA
    double corr1, corr2, corr3, corr4, corr5, corr6, corr7, corr8;
    corr1 = 0.0;
    corr2 = 0.0;
    corr3 = 0.0;
    corr4 = 0.0;
    corr5 = 0.0;
    corr6 = 0.0;
    corr7 = 0.0;
    corr8 = 0.0;
    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            // 1 2 7      1 8 7
            // 8 C 4      2 C 6
            // 3 6 5      3 4 5
            double max_corr = 0.0;
            double this_var = curr_variance_double_in_window->at(j)->at(i);
            if (this_var <= 0.001)
            {this_var = 1.0;}
            double other_var = 0.0;
            if (i != 0) {
                other_var = curr_variance_double_in_window->at(j)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr8 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(7) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr8 > max_corr) {
                    max_corr = corr8;
                }
            }
            if (j != 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr2 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(1) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr2 > max_corr) {
                    max_corr = corr2;
                }
            }
            if (i != w-1) {
                other_var = curr_variance_double_in_window->at(j)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr4 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(3) / (WINDOW_SIZE* qSqrt(other_var * this_var )));
                if (corr4 > max_corr) {
                    max_corr = corr4;
                }
            }
            if (j != h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr6 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(5) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
               if (corr6 > max_corr) {
                   max_corr = corr6;
                }
            }
            if (i != 0 && j != 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr1 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(0) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr1 > max_corr) {
                    max_corr = corr1;
                }
            }
            if (i != 0 && j!= h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i-1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr3 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(2) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr3 > max_corr) {
                    max_corr = corr3;
                }
            }
            if ( i != w-1 && j!= 0) {
                other_var = curr_variance_double_in_window->at(j-1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr7 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(6) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr7 > max_corr) {
                    max_corr = corr7;
                }
            }
            if (i!= w-1 && j!= h-1) {
                other_var = curr_variance_double_in_window->at(j+1)->at(i+1);
                if (other_var <= 0.001) {
                    other_var = 1.0;
                }
                corr5 = fabs(curr_sum_for_correlation_double_in_window->at(j)->at(i)->at(4) / (WINDOW_SIZE* qSqrt(other_var * this_var)));
                if (corr5 > max_corr) {
                    max_corr = corr5;
                }
            }
            this->curr_correlation_double_in_window->at(j)->replace(i, max_corr);
        }
    }
}


void dr_group_image_buffer::set_do_sliding_windowing(bool doit) {
    if (doit) {
        this->do_sliding_windowing = true;
        this->buffer_for_sliding_windowing = new QVector<QImage *>();
    }
}


void dr_group_image_buffer::set_sliding_window_size(int size) {
    this->sliding_window_buffer_size = size;
}


int dr_group_image_buffer::get_window_buffer_size() {
    return this->buffer_for_sliding_windowing->size();
}


