/***********************************************************************************************
 * Author Chong Shao (cshao@cs.unc.edu), Alfred Zhong
 *
 * This file is part of data reduction project.
 *
 * UNC has filed for patent protection on the code and intends to license it for commercial use.
 ***********************************************************************************************/


#include <iostream>
#include <dr_data.h>

using namespace std;

int RANSAC_METRIC = 7;
video_analysis::video_analysis(dr_group_image_buffer *img_buffer) 
{
    // first implement the multiple plots
    // then change it into one plot, with interactive stuff to switch between segments
    current_x = 10;
    current_y = 30;

    show_time_series_mean = false;
    show_time_series_stdev = false;

    ///TODO: modify code there to change everything turned on/off by boolean variables

    show_mean_images = false;
    show_bin_diff_images = false;
    show_one_mean_image = false;
    show_variance_image = false;
    show_correlation_image = false;
    show_time_series = false;
    show_pixel_hist = false;
    show_least_square_fitted_line = false;
    show_least_square_fitted_lines = false;
    show_windowed_correlation_image = false;
    show_pixel_cumulative_hist = false;
    show_normalized_pixel_cumulative_hist = false;
    show_variance_vs_mean_intensity_plot = false;
    show_least_square_fitted_line_on_variance_vs_mean_intensity_plot = false;

    cout << "before video initialization" << endl;
    video = new GLWidget(img_buffer, NULL);
    cout << "after video initialization" << endl;
    //plot a pixel's time series
    char plot_name[256];

    // one variance label
    int h = 0; int w= 0;

    if ( show_variance_image ) {
        video->img_buffer->get_group_stat(&groupStat);
        cout << "get group stat in video analysis" << endl;
        h = groupStat.var_images->at(0)->height();
        w = groupStat.var_images->at(0)->width();
        one_variance_index = 0;
    }

    // plot a new variance versus mean intensity scatter plot
    // here add a for loop to show multiple variance vs mean intensity plots

    if ( show_variance_vs_mean_intensity_plot ) {
        sprintf(plot_name,"Variance vs Mean Intensity");
        if (USE_LEAST_SQUARES || USE_MAXMIN_AND_MINMAX || USE_RANSAC || USE_THRESHOLD_BASED_TRIAL_NUMBER) {
            for (int i = 0; i < NUMBER_OF_SEGS; i++) {
                // use var_doubles
                // plot multiple variance vs mean intenstiy scatter plots
                // and the least square line fit plots (one day i will change this name into fitted line plots)
                // need to add the following variable into a QVector..
                double slope, offset;
                QVector<double> * const ys = new QVector<double>();
                QVector<double> * const xs = new QVector<double>();

                plot_pixel_variance_vs_mean_intensity_scatter_plot(slope, offset, ys, xs, i); // WORKING:: need to change this function to add ifelse branches
                printf("in main, slope and offset for current no. %d fragment: %f %f\n", i, slope, offset);

                plot_line_fit_plot(slope, offset, ys, xs, i);
            }
        }
    }

    if (USE_LEAST_SQUARES || USE_MAXMIN_AND_MINMAX || USE_RANSAC || USE_THRESHOLD_BASED_TRIAL_NUMBER) {
        /// TODO: change these handle button function names;;;
        QObject::connect(video, SIGNAL(mousePressed()), this, SLOT(handle_button4()));
    }

    QObject::connect(video, SIGNAL(mousePressed()), this, SLOT(handle_button3()));
    QObject::connect(video, SIGNAL(mousePressed()), this, SLOT(handle_button()));
    QObject::connect(video, SIGNAL(mousePressed()), this, SLOT(handle_button2()));

    QObject::connect(video, SIGNAL(mouseMoved(QString)), this, SLOT(pixel_timeseries_update(QString)));
    QObject::connect(video, SIGNAL(mouseMoved(QString)), this, SLOT(pixel_histogram_update(QString)));

    QObject::connect(video, SIGNAL(mouseMoved(QString)), this, SLOT(pixel_cumulative_histogram_update(QString)));
    QObject::connect(video, SIGNAL(mouseMoved(QString)), this, SLOT(pixel_normalized_histogram_update(QString)));
    /// TODO: change this if later
    if (show_time_series) {
        video->setGeometry(410,100,w,h);
        video->show();
    }
}


void video_analysis::plot_correlation_r_value_image() {
}

void video_analysis::plot_time_series() {
    int buffer_filled_size = video->img_buffer->get_size();
    double data_out[buffer_filled_size];

    int time_steps =  buffer_filled_size;
    //generate index array
    double indices[time_steps];
    for (int i=0; i<time_steps; i++) {
        indices[i] = i;
    }

    QVector<QPointF> points(time_steps);
    for (int i = 0; i < time_steps; i++)
    {
        QPointF pt(indices[i], data_out[i]);
        points[i] = pt;
    }
    char plot_name[256];
    sprintf(plot_name," Time Series, Pixel(%d, %d)",current_x,current_y);
}


void video_analysis::set_hist_bins(unsigned n) {
    h_bins = n;
}


void video_analysis::plot_histogram() {
    // class variables used:
    //       video <- this is a GLWidget object, which contains image buffer
    //       h_bins: number of histogram bins (?)
    //       current_x / current_y: somehow these two will automatically update when the GUI is running
    //
    //       All this small window stuff are initialized at the beginning when the class was constructed.
    //           in order to make a new one, go change the initializer (constructor).

    int time_steps = video->img_buffer->get_size(); // retrieve time steps

    int bins = h_bins;  // get number of bins

    QImage *h_max_image, *h_min_image;
    dr_group_stat g_stat;
    video->img_buffer->get_group_stat(&g_stat);   // callback function returns g_stat object
    h_max_image = g_stat.max_image;             // from g_stat get max and min images
    h_min_image = g_stat.min_image;

    //get max/min value
    double max = (double)qGray(h_max_image->pixel(current_x,current_y));
    double min = (double)qGray(h_min_image->pixel(current_x,current_y));
  //  double step_size = (max - min) / bins;   // step_size depends on the difference btwn max and min.
     if( min >= max ) {
        cout<<"Error: in dr_group_image_buffer::trace_one_pixelwithDistributionPlot, a pixel's minimum value must be smaller maximum value."<<endl;
        return;
    }

    QImage temp_img;
    //add data to the histogram
    for ( int k = 0; k < time_steps; k++ ) {
        video->img_buffer->get_frame(k, &temp_img);   // call back to fill up temp_img
   //     gsl_histogram_increment(h,  (double)qGray( temp_img.pixel(current_x, current_y)));
    }
    char plot_name[256];
    sprintf(plot_name,"Histogram, Pixel(%d, %d)",current_x,current_y);
}


void video_analysis::plot_line_fit_plot(double slope, double offset, const QVector<double> * y, const QVector<double> * x, int plot_index) {
    //line_fit_plots.at(plot_index)->setGeometry(0,0,640,400);
    QVector<double>  xs;
    QVector<double>  ys;
    QVector<double>  fitted_xs;
    QVector<double>  fitted_ys;
    int i;
    int len = y->size();
    // even though both x and y has first dimension >1 and second dimension 1, there the [] references are inverted...
    for (i = 0; i < len; i++) {
        xs.append((*x)[i]);
        ys.append((*y)[i]);
    }
    for (i = -10; i < 256; i++) {
        fitted_xs.push_back(static_cast<double>(i));
        fitted_ys.push_back(slope * static_cast<double>(i) + offset);
    }
}


void video_analysis::plot_pixel_variance_vs_mean_intensity_scatter_plot(double & slope, double & offset, QVector<double> * ys_out, QVector<double> * xs_out, int plot_index) {

    //   variance_vs_mean_intensity_scatter_plots.at(plot_index)->setGeometry(0,0,640,400);
    std::vector<double> xs;
    std::vector<double> ys;
    uint w=0, h=0;
    int meanIntensity;
    double varianceIntensity;
    uint i, j;

    dr_group_stat groupStat;
    video->img_buffer->get_group_stat(&groupStat);
    QVector<double> * variance = new QVector<double>(); // an array that keeps variance values in double, if it does not have terrible bugs, then meanImage-> pixel(j, i) is the same pixel location
    // in variance[i* w + j]
    QImage * meanImage = new QImage(*(groupStat.mean_images->at(0)));
    w = groupStat.mean_images->at(0)->width();
    h = groupStat.mean_images->at(0)->height();
    for (i = 0; i < h*w; i++) {
        variance->append(groupStat.var_doubles->at(plot_index)->at(i));
    }
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            int curr_pixel = groupStat.mean_images->at(plot_index)->pixel(j,i);
            meanImage->setPixel(j,i,qRgb(curr_pixel,curr_pixel,curr_pixel));
        }
    }
    // do x and y ordering
    // first get var image size

    TNT::Array2D<double> y_for_least_sq(236,1,0.0); // this only keeps min, for ransac we need the pts before that
    QVector<int> bins_count(236,0);  // bins count is for counting the bins with 10 more values, we only use this values in least squares fit
    int highCount = 0;
    int lowCount = 0;
    QVector<double> sumVarianceOfVarianceForEachMeanIntensity(236,0.0);
    for (i = 0; i < h; i++) {
        for (j = 0; j < w; j++) {
            // get the intensity value at the mean image, which is the x value
            //  then get the corresponding value at the variance image, which is the y value
            meanIntensity = qGray(meanImage->pixel(j,i));
            varianceIntensity= variance->at(i*w+j);
            //        if (varianceIntensity < 0.000001) {
            //         printf("zero variance at pixel location: %d %d \n", j , i );
            //      }
            if (varianceIntensity > 0.1) {
                highCount += 1;
            }
            else {
                lowCount += 1;
            }

            //   cout << "HERE GG\n\n";
            xs.push_back(static_cast<double>(meanIntensity));
            ys.push_back(varianceIntensity);

            // MOVE or change the dim of this thing
            // record the minimum y for each x
            // here use a array to keep number of balls in each bin, if it exceeds 10, record the min value.
            if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                sumVarianceOfVarianceForEachMeanIntensity[meanIntensity-10] += varianceIntensity;
                bins_count[meanIntensity-10] += 1;
                if (bins_count[meanIntensity-10] == 10) {
                    y_for_least_sq[meanIntensity-10][0] = varianceIntensity;
                }
                else {
                    if (bins_count[meanIntensity-10] >= 10 && y_for_least_sq[meanIntensity-10][0] > varianceIntensity) {
                        y_for_least_sq[meanIntensity-10][0] = varianceIntensity;
                    }
                }
            }
        }
    }

    // output the x and y into a txt file
    QVector<double> outputX(h*w, 0.0);
    QVector<double> outputY(h*w, 0.0);
    // fill these two arrays with xs and ys
    for (uint i = 0; i < xs.size(); i++) {
        outputX[i] = xs[i];
        outputY[i] = ys[i];
    }
    // for (int i = 0; i < outputX.size(); i++) {
    //   cout << outputX[i] << "," << outputY[i] << "\n";
    // }
    if (WRITE_XANDY_INTO_TXT) {
        output_x_and_y_into_txt(&outputX, &outputY);
    }

    printf("========================\nhigh count vs low count: %d %d \n========================\n", highCount, lowCount);

    // fit a line to the minimum value in each slot
    // before we do that,
    //  we should remove all the y values that has zero values.. and we don't use corresponding x values
    // use a vector

    // we can make a branch here so that we can choose doing a least square fit or a ransac
    // in any case we return the same thing: the points we peeked, the slope and the offset of the fitted line.
    if (USE_LEAST_SQUARES) {
        // two vectors storing the values used for fitting ( the bins with 10 or more values)
        std::vector<double> ys_for_least_sq_without_zero;
        std::vector<double> xs_for_least_sq_without_zero;

        // we record the length of the vector... [REFACTOR: change it into just calculate the length of the vector later]
        int len;
        len = 0;
        for (i = 0; i < 236; i++) {
            if (bins_count[i] >= 10) {
                ys_for_least_sq_without_zero.push_back(y_for_least_sq[i][0]);
                xs_for_least_sq_without_zero.push_back(static_cast<double>(i)+10.0);
                len+=1;
            }
        }

        // in order to do the least squares fit, we need to fill in the TNT::Array2D
        TNT::Array2D<double> A_for_least_sq(len,2,1.0);
        TNT::Array1D<double> x_for_least_sq(len,0.0);
        TNT::Array2D<double> y_tnt_for_least_sq_without_zero(len,1,0.0);
        TNT::Array1D<double> y_tnt_for_least_sq_without_zero_1d(len,0.0);
        for (i = 0; i < len; i++) {
            A_for_least_sq[i][0] = xs_for_least_sq_without_zero[i];
            x_for_least_sq[i] = xs_for_least_sq_without_zero[i]; // [REFACTOR: here we convert vector into TNTarray, and in another function we convert it back to vector. consider refactor]
            y_tnt_for_least_sq_without_zero[i][0] = ys_for_least_sq_without_zero[i];
            y_tnt_for_least_sq_without_zero_1d[i] = ys_for_least_sq_without_zero[i];
        }

        JAMA::QR<double> * qr_A_for_least_sq = new JAMA::QR<double>(A_for_least_sq);
        TNT::Array2D<double> x;
        x = qr_A_for_least_sq->solve(y_tnt_for_least_sq_without_zero);
        slope = x[0][0];
        offset = x[1][0];
        QVector<double> y_tnt_for_least_sq_without_zero_1d_qvector(len,0.0);
        for (i = 0; i < len; i++) {
            y_tnt_for_least_sq_without_zero_1d_qvector[i] = y_tnt_for_least_sq_without_zero_1d[i];
        }
        *ys_out = y_tnt_for_least_sq_without_zero_1d_qvector;
        QVector<double> xs_for_least_sq_qvector(len,0.0);
        for (i = 0; i < len; i++) {
            xs_for_least_sq_qvector[i] = x_for_least_sq[i];
        }
        *xs_out = xs_for_least_sq_qvector;
        printf("================\nResult:\nslope of the fitted line: %f\noffset of the fitted line: %f and x's dimension %d %d \n================\n\n",x[0][0],x[1][0],x.dim1(),x.dim2());
        qr_A_for_least_sq = NULL;
        cout << "HERE DD\n";
    }
    // if choose to use ransac, then it goes here.

    else if (USE_RANSAC) {
        if (! USE_THRESHOLD_BASED_TRIAL_NUMBER ) {
            // if use ransac,
            // we may want to obtain the first 30% stable region
            // calculate it here.
            QVector<double> thingsToPutInXsout;
            QVector<double> thingsToPutInYsout;
            int smallVarianceLength = 70; // 236 * 0.3
            QVector<double> * lineParams = new QVector<double>(2,0.0);
            QVector<double> weights(236,1.0);

            if (RANSAC_METRIC == 4) {
                thingsToPutInXsout = *(new QVector<double>(smallVarianceLength, 0.0));
                thingsToPutInYsout = *(new QVector<double>(smallVarianceLength, 0.0));
            }
            else {
                thingsToPutInXsout = *(new QVector<double>(h*w, 0.0));
                thingsToPutInYsout = *(new QVector<double>(h*w, 0.0));
            }

            if (RANSAC_METRIC == 1){
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 3) {
                // parameters needed will be:
                //  a list of weight for each bin
                //  x and y
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // calculate weight here
                // the weights are (squared) inverse proportional to the number of pts in each bin
                for (int i = 0; i < 236; i++)  {
                    weights[i] = (1.0 / ((1.0+ static_cast<double>(bins_count[i])))) * (1.0 / (1.0+ static_cast<double>(bins_count[i])));
                }

                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 4) {
                // the following parameters might be used in ransac:
                QVector<double> varianceOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            varianceOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    varianceOfVarianceForEachMeanIntensity[i] = varianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                // get the first 30% variance indices
                QVector<int> smallVarianceIndices(smallVarianceLength,-1);
                // do an intertion sort, may change it into merge sort later
                TNT::Array2D<double> sortedVariances(236,2,0.0);
                for (i = 0; i < 236; i++) {
                    for (j = 0; j < i; j++) {
                        if (sortedVariances[j][0] != 0.0 && sortedVariances[j][0] > varianceOfVarianceForEachMeanIntensity[i]) {
                            // insert, push all things after more step more
                            for (int k = i; k > j; k--) {
                                sortedVariances[k][0] = sortedVariances[k-1][0];
                                sortedVariances[k][1] = sortedVariances[k-1][1];
                            }
                            // insert at that place
                            sortedVariances[j][0] = varianceOfVarianceForEachMeanIntensity[i];
                            sortedVariances[j][1] = static_cast<double>(i);
                        }
                    }
                }
                // fill in the 30% first indices
                for (i = 0; i < smallVarianceLength; i++) {
                    smallVarianceIndices[i] = sortedVariances[i][1];
                }

                // do ransac
                // still use the binned points
                //  then what are the points needed for display?
                //  all the binned points?
                // we call ransac(xs, ys) to return the line
                // then we just use the original xs , ys as xs_out, ys_out

                // change the xs_out and ys_out
                // how they are the variances points we want
                // first use a vector to record the smallVariances xs and ys
                std::vector<double> smallVarianceXs;
                std::vector<double> smallVarianceYs;
                // check that there are actually things stored in smallVarianceIndices
                printf("small variance length: %d\n", smallVarianceLength);
                printf("small variance indices lnegth: %d\n", smallVarianceIndices.size());
                for (i = 0; i < smallVarianceIndices.size(); i++) {
                    printf("smallVarianceIndice[i]: %d\n", smallVarianceIndices[i]);
                }
                for (i = 0; i < xs.size(); i++) {
                    for (j = 0; j < smallVarianceLength; j++) {
                        if (xs[i] == smallVarianceIndices[j]) {
                            smallVarianceXs.push_back(smallVarianceIndices[j]);
                            smallVarianceYs.push_back(ys[i]);
                            cout << "======================\n pushed one \n =======================\n";
                            break;
                        }
                    }
                }
                // check that there are actually things stored in smallVarianceXs
                for (int i = 0; i < smallVarianceXs.size(); i++) {
                    printf("smallVarianceXs[i]: %f\n", smallVarianceXs[i]);
                }
                for (i = 0; i < smallVarianceXs.size(); i++) {
                    thingsToPutInXsout[i] = smallVarianceXs[i];
                    thingsToPutInYsout[i] = smallVarianceYs[i];
                }
                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 5) {
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 6) {
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // fill in the array of weights
                // the weights are inverse proportional to the standard deviation
                // first compute the standard deviation

                QVector<double> stDevOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            stDevOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    if (bins_count[i] >= 2) {
                        stDevOfVarianceForEachMeanIntensity[i] = sqrt(stDevOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]));
                    }
                    else {
                        stDevOfVarianceForEachMeanIntensity[i] = 0.0;
                    }
                }

                // fill in the weight array
                for (int i = 0; i < 236; i++) {
                    if (stDevOfVarianceForEachMeanIntensity[i] < 0.000001) {  //stdev is zero means that it is not defined (number of pts in a bin is 0) or it has only one pt
                        weights[i] = 0.0;
                    }
                    else {
                        weights[i] = 1.0 / (1.0 + stDevOfVarianceForEachMeanIntensity[i]);
                    }
                }

                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 7) {
                // make the weights to consider both std deviation and number of pts in the bin

                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // fill in the array of weights
                // the weights are inverse proportional to the standard deviation
                // first compute the standard deviation

                QVector<double> stDevOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            stDevOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    if (bins_count[i] >= 2) {
                        stDevOfVarianceForEachMeanIntensity[i] = sqrt(stDevOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]));
                    }
                    else {
                        stDevOfVarianceForEachMeanIntensity[i] = 0.0; //
                    }
                }
                // fill in the weight array
                for (int i = 0; i < 236; i++) {
                    if (stDevOfVarianceForEachMeanIntensity[i] < 0.000001) {  //stdev is zero means that it is not defined (number of pts in a bin is 0) or it has only one pt
                        weights[i] = 0.0;
                    }
                    else {
                        weights[i] = 1.0 / ( 1.0 + (stDevOfVarianceForEachMeanIntensity[i] * static_cast<double>(bins_count[i]) ) );
                    }
                }
                ransac(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }
            // check the line's parameters are actually passed:
            slope = (*lineParams)[0];
            offset = (*lineParams)[1];

            printf("\nin scatter plot(), after ransac, the line parameters: %f, %f\n", slope, offset);

            // change the things that xs_out and ys_out are pointing to
            *xs_out = thingsToPutInXsout;
            *ys_out = thingsToPutInYsout;
        }

        else if (USE_THRESHOLD_BASED_TRIAL_NUMBER) { // number of trials depend on whether you get a better result after, say 1000 trials
            // if use ransac,
            // we may want to obtain the first 30% stable region
            // calculate it here.
            QVector<double> thingsToPutInXsout;
            QVector<double> thingsToPutInYsout;
            int smallVarianceLength = 70; // 236 * 0.3
            QVector<double> * lineParams = new QVector<double>(2,0.0);
            QVector<double> weights(236,1.0);

            if (RANSAC_METRIC == 4) {
                thingsToPutInXsout = *(new QVector<double>(smallVarianceLength, 0.0));
                thingsToPutInYsout = *(new QVector<double>(smallVarianceLength, 0.0));
            }
            else {
                thingsToPutInXsout = *(new QVector<double>(h*w, 0.0));
                thingsToPutInYsout = *(new QVector<double>(h*w, 0.0));
            }

            if (RANSAC_METRIC == 1){
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 3) {
                // parameters needed will be:
                //  a list of weight for each bin
                //  x and y
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // calculate weight here
                // the weights are (squared) inverse proportional to the number of pts in each bin
                for (int i = 0; i < 236; i++)  {
                    weights[i] = (1.0 / ((1.0+ static_cast<double>(bins_count[i])))) * (1.0 / (1.0+ static_cast<double>(bins_count[i])));
                }

                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 4) {
                // the following parameters might be used in ransac:
                QVector<double> varianceOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            varianceOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    varianceOfVarianceForEachMeanIntensity[i] = varianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                // get the first 30% variance indices
                QVector<int> smallVarianceIndices(smallVarianceLength,-1);
                // do an intertion sort, may change it into merge sort later
                TNT::Array2D<double> sortedVariances(236,2,0.0);
                for (i = 0; i < 236; i++) {
                    for (j = 0; j < i; j++) {
                        if (sortedVariances[j][0] != 0.0 && sortedVariances[j][0] > varianceOfVarianceForEachMeanIntensity[i]) {
                            // insert, push all things after more step more
                            for (int k = i; k > j; k--) {
                                sortedVariances[k][0] = sortedVariances[k-1][0];
                                sortedVariances[k][1] = sortedVariances[k-1][1];
                            }
                            // insert at that place
                            sortedVariances[j][0] = varianceOfVarianceForEachMeanIntensity[i];
                            sortedVariances[j][1] = static_cast<double>(i);
                        }
                    }
                }
                // fill in the 30% first indices
                for (i = 0; i < smallVarianceLength; i++) {
                    smallVarianceIndices[i] = sortedVariances[i][1];
                }

                // do ransac
                // still use the binned points
                //  then what are the points needed for display?
                //  all the binned points?
                // we call ransac(xs, ys) to return the line
                // then we just use the original xs , ys as xs_out, ys_out

                // change the xs_out and ys_out
                // how they are the variances points we want
                // first use a vector to record the smallVariances xs and ys
                std::vector<double> smallVarianceXs;
                std::vector<double> smallVarianceYs;
                // check that there are actually things stored in smallVarianceIndices
                printf("small variance length: %d\n", smallVarianceLength);
                printf("small variance indices lnegth: %d\n", smallVarianceIndices.size());
                for (i = 0; i < smallVarianceIndices.size(); i++) {
                    printf("smallVarianceIndice[i]: %d\n", smallVarianceIndices[i]);
                }
                for (i = 0; i < xs.size(); i++) {
                    for (j = 0; j < smallVarianceLength; j++) {
                        if (xs[i] == smallVarianceIndices[j]) {
                            smallVarianceXs.push_back(smallVarianceIndices[j]);
                            smallVarianceYs.push_back(ys[i]);
                            cout << "======================\n pushed one \n =======================\n";
                            break;
                        }
                    }
                }
                // check that there are actually things stored in smallVarianceXs
                for (int i = 0; i < smallVarianceXs.size(); i++) {
                    printf("smallVarianceXs[i]: %f\n", smallVarianceXs[i]);
                }
                for (i = 0; i < smallVarianceXs.size(); i++) {
                    thingsToPutInXsout[i] = smallVarianceXs[i];
                    thingsToPutInYsout[i] = smallVarianceYs[i];
                }
                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 5) {
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 6) {
                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // fill in the array of weights
                // the weights are inverse proportional to the standard deviation
                // first compute the standard deviation
                QVector<double> stDevOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            stDevOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    if (bins_count[i] >= 2) {
                        stDevOfVarianceForEachMeanIntensity[i] = sqrt(stDevOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]));
                    }
                    else {
                        stDevOfVarianceForEachMeanIntensity[i] = 0.0;
                    }
                }
                // fill in the weight array
                for (int i = 0; i < 236; i++) {
                    if (stDevOfVarianceForEachMeanIntensity[i] < 0.000001) {  //stdev is zero means that it is not defined (number of pts in a bin is 0) or it has only one pt
                        weights[i] = 0.0;
                    }
                    else {
                        weights[i] = 1.0 / (1.0 + stDevOfVarianceForEachMeanIntensity[i]);
                    }
                }

                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }

            else if (RANSAC_METRIC == 7) {
                // make the weights to consider both std deviation and number of pts in the bin

                // fill these two arrays with xs and ys
                for (int i = 0; i < xs.size(); i++) {
                    thingsToPutInXsout[i] = xs[i];
                    thingsToPutInYsout[i] = ys[i];
                }
                // fill in the array of weights
                // the weights are inverse proportional to the standard deviation
                // first compute the standard deviation

                QVector<double> stDevOfVarianceForEachMeanIntensity(236,0.0); // an array that stores the variance of all the bins
                // To fill in the above array, we need to first calculate the mean of all bins
                QVector<double> meanVarianceOfVarianceForEachMeanIntensity(236,0.0);

                // do a for loop for it
                for (i = 0; i < 236; i++) {
                    // calculate the mean, fill in the second array
                    meanVarianceOfVarianceForEachMeanIntensity[i] = sumVarianceOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]);
                }
                double diffFromMeanVariance;
                for (i = 0; i < h; i++) {
                    for (j = 0; j < w; j++) {
                        meanIntensity = qGray(meanImage->pixel(j,i));
                        varianceIntensity = variance->at(i*w+j);
                        if (meanIntensity-10 < 236 && meanIntensity-10 > -1) {
                            // retrieve mean intensity again for index
                            // sum up the diff from mean for caclulating the variance of variance
                            diffFromMeanVariance = varianceIntensity - meanVarianceOfVarianceForEachMeanIntensity[meanIntensity-10];
                            stDevOfVarianceForEachMeanIntensity[meanIntensity-10] += diffFromMeanVariance * diffFromMeanVariance;
                        }
                    }
                }
                // divide by n to get the variance
                for (i = 0; i < 236; i++) {
                    if (bins_count[i] >= 2) {
                        stDevOfVarianceForEachMeanIntensity[i] = sqrt(stDevOfVarianceForEachMeanIntensity[i] / static_cast<double>(bins_count[i]));
                    }
                    else {
                        stDevOfVarianceForEachMeanIntensity[i] = 0.0; //
                    }
                }
                // fill in the weight array
                for (int i = 0; i < 236; i++) {
                    if (stDevOfVarianceForEachMeanIntensity[i] < 0.000001) {  //stdev is zero means that it is not defined (number of pts in a bin is 0) or it has only one pt
                        weights[i] = 0.0;
                    }
                    else {
                        weights[i] = 1.0 / ( 1.0 + (stDevOfVarianceForEachMeanIntensity[i] * static_cast<double>(bins_count[i]) ) );
                        //      cout << weights[i];
                    }
                }
                ransac_with_dynamic_threshold(lineParams, &thingsToPutInXsout, &thingsToPutInYsout, &weights);
            }
            // check the line's parameters are actually passed:
            slope = (*lineParams)[0];
            offset = (*lineParams)[1];

            printf("\nin scatter plot(), after ransac, the line parameters: %f, %f\n", slope, offset);

            // change the things that xs_out and ys_out are pointing to
            *xs_out = thingsToPutInXsout;
            *ys_out = thingsToPutInYsout;
        }
    }

    else if (USE_MAXMIN_AND_MINMAX) {
        // do the method I proposed
        // minmax and maxmin,
        // and get the bisector
        // and inversely use extereme value theory, to include the points, and do ransac? or least squares?
        //   > right now use least squares.

        // parameters needed: variance (double) array as y values, x values are obtained from mean image on the corresponding pixel location.
        slope = 0.0;
        offset = 0.0;
    }
    else { // USE TWO STABLE REGIONS...

    }
    video->img_buffer->set_slope_and_offset(slope, offset);

    // delete variance;
    delete(meanImage);
    meanImage = NULL;
    return;
}


void video_analysis::plot_pixel_cumulative_distribution_histogram() {
    // class variables used:
    //       video <- this is a GLWidget object, which contains image buffer
    //       h_bins: number of histogram bins (?)
    //       current_x / current_y: somehow these two will automatically update when the GUI is running
    //
    //       All this small window stuff are initialized at the beginning when the class was constructed.
    //           in order to make a new one, go change the initializer (constructor).

    int time_steps = video->img_buffer->get_size(); // retrieve time steps

    int bins = 256;  // get number of bins

    // assume this histogram is from -6 to 6, 601 bins
    QVector<double> pixel_distribution(256,0.0);
    // fill in intensity over time.
    QVector<double> intensity_over_time(time_steps,0.0);
    QImage temp_img;
    for (int k = 0; k < time_steps; k++) {
        video->img_buffer->get_frame(k, &temp_img);   // call back to fill up temp_img
        intensity_over_time[k] = static_cast<double>(qGray(temp_img.pixel(current_x, current_y)));
        //     printf("intensity: %f\n", intensity_over_time[k]);
    }
    // need to do a normalization proceses here...
    // fill in pixel distribution.
    for (int k = 0; k < time_steps; k++) {
        int bin = static_cast<int>(intensity_over_time[k]);
        if (bin < 0) {
            pixel_distribution[0] += 1;
        }
        else if (bin > 255) {
            pixel_distribution[255] += 1;
        }
        else {
            pixel_distribution[bin] += 1;
        }
    }
    // normalize the pixel distribution so that it integrates up to 1
    for (int k = 0; k < 256; k++) {
        // printf("pixel distribution k: %f\n", pixel_distribution[k]);
        pixel_distribution[k] = pixel_distribution[k] / static_cast<double>(time_steps);

    }

    double current_cumulative_val = 0.0;
    // calculate the cumulative distribution

    //calculate the histogram using gsl
    //add data to the histogram
    for ( int k = 0; k < 256; k++ ) {
        current_cumulative_val += pixel_distribution[k];
    }
}


void video_analysis::plot_normalized_pixel_cumulative_distribution_histogram() {
    // class variables used:
    //       video <- this is a GLWidget object, which contains image buffer
    //       h_bins: number of histogram bins (?)
    //       current_x / current_y: somehow these two will automatically update when the GUI is running
    //
    //       All this small window stuff are initialized at the beginning when the class was constructed.
    //           in order to make a new one, go change the initializer (constructor).

    int time_steps = video->img_buffer->get_size(); // retrieve time steps

    int bins = 601;  // get number of bins
    // assume this histogram is from -6 to 6, 601 bins
    QVector<double> pixel_distribution(601,0.0);
    QVector<double> intensity_over_time(time_steps,0.0);
    QImage temp_img;
    // fill in intensity over time.
    double intensity_sum = 0.0;
    for (int k = 0; k < time_steps; k++) {
        video->img_buffer->get_frame(k, &temp_img);   // call back to fill up temp_img
        intensity_over_time[k] = static_cast<double>(qGray(temp_img.pixel(current_x, current_y)));
        intensity_sum += intensity_over_time[k];
        //     printf("intensity: %f\n", intensity_over_time[k]);
    }

    double intensity_mean = intensity_sum / static_cast<double>(time_steps);
    double sum_for_standard_deviation = 0.0;
    // calculate the standard deviation

    for (int k = 0; k < time_steps; k++) {
        double diff_from_mean = intensity_over_time[k] - intensity_mean;
        sum_for_standard_deviation += diff_from_mean * diff_from_mean;
    }

    double std_dev = sqrt( sum_for_standard_deviation / static_cast<double>(time_steps) );

    // normalize the intensity value
    for (int k = 0; k < time_steps; k++) {
        intensity_over_time[k] = intensity_over_time[k] - intensity_mean;
        intensity_over_time[k] = intensity_over_time[k] / std_dev;
    }

    // fill in pixel distribution.
    for (int k = 0; k < time_steps; k++) {
        int bin = static_cast<int>(intensity_over_time[k] * 100.0 );
        //    printf("bin: %d\n", bin);
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
    for (int k = 0; k < 601; k++) {
        // printf("pixel distribution k: %f\n", pixel_distribution[k]);
        pixel_distribution[k] = pixel_distribution[k] / static_cast<double>(time_steps);
    }

    double current_cumulative_val = 0.0;

    //add data to the histogram
    for ( int k = 0; k < 601; k++ ) {
        for ( int l = 0; l < static_cast<int>(current_cumulative_val * 100); l++) {
   //         gsl_histogram_increment(h, k-300);
        }
   //     current_cumulative_val += pixel_distribution[k];
    }
}


void video_analysis::pixel_timeseries_update(QString str) {
    //cout<<str.toUtf8().constData()<<endl;
    //parse to get the x,y coordinate out
    int i,n;
    n = str.size();
    for(i=0; i<n; i++)
        if (str[i] == QChar(','))
            break;
    current_x = str.left(i).toDouble();
    current_y = str.mid(i+2, 0xffffffff).toDouble();

    if ( current_x>=0 && current_x<video->img_buffer->get_width() && current_y>=0 && current_y<video->img_buffer->get_height() )
        plot_time_series();
}


void video_analysis::pixel_histogram_update(QString str) {
    //cout<<str.toUtf8().constData()<<endl;
    //parse to get the x,y coordinate out
    int i,n;
    n = str.size();
    for(i=0; i<n; i++)
        if (str[i] == QChar(','))
            break;
    current_x = str.left(i).toDouble();
    current_y = str.mid(i+2, 0xffffffff).toDouble();

    if ( current_x>=0 && current_x<video->img_buffer->get_width() && current_y>=0 && current_y<video->img_buffer->get_height() )
        plot_histogram();
}


void video_analysis::pixel_cumulative_histogram_update(QString str) {
    //cout<<str.toUtf8().constData()<<endl;
    //parse to get the x,y coordinate out
    int i,n;
    n = str.size();
    for(i=0; i<n; i++)
        if (str[i] == QChar(','))
            break;
    current_x = str.left(i).toDouble();
    current_y = str.mid(i+2, 0xffffffff).toDouble();

    if ( current_x>=0 && current_x<video->img_buffer->get_width() && current_y>=0 && current_y<video->img_buffer->get_height() )
        plot_pixel_cumulative_distribution_histogram();
}


void video_analysis::pixel_normalized_histogram_update(QString str) {
    //parse to get the x,y coordinate out
    int i,n;
    n = str.size();
    for(i=0; i<n; i++)
        if (str[i] == QChar(','))
            break;
    current_x = str.left(i).toDouble();
    current_y = str.mid(i+2, 0xffffffff).toDouble();

    if ( current_x>=0 && current_x<video->img_buffer->get_width() && current_y>=0 && current_y<video->img_buffer->get_height() )
        plot_normalized_pixel_cumulative_distribution_histogram();
}


void video_analysis::windowed_corr_pixel_update(QString str) {
    int i, n;
    n = str.size();
    for(i=0; i<n; i++)
        if (str[i] == QChar(','))
            break;
    current_x = str.left(i).toDouble();
    current_y = str.mid(i+2, 0xffffffff).toDouble();

    if ( current_x>=0 && current_x<video->img_buffer->get_width() && current_y>=0 && current_y<video->img_buffer->get_height() )
        plot_windowed_corr_timeseries_plot();
}


void video_analysis::set_show_time_series_mean(bool b) {
    show_time_series_mean = b;
}


void video_analysis::set_show_time_series_stdev(bool b) {
    show_time_series_stdev = b;
}


void video_analysis::ransac(QVector<double> * lineParams, QVector<double> * x, QVector<double> * y, QVector<double> * weights) {
    // peak random two points, and fit a line to it, check how many points agree with that line
    int numberOfTrials = 10000;
    // the record array records the score, slope and offset for every trial
    TNT::Array2D<double> records(numberOfTrials, 3, 0.0);

    int i, j;
    int index1, index2;
    double curr_slope, curr_offset, a, b, c;
    double x1, y1, x2, y2;
    double minimum_overall_score = 0.0;
    int minimum_overall_score_index = 0;

    for (i = 0; i < numberOfTrials; i++) {
        // peak random two points
        index1 = rand_int(0, x->size());
        x1 = ((*x)[index1]);
        y1 = ((*y)[index1]);
        do {
            index2 = rand_int(0, x->size());
            x2 = ((*x)[index2]);
            y2 = ((*y)[index2]);

        } while (x2 == x1);

        // fit a line to it
        curr_slope = (y2 - y1) / (x2 - x1);
        curr_offset = y1 - curr_slope * x1;
        // convert it into ax + by + c = 0 form
        a = curr_slope;
        b = -1.0;
        c = curr_offset;
        double divisor = sqrt(pow(a,2.0) + pow(b,2.0));
        //  printf("ransac, divisor: %f\n", divisor);

        // calculate the overall score of this fit
        // there are multiple ways of calculating scores...
        // first approach: sum of distancese over all pts (wrong!)
        // second approach: counting the number of pts in a band, band width depends on variance
        // third approach: counting the weighted distance over all pts, the weight is inverse porportional to the variance
        // fourth approach: only calculate the sumed distance score on the stable regions ( first 30 % percent variance ),
        //   NOTE that the variance means the variance of variance vs mean intensity, i.e., the place that variance does not change much
        // -----
        // the first one is implemented
        //  now implementing fourth one.
        // -----
        // this might be slow...
        // print something to see how slow it is
        double curr_dist, curr_x, curr_y;
        if (RANSAC_METRIC == 1) {
            double overall_score = 0.0; // the smaller the score is, the better.
            for (j = 0; j < x->size(); j++) {
                // calculate how far it is from that line
                curr_x = (*x)[j];
                curr_y = (*y)[j];
                if (curr_x >= 10 && curr_x < 246) {
                    curr_dist = abs(a * curr_x + b * curr_y + c) / divisor;
                    overall_score += curr_dist;
                }
            }
            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score < 0.00001) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
        }

        else if (RANSAC_METRIC == 3) {
            // stable region: calculate the weighted sum of distances according to number
            double overall_score = 0.0;
            for (j = 0; j < x->size(); j++) {
                curr_x = (*x)[j];
                curr_y = (*y)[j];
                //   printf("current x: %d\n" ,static_cast<int>(curr_x - 10));
                if (curr_x >= 10 && curr_x < 246) {
                    double weight = (*weights)[static_cast<int>(curr_x - 10)];

                    //  printf("weight times divisor: %f\n",(weight * divisor));
                    curr_dist = abs(a * curr_x + b * curr_y + c) / (weight * divisor) ;
                    //   printf("current distance: %f\n", curr_dist);
                }
                overall_score += curr_dist;
            }
            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score == 0.0) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
        }
        else if (RANSAC_METRIC == 4) {
            //fourth approach: only calculate the sumed distance score on the stable regions ( first 30 % percent variance )
            double overall_score = 0.0;
            // if the metric 4 is used, then the x and y passed should
            int lenOfSmallVariance = x->size();
            for (j = 0; j < lenOfSmallVariance; j++) {
                curr_x = (*x)[j];
                curr_y = (*y)[j];
                curr_dist = std::abs(a * curr_x + b * curr_y +c) / divisor;
                overall_score += curr_dist;
            }

            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score == 0.0) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
        }

        else if (RANSAC_METRIC == 5) {
            // metric 5: using sum of vertical distance
            double overall_score = 0.0; // the smaller the score is, the better.

            for (j = 0; j < x->size(); j++) {
                // calculate how far it is from that line
                curr_x = (*x)[j];
                curr_y = (*y)[j];
                curr_dist = abs(curr_y - curr_slope * curr_x - curr_offset);
                overall_score += curr_dist;
            }
            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score < 0.0000001) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
        }

        else if (RANSAC_METRIC == 6) {
            // metric 6: using sum of weighted vertical distance
            double overall_score = 0.0; // the smaller the score is, the better.

            for (j = 0; j < x->size(); j++) {
                // calculate how far it is from that line, vertically, with weight that is inverse proportional to the standard deviation
                curr_x = (*x)[j];
                if (curr_x >= 10 && curr_x < 246) {
                    curr_y = (*y)[j];
                    curr_dist = abs(curr_y - curr_slope * curr_x - curr_offset) * (*weights)[static_cast<int>(curr_x) - 10];
                    overall_score += curr_dist;
                }
            }
            printf("ransac considering vertical distance only, calculating the score for current line: %f\n", overall_score);
            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score < 0.0000001) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
        }

        else if (RANSAC_METRIC == 7) {
            double overall_score = 0.0;
            //  QVector<double> outputX;
            //  QVector<double> outputY;
            for (j = 0; j < x->size(); j++) {
                curr_x = (*x)[j];
                if (curr_x >= 10 && curr_x < 246) {
                    //         outputX.push_back(curr_x);
                    //         outputY.push_back(curr_y);
                    curr_y = (*y)[j];
                    curr_dist = abs(curr_y - curr_slope * curr_x - curr_offset) * weights->at(static_cast<int>(curr_x) - 10);
                    overall_score += curr_dist;
                }
            }

            records[i][0] = overall_score;
            records[i][1] = curr_slope;
            records[i][2] = curr_offset;
            if (minimum_overall_score < 0.0000001) {
                minimum_overall_score = overall_score;
            }
            else if ( overall_score < minimum_overall_score ) {
                minimum_overall_score = overall_score;
                minimum_overall_score_index = i;
            }
            //       *x = outputX;
            //       *y = outputY;
        }
    }
    printf("==============\nransac, overall minimum index: %d, overall minimum score %f slope %f offset %f\n================n", minimum_overall_score_index, records[minimum_overall_score_index][0],records[minimum_overall_score_index][1],records[minimum_overall_score_index][2]);

    (*lineParams)[0] = (records[minimum_overall_score_index][1]);
    (*lineParams)[1] = (records[minimum_overall_score_index][2]);
    return;
}


int video_analysis::rand_int(int low, int high) {
    qsrand(QTime::currentTime().msec());
    return qrand () % ((high + 1) - low) + low;
}


void video_analysis::output_x_and_y_into_txt(QVector<double> * xs, QVector<double> * ys) {
    if (xs->size() != ys->size()) {
        printf("x and y does not have the same size in output_x_and_y_into_txt!\n");
        return;
    }
    QFile file("outputXsandYs.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    for (int i = 0; i < xs->size(); i++) {
        out << (*xs)[i] << "," << (*ys)[i] << "\n";
    }
    file.close();
    return;
}


int video_analysis::calculate_current_number_of_trails(double threshold, double inlierRatio, int m) {
    int k;
    k = ceil(log(1.0-threshold) / log(1.0-pow(inlierRatio, static_cast<double>(m))));
    return k;
}


void video_analysis::set_bin_diff_imgs(QVector<QImage*> * bin_diff_imgs) {
    this->bin_diff_imgs = bin_diff_imgs;
    return;
}


QRgb video_analysis::gray_scale_to_blackbody_rgb(int gray_scale_intensity) {
    double red_intensity = 0;
    double green_intensity = 0;
    double blue_intensity = 0;
    if (gray_scale_intensity < 255 / 3) {
        red_intensity = gray_scale_intensity * 3;
    }
    else {
        if (gray_scale_intensity < 510 / 3) {
            red_intensity = 255;
            green_intensity = (gray_scale_intensity - (255 / 3)) * 3 ;
        }
        else {
            red_intensity = 255;
            green_intensity = 255;
            blue_intensity = (gray_scale_intensity - (510 / 3)) * 3 ;
        }
    }
    return qRgb(red_intensity, green_intensity, blue_intensity);
}


void video_analysis::plot_windowed_corr_timeseries_plot() {
    dr_group_stat groupStat;
    video->img_buffer->get_group_stat(&groupStat);
    QVector<QImage *> * windowed_corr_images  = groupStat.window_corr_images;
    QVector<QPointF> points(windowed_corr_images->size());

    for (int i = 0; i < windowed_corr_images->size(); i++) {
        QPointF pt(i,windowed_corr_images->at(i)->pixel(current_x, current_y));
        points[i] = pt;
    }
    char plot_name[256];
    sprintf(plot_name,"Corr Time Series, Pixel(%d, %d)", current_x, current_y);
}
