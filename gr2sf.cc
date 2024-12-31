#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <unistd.h>

using namespace std;
double r_max, q_max, r_step, q_step;
int r_size, q_size;
float *bess_mat;
int r_division = 10;
int box_size;
int init_bess_mat(int box_size_, int r_division_, double q_max)
{
    box_size = box_size_;
    r_division = r_division_;
    r_max = box_size;
//    q_max = r_division*2;
    r_step = 1./r_division;
    r_size = r_division * box_size;// = r_max/r_step;
    q_size = r_division * box_size;//q_max/q_step;
    q_step = q_max/q_size;

    bess_mat = new float[q_size * r_size];
    float *ptr = bess_mat;
    double r = 0.5 * r_step;
    for(int j = 0; j < r_size; j++, r+= r_step){
        double q = q_step;
        for(int i = 0; i < q_size; i++, q+= q_step, ptr++){
            double qr = q*r;
            *ptr = r*gsl_sf_bessel_J0(qr);
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{
    int box_size_ = 500;
    int r_div_ = 10;
    double qmax = r_div_*2;
    char *outfile;
    char filename[] = "gr2sf.dat";
    outfile = filename;
    char c;
    while ((c = getopt (argc, argv, "l:d:f:q:")) != -1){
        switch (c)
        {
            case 'l':
                box_size_ = atoi(optarg); // box_size
                break;
            case 'd':
                r_div_ = atoi(optarg); // division of radius
                break;
            case 'f':
                outfile = optarg; // output file
                break;
            case 'q':
                qmax= atof(optarg);
                break;
        }
    }
    cout << "box_size "<< box_size_<<endl;
    cout << "r_div "<< r_div_<<endl;
    cout << "q_max "<< qmax<<endl;
    cout << "output file "<< outfile<<endl;
    init_bess_mat(box_size_, r_div_, qmax);
    ofstream fo(outfile);
    fo.write(reinterpret_cast<char *>(bess_mat), r_size *q_size*sizeof(float));
}
