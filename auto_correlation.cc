#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <sys/time.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <cblas.h>
#include <cstring>
#endif
/*
 */

using namespace std;
int r_size; // = 0.5* box_size /r_step
static float *gr =  0;
static float *corr = 0;
int seg_len;
int n_corr_seg;
void calc_nrm(float *gr, int seg_len, int n_corr_seg, float *nrm){
    float *gr_ptr = gr;
    for(int i = 0; i < n_corr_seg; i++, gr_ptr += seg_len){
        nrm[i] = cblas_snrm2(seg_len, gr_ptr, 1);
    }
    cout << "nrm 0 "<<nrm[0]<< endl;
}
void calc_corr(float *c, float *gr_0, float *gr, float *nrm_0, float *nrm, int seg_len, int n_corr_seg)
{
    float *grptr = gr;
    float *grptr0 = gr_0;
    for(int i = 0; i < n_corr_seg; i++, grptr+= seg_len, grptr0+=seg_len){
        c[i] = cblas_sdot(seg_len, grptr, 1,  grptr0, 1)/nrm[i]/nrm_0[i];
    }
}

void add_corr(float *c, float *gr, float *gr1, int seg_len, int n_corr_seg)
{
    float *grptr = gr;
    float *grptr0 = gr1;
    for(int i = 0; i < n_corr_seg; i++, grptr+= seg_len, grptr0+=seg_len){
        c[i] += cblas_sdot(seg_len, grptr, 1,  grptr0, 1);
    }
}

void calc_t_corr(float *c, int dt, int nt,  float *gr, int r_size, int seg_len, int n_corr_seg)
{
    for(int i = 0; i < nt-dt; i++){
        add_corr(c, gr+i*r_size, gr+(i+dt)*r_size, seg_len, n_corr_seg );
    }
    double scale = 1./(nt-dt);
    for(int i = 0; i < n_corr_seg; i++) c[i] *= scale;
}


int main(int argc, char *argv[])
{
    int r_size = 5000;
    int nt = 1000;
    int seg_len = 100;
    int n_corr_seg = 50;
    float *gr = new float[r_size*nt];
    ifstream fi(argv[1]);
    fi.read(reinterpret_cast <char *>(gr), r_size *nt*sizeof(float));
    float *corr = new float[500*n_corr_seg];
    memset(corr, 0, 500*n_corr_seg*sizeof(float));
    for(int i = 0; i <500; i++){
        calc_t_corr(corr+i*n_corr_seg, i+1, 1000, gr, r_size, seg_len, n_corr_seg);
    }
    ofstream fo("corr");
    fo.write(reinterpret_cast<char *>(corr), 500*n_corr_seg*sizeof(float));
    return 0;
}

