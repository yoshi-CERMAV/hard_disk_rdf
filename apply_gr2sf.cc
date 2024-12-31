#include <iostream>
#include <fstream>
#include <Accelerate/Accelerate.h>
#include <sys/stat.h>
#include <gsl/gsl_sf_bessel.h>
using namespace std;
double r_max, q_max, r_step, q_step;
int r_size, q_size;
float *bess_mat=0;
float *gr = 0;
float *sf = 0;
float *scale_array = 0;
int r_division = 10;
int n;

int get_number(char *posfile)
{
    int n;
    ifstream fi(posfile);
    if(!fi) cerr <<"could not open "<< posfile<<endl;
    fi.read(reinterpret_cast <char *> (&n), sizeof(int));
    return n;
}


int apply_gr2sf(char *mat_file, int size, int r_size_, int q_size_, int r_end, char *gr_file, char *pos_file, char *out_file)
{
    r_size = r_size_;
    q_size = q_size_;
    struct stat stat_buf;
    stat(mat_file, &stat_buf);
    int len =r_size*q_size*sizeof(float);
    assert(( void("size does not match"), stat_buf.st_size == len) );
    bess_mat = new float[r_size *q_size];
    ifstream fi(mat_file);
    fi.read(reinterpret_cast<char *>(bess_mat), r_size *q_size*sizeof(float));
    cout <<"r_size "<<r_size <<endl;
    stat(gr_file, &stat_buf);
    len =r_size*sizeof(float);
    cout << "stat_buf st_size"<<stat_buf.st_size<<endl;
    cout <<!(stat_buf.st_size%len)<<endl;
    assert((void("size does not match"), !(stat_buf.st_size%len))) ;
    int nframe = stat_buf.st_size / len;
    if(gr) delete [] gr;
    gr = new float[nframe * r_size];
    ifstream fi_gr(gr_file);
    fi_gr.read(reinterpret_cast <char *>(gr), len*nframe);
    
    if(pos_file){
        n = get_number(pos_file);
    }
    cout <<"n = "  <<n<<endl;
    
    double phi = M_PI*n/size/size;
    cout <<"phi = "<< phi<<endl;
    r_step = 1./r_division;
    double number_density = 1.0*n/size/size;
    float scale = M_PI*2*r_step * number_density;
    cout << "scale = "<< scale<<endl;

//    scale_array = new float[r_size];
//    std::fill_n(scale_array, r_size, 1);
    int transition = r_end/2;
    double step = 1./transition;
    int offset = r_end-transition;
    for(int j = 0; j < nframe; j++){
        for(int i = 0; i !=transition; i++ ){
            gr[r_size*j + i+offset] *= (1- step*i);
        }
    }
    
    if(sf) delete [] sf;
    sf = new float[nframe * r_size];
    std::fill_n(sf, q_size * nframe, 1);
    cblas_sgemm(CblasColMajor, CblasNoTrans,CblasNoTrans, q_size, nframe, r_end, scale, bess_mat, q_size, gr, r_size, 1,  sf, q_size);
    ofstream fo(out_file);
    fo.write(reinterpret_cast<char *>(sf), q_size * nframe * sizeof(float));
    return 0;
}


int main(int argc, char *argv[])
{
    int box_size = 500;
    int r_size = 5000;
    int q_size = 5000;
    char *gr2sf_file, *gr_file, *sf_file, *pos_file, *out_file ;
    char filename[] = "gr2sf.dat";
    gr2sf_file = filename;
    char gr_filename[] = "gr.dat";
    gr_file = gr_filename;
    char sf_filename[] = "sf.dat";
    sf_file = sf_filename;
    char pos_filename[] = "pos.dat";
    pos_file = 0;
    char c;
    float cutoff = 0.5;
    while ((c = getopt (argc, argv, "l:m:g:f:p:n:q:r:c:")) != -1){
        switch (c)
        {
            case 'l':
                box_size = atoi(optarg); // box_size
                cout << "box_size "<<box_size<<endl;
                break;
            case 'm':
                gr2sf_file = optarg;
                cout << "gs2sf_file "<<gr2sf_file<<endl;
                break;
            case 'g':
                gr_file = optarg;
                cout <<"gr_file "<<gr_file<<endl;
                break;
            case 'f':
                sf_file = optarg;
                cout <<"sf_file " <<sf_file<<endl;
                break;
            case 'p':
                pos_file = optarg;
                cout <<"pos_file "<<pos_file<<endl;
                break;
            case 'n':
                n = atoi(optarg);
                cout << "n "<<n <<endl;
                break;
            case 'q':
                q_size = atoi(optarg);
                cout <<"q_size "<< q_size<<endl;
                break;
            case 'r':
                r_size = atoi(optarg);
                cout << "r_size "<< r_size <<endl;
                break;
            case 'c':
                cout << optarg<<endl;
                cutoff = atof(optarg);
                cout <<"cutoff " <<cutoff<<endl;
                break;
        }
    }
                        
    int r_end = r_size * cutoff;
    cout << "r_end "<<r_end<<endl;
    apply_gr2sf(gr2sf_file, box_size, r_size, q_size, r_end, gr_file, pos_file, sf_file);
    return 0;
}
