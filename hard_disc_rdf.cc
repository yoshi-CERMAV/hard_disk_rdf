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
 1. calculation of line-picking probability density in the sampling area of radius R.
 2.
 */

using namespace std;

int box_size; //lateral size of the sampling square
double radius2;
double phi; // volume (area) fraction
int area; // area of the sampling square
double circle_area;//
int n; // number of particles
int n_iteration; // basic unit of iteration
double number_density; // number density in the square box
double number_density1; // number density in the sampling circle
double a; //starting center_to_center distance
int ix, iy;
double x_step;
double y_step;
int n1; // real number of particles

int r_division = 0; //division of a unit length
double r_step; //step in the radial distance
double r_max;
double q_step;
double q_max;
int r_size; // = 0.5* box_size /r_step
int q_size; // =

int hist_size = 0;
int v_hist_size;
int v_hist_bin;
static int *v_hist = 0;
int h_hist_size = 0;
int h_hist_bin;
static int *h_hist = 0;
static float *bess_mat = 0;
static int *hist; // counting the pairs at given distances
static float *normalized_hist = 0; // normalized by the total number of pairs
static float *normalized_hist0 = 0;
//float *gr0 = 0;
static float *gr = 0;


static float *line_picking = 0; //line_picking probability
static float *r_line_picking = 0; // inverse of line_picking probability for normalization
static float *sf = 0;

char filename[255];

static int *address_keeper = 0;  //The whole box is divided into aquare tiles of size 2/sqrt(2) so that only one disc can be present on a tile
int address_lda, address_lda1; //number of tiles in one direction
int address_size;// the total number of tiles

double sqrt2_2 = sqrt(2.)*0.5;
double sqrt3_2 = sqrt(3)*0.5;

double center_x;//
double center_y;

double new_x, new_y; //potential new position
double jump_x, jump_y; // for periodic boundary conditions
typedef struct pos{
    double x, y;
    int address; // the tile position
} position_structure;
static position_structure *disc_centers = 0;

default_random_engine generator;
normal_distribution<double> *distribution = 0;
uniform_int_distribution<std::mt19937::result_type> *flat_distribution = 0;
mt19937 *rng  = 0;
double step_factor = 1;
int reject_count;
int accept_count;

int allocate_float(float **ptr, int n)
{
    if(*ptr) delete[] *ptr;
    *ptr = new float[n];
    return 0;
}

int allocate_int(int **ptr, int n)
{
    if(*ptr) delete[] *ptr;
    *ptr = new int[n];
    return 0;
}

int distance_check(int k)
{
    double x = disc_centers[k].x + jump_x - new_x;
    double y = disc_centers[k].y + jump_y - new_y;
    if(x < -10) {
        cerr <<x  << " "<<jump_x <<" "<< new_x <<" "<< k;
        exit(0);
    }
    if(x > 10) {
        cerr << x << " "<<jump_x <<" "<< new_x <<" "<< k;
        exit(0);
    }
    if(y < -10 || y > 10) {
        cerr << "y "<<y << " "<<jump_y <<" "<< new_y <<" "<< k;
        exit(0);
    }
    
    double D = x*x+y*y;
    return x*x+y*y < 4; // true if the distance is closer than 2.
}

int pbc(int i, double &shift) //periodic boundary condition
{
    if (i < 0){
        i += address_lda;
        shift = -box_size;
    }
    else if (i > address_lda1){
        i -= address_lda;
        shift = box_size;
    }
    else shift = 0;
    return i;
}

void vertical_histogram()
{
    memset(v_hist, 0, box_size * sizeof(int));
    for(int i = 0; i< n; i++){
        v_hist[int((disc_centers[i].y)/v_hist_bin)]++;
    }
}

void horizontal_histogram()
{
    memset(h_hist, 0, box_size * sizeof(int));
    for(int i = 0; i< n; i++){
        h_hist[int((disc_centers[i].x)/h_hist_bin)]++;
    }
}
void horizontal_histogram(ofstream &fo)
{
    horizontal_histogram();
    fo.write(reinterpret_cast<char *>(h_hist), sizeof(int)*h_hist_size);
}

int test_move(position_structure &pos){
    new_x = pos.x + (*distribution)(generator);
    new_y = pos.y + (*distribution)(generator);
    while (new_x > box_size) new_x -= box_size;
    while (new_y > box_size) new_y -= box_size;
    while (new_x<0) new_x += box_size;
    while (new_y<0) new_y += box_size;
    // now the new_x and new_y is between 0 and box_size
    int ix = int(new_x* sqrt2_2);
    int iy = int(new_y* sqrt2_2);
    // the tile size is sqrt(2)/2
    int new_address = ix + address_lda * iy;
    int i_old = pos.address%address_lda;
    int j_old = pos.address/address_lda;
    for(int j = -2; j != 3; j++){
        int j1 = pbc(iy+j, jump_y);
        j1 *= address_lda;
        for(int i = -2; i !=3; i++){
            int i1 = pbc(ix + i, jump_x);
            int k =address_keeper[i1 + j1];
            if( k != -1 &&  (i1 + j1) !=pos.address) {
                if(distance_check(k))
                    return 0;
            }
        }
    }
    accept_count++;
    pos.x = new_x;
    pos.y = new_y;
    pos.address = new_address;
    return 1;
}

int fill_initial_position()
{
    x_step = box_size/(double)ix;
    y_step = box_size/(double)iy;
    double x_step05 = x_step*0.5;
    int count = 0;
    for(int j = 0; j < iy; j++){
        for (int i = 0; i < ix; i++, count++){
            double x = i*x_step;
            if (j%2) x += x_step05;
            disc_centers[count].x = x;
            x *= sqrt2_2;
            double y = j * y_step;
            disc_centers[count].y = y;
            y *= sqrt2_2;
            int address = int(x) + address_lda * int(y);
            disc_centers[count].address = address;
            address_keeper[address] = count;
        }
    }
    return 0;
}

void calc_center_to_center()
{
    area = box_size*box_size;
    circle_area = area*M_PI*0.25;
    n = area/M_PI*phi;
    a = box_size / sqrt(n*sqrt3_2);
}

int rounding_number_of_discs()
{
    calc_center_to_center();
    cout << "initial number of cylinders "<< n<< endl;
    ix = box_size / a;
    iy = box_size / a / sqrt3_2;
    n = ix * iy;
    cout << "update number of cylinders "<< n<< endl;
    phi = n*M_PI/area;
    cout << "update phi "<< phi<< endl;
    number_density = (1.0*n)/area;
    return n;
}

int rounding_number_of_discs(double p, int s)
{
    phi = p;
    box_size = s;
    return rounding_number_of_discs();
}

void histogram_alloc()
{
    hist_size = box_size*r_division;
    v_hist_size = box_size/v_hist_bin;
    allocate_int(&v_hist, v_hist_size);
    h_hist_size = box_size/h_hist_bin;
    allocate_int(&h_hist, h_hist_size);
    allocate_int(&hist,hist_size);
    allocate_float(&normalized_hist, hist_size);
    allocate_float(&line_picking, hist_size);
    allocate_float(&r_line_picking, hist_size);
}

int allocate_centers()
{
    if(flat_distribution) delete flat_distribution;
    flat_distribution = new uniform_int_distribution<std::mt19937::result_type> (0,n);
    if (disc_centers) delete[] disc_centers;
    disc_centers = new position_structure[n];
    return 0;
}

void allocate_normal_distribution (double x)
{
    if(distribution) delete distribution;
    distribution = new normal_distribution<double> (0, x*(0.2*a-1)); // half the space between cylinders
}

void init_box()
{
    radius2 = box_size * 0.5;
    center_x = radius2;
    center_y = radius2;
    radius2 *= radius2;
    address_lda1 = box_size * sqrt2_2;
    address_lda = address_lda1 + 1;
    address_size = address_lda * address_lda;
    allocate_int(&address_keeper, address_size);
    for(int i = 0; i < address_size;i++) address_keeper[i]= -1;
}

int init(double p, int s, int r_div)
{
    box_size = s;
    phi = p;
    r_division = r_div;
    r_max = box_size;
    q_max = r_division;
    r_step = 1./q_max;
    r_size = r_division * box_size;// = r_max/r_step;
    q_step = 1./r_max;
    q_size = r_division * box_size;//q_max/q_step;
    
    init_box();
    rounding_number_of_discs();
    cout << r_size<<" "<<q_size<<endl;
    histogram_alloc();
    allocate_centers();
    cout << "number of disc "<<n << endl;
    cout <<" center to center dist "<< a <<endl;
    fill_initial_position();
    
    return 0;
}

int move(int count)
{
    int address = disc_centers[count].address;
    if(test_move(disc_centers[count])){
        address_keeper[address] = -1;
        address_keeper[disc_centers[count].address] = count;
    }else{
        reject_count++;
    }
    return 0;
}

int out_of_radius(int i)
{
    double x = disc_centers[i].x-center_x;
    double y = disc_centers[i].y-center_y;
    return (x*x + y*y) > radius2;
}

inline double distance(int i, int j)
{
    double x = disc_centers[i].x
    - disc_centers[j].x ;
    double y = disc_centers[i].y
    - disc_centers[j].y ;
    return sqrt(x*x+y*y);//hypot(x, y);//
}

long int accumulate_hist()
{
    double dist;
    long int count=0;
    int *in= new int[n];
    for(int i = 0; i < n; i++){
        if(out_of_radius(i)) continue;
        else in[count++] = i;
    }
    long int m = count;
    for(int i = 0; i < m; i++){
        int i4 = (i/4)*4;
        double distk[8];
        int posk[8];
        for(int j = 0; j < i4; j+=4){
            for(int k = 0; k< 4; k++){
                distk[k] = distance(in[i], in[j+k]);
                posk[k] = int(distk[k]*r_division);
                hist[posk[k]] ++;
            }
        }
        for(int j = i4; j < i; j++){
            dist = distance(in[i], in[j]);
            int pos = int(dist*r_division);
            hist[pos] ++;
        }
        
    }
    delete[] in;
    number_density1 = m/circle_area;
    return m;
}

int dump_hist(const char filename[])
{
    ofstream fo(filename);
    for(int i = 0; i < hist_size; i++){
        fo << i << " "<< normalized_hist[i]<<" "<< line_picking[i]<<endl;
    }
    return 0;
}

void dump(float *hist, const char filename[])
{
    ofstream fo(filename);
    fo.write(reinterpret_cast<char *>(hist), sizeof(float)*hist_size);
}


int hist_normalize(float *n_hist, double scale)
{
    for(int i = 0; i < hist_size; i++){
        n_hist[i] = hist[i] * scale;
    }
    return 0;
}

int hist_normalize(double scale)
{
    return hist_normalize(normalized_hist, scale);
}

int line_picking_init()
{
    double hist_step = 1./hist_size;
    double r = 0.5*hist_step;
    double A = 4/M_PI;
    for(int i = 0; i < hist_size; i++, r+= hist_step){
        line_picking[i] = A*r*(acos(r) - r *sqrt(1-r*r))*hist_step*4;
        r_line_picking[i] = 1./line_picking[i];
    }
    return 0;
}

int move_ntimes(int n_times)
{
    for(int i = 0; i< n_times; i++){
        int k =(*flat_distribution)(*rng);
        move(k);
    }
    return 0;
}

int calc_gr(float *p, float *gr)
{
    for(int i = 0; i < hist_size; i++){
        gr[i] = p[i] * r_line_picking[i];
        gr[i] -=1;
    }
    return 0;
}

void calc_sf(float  *gr)
{
    double scale = r_step *2*M_PI* number_density;
    cout << "scale "<< scale <<endl;
    std::fill_n(sf, q_size, 1);
    cblas_sgemv(CblasColMajor, CblasNoTrans, q_size, r_size/2, scale, bess_mat, q_size, gr, 1, 1, sf, 1 );
 }

void dump_positions(char filename[])
{
    ofstream fo(filename);
    fo.write(reinterpret_cast<char *> (&n), sizeof(int));
    for(int i = 0; i < n; i++){
        fo.write(reinterpret_cast<char *>(&(disc_centers[i])), sizeof(position_structure));
    }
    
}
void read_positions(char filename[])
{
    ifstream fi(filename);
    fi.read(reinterpret_cast<char *>(&n), sizeof(int));
    if(disc_centers) delete [] disc_centers;
    disc_centers = new position_structure[n];
    for(int i = 0; i < n; i++){
        fi.read(reinterpret_cast<char *>(&(disc_centers[i])), sizeof(position_structure));
    }
}

int run(float *gr)
{
    
    struct timespec start, end;
//    ios_base::sync_with_stdio(false);

    memset(hist, 0, hist_size*sizeof(int));
    int sum;
    
//    clock_gettime(CLOCK_MONOTONIC, &start);
    cout << n_iteration<<endl;
    move_ntimes(n_iteration);
//    clock_gettime(CLOCK_MONOTONIC, &end);
//    double time_taken = (end.tv_sec - start.tv_sec);
//    cout << "shuffle : " <<  time_taken << endl;
  
    clock_gettime(CLOCK_MONOTONIC, &start);
    sum = accumulate_hist();
    cout <<"sum "<< sum <<endl;
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time_taken = (end.tv_sec - start.tv_sec);
    cout << "accumulate hist : " <<  time_taken << endl;
    double sum2 = sum*((sum-1)*0.5);
    cout <<"sum2 "<< sum2<<endl;
    hist_normalize(normalized_hist, 1./sum2);
//    clock_gettime(CLOCK_MONOTONIC, &start);
    cout << "normalized"<<endl;
    cout<< "hist size "<<hist_size<<endl;
    cout <<gr <<endl;
    calc_gr(normalized_hist, gr);
    cout << "gr_calculated"<< endl;
     
    return 0;
}

double rejection_ratio()
{
    reject_count = 0;
    accept_count = 0;
    move_ntimes(n_iteration);
    double r_ratio =reject_count*1.0/(reject_count + accept_count);
    cout << "rejection ratio "<< r_ratio <<endl;
    return r_ratio;
}

int tune()
{
    n_iteration = n*10;
    while(rejection_ratio() > 0.5) {
        step_factor *=0.8;
        allocate_normal_distribution(step_factor);
    }
    return 0;
}


int main(int argc, char *argv[])
{
    v_hist_bin = 10; //for diagnosis horizontal fluctuation
    h_hist_bin = 10; //for diagnosis vertical fluctuation
    
    init(0.6, 500, 10); // phi, box)size,  r_div = 10; q_step = 0.5/box_size, q_size = box_size*rdiv;
    line_picking_init();
    ofstream fo1("line_picking.dat");
    fo1.write(reinterpret_cast<char *>(line_picking), sizeof(float)*hist_size);
    fo1.close();
    cout << "line picked"<<endl;
    
    random_device dev;
    rng = new mt19937(dev());
    n_iteration = n*10;
    allocate_float(&gr, r_size);
    allocate_normal_distribution(1.);
    tune();
    cout << "step_factor : "<< step_factor <<endl;
    n_iteration = n*10;
    struct timespec start, end;
    for(int k = 0; k < 10; k++){
        clock_gettime(CLOCK_MONOTONIC, &start);
        for(int i = 0; i < 100; i++){
            move_ntimes(n_iteration);
        }
        clock_gettime(CLOCK_MONOTONIC, &end);
        double time_taken = (end.tv_sec - start.tv_sec);
         cout << k <<"/10  "<< time_taken<<endl;
   }
    ofstream fo("gr6");
    for(int i = 0; i < 1000; i++){
        run(gr);
        fo.write(reinterpret_cast<char *>(gr), sizeof(float)*r_size);
    }
    dump_positions("pos6");
    return 0;
}

