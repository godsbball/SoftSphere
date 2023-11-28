//######################################################################################
//declaration of constants
int N;
int L_mesh, L2_mesh, cell_num;
int t_rlx, tobs, frame; //input

double h;

//######################################################################################
//Structure definition
typedef struct {
    double q0;
    double q1;
} Time;

typedef struct {
    double x0;
    double y0;
    double z0;


    double x;
    double y;
    double z;

    double vx;
    double vy;
    double vz;

    double _2r;
    double Fx;
    double Fy;
    double Fz;


} Particle;

typedef struct {
    Particle* particles;
    Time* time;

    int* cell_ind;
    int** cell;

    double V_hs;
    double amp_fac;

    double e;
    double K;
    double K_prime;
    double V;
    double P;
    double* delta;


} Sample;

typedef struct {
    unsigned int ux, uy;
} Seed;
