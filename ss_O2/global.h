//par
int label;

int N;
float pf0;
float mesh_fac;
int tobs;

//initial conditions
float V_hs;
float amp_fac;

float pf;

//host
float *x_h, *y_h;
float *px_h, *py_h;
float *_2r_h;

int L_mesh, L2_mesh, cell_num;
int num_particle_per_cell;

//########################################################

int block_size, NumBlock; //WARNING, device parameters (host)
int BSF;
int ct, ct2;

__constant__ int dct1; //WARNING, device parameters (device)
__constant__ int dct2;

__device__ float d_h; //time step (device)
__device__ float d_h2;

__constant__ float mesh_d;
__constant__ int L_mesh_d;
__constant__ int L2_mesh_d;
__constant__ int cell_num_d;


__constant__ int num_particle_per_cell_d;
__constant__ int L_d;

__constant__ float amp_fac_d;


__constant__ float f_inc;
__constant__ float f_dec;
__constant__ float f_alpha;
__device__ float alpha;
__device__ int ii_bac;

__constant__ float hmax;
__constant__ float alpha_start;
