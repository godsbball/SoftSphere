#include <cuda_runtime.h>
#include <cooperative_groups.h>

#include <stdio.h>

#include <stdlib.h>
#include <math.h>
#include <time.h>

// for memset
#include  <cstring>

#define SCALE 2.328306436e-10
static unsigned int ux, uy;
double xor64( )
{
    unsigned int t = (ux^(ux<<8));
    ux = uy;
    uy = (uy^(uy>>22))^(t^(t>>9));
    return ((double) uy) * SCALE; 
}        

float h, h2;
static int L = 10;


#include "global.h"
#include "host.c"
#include "ovito_plot.c"

#include "cell.c"
#include "device.c"


int main(int argc, char** argv)
{

    int label;
    sscanf (argv[1],"%d",&label);

    char file_name[50];

    FILE *fin0;
    sprintf(file_name,"par.dat");
    fin0 = fopen(file_name, "r");

    int N_p;
    fscanf(fin0,"%d\n", &N_p);
    fscanf(fin0,"%f\n", &pf0);
    fscanf(fin0,"%f\n", &mesh_fac);
    fscanf(fin0,"%d\n", &tobs);
    float frac_total_part;
    fscanf(fin0,"%f\n", &frac_total_part);
    fscanf(fin0,"%f\n", &h);

    fclose(fin0);


    N = 1 << N_p;

// 	srand( 300101985 );
    int seed = (int)time(NULL);
//     int seed = 124095837;

    srand( seed );
    int ii;
    for (ii=0; ii<100; ++ii) ux = rand();
    for (ii=0; ii<100; ++ii) uy = rand();
    


//#################################################################################
//host_variables_initialization
    int NO_BYTES = N*sizeof(float);
    
    x_h = (float*)malloc(NO_BYTES);
    y_h = (float*)malloc(NO_BYTES);

    px_h = (float*)malloc(NO_BYTES);
    py_h = (float*)malloc(NO_BYTES);

    _2r_h = (float*)malloc(NO_BYTES);

//cell
    L_mesh = L/(float)mesh_fac;

    L2_mesh = L_mesh*L_mesh;
    cell_num = L2_mesh;
    num_particle_per_cell = (int)(frac_total_part*N);

    int NO_BYTES_cell = cell_num*num_particle_per_cell*sizeof(int);
    int NO_BYTES_cell_ind = (cell_num)*sizeof(int);

//     int *cell_ind_h;
//     int *cell_h;
//
//     cell_h = (int*)malloc(NO_BYTES_cell);
//     cell_ind_h = (int*)malloc(NO_BYTES_cell_ind);
//
//     for (ii=0; ii < cell_num; ++ii)
//         cell_ind_h[ii] = 0;




//#################################################################################
//initial conditions importation (host)

    sprintf(file_name,"./input/c_000.dat");
    FILE *fin2;
    fin2 = fopen(file_name, "r");

    for (ii=0; ii < N; ++ii)
        fscanf(fin2,"%f %f %f\n", &x_h[ii], &y_h[ii], &_2r_h[ii]);

    fclose(fin2);


    sprintf(file_name,"./input/par_000.dat");
    FILE *fin3;
    fin3 = fopen(file_name, "r");

    fscanf(fin3,"%d %f", &ii, &amp_fac);

    fclose(fin3);

//#################################################################################
//initial conditions generation (host)

//     r_dist();
//
//     amp_fac = 1.e-3;
//     V_hard_sphere();
//
//
//     while( V_hs/( (float)(L*L)/pow(amp_fac,2) ) < pf0 ){
//         amp_fac += 0.0001;
//     }
//
//     pf =V_hs/( (float)(L*L)/pow(amp_fac,2));
    for (ii=0; ii<N; ++ii){

//         x_h[ii] = (float)(L)*xor64();
//         y_h[ii] = (float)(L)*xor64();

        px_h[ii] =  1e-2*(2.*xor64() - 1);
        py_h[ii] =  1e-2*(2.*xor64() - 1);

    }

    V_hard_sphere();
    pf = V_hs/( (float)(L*L)/pow(amp_fac,2));

    h2 = 0.5*h;

    printf("\n#############################\n");
    printf("N = %d \n",N);
    printf("tobs = %d \n", tobs);
    printf("h = %g \n", h);
    printf("seed = %d \n", seed);
    printf("#############################\n");

    printf("\ncell_num: %d   cell_num*num_particle_per_cell: %d L_mesh: %d\n", cell_num, cell_num*num_particle_per_cell, L_mesh);
    printf("num_particle_per_cell: %d mesh_fac: %g\n\n", num_particle_per_cell, mesh_fac);

    printf("amp_fac: %g  pf: %g\n\n", amp_fac, pf);
//#################################################################################
//device_variables_initialization

    float *x_d, *y_d;
    float *_2r_d;

    int *cell_ind_d;
    int *cell_d;

    float *px_d, *py_d;
    float *fx_d, *fy_d;

    float *u_d;


    cudaMalloc((float **)&x_d, NO_BYTES);
    cudaMemcpy(x_d,x_h, NO_BYTES, cudaMemcpyHostToDevice);

    cudaMalloc((float **)&y_d, NO_BYTES);
    cudaMemcpy(y_d,y_h, NO_BYTES, cudaMemcpyHostToDevice);

    cudaMalloc((float **)&_2r_d, NO_BYTES);
    cudaMemcpy(_2r_d,_2r_h,NO_BYTES, cudaMemcpyHostToDevice);

    cudaMalloc((float **)&px_d, NO_BYTES);
    cudaMemcpy(px_d,px_h, NO_BYTES, cudaMemcpyHostToDevice);

    cudaMalloc((float **)&py_d, NO_BYTES);
    cudaMemcpy(py_d,py_h, NO_BYTES, cudaMemcpyHostToDevice);

//     cudaMalloc((float **)&px_d, NO_BYTES);
//     cudaMemset(px_d, 0, NO_BYTES);
//
//     cudaMalloc((float **)&py_d, NO_BYTES);
//     cudaMemset(py_d, 0, NO_BYTES);



    cudaMalloc((int **)&cell_d, NO_BYTES_cell);
    cudaMalloc((int **)&cell_ind_d, NO_BYTES_cell_ind);

    cudaMemset(cell_d, 0, NO_BYTES_cell);
    cudaMemset(cell_ind_d, 0, NO_BYTES_cell_ind);

    cudaMalloc((float **)&fx_d, NO_BYTES);
    cudaMalloc((float **)&fy_d, NO_BYTES);
    cudaMemset(fx_d, 0, NO_BYTES);
    cudaMemset(fy_d, 0, NO_BYTES);

    cudaMalloc((float **)&u_d, NO_BYTES);
    cudaMemset(u_d, 0, NO_BYTES);


//#################################################################################
//device_structure_definition
    
    device_structure_definition();

    printf("\ngrid: %d,  block: %d,  c1*grid*block = %d,  ct,ct2 = %d,%d,  BSF = %d\n", NumBlock,block_size,ct*NumBlock*block_size,ct,ct2,BSF);



    cudaMemcpyToSymbol(dct1, &ct, sizeof(int));//device parameter
    cudaMemcpyToSymbol(dct2, &ct2, sizeof(int));

    cudaMemcpyToSymbol(mesh_d, &mesh_fac, sizeof(float));
    cudaMemcpyToSymbol(L_mesh_d, &L_mesh, sizeof(int));
    cudaMemcpyToSymbol(num_particle_per_cell_d, &num_particle_per_cell, sizeof(int));

    cudaMemcpyToSymbol(d_h, &h, sizeof(float));
    cudaMemcpyToSymbol(d_h2, &h2, sizeof(float));

    cudaMemcpyToSymbol(L_d, &L, sizeof(int));
    cudaMemcpyToSymbol(L_mesh_d, &L_mesh, sizeof(int));
    cudaMemcpyToSymbol(L2_mesh_d, &L2_mesh, sizeof(int));
    cudaMemcpyToSymbol(cell_num_d, &cell_num, sizeof(int));

    cudaMemcpyToSymbol(amp_fac_d, &amp_fac, sizeof(float));

//     float temp;
//
//     temp = 1.1;
//     cudaMemcpyToSymbol(f_inc, &temp, sizeof(float));
//     temp = 0.5;
//     cudaMemcpyToSymbol(f_dec, &temp, sizeof(float));
//     temp = 0.99;
//     cudaMemcpyToSymbol(f_alpha, &temp, sizeof(float));
//     temp = 0.2;
//     cudaMemcpyToSymbol(alpha, &temp, sizeof(float));
//     cudaMemcpyToSymbol(alpha_start, &temp, sizeof(float));
//
//     ii = 0;
//     cudaMemcpyToSymbol(ii_bac, &ii, sizeof(int));
//
//     temp = 30*h;
//     cudaMemcpyToSymbol(hmax, &temp, sizeof(float));



//#################################################################################
//backup_variables_definition

    int temp_array_byte_size = sizeof(float)*NumBlock;

    float *d_temp;
    cudaMalloc((void**)&d_temp,temp_array_byte_size);
    cudaMemset(d_temp, 0., temp_array_byte_size);

    float *d_temp2;
    cudaMalloc((void**)&d_temp2,temp_array_byte_size);
    cudaMemset(d_temp2, 0., temp_array_byte_size);

    float *d_temp3;
    cudaMalloc((void**)&d_temp3,temp_array_byte_size);
    cudaMemset(d_temp3, 0., temp_array_byte_size);


    float *dE;
    cudaMalloc((void**)&dE, temp_array_byte_size);
    cudaMemset(dE, 0., temp_array_byte_size);

    float *dVec;
    cudaMalloc((void**)&dVec, temp_array_byte_size);
    cudaMemset(dVec, 0., temp_array_byte_size);


//#################################################################################
//#################################################################################
//run_device
//     cudaError_t cudaStatus;

//     int jj;
    float energy_bac0, energy_bac1, energy_bac2;
    energy_bac0 = energy_bac1 = energy_bac2 = 1.;


    for (ii=0; ii<tobs; ++ii){

        EvolveP<<<NumBlock,block_size>>>(x_d, y_d, px_d, py_d, fx_d, fy_d, cell_ind_d, cell_d, _2r_d, u_d);
        cudaDeviceSynchronize();

        EvolveQ<<<NumBlock, block_size>>>(x_d, y_d, px_d, py_d);
        cudaDeviceSynchronize();

        initializeCellInd<<<NumBlock, block_size>>>(cell_ind_d);
        cudaDeviceSynchronize();
        initializeCells<<<NumBlock,block_size>>>(x_d, y_d, cell_ind_d, cell_d);
        cudaDeviceSynchronize();

        EvolveP<<<NumBlock,block_size>>>(x_d, y_d, px_d, py_d, fx_d, fy_d, cell_ind_d, cell_d, _2r_d, u_d);
        cudaDeviceSynchronize();


        if(ii%1000==0){
            Energy_reduction<<<NumBlock,block_size>>>(px_d, py_d, u_d, d_temp, d_temp2);
            Energy_FinalReduction<<<1,BSF>>>(d_temp, d_temp2, dE);
            cudaDeviceSynchronize();

            cudaMemcpy(&energy_bac0, &dE[0], sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&energy_bac1, &dE[1], sizeof(float), cudaMemcpyDeviceToHost);
            cudaMemcpy(&energy_bac2, &dE[2], sizeof(float), cudaMemcpyDeviceToHost);

            printf("\ntemperature: %.3g\t potential_e/N: %.3g\t tot_e: %.3g\t ii: %d\n", energy_bac0/(float)N, energy_bac1/(float)N, energy_bac2, ii);
        }

    }

    Energy_reduction<<<NumBlock,block_size>>>(px_d, py_d, u_d, d_temp, d_temp2);
    Energy_FinalReduction<<<1,BSF>>>(d_temp, d_temp2, dE);
    cudaDeviceSynchronize();

//#################################################################################
//Copy from device to host

    int *cell_ind_h_bac = (int*)malloc(NO_BYTES_cell_ind);
    cudaMemcpy(cell_ind_h_bac, cell_ind_d, NO_BYTES_cell_ind, cudaMemcpyDeviceToHost);

    int *cell_h_bac = (int*)malloc(NO_BYTES_cell);
    cudaMemcpy(cell_h_bac, cell_d, NO_BYTES_cell, cudaMemcpyDeviceToHost);

    cudaMemcpy(x_h, x_d, NO_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(y_h, y_d, NO_BYTES, cudaMemcpyDeviceToHost);

    cudaMemcpy(px_h, px_d, NO_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(py_h, py_d, NO_BYTES, cudaMemcpyDeviceToHost);

    float *u_h = (float*)malloc(NO_BYTES);
    cudaMemcpy(u_h, u_d, NO_BYTES, cudaMemcpyDeviceToHost);

    float *dh_cpy = (float*)malloc(temp_array_byte_size);
    cudaMemcpy(dh_cpy, dE, temp_array_byte_size, cudaMemcpyDeviceToHost);
    printf("\nenergy_device: %g\n", dh_cpy[2]);


    cudaFree(cell_ind_d);
    cudaFree(cell_d);

    cudaFree(x_d);
    cudaFree(y_d);

    cudaFree(px_d);
    cudaFree(py_d);

    cudaFree(fx_d);
    cudaFree(fy_d);

    cudaFree(dE);
    cudaFree(u_d);

    cudaFree(dVec);

    float E_check = 0;
    for (ii=0; ii < N; ++ii){


        E_check += 0.5*(px_h[ii]*px_h[ii] + py_h[ii]*py_h[ii]);
        E_check += u_h[ii];
    }

    printf("energy_host: %g\n", E_check);


//#################################################################################
    printf("\n\n############### device task concluded\n\n");
//plot and check


    sprintf(file_name,"test_h.dat");
    FILE *fout0;
    fout0 = fopen(file_name, "w");

    for (ii=0; ii < N; ++ii)
        fprintf(fout0,"%g %g %g\n", px_h[ii], py_h[ii], u_h[ii]);

    fclose(fout0);
    ovito_plot(0);

    return 0;
}
