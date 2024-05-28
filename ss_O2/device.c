__global__ void initializeCellInd(int *cell_ind) {


    int tid = threadIdx.x;
    int ii = (blockIdx.x * blockDim.x)*dct1;

    int jj;

    int idC;

    for (jj=0; jj<dct1; jj++)
    {
        idC = ii+jj*blockDim.x+tid;

        if(idC < cell_num_d)
            cell_ind[idC] = 0;

    }

}

__global__ void initializeCells(float *x, float *y, int *cell_ind, int *cell) {


    int tid = threadIdx.x;
    int ii = (blockIdx.x * blockDim.x)*dct1;

    int jj;

    int idP, idC, index_in_cell;

    for (jj=0; jj<dct1; jj++)
    {
        idP = ii+jj*blockDim.x+tid;
        idC = (int)(x[idP] / mesh_d) + L_mesh_d * (int)(y[idP] / mesh_d);

        index_in_cell = atomicAdd(&cell_ind[idC], 1);

        cell[idC * num_particle_per_cell_d + index_in_cell] = idP;

    }

}


__global__ void EvolveQ(float *x, float *y, float *px, float *py)
{

    int tid = threadIdx.x;
    int ii = (blockIdx.x * blockDim.x)*dct1;

    int jj, idP;

    float eps_bac = 1.e-6;

    for (jj=0; jj < dct1; jj++)
    {

        idP = ii+jj*blockDim.x+tid;

        x[idP] += d_h*px[idP];
        y[idP] += d_h*py[idP];

        while(x[idP] >= (float)L_d)
            x[idP] -= (float)L_d;

        while(x[idP] < 0. - eps_bac)
            x[idP] += (float)L_d;

        while(y[idP] >= (float)L_d)
            y[idP] -= (float)L_d;

        while(y[idP] < 0. - eps_bac)
            y[idP] += (float)L_d;
    }
}




__global__ void EvolveP(float *x, float *y, float *px, float *py, float *Fx, float *Fy, int *cell_ind, int *cell, float *_2r, float *U) {
    int tid = threadIdx.x;
    int ii = (blockIdx.x * blockDim.x) * dct1;

    int jj, kk, ll, idP, idC, II, list[4];

    float d, D, fx, fy, delta_x, delta_y, u;

    for (jj = 0; jj < dct1; jj++) {
        idP = ii + jj * blockDim.x + tid;
        idC = (int)(x[idP] / mesh_d) + L_mesh_d * (int)(y[idP] / mesh_d);

        fx = fy = u = 0.0f;

        for (kk = 0; kk < 9; ++kk) {
            neighbor_cell_list(kk, idC, list);

            for (ll = 0; ll < cell_ind[list[0]]; ll++) {
                II = cell[list[0] * num_particle_per_cell_d + ll];

                //list[1] and list[2]: spatial shift for PBC.
                delta_x = x[idP] + (float)list[1] - x[II];
                delta_y = y[idP] + (float)list[2] - y[II];

                d = sqrt(delta_x * delta_x + delta_y * delta_y);
                D = 0.5 * amp_fac_d * (_2r[idP] + _2r[II]);

                if (d < D && II != idP) {
                    fx += -(1.0f - d / D) * (-delta_x) / (d * D);
                    fy += -(1.0f - d / D) * (-delta_y) / (d * D);
                    u += 0.5 * (1.0f - d / D) * (1.0f - d / D);
                }
            }
        }

        // Store total force in global memory
        Fx[idP] = fx;
        Fy[idP] = fy;

        // Other computations and momentum updates can be done here...

        U[idP] = u;

        px[idP] += d_h2*fx;
        py[idP] += d_h2*fy;
    }
}

__global__ void Energy_reduction(float *input_x, float *input_y, float *input_z,  float *output_1, float *output_2)
{
    int tid = threadIdx.x;
    int ii = (blockIdx.x * blockDim.x)*dct1;

    __shared__ float i_data_x[512];
    __shared__ float i_data_y[512];

    i_data_x[tid] = 0;
    i_data_y[tid] = 0;

    int idP, jj;

    for (jj = 0; jj < dct1; ++jj){

        idP = ii + jj*blockDim.x + tid;

        i_data_x[tid] += input_x[idP]*input_x[idP] + input_y[idP]*input_y[idP];
        i_data_y[tid] += input_z[idP];
    }

    __syncthreads();


    if(tid < 256){
        i_data_x[tid] += i_data_x[tid + 256];
        i_data_y[tid] += i_data_y[tid + 256];
    }
    __syncthreads();

    if(tid < 128){
        i_data_x[tid] += i_data_x[tid + 128];
        i_data_y[tid] += i_data_y[tid + 128];
    }
        __syncthreads();

    if(tid < 64){
        i_data_x[tid] += i_data_x[tid + 64];
        i_data_y[tid] += i_data_y[tid + 64];
    }
    __syncthreads();


    if(tid < 32)
    {
        volatile float *vsmem_x = i_data_x;
        volatile float *vsmem_y = i_data_y;

        vsmem_x[tid] += vsmem_x[tid+32];
        vsmem_x[tid] += vsmem_x[tid+16];
        vsmem_x[tid] += vsmem_x[tid+8];
        vsmem_x[tid] += vsmem_x[tid+4];
        vsmem_x[tid] += vsmem_x[tid+2];
        vsmem_x[tid] += vsmem_x[tid+1];

        vsmem_y[tid] += vsmem_y[tid+32];
        vsmem_y[tid] += vsmem_y[tid+16];
        vsmem_y[tid] += vsmem_y[tid+8];
        vsmem_y[tid] += vsmem_y[tid+4];
        vsmem_y[tid] += vsmem_y[tid+2];
        vsmem_y[tid] += vsmem_y[tid+1];

    }

    if( tid == 0 ){
        output_1[blockIdx.x] = i_data_x[0];
        output_2[blockIdx.x] = i_data_y[0];
    }
}

__global__ void Energy_FinalReduction(float *x, float *y, float *dmag)
{
    __shared__ float mx[512];
    __shared__ float my[512];

    int jj;
    int tid = threadIdx.x;

    mx[tid] = 0;
    my[tid] = 0;

    for (jj=0; jj < dct2; ++jj){
        mx[tid] += x[jj*blockDim.x+tid];
        my[tid] += y[jj*blockDim.x+tid];
    }
    __syncthreads();


    if (blockDim.x > 256 && tid < 256){
        mx[tid] += mx[tid+256];
        my[tid] += my[tid+256];
    }
    __syncthreads();

    if (blockDim.x > 128 && tid < 128){
        mx[tid] += mx[tid+128];
        my[tid] += my[tid+128];
    }
    __syncthreads();

    if (blockDim.x > 64 && tid < 64){
        mx[tid] += mx[tid+64];
        my[tid] += my[tid+64];
    }
    __syncthreads();

    if (tid < 32)
    {
        volatile float *vx = mx;
        volatile float *vy = my;

        if (blockDim.x > 32){
            vx[tid] += vx[tid + 32];
            vy[tid] += vy[tid + 32];
        }

        if (blockDim.x > 16){
            vx[tid] += vx[tid + 16];
            vy[tid] += vy[tid + 16];
        }

        if (blockDim.x > 8){
            vx[tid] += vx[tid +  8];
            vy[tid] += vy[tid +  8];
        }

        if (blockDim.x > 4){
            vx[tid] += vx[tid +  4];
            vy[tid] += vy[tid +  4];
        }


        if (blockDim.x > 2){
            vx[tid] += vx[tid +  2];
            vy[tid] += vy[tid +  2];
        }
        vx[tid] += vx[tid +  1];
        vy[tid] += vy[tid +  1];
    }

    if (tid == 0)
    {
        dmag[0] = 0.5*mx[0];
        dmag[1] = my[0];
        dmag[2] = 0.5*mx[0]+my[0];

    }
}









