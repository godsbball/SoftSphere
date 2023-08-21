//Hard Sphere System, Monte Carlo
//L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023
//**************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#define SCALE 2.328306436e-10

static int t_rlx = (int)0;
static int tobs = (int)1.5e6;

static int num_threads = 1;


//######################################################################################
//declaration of constants
// int N = 1000;
// int L = 10;
// double mesh_fac = 2.5;
// double amp_fac0 = 0.2;

int N;
int L;
double mesh_fac;
double amp_fac0;

int L2;
int L_mesh, L2_mesh, cell_num;

double k = 5.25;


double h;
double s2;

double alpha1, alpha2, alpha3, alpha4;

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
    double V;

} Sample;

typedef struct {
    unsigned int ux, uy;
} Seed;
//######################################################################################

//random number generator
double xor64(int ii, Seed* seeds ){
    unsigned int t = (seeds[ii].ux^(seeds[ii].ux<<8));
    seeds[ii].ux = seeds[ii].uy;
    seeds[ii].uy = (seeds[ii].uy^(seeds[ii].uy>>22))^(t^(t>>9));
    return ((double) seeds[ii].uy) * SCALE;
}
//######################################################################################

//Paralellized functions//######################################################################################


void neighbor_cell_list(int label_bac, int index, int list[]){

    int bac, sx, sy;

    switch (label_bac%9){

        case 0:
            list[0] = index;
            list[1] = 0;
            list[2] = 0;
            list[3] = 0;

            if( (int)(label_bac/9 )==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }
            break;

        case 1:
            bac = index - 1; //left
            sx = sy = 0;
            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 2:
            bac = index + 1; //right
            sx = sy = 0;
            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 3:
            bac = index - L_mesh; //bottom
            sx = sy = 0;
            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 4:
            bac = index + L_mesh; //top
            sx = sy = 0;
            if( index%L2_mesh >= L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 5:
            bac = index + L_mesh - 1; //top/left
            sx = sy = 0;

            if(index%L2_mesh>=L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }


            break;

        case 6:
            bac = index + L_mesh + 1; //top/right
            sx = sy = 0;

            if(index%L2_mesh>=L_mesh*(L_mesh-1)){
                bac -= L2_mesh;
                sy = -L;
            }

            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 7:
            bac = index - L_mesh - 1; //bottom/left
            sx = sy = 0;

            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            if(index%L_mesh==0){
                bac += L_mesh;
                sx = L;
            }


            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;

        case 8:
            bac = index - L_mesh + 1; //bottom/rigth
            sx = sy = 0;


            if(index%L2_mesh<L_mesh){
                bac += L2_mesh;
                sy = L;
            }

            if( (index+1)%L_mesh==0 ){
                bac -= L_mesh;
                sx = -L;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh;

                if( list[0] >= cell_num ){

                    list[0] -= cell_num;
                    list[3] = -L;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh;

                if( list[0] < 0 ){

                    list[0] += cell_num;
                    list[3] = L;
                }
            }

            break;
    }

    return;
}


void energy(int thread_id, Sample* samples){

    int ii,jj,kk,ll,II;
    double d, delta_x, delta_y, delta_z, D;
    double K_bac, V_bac;

    int list[4];

    double x_bac, y_bac, z_bac;


    K_bac = V_bac = 0;

    for(ii=0; ii<N; ++ii){


        K_bac += samples[thread_id].particles[ii].vx*samples[thread_id].particles[ii].vx + samples[thread_id].particles[ii].vy*samples[thread_id].particles[ii].vy + samples[thread_id].particles[ii].vz*samples[thread_id].particles[ii].vz;

        x_bac = samples[thread_id].particles[ii].x;
        y_bac = samples[thread_id].particles[ii].y;
        z_bac = samples[thread_id].particles[ii].z;


        ll = (int)( x_bac/mesh_fac ) + L_mesh*(int)( y_bac/mesh_fac ) + L2_mesh*(int)( z_bac/mesh_fac );

        for(kk=0; kk<27; ++kk){

            neighbor_cell_list(kk, ll, list);

            for(jj=0; jj < samples[thread_id].cell_ind[list[0]]; ++jj){

                II = samples[thread_id].cell[ list[0] ][jj];

                delta_x = x_bac + list[1] - samples[thread_id].particles[II].x;
                delta_y = y_bac + list[2] - samples[thread_id].particles[II].y;
                delta_z = z_bac + list[3] - samples[thread_id].particles[II].z;


                d = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z  );
                D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[II]._2r);

                if( d < D && II!=ii )
                    V_bac+= 0.5*(1. - d/D)*(1. - d/D);
            }
        }
    }

    K_bac *= 0.5;
//     V_bac /= (double)N;

    samples[thread_id].K = K_bac;
    samples[thread_id].V = V_bac;
    samples[thread_id].e = K_bac+V_bac;

    return;
}

void evolve_Q(int thread_id, Sample* samples, double ALPHA){

    int jj,kk,ll,mm;

    for(jj=0; jj < N; ++jj){


        kk = (int)(samples[thread_id].particles[jj].x/mesh_fac)+L_mesh*(int)(samples[thread_id].particles[jj].y/mesh_fac)+L2_mesh*(int)(samples[thread_id].particles[jj].z/mesh_fac);

//         printf("%g %g\n", ALPHA*samples[thread_id].particles[jj].vx, samples[thread_id].particles[jj].vx);


        samples[thread_id].particles[jj].x += ALPHA*samples[thread_id].particles[jj].vx;
        samples[thread_id].particles[jj].y += ALPHA*samples[thread_id].particles[jj].vy;
        samples[thread_id].particles[jj].z += ALPHA*samples[thread_id].particles[jj].vz;


        while(samples[thread_id].particles[jj].x < 0.)
            samples[thread_id].particles[jj].x += L;

        while(samples[thread_id].particles[jj].x >= L)
            samples[thread_id].particles[jj].x -= L;

        while(samples[thread_id].particles[jj].y < 0.)
            samples[thread_id].particles[jj].y += L;

        while(samples[thread_id].particles[jj].y >= L)
            samples[thread_id].particles[jj].y -= L;

        while(samples[thread_id].particles[jj].z < 0.)
            samples[thread_id].particles[jj].z += L;

        while(samples[thread_id].particles[jj].z >= L)
            samples[thread_id].particles[jj].z -= L;

        ll = (int)(samples[thread_id].particles[jj].x/mesh_fac)+L_mesh*(int)(samples[thread_id].particles[jj].y/mesh_fac)+L2_mesh*(int)(samples[thread_id].particles[jj].z/mesh_fac);

        if( kk != ll ){

            for(mm=0; mm < samples[thread_id].cell_ind[kk]; ++mm)
                if( samples[thread_id].cell[kk][mm] == jj )
                    break;

            samples[thread_id].cell_ind[kk] -= 1;
            samples[thread_id].cell[kk][mm] = samples[thread_id].cell[kk][ samples[thread_id].cell_ind[kk] ];

            samples[thread_id].cell[ll][ samples[thread_id].cell_ind[ll] ] = jj;
            samples[thread_id].cell_ind[ll] += 1;
        }
    }

    return;
}

double findMin(double num1, double num2) {

    if (num1 < num2)
        return num1;

    else
        return num2;
}

void evolve_P(int thread_id, Sample* samples, double ALPHA){


    double fx, fy, fz;
    int ii, jj,kk,ll, II;

    int list[4];
    double x_bac, y_bac, z_bac;
    double d, delta_x, delta_y, delta_z, D, d_prime;

    for(ii=0; ii<N; ++ii){

        x_bac = samples[thread_id].particles[ii].x;
        y_bac = samples[thread_id].particles[ii].y;
        z_bac = samples[thread_id].particles[ii].z;

        ll = (int)( x_bac/mesh_fac ) + L_mesh*(int)( y_bac/mesh_fac ) + L2_mesh*(int)( z_bac/mesh_fac );

        fx = fy = fz = 0;
        for(kk=0; kk<27; ++kk){

            neighbor_cell_list(kk, ll, list);

            for(jj=0; jj < samples[thread_id].cell_ind[list[0]]; ++jj){

                II = samples[thread_id].cell[ list[0] ][jj];

                delta_x = x_bac + list[1] - samples[thread_id].particles[II].x;
                delta_y = y_bac + list[2] - samples[thread_id].particles[II].y;
                delta_z = z_bac + list[3] - samples[thread_id].particles[II].z;

                delta_x = findMin( (double)L-delta_x, delta_x );
                delta_y = findMin( (double)L-delta_y, delta_y );
                delta_z = findMin( (double)L-delta_z, delta_z );

                d = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z  );
                D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[II]._2r);

                if( d < D && II!=ii ){

                    d_prime = sqrt( delta_x*delta_x + delta_y*delta_y );

                    fx += -( 1. - d/D )*(-delta_x)/d;
                    fy += -( 1. - d/D )*(-delta_y)/d;
                    fz += -( 1. - d/D )*(-delta_z)/d;
                }
            }
        }

        samples[thread_id].particles[ii].vx += ALPHA*fx;
        samples[thread_id].particles[ii].vy += ALPHA*fy;
        samples[thread_id].particles[ii].vz += ALPHA*fz;
    }

    return;
}

void  time_step(int thread_id, Sample* samples){

    evolve_Q(thread_id, samples, alpha1);

    evolve_P(thread_id, samples, alpha2);
    evolve_Q(thread_id, samples, alpha2);
    evolve_P(thread_id, samples, alpha2);

    evolve_Q(thread_id, samples, alpha3);
    evolve_P(thread_id, samples, alpha4);
    evolve_Q(thread_id, samples, alpha3);

    evolve_P(thread_id, samples, alpha2);
    evolve_Q(thread_id, samples, alpha2);
    evolve_P(thread_id, samples, alpha2);

    evolve_Q(thread_id, samples, alpha1);

    return;
}

//######################################################################################
//Parallel Tempering

//######################################################################################
//Utilities

void testz(int thread_id, Sample* samples){

    double d, delta_x, delta_y, delta_z, D;
    int ii,jj;

    for(ii=0; ii < N; ++ii)
        for(jj=ii+1; jj < N; ++jj){

            delta_x = samples[thread_id].particles[ii].x - samples[thread_id].particles[jj].x;
            delta_y = samples[thread_id].particles[ii].y - samples[thread_id].particles[jj].y;
            delta_z = samples[thread_id].particles[ii].z - samples[thread_id].particles[jj].z;

            d = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
            d = sqrt(d);

            D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[jj]._2r);

            if(d < D){
                printf("Overlap Warning: d=%lf, D = %lf, ii = %d jj = %d sample = %d\n", d, D, ii, jj, thread_id);
                return;

            }
        }

    printf("No Overlap detected: sample = %d\n", thread_id);
    return;
}


void ovito_plot(int thread_id, Sample* samples, int label, int tt){

    int ii;

    char file_name[50];


    FILE *fout_conf;
    sprintf(file_name,"./output/conf/c_%02d_%02d.xyz", label, thread_id);
    fout_conf = fopen(file_name, "a");

    fprintf( fout_conf,"ITEM: TIMESTEP\n");
    fprintf( fout_conf,"%d\n", tt);

    fprintf( fout_conf,"ITEM: NUMBER OF ATOMS\n");
    fprintf( fout_conf,"%d\n", N);

    fprintf( fout_conf,"ITEM: BOX BOUNDS pp pp pp\n");
    fprintf( fout_conf,"%d %lf\n", 0,(double)L/samples[thread_id].amp_fac);
    fprintf( fout_conf,"%d %lf\n", 0,(double)L/samples[thread_id].amp_fac);
    fprintf( fout_conf,"%d %lf\n", 0,(double)L/samples[thread_id].amp_fac);

    fprintf(fout_conf, "ITEM: ATOMS id type x y z radius Transparency\n");

    for(ii=0; ii < N; ++ii)
        fprintf( fout_conf,"%d %d %lf %lf %lf %lf %lf\n", ii, 1, samples[thread_id].particles[ii].x/samples[thread_id].amp_fac, samples[thread_id].particles[ii].y/samples[thread_id].amp_fac, samples[thread_id].particles[ii].z/samples[thread_id].amp_fac, 0.5*samples[thread_id].particles[ii]._2r, 0.2);

    fclose(fout_conf);

    FILE *fout_conf2;
    sprintf(file_name,"./output/conf/c_%02d_%02d.dat", label, thread_id);
    fout_conf2 = fopen(file_name, "w");

    for(ii=0; ii < N; ++ii)
        fprintf( fout_conf2,"%lf %lf %lf %lf\n", samples[thread_id].particles[ii].x, samples[thread_id].particles[ii].y, samples[thread_id].particles[ii].z, samples[thread_id].particles[ii]._2r);

    fclose(fout_conf2);


    FILE *fout_conf3;
    sprintf(file_name,"./output/conf/par_%02d_%02d.dat", label, thread_id);
    fout_conf3 = fopen(file_name, "w");

    fprintf( fout_conf3,"%d %d %lf %lf", N, L, samples[thread_id].amp_fac, mesh_fac);

    fclose(fout_conf3);

}

double EA_order_parameter(int ll, int mm, Sample* samples){

    int ii, idA, idB;
    double delta_x, delta_y, delta_z;

    double delta = 0.;;

    for(ii=0; ii < N; ++ii){

        delta_x = samples[ll].particles[ii].x - samples[mm].particles[ii].x;
        delta_y = samples[ll].particles[ii].y - samples[mm].particles[ii].y;
        delta_z = samples[ll].particles[ii].z - samples[mm].particles[ii].z;

        delta += delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;

    }

    return delta/(double)N;
}

void scatt_f(int thread_id, Sample* samples, int time_bac){

    double dr_norm, dx, dy, dz;
    int jj, kk;

    samples[thread_id].time[time_bac].q1 = 0.;


    for (int jj = 0; jj < N; jj++){

        dx = samples[thread_id].particles[jj].x - samples[thread_id].particles[jj].x0;
        dy = samples[thread_id].particles[jj].y - samples[thread_id].particles[jj].y0;
        dz = samples[thread_id].particles[jj].z - samples[thread_id].particles[jj].z0;

        dr_norm  = dx*dx + dy*dy + dz*dz;
        dr_norm = sqrt(dr_norm);

        if (dr_norm > 0)
            samples[thread_id].time[time_bac].q1 += sin(k*dr_norm/samples[thread_id].particles[jj]._2r)/(k*dr_norm/samples[thread_id].particles[jj]._2r);

        else
            samples[thread_id].time[time_bac].q1 += 1;
    }

    samples[thread_id].time[time_bac].q1 /= (double)N;
}

void V_hard_sphere(int thread_id, Sample* samples){

    int ii;

    samples[thread_id].V_hs = 0.;
    for(ii=0; ii < N; ++ii)
        samples[thread_id].V_hs += pow(samples[thread_id].particles[ii]._2r,3);
    samples[thread_id].V_hs *= (4./3.)*M_PI/8.;

    printf("sample: %d,\t\tPF = %lf,\t\tamp_fac = %lf\n", thread_id, samples[thread_id].V_hs/(L2*L/pow(samples[thread_id].amp_fac,3)), samples[thread_id].amp_fac);

}

//###########################################################################

int main(int argc, char** argv) {

    int label;
    sscanf (argv[1],"%d",&label);

    h = 1.e-5;
    s2 = 1.0/(4.0-pow(4.0,1.0/3.0));

    alpha1 = 0.5*h*s2;
    alpha2 = s2*h;
    alpha3 = 0.5*(1.0-3.0*s2)*h;
    alpha4 = (1.0-4.0*s2)*h;

    char file_name[50];

    double bac_double;

    FILE *fin;
    sprintf(file_name,"par.dat");
    fin = fopen(file_name, "r");
    fscanf(fin,"%d %d %lf %lf %lf %lf\n", &N, &L, &amp_fac0, &mesh_fac, &bac_double, &bac_double);

    fclose(fin);;

    L2 = L*L;

    L_mesh = L/mesh_fac;
    L2_mesh = L_mesh*L_mesh;

    cell_num = L2_mesh*L_mesh;

//###########################################################################
//initial conditions definitions

//loop counters
    int ii,jj,kk;

//structures initialization
    Sample* samples = (Sample*)malloc(num_threads * sizeof(Sample));
    Seed* seeds = (Seed*)malloc( (num_threads+1) * sizeof(Seed));


    for (ii = 0; ii < num_threads; ii++) {
        samples[ii].particles = (Particle*)malloc(N * sizeof(Particle));
        samples[ii].time = (Time*)malloc(tobs * sizeof(Time));

        samples[ii].amp_fac = amp_fac0;

        samples[ii].V = L2*L;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }
    }

//seeds random numbers settings
        srand( (int)time(NULL) );
//     srand( 27278911 );

    for (jj = 0; jj < num_threads+1; jj++) {
        for (ii=0; ii<100; ii++) seeds[jj].ux = rand();
        for (ii=0; ii<100; ii++) seeds[jj].uy = rand();
    }

//initial configuration setting,

    FILE *fin2;
    sprintf(file_name,"xyz.dat");
    fin2 = fopen(file_name, "r");

    for (ii = 0; ii < N; ii++)
        fscanf(fin2,"%lf %lf %lf %lf\n", &samples[0].particles[ii].x0, &samples[0].particles[ii].y0, &samples[0].particles[ii].z0, &samples[0].particles[ii]._2r);

    fclose(fin2);


    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x = samples[0].particles[jj].x0;
            samples[ii].particles[jj].y = samples[0].particles[jj].y0;
            samples[ii].particles[jj].z = samples[0].particles[jj].z0;

            samples[ii].particles[jj]._2r = samples[0].particles[jj]._2r;
        }

    for(ii=0; ii < N; ++ii){

        jj = (int)(samples[0].particles[ii].x/mesh_fac)+L_mesh*(int)(samples[0].particles[ii].y/mesh_fac)+L2_mesh*(int)(samples[0].particles[ii].z/mesh_fac);

        samples[0].cell[ jj ][ samples[0].cell_ind[jj] ] = ii;
        samples[0].cell_ind[jj] += 1;
    }


    for (ii = 1; ii < num_threads; ii++)
        for (jj = 0; jj < cell_num; jj++)
            samples[ii].cell_ind[jj] = samples[0].cell_ind[jj];


    for (ii = 1; ii < num_threads; ii++)
        for (jj = 0; jj < cell_num; jj++)
            for (kk = 0; kk < N; kk++)
                samples[ii].cell[jj][kk] = samples[0].cell[jj][kk];

//seting velocity distribution
    bac_double = 0.;
    double A = 1.;
    for(ii=0; ii < N; ++ii){

//         samples[0].particles[ii].vx = samples[0].particles[ii].vy = samples[0].particles[ii].vz = 0.;

        samples[0].particles[ii].vx = A*(2.*xor64(0,seeds) - 1.);
        samples[0].particles[ii].vy = 0.;
        samples[0].particles[ii].vz = 0.;

        bac_double += pow( samples[0].particles[ii].vx, 2. ) + pow( samples[0].particles[ii].vy, 2. ) + pow( samples[0].particles[ii].vz, 2. );
    }

//reseting counters


//reseting counters

//##################################
//#Dynamics

    printf("\n#############################\n");
    printf("label = %d, N = %d\n\n", label, N);
    printf("V = %d, mesh_fac = %g, cell_num = %d\n", L2*L, mesh_fac,cell_num);

    printf("num_threads = %d\n\n", num_threads);

    printf("tobs = %.1e, t_rlx = %d, \n", (double)tobs, t_rlx);
    printf("#############################\n\n");

    jj = 0;
    for (ii = 0; ii < tobs; ii++){

        time_step( 0., samples);

        if(ii%5000==0){
            energy(0,samples);
            samples[0].time[jj].q0 = samples[0].e;
            ovito_plot(0, samples, label, jj);
            jj += 1;
        }
    }


//##################################
//plot

    printf("\n");

    FILE *fout0;
    sprintf(file_name,"./output/energy_%.2d.dat",label);
    fout0 = fopen(file_name, "w");

    for (ii = 0; ii < jj; ii++)
        fprintf(fout0,"%lf\n", samples[0].time[ii].q0);

    fclose(fout0);

    return 0;
}

