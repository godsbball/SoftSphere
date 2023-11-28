//Soft Sphere System, MD
//L.H.Miranda-Filho., ITP-CAS, lucmiranda@gmail.com, 2023
//**************************************************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//######################################################################################
//Variables definitio//#################################################################


static int num_threads = 1;

static int L = 10;
static int L2 = 100;

// static double h = 5.e-3;

static double mesh_fac = 2.5;
double eps;

#include "global.h"

#include "cell.c"
#include "measurements.c"
#include "integrator.c"
//######################################################################################

//random number generator
#define SCALE 2.328306436e-10

double xor64(int ii, Seed* seeds ){
    unsigned int t = (seeds[ii].ux^(seeds[ii].ux<<8));
    seeds[ii].ux = seeds[ii].uy;
    seeds[ii].uy = (seeds[ii].uy^(seeds[ii].uy>>22))^(t^(t>>9));
    return ((double) seeds[ii].uy) * SCALE;
}
//######################################################################################

//Paralellized functions//######################################################################################

void  time_step(int thread_id, Sample* samples){

    //non-linear eq. integration
    evolve_Q(thread_id, samples);
    evolve_P(thread_id, samples);
    evolve_Q(thread_id, samples);

    //linearized eq. integration
    evolve_DELTA_Q(thread_id, samples);
    evolve_DELTA_P(thread_id, samples);
    evolve_DELTA_Q(thread_id, samples);

    return;
}


void  time_wait(int thread_id, Sample* samples){

    //non-linear eq. integration
    evolve_Q(thread_id, samples);
    evolve_P(thread_id, samples);
    evolve_Q(thread_id, samples);

    return;
}

void  time_wait2(int thread_id, Sample* samples){

    //non-linear eq. integration
    evolve_Q(thread_id, samples);
    evolve_P_relax(thread_id, samples);
    evolve_Q(thread_id, samples);

    return;
}

#include "utilities.c"


//###########################################################################

int main(int argc, char** argv) {

    int label;
    sscanf (argv[1],"%d",&label);


    char file_name[50];

    FILE *fin;
    sprintf(file_name,"par.dat");
    fin = fopen(file_name, "r");

    double pf0;
    fscanf(fin,"%d\n", &N);
    fscanf(fin,"%lf\n", &pf0);
    fscanf(fin,"%d\n", &tobs);
    int tw;
    fscanf(fin,"%d\n", &tw);
    fscanf(fin,"%d\n", &frame);
    fscanf(fin,"%lf\n", &eps);
    double h_bac;
    fscanf(fin,"%lf\n", &h_bac);
    double T0;
    fscanf(fin,"%lf\n", &T0);


    fclose(fin);

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
//         samples[ii].time = (Time*)malloc(t_rlx/frame * sizeof(Time));

        samples[ii].amp_fac = 1.e-4;

        samples[ii].V = L2*L;

        samples[ii].cell_ind = (int*)malloc(cell_num * sizeof(int));
        samples[ii].cell = (int**)malloc(cell_num * sizeof(int*));
        for (int jj = 0; jj < cell_num; jj++) {
            samples[ii].cell[jj] = (int*)malloc(N * sizeof(int));
        }

        samples[ii].delta = (double*)malloc(6*N * sizeof(double));

    }

//seeds random numbers settings

    int seed_int = (int)time(NULL);
//         srand( seed_int );
//     srand( 1700636963 );

    for (jj = 0; jj < num_threads+1; jj++) {
        for (ii=0; ii<100; ii++) seeds[jj].ux = rand();
        for (ii=0; ii<100; ii++) seeds[jj].uy = rand();
    }

//initial configuration setting,


    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < N; jj++){
            samples[ii].particles[jj].x = L*xor64(0,seeds);
            samples[ii].particles[jj].y = L*xor64(0,seeds);
            samples[ii].particles[jj].z = 0.;

            samples[ii].particles[jj].vx = 0.;
            samples[ii].particles[jj].vy = 0.;
            samples[ii].particles[jj].vz = 0.;

            if(xor64(0,seeds)>0.5)
                samples[ii].particles[jj]._2r = 1.6666666666666667;

            else
                samples[ii].particles[jj]._2r = 0.8333333333333333;

        }


    for(ii=0; ii < N; ++ii){

        jj = (int)(samples[0].particles[ii].x/mesh_fac)+L_mesh*(int)(samples[0].particles[ii].y/mesh_fac)+L2_mesh*(int)(samples[0].particles[ii].z/mesh_fac);

        samples[0].cell[ jj ][ samples[0].cell_ind[jj] ] = ii;
        samples[0].cell_ind[jj] += 1;
    }


    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < cell_num; jj++)
            samples[ii].cell_ind[jj] = samples[0].cell_ind[jj];


    for (ii = 0; ii < num_threads; ii++)
        for (jj = 0; jj < cell_num; jj++)
            for (kk = 0; kk < N; kk++)
                samples[ii].cell[jj][kk] = samples[0].cell[jj][kk];


    V_hard_sphere(0, samples);

    for (ii = 0; ii < num_threads; ii++)
        while( samples[ii].V_hs/(L2/pow(samples[ii].amp_fac,2)) < pf0 ){

            V_hard_sphere(0, samples);
            samples[ii].amp_fac += 0.0001;
        }

//seting velocity distribution

    double x_bac, y_bac, z_bac;

//seting tangent space
    x_bac = y_bac = z_bac = 0.;

    for(ii=0; ii < N; ++ii){

        //delta position
        samples[0].delta[ii] = 2.*xor64(0,seeds) - 1.;
        samples[0].delta[ii+N] = 2.*xor64(0,seeds) - 1.;
        samples[0].delta[ii+2*N] = 0.;


        //delta velocity
        samples[0].delta[ii+3*N] = 2.*xor64(0,seeds) - 1.;
        samples[0].delta[ii+4*N] = 2.*xor64(0,seeds) - 1.;
        samples[0].delta[ii+5*N] = 0.;

        x_bac += samples[0].delta[ii]*samples[0].delta[ii]         + samples[0].delta[ii+3*N]*samples[0].delta[ii+3*N];
        x_bac += samples[0].delta[ii+N]*samples[0].delta[ii+N]     + samples[0].delta[ii+4*N]*samples[0].delta[ii+4*N];
        x_bac += samples[0].delta[ii+2*N]*samples[0].delta[ii+2*N] + samples[0].delta[ii+5*N]*samples[0].delta[ii+5*N];

    }

    x_bac = sqrt(x_bac);

    for(ii=0; ii < N; ++ii){

        //delta position
        samples[0].delta[ii] /= x_bac;
        samples[0].delta[ii+N] /= x_bac;
//         samples[0].delta[ii+2*N] /= x_bac;


        //delta velocity
        samples[0].delta[ii+3*N] /= x_bac;
        samples[0].delta[ii+4*N] /= x_bac;
//         samples[0].delta[ii+5*N] /= x_bac;
    }


//reseting counters

//##################################
//#Dynamics

    printf("\n#############################\n");
    printf("label = %d, N = %d\n\n", label, N);
    printf("V = %d, mesh_fac = %g, cell_num = %d, seed = %d\n", L2*L, mesh_fac,cell_num, seed_int);

    printf("num_threads = %d\n\n", num_threads);

    printf("tobs = %.1e,\n", (double)tobs);


    energy(0,samples);
    double e0 = samples[0].e;
    double lambda = 0.;


    printf("P0 = %g, h = %g, eps = %g, pf = %.3g\n", samples[0].P, h_bac, eps, samples[0].V_hs/(L2/pow(samples[0].amp_fac,2)));
    printf("T0 = %g\n", T0);
    printf("#############################\n\n");

    V_hard_sphere(0, samples);
    printf("\nK=%g, V=%g, e=%g, P=%g", samples[0].K,  samples[0].V, samples[0].e,  samples[0].P);



    testz(0,samples);
    printf("\ninitial condition preparation\n");

//##################################
//#initial condition generation
    h=h_bac;
    init_cond(0,samples);
    h=h_bac;

    printf("\n");

    printf("\ntemperature adjustments\n");

    if(T0 > 0)
        while( temperature(0,samples) < T0 ){
            time_wait2(0,samples);
            printf("\rT = %.2e ", temperature(0,samples));
            fflush(stdout);
        }

    for(ii=0; ii < tw; ++ii)
        time_wait(0,samples);


//##################################

    printf("\n");

    energy(0,samples);
    printf("\nK=%g, V=%g, e=%g, P=%g, T = %g", samples[0].K,  samples[0].V, samples[0].e,  samples[0].P, samples[0].K/(double)N);
    printf("\n");

    for (jj = 0; jj < N; jj++){
        samples[0].particles[jj].x0 = samples[0].particles[jj].x;
        samples[0].particles[jj].y0 = samples[0].particles[jj].y;
        samples[0].particles[jj].z0 = samples[0].particles[jj].z;
    }
    e0 = samples[0].e;

    FILE *fout0;
    sprintf(file_name,"./output/energy_%.2d.dat",label);
    fout0 = fopen(file_name, "w");

    jj = 0;
    printf("\n\ntobs\n");

    for (ii = 1; ii < tobs; ii++){

        time_step( 0., samples);
        lambda += log( mod_renorm(0., samples) );

        if(ii%frame==0){

//             printScrollingValue(ii, tobs);
            energy(0,samples);
            printf("\rtime = %.2f error_energy = %.2e energy = %.2e  ", (double)ii * h, fabs(samples[0].e - e0) / e0, samples[0].e);
            fflush(stdout);

            fprintf(fout0,"%g %g %g %g %g %g\n", (double)ii*h, fabs(samples[0].e - e0)/e0, samples[0].K, samples[0].V, samples[0].P, lambda/h/(double)ii);
            ovito_plot(0, samples, label, jj);
            jj += 1;

        }
    }

     printf("\n");
//     ovito_plot(0, samples, label, 0);


    energy(0,samples);
    printf("\nK=%g, V=%g, e=%g, P=%g", samples[0].K,  samples[0].V, samples[0].e,  samples[0].P);

    testz(0,samples);

//##################################
//plot


    return 0;

}

