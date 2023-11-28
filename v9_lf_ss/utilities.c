void testz(int II, Sample* samples){

    double d, delta_x, delta_y, delta_z, D;
    int ii,jj;

    double K_bac, V_bac, vx_bac, vy_bac, vz_bac;
    K_bac = V_bac = vx_bac = vy_bac = vz_bac = 0.;


    for(ii=0; ii < N; ++ii){

        vx_bac = samples[II].particles[ii].vx;
        vy_bac = samples[II].particles[ii].vy;
        vz_bac = samples[II].particles[ii].vz;

        vx_bac *= vx_bac;
        vy_bac *= vy_bac;
        vz_bac *= vz_bac;

        K_bac += vx_bac + vy_bac + vz_bac;

        for(jj=ii+1; jj < N; ++jj){

            delta_x = samples[II].particles[ii].x - samples[II].particles[jj].x;
            delta_y = samples[II].particles[ii].y - samples[II].particles[jj].y;
            delta_z = samples[II].particles[ii].z - samples[II].particles[jj].z;

            delta_x = fabs(delta_x);
            delta_y = fabs(delta_y);
            delta_z = fabs(delta_z);

            delta_x = fmin( (double)L-delta_x, delta_x );
            delta_y = fmin( (double)L-delta_y, delta_y );
            delta_z = fmin( (double)L-delta_z, delta_z );

            d = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
            d = sqrt(d);

            D = 0.5*samples[II].amp_fac*(samples[II].particles[ii]._2r + samples[II].particles[jj]._2r);


            if(d < D){
                V_bac += 0.5*eps*(1. - d/D)*(1. - d/D); //warning, single counting

            }
        }
    }

    K_bac *= 0.5;
    printf("\nperformDoubleCheck: K = %g, V = %g, e = %g,     sample = %d\n", K_bac, V_bac, K_bac+V_bac, II);


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
    fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[thread_id].amp_fac);
    fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[thread_id].amp_fac);
    fprintf( fout_conf,"%d %g\n", 0,(double)L/samples[thread_id].amp_fac);

    fprintf(fout_conf, "ITEM: ATOMS id type x y z color color color radius Transparency\n");

    double rr,gg,bb, bac_double;

    rr = 128/255;
    gg = 128/255;
    bb = 128/255;

    for(ii=0; ii < N; ++ii){
        fprintf( fout_conf,"%d %d %g %g %g %g %g %g %g %g\n", ii, 1, samples[thread_id].particles[ii].x/samples[thread_id].amp_fac, samples[thread_id].particles[ii].y/samples[thread_id].amp_fac, samples[thread_id].particles[ii].z/samples[thread_id].amp_fac, rr, gg, bb, 0.5*samples[thread_id].particles[ii]._2r, 0.15);
    }
    fclose(fout_conf);

    FILE *fout_conf2;
    sprintf(file_name,"./output/conf/c_%02d_%02d_%04d.dat", label, thread_id, tt);
    fout_conf2 = fopen(file_name, "w");

    for(ii=0; ii < N; ++ii){
        fprintf( fout_conf2,"%g %g %g ", samples[thread_id].particles[ii].x, samples[thread_id].particles[ii].y, samples[thread_id].particles[ii].z);
        fprintf( fout_conf2,"%g %g %g %g\n", samples[thread_id].particles[ii].vx, samples[thread_id].particles[ii].vy, samples[thread_id].particles[ii].vz, samples[thread_id].particles[ii]._2r);

    }

    fclose(fout_conf2);


    FILE *fout_conf3;
    sprintf(file_name,"./output/conf/par_%02d_%02d_%04d.dat", label, thread_id, tt);
    fout_conf3 = fopen(file_name, "w");

    fprintf( fout_conf3,"%d %g %g", N, samples[thread_id].amp_fac, mesh_fac);

    fclose(fout_conf3);

    FILE *fout_conf4;
    sprintf(file_name,"./output/conf/delta_%02d_%02d_%04d.dat", label, thread_id, tt);
    fout_conf4 = fopen(file_name, "w");

    for(ii=0; ii < 2*N; ++ii)
        fprintf(fout_conf4,"%g\n", samples[thread_id].delta[ii]);

    for(ii=3*N; ii < 5*N; ++ii)
        fprintf(fout_conf4,"%g\n", samples[thread_id].delta[ii]);

    fclose(fout_conf4);

    return;
}


void V_hard_sphere(int thread_id, Sample* samples){

    int ii;

    samples[thread_id].V_hs = 0.;
    for(ii=0; ii < N; ++ii)
        samples[thread_id].V_hs += pow(samples[thread_id].particles[ii]._2r,2);
//         samples[thread_id].V_hs += pow(samples[thread_id].particles[ii]._2r,3);

//     samples[thread_id].V_hs *= (4./3.)*M_PI/8.;
    samples[thread_id].V_hs *= M_PI/4.;


//     printf("sample: %d,\t\tPF = %lf,\t\tamp_fac = %lf\n", thread_id, samples[thread_id].V_hs/(L2*L/pow(samples[thread_id].amp_fac,3)), samples[thread_id].amp_fac);
//     printf("sample: %d,\t\tPF = %lf,\t\tamp_fac = %lf\n", thread_id, samples[thread_id].V_hs/(L2/pow(samples[thread_id].amp_fac,2)), samples[thread_id].amp_fac);

    return;

}

void printScrollingValue(double value, double targetValue) {
    printf("\rValue: %.4g/%g ", value, targetValue);
    fflush(stdout);

    return;
}


void init_cond(int thread_id, Sample* samples){

    double h_tmax = h;
    double P_par, v_mod, F_mod, f_inc, f_dec, f_alpha, alpha_start, alpha;

    f_inc = 1.1;
    f_dec = 0.5;
    f_alpha = 0.99;

    alpha = alpha_start = 0.1;

    f_inc = 1.1;
    f_dec = 0.5;
    f_alpha = 0.99;

    alpha = alpha_start = 0.1;

    int ii, jj, kk, int_bac;

    ii = 0;
    while(samples[0].e > 1.e-3 && ii < 7e2){

        int_bac = 0;
        for (jj = 0; jj < 100; jj++){
            time_wait( 0., samples);

            if(ii == 6){
                P_par = 0;
                for (kk = 0; kk < N; kk++)
                    P_par += samples[0].particles[kk].vx*samples[0].particles[kk].Fx + samples[0].particles[kk].vy*samples[0].particles[kk].Fy;

                if(P_par > 0)
                    int_bac = 1;
            }
        }

        if(ii%1==0){
            energy(0,samples);
//             printf("\nK = %.3g, V = %.3g, e/N = %.3g, h = %.3g, ii = %d ", samples[0].K, samples[0].V, samples[0].e/N,  h, ii);
            printf("\rK = %.2e, V = %.2e, e = %.3e, h = %.2e, ii = %d   ", samples[0].K, samples[0].V, samples[0].e,  h, ii);
            fflush(stdout);
        }

        P_par = v_mod = F_mod = 0;
        for (jj = 0; jj < N; jj++){

            P_par += samples[0].particles[jj].vx*samples[0].particles[jj].Fx + samples[0].particles[jj].vy*samples[0].particles[jj].Fy;
            v_mod += samples[0].particles[jj].vx*samples[0].particles[jj].vx + samples[0].particles[jj].vy*samples[0].particles[jj].vy;

            F_mod += samples[0].particles[jj].Fx*samples[0].particles[jj].Fx + samples[0].particles[jj].Fy*samples[0].particles[jj].Fy;
        }

        v_mod = sqrt(v_mod);
        F_mod = sqrt(F_mod);


        for (jj = 0; jj < N; jj++){

            samples[0].particles[jj].vx = (1-alpha)*samples[0].particles[jj].vx + alpha*v_mod*samples[0].particles[jj].Fx/F_mod;
            samples[0].particles[jj].vy = (1-alpha)*samples[0].particles[jj].vy + alpha*v_mod*samples[0].particles[jj].Fy/F_mod;
        }


        if(P_par>0 || int_bac == 1){

            h *= f_inc;
            if(h>h_tmax)
                h = h_tmax;

            alpha *= f_alpha;

        }

        else{

            h *= f_dec;

            for (jj = 0; jj < N; jj++){
                samples[0].particles[jj].vx = 0.;
                samples[0].particles[jj].vy = 0.;
            }

            alpha = alpha_start;
        }

        ii += 1;
    }

    return;

}
