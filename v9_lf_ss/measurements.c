void energy(int thread_id, Sample* samples){

    int ii,jj,kk,ll,II;
    double d, delta_x, delta_y, delta_z, D;
    double K_bac, V_bac, vx_bac, vy_bac, vz_bac;

    int list[4];

    double x_bac, y_bac, z_bac;
    K_bac = V_bac = vx_bac = vy_bac = vz_bac = 0.;


    for(ii=0; ii<N; ++ii){

        vx_bac += samples[thread_id].particles[ii].vx;
        vy_bac += samples[thread_id].particles[ii].vy;
        vz_bac += samples[thread_id].particles[ii].vz;

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
                    V_bac+= 0.5*eps*(1. - d/D)*(1. - d/D); //warning, double counting
            }
        }
    }

    K_bac *= 0.5;
    V_bac *= 0.5;  //warning, double counting

    vx_bac *= vx_bac;
    vy_bac *= vy_bac;
    vz_bac *= vz_bac;

    samples[thread_id].K = K_bac;
    samples[thread_id].V = V_bac;
    samples[thread_id].e = K_bac+V_bac;
    samples[thread_id].P = sqrt(vx_bac + vy_bac + vz_bac);

    return;
}

double temperature(int thread_id, Sample* samples){

    int ii;
    double bac = 0.;

    for(ii=0; ii<N; ++ii)
        bac += samples[thread_id].particles[ii].vx*samples[thread_id].particles[ii].vx + samples[thread_id].particles[ii].vy*samples[thread_id].particles[ii].vy;

    return 0.5*bac/(double)N;
}

double mod_renorm(int thread_id, Sample* samples){

    int ii;
    double mod = 0.;
    for(ii=0; ii<6*N; ++ii)
        mod += samples[thread_id].delta[ii]*samples[thread_id].delta[ii];

    mod = sqrt(mod);

    for(ii=0; ii<6*N; ++ii)
        samples[thread_id].delta[ii] /= mod;

    return mod;
}
