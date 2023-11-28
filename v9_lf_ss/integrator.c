void evolve_Q(int thread_id, Sample* samples){

    int jj,kk,ll,mm;

    for(jj=0; jj < N; ++jj){


        kk = (int)(samples[thread_id].particles[jj].x/mesh_fac)+L_mesh*(int)(samples[thread_id].particles[jj].y/mesh_fac)+L2_mesh*(int)(samples[thread_id].particles[jj].z/mesh_fac);
//         printf("%g %g\n", ALPHA*samples[thread_id].particles[jj].vx, samples[thread_id].particles[jj].vx);


        samples[thread_id].particles[jj].x += 0.5*h*samples[thread_id].particles[jj].vx;
        samples[thread_id].particles[jj].y += 0.5*h*samples[thread_id].particles[jj].vy;
//         samples[thread_id].particles[jj].z += 0.5*h*samples[thread_id].particles[jj].vz;


        while(samples[thread_id].particles[jj].x < 0.)
            samples[thread_id].particles[jj].x += L;

        while(samples[thread_id].particles[jj].x >= L)
            samples[thread_id].particles[jj].x -= L;

        while(samples[thread_id].particles[jj].y < 0.)
            samples[thread_id].particles[jj].y += L;

        while(samples[thread_id].particles[jj].y >= L)
            samples[thread_id].particles[jj].y -= L;

//         while(samples[thread_id].particles[jj].z < 0.)
//             samples[thread_id].particles[jj].z += L;
//
//         while(samples[thread_id].particles[jj].z >= L)
//             samples[thread_id].particles[jj].z -= L;

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

void evolve_P(int thread_id, Sample* samples){


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

                d = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z  );
                D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[II]._2r);

                if( d < D && II!=ii ){

                    fx += -eps*( 1. - d/D )*(-delta_x)/d/D;
                    fy += -eps*( 1. - d/D )*(-delta_y)/d/D;
//                     fz += -eps*( 1. - d/D )*(-delta_z)/d/D;
                }
            }
        }

        samples[thread_id].particles[ii].Fx = fx;
        samples[thread_id].particles[ii].Fy = fy;

        samples[thread_id].particles[ii].vx += h*fx;
        samples[thread_id].particles[ii].vy += h*fy;
//         samples[thread_id].particles[ii].vz += h*fz;
    }

    return;
}


void evolve_P_relax(int thread_id, Sample* samples){


    double fx, fy, fz;
    int ii, jj,kk,ll, II;

    int list[4];
    double x_bac, y_bac, z_bac;
    double d, delta_x, delta_y, delta_z, D, d_prime;

    double rlx_alpha = 0.05;
    double T0 = 1.e-4;

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

                d = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z  );
                D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[II]._2r);

                if( d < D && II!=ii ){

                    fx += -eps*( 1. - d/D )*(-delta_x)/d/D;
                    fy += -eps*( 1. - d/D )*(-delta_y)/d/D;
//                     fz += -eps*( 1. - d/D )*(-delta_z)/d/D;
                }
            }
        }

        samples[thread_id].particles[ii].Fx = fx;
        samples[thread_id].particles[ii].Fy = fy;

        samples[thread_id].particles[ii].vx += h*(fx - rlx_alpha*(temperature(thread_id, samples)/T0 - 1.)*samples[thread_id].particles[ii].vx );
        samples[thread_id].particles[ii].vy += h*(fy - rlx_alpha*(temperature(thread_id, samples)/T0 - 1.)*samples[thread_id].particles[ii].vy );
//         samples[thread_id].particles[ii].vz += h*fz;
    }

    return;
}

void evolve_DELTA_Q(int thread_id, Sample* samples){

    int jj;

    for(jj=0; jj < N; ++jj){

        samples[thread_id].delta[jj]     += 0.5*h*samples[thread_id].delta[jj+3*N];
        samples[thread_id].delta[jj+N]   += 0.5*h*samples[thread_id].delta[jj+4*N];
//         samples[thread_id].delta[jj+2*N] += 0.5*h*samples[thread_id].delta[jj+5*N];
    }

    return;
}

void evolve_DELTA_P(int thread_id, Sample* samples){


    double fx, fy, fz, fxx, fyy, fzz;
    int ii, jj,kk,ll, II;

    int list[4];
    double x_bac, y_bac, z_bac;
    double d, delta_x, delta_y, delta_z, D, d_prime;

    double Dx, Dy, Dz, Dd;

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

                d = sqrt( delta_x*delta_x + delta_y*delta_y + delta_z*delta_z  );
                D = 0.5*samples[thread_id].amp_fac*(samples[thread_id].particles[ii]._2r + samples[thread_id].particles[II]._2r);

                if( d < D && II!=ii ){

                    delta_x /= samples[thread_id].amp_fac;
                    delta_y /= samples[thread_id].amp_fac;
                    delta_z /= samples[thread_id].amp_fac;

                    D /= samples[thread_id].amp_fac;
                    d /= samples[thread_id].amp_fac;

                    Dx = samples[thread_id].delta[ii]      - samples[thread_id].delta[II];
                    Dy = samples[thread_id].delta[ii+N]    - samples[thread_id].delta[II+N];
                    Dz = samples[thread_id].delta[ii+2*N]  - samples[thread_id].delta[II+2*N];

                    Dd = delta_x*Dx + delta_y*Dy + delta_z*Dz;
                    Dd /= d;

                    fxx = Dx*d - delta_x*Dd;
                    fxx /= (d*d);
                    fxx -= Dx/D;
                    fxx *= (eps/D);

                    fyy = Dy*d - delta_y*Dd;
                    fyy /= (d*d);
                    fyy -= Dy/D;
                    fyy *= (eps/D);

                    fzz = Dz*d - delta_z*Dd;
                    fzz /= (d*d);
                    fzz -= Dz/D;
                    fzz *= (eps/D);

                    fx += fxx;
                    fy += fyy;
                    fz += fzz;

                }
            }
        }

        samples[thread_id].delta[ii+3*N] += h*fx;
        samples[thread_id].delta[ii+4*N] += h*fy;
//         samples[thread_id].delta[ii+5*N] += h*fz;
    }

    return;
}
