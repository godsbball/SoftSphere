void ovito_plot(int time){

    int ii;

    char file_name[50];

    FILE *fout_conf;
    sprintf(file_name,"./output/conf/c_%03d.xyz", label);
    fout_conf = fopen(file_name, "w");
//     fout_conf = fopen(file_name, "a");

    fprintf( fout_conf,"ITEM: TIMESTEP\n");
    fprintf( fout_conf,"%d\n", time);

    fprintf( fout_conf,"ITEM: NUMBER OF ATOMS\n");
    fprintf( fout_conf,"%d\n", N);

    fprintf( fout_conf,"ITEM: BOX BOUNDS pp pp pp\n");
    fprintf( fout_conf,"%d %g\n", 0,(float)L/amp_fac);
    fprintf( fout_conf,"%d %g\n", 0,(float)L/amp_fac);
    fprintf( fout_conf,"%d %g\n", 0, (float)L/amp_fac);

    fprintf(fout_conf, "ITEM: ATOMS id type x y z Color Color Color radius Transparency\n");

    float rr,gg,bb;
    float bac_color;

    rr = 128/255;
    gg = 128/255;
    bb = 128/255;

    for(ii=0; ii < N; ++ii){

        bac_color = sqrt((x_h[ii] - 0.5*(float)L)*(x_h[ii] - 0.5*(float)L) + (y_h[ii] - 0.5*(float)L)*(y_h[ii] - 0.5*(float)L));

        if(bac_color > 0.97*0.5*(float)L)
            rr=gg=bb=0;

        else{

            if(_2r_h[ii] > 1.){
                rr = 65./255.;
                gg = 105./255.;
                bb = 225./255.;
            }

            else{
                rr = 225./255.;
                gg = 165./255.;
                bb = 0/255.;
            }

        }

        fprintf( fout_conf,"%d %d %g %g %g %g %g %g %g %g\n", ii, 1, x_h[ii]/amp_fac,y_h[ii]/amp_fac, 0., rr, gg, bb, 0.5*_2r_h[ii], 0.15);

    }
    fclose(fout_conf);

//     FILE *fout_conf2;
//     sprintf(file_name,"./output/conf/c_%03d.dat", label);
//     fout_conf2 = fopen(file_name, "w");
//
//     for(ii=0; ii < N; ++ii){
//         fprintf( fout_conf2,"%g %g\n",x,y);
//         fprintf( fout_conf2,"%g",_2r);
//
//     }
//
//     fclose(fout_conf2);
//
//
//     FILE *fout_conf3;
//     sprintf(file_name,"./output/conf/par_%03d.dat", label);
//     fout_conf3 = fopen(file_name, "w");
//
//     fprintf( fout_conf3,"%d %g", N, amp_fac);
//
//     fclose(fout_conf3);

    return;
}

