void r_dist(){

    int ii;

    double mean = 0.;

    for(ii=0; ii < N; ++ii){

        if(xor64() > 0.5)
            _2r_h[ii] = 1.;
        else
            _2r_h[ii] = 1.4;

        mean += _2r_h[ii];
    }

    mean /= (double)N;

    for(ii=0; ii < N; ++ii)
        _2r_h[ii] /= mean;

    return;
}

void V_hard_sphere(){

    int ii;

    V_hs = 0.;
    for(ii=0; ii < N; ++ii)
        V_hs += pow(0.5*_2r_h[ii],2);

    V_hs *= M_PI;

    return;
}

void device_structure_definition()
{

    switch (N) {

    case 1 << 25:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 512;
        break;

    case 1 << 24:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 256;
        break;

    case 1 << 23:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 128;
        break;


    case 1 << 22:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 64;
        break;

    case 1 << 21:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 32;
        break;

    case 1 << 20:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 16;
        break;

    case 1 << 19:
        block_size = 512;
        NumBlock = (1<<7);
        ct = 8;
        break;

    case 1 << 18:
	   block_size = 512;
	   NumBlock = (1<<7);
	   ct = 4;
	   break;

    case 1 << 17:
	   block_size = 512;
	   NumBlock = (1<<7);
	   ct = 2;
	   break;

	default :
	   block_size = 512;
	   NumBlock = N/block_size;
	   ct = 1;
    }

    BSF = NumBlock < block_size ? NumBlock : block_size;
    ct2 = NumBlock > block_size ? NumBlock/block_size : 1;

    return;
}
