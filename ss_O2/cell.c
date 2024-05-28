__device__ void neighbor_cell_list(int label_bac, int index, int list[]){

    int bac, sx, sy;

    switch (label_bac%9){

        case 0:
            list[0] = index;
            list[1] = 0;
            list[2] = 0;
            list[3] = 0;

            if( (int)(label_bac/9 )==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }
            break;

        case 1:
            bac = index - 1; //left
            sx = sy = 0;
            if(index%L_mesh_d==0){
                bac += L_mesh_d;
                sx = L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 2:
            bac = index + 1; //right
            sx = sy = 0;
            if( (index+1)%L_mesh_d==0 ){
                bac -= L_mesh_d;
                sx = -L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 3:
            bac = index - L_mesh_d; //bottom
            sx = sy = 0;
            if(index%L2_mesh_d<L_mesh_d){
                bac += L2_mesh_d;
                sy = L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 4:
            bac = index + L_mesh_d; //top
            sx = sy = 0;
            if( index%L2_mesh_d >= L_mesh_d*(L_mesh_d-1)){
                bac -= L2_mesh_d;
                sy = -L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 5:
            bac = index + L_mesh_d - 1; //top/left
            sx = sy = 0;

            if(index%L2_mesh_d>=L_mesh_d*(L_mesh_d-1)){
                bac -= L2_mesh_d;
                sy = -L_d;
            }

            if(index%L_mesh_d==0){
                bac += L_mesh_d;
                sx = L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }


            break;

        case 6:
            bac = index + L_mesh_d + 1; //top/right
            sx = sy = 0;

            if(index%L2_mesh_d>=L_mesh_d*(L_mesh_d-1)){
                bac -= L2_mesh_d;
                sy = -L_d;
            }

            if( (index+1)%L_mesh_d==0 ){
                bac -= L_mesh_d;
                sx = -L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 7:
            bac = index - L_mesh_d - 1; //bottom/left
            sx = sy = 0;

            if(index%L2_mesh_d<L_mesh_d){
                bac += L2_mesh_d;
                sy = L_d;
            }

            if(index%L_mesh_d==0){
                bac += L_mesh_d;
                sx = L_d;
            }


            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;

        case 8:
            bac = index - L_mesh_d + 1; //bottom/rigth
            sx = sy = 0;


            if(index%L2_mesh_d<L_mesh_d){
                bac += L2_mesh_d;
                sy = L_d;
            }

            if( (index+1)%L_mesh_d==0 ){
                bac -= L_mesh_d;
                sx = -L_d;
            }

            list[0] = bac;
            list[1] = sx;
            list[2] = sy;
            list[3] = 0;

            if( (int)(label_bac/9)==1 ){

                list[0] += L2_mesh_d;

                if( list[0] >= cell_num_d ){

                    list[0] -= cell_num_d;
                    list[3] = -L_d;
                }
            }

            if( (int)(label_bac/9)==2 ){

                list[0] -= L2_mesh_d;

                if( list[0] < 0 ){

                    list[0] += cell_num_d;
                    list[3] = L_d;
                }
            }

            break;
    }

    return;
}
