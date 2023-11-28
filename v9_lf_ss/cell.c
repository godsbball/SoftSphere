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
