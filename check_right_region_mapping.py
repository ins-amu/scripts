import os
PRD = os.environ['PRD']
region_mapping_corr = os.environ['region_mapping_corr']
os.chdir(os.path.join(PRD, 'surface'))
from copy import deepcopy
from pylab import *
from mpl_toolkits.mplot3d import Axes3D


def check_region_mapping(texture, region_mapping_corr):
    vert = loadtxt('rh_vertices_low.txt')
    trian = loadtxt('rh_triangles_low.txt')
    labels = np.unique(texture)
    new_texture = deepcopy(texture)


    for indx_lab, i in enumerate(labels):
        fig = figure(figsize=(15,15))
        print indx_lab+1
        #ax = fig.add_subplot(ceil(len(labels)/5.), 5, indx_lab+1, projection='3d')
        ax = fig.gca(projection='3d')
        list_vert = (nonzero(np.round(texture)==i)[0]).tolist()
        if len(list_vert)>0:
            print '#####################'
            print i
            print '#####################'

            #print x_pos

            for curr_ind in reversed(list_vert):
                #print curr_ind
                x_pos, y_pos, z_pos = [], [], []
                (list_pos, _) = np.nonzero(trian == curr_ind)
                #print i
                #print str(vert[curr_ind]) +' : ' + str(curr_ind)
                for int_curr in range(len(list_pos)):

                    t1 = trian[list_pos[int_curr], 0]
                    t2 = trian[list_pos[int_curr], 1]
                    t3 = trian[list_pos[int_curr], 2]
                    #print(str(t1)+',_,'+str(t2)+',_,'+str(t3))
                    if round(texture[t1]) == i and round(texture[t2])==i:
                        x_pos.extend([vert[t1,0], vert[t2,0]])
                    if round(texture[t1]) == i and round(texture[t3])==i:
                        x_pos.extend([vert[t1,0], vert[t3,0]])
                    if round(texture[t2]) == i and round(texture[t2])==i:
                        x_pos.extend([vert[t2,0], vert[t2,0]])
                    if round(texture[t1]) == i and round(texture[t2])==i:
                        y_pos.extend([vert[t1,1], vert[t2,1]])
                    if round(texture[t1]) == i and round(texture[t3])==i:
                        y_pos.extend([vert[t1,1], vert[t3,1]])
                    if round(texture[t2]) == i and round(texture[t2])==i:
                        y_pos.extend([vert[t2,1], vert[t2,1]])
                    if round(texture[t1]) == i and round(texture[t2])==i:
                        z_pos.extend([vert[t1,2], vert[t2,2]])
                    if round(texture[t1]) == i and round(texture[t3])==i:
                        z_pos.extend([vert[t1,2], vert[t3,2]])
                    if round(texture[t2]) == i and round(texture[t2])==i:
                        z_pos.extend([vert[t2,2], vert[t2,2]])
                ax.plot(x_pos, y_pos, z_pos)
                xlabel('x')
                ylabel('y')
            c_withdraw = []
            change_val = []
            for indx_curr, vert_curr in enumerate(list_vert):
                list_pos0, list_pos1 = nonzero(trian==vert_curr)
                res_curr = []
                for int_curr in range(len(list_pos0)):

                    #print int_curr
                    if list_pos1[int_curr] ==0:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],1]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],2]]))
                    if list_pos1[int_curr] ==1:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],0]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],2]]))
                    if list_pos1[int_curr] ==2:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],0]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr],1]]))
                if len([x for x in res_curr if x==i])<region_mapping_corr*len(res_curr):
                    print res_curr
                    print  vert[vert_curr]
                    c_withdraw.append(vert[vert_curr])
                    counter = 0
                    good_val = 0.0
                    for curr_val in set(res_curr):
                        if len([x for x in res_curr if x ==curr_val])>counter:
                            counter = len([x for x in res_curr if x ==curr_val])
                            good_val = curr_val
                    change_val.append(good_val)
                    print good_val
                    #print list_vert[indx_curr]
                    new_texture[np.floor(list_vert[indx_curr])] = good_val
                    #import pdb; pdb.set_trace()

        if len(c_withdraw)>0:
            #print c_withdraw
            c_withdraw = np.array(c_withdraw)
            ax.scatter(c_withdraw[:,0],c_withdraw[:,1],c_withdraw[:,2], color='b', s=80)

        show()

    return new_texture

if __name__ == '__main__':
    texture = loadtxt('rh_region_mapping_low_not_corrected.txt')
    while True:
        new_texture = check_region_mapping(texture, region_mapping_corr)
        choice_user = raw_input("Do you want:\n"\
                            "1) to continue the pipeline with this correction\n"\
                            "2) rerun this correction with another value\n"\
                            "3) run another time with the same parameter value\n"\
                            "4) run another time with another parameter value\n"\
                            "(answer: 1, 2, 3 or 4)? \n")
        print "you chose " + choice_user
        if int(choice_user)==1:
            print "keep that correction"
            savetxt('rh_region_mapping_low.txt', new_texture)
            quit()
        elif int(choice_user)==2:
            print('rerun with another value')
            param_corr = np.float(raw_input('enter new value for the correction parameter: \n'))
        elif  int(choice_user)==3:
            print('run another time with same correction value')
            texture = deepcopy(new_texture)
        elif int(choice_user)==4:
            print('run another time with new correction value')
            param_corr = np.float(raw_input('enter new value for the correction parameter: \n'))
            texture = deepcopy(new_texture)


