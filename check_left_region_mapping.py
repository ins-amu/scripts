import os
PRD = os.environ['PRD']
os.chdir(os.path.join(PRD, 'surface'))
from copy import deepcopy
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

vert = loadtxt('lh_vertices_low.txt')
trian = loadtxt('lh_triangles_low.txt') 
texture = loadtxt('lh_region_mapping_low.txt')
new_texture = deepcopy(texture)
labels = np.unique(texture)

for i in labels:
    list_vert = (nonzero(np.round(texture)==i)[0]).tolist()
    if len(list_vert)>0:
        print '#####################'
        print i
        print '#####################'

        #print x_pos

        fig = figure(figsize=(15,15))
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
            ax = fig.gca(projection='3d')
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
            if len([x for x in res_curr if x==i])<3*len(res_curr)/7.:
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