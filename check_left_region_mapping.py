import os
PRD = os.environ['PRD']
region_mapping_corr = float(os.environ['region_mapping_corr'])
os.chdir(os.path.join(PRD, 'surface'))
from copy import deepcopy
from pylab import *
from mpl_toolkits.mplot3d import Axes3D

def calculate_connected(texture):
    "find if the regions are connected components using Breadth-first seach"
    vert = loadtxt('lh_vertices_low.txt')
    trian = loadtxt('lh_triangles_low.txt')
    labels = np.unique(texture) 
    res = []
    for ilab in labels:
        ivert = vert[np.nonzero(texture==ilab)]
        itrian=[]
        for itri in np.nonzero(texture==ilab)[0].tolist():
            itrian.extend(trian[np.nonzero(trian==itri)[0]])
        itrian = np.array(itrian).astype('int')
        queue = [itrian[0, 0]]
        V = [itrian[0, 0]]
        while len(queue) > 0:
            iQ = queue.pop()
            iedges = list(itrian[np.argwhere(itrian==iQ)[:,0]].flatten())
            while len(iedges)>0:
                ineig = iedges.pop()
                if ineig not in V and texture[ineig]==ilab:
                    V.append(ineig)
                    queue.append(ineig)
        res.append((ilab, ivert.shape[0]-len(V)))
    return res

def find_both_components(texture, labels):
    " find the two subgraphs"
    vert = loadtxt('lh_vertices_low.txt')
    trian = loadtxt('lh_triangles_low.txt')
    ivert = vert[np.nonzero(texture==ilab)]

    itrian=[]
    for itri in np.nonzero(texture==ilab)[0].tolist():
        itrian.extend(trian[np.nonzero(trian==itri)[0]])
    itrian = np.array(itrian).astype('int')

    # first region
    queue = [itrian[0, 0]]
    V1 = [itrian[0, 0]]
    while len(queue) > 0:
        iQ = queue.pop()
        iedges = list(itrian[np.argwhere(itrian==iQ)[:,0]].flatten())
        while len(iedges)>0:
            ineig = iedges.pop()
            if ineig not in V1 and texture[ineig]==ilab:
                V1.append(ineig)
                queue.append(ineig)

    # second region
    istart = 1
    while itrian[istart, 0] in V1
        istart += 1   
    queue = [itrian[istart, 0]]
    V2 = [itrian[istart, 0]]
    while len(queue) > 0:
        iQ = queue.pop()
        iedges = list(itrian[np.argwhere(itrian==iQ)[:,0]].flatten())
        while len(iedges)>0:
            ineig = iedges.pop()
            if ineig not in V2 and texture[ineig]==ilab:
                V2.append(ineig)
                queue.append(ineig)

    return (V1, V2)

def correct_sub_region(wrong_comp):

def check_region_mapping(texture, param_corr, labels):
    vert = loadtxt('lh_vertices_low.txt')
    trian = loadtxt('lh_triangles_low.txt')
    new_texture = deepcopy(texture)
    # labels = np.unique(texture)
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
                    if list_pos1[int_curr]==0:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 1]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 2]]))
                    if list_pos1[int_curr] ==1:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 0]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 2]]))
                    if list_pos1[int_curr] ==2:
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 0]]))
                        res_curr.append(np.round(texture[trian[list_pos0[int_curr], 1]]))
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
            if len(c_withdraw)>0:
                #print c_withdraw
                c_withdraw = np.array(c_withdraw)
                ax.scatter(c_withdraw[:,0],c_withdraw[:,1],c_withdraw[:,2], color='b', s=80)
            show()
    return new_texture

if __name__ == '__main__':

    texture = loadtxt('lh_region_mapping_low.txt')

    while True:
        res = np.array(calculate_connected(texture))
        wrong_labels = res[np.where(res[:,1]>0.),0][0]
        new_texture = check_region_mapping(texture, region_mapping_corr, wrong_labels)
        
        choice_user = raw_input("Do you want to get rid of region with:\n"\
                            "1) %{1} nodes\n"\
                            "2) %{2} nodes\n"\
                            "3) continue the pipeline\n".format({1:str(V1), 2:str(v2)}))\
                            "(answer: 1, 2, 3 or 4)? \n")
        print "you chose " + choice_user
        if int(choice_user)==3:
            print "keep that correction"
            savetxt('lh_region_mapping_low.txt', new_texture)
            quit()
        elif int(choice_user)==1:
            print('rerun with another value')
            param_corr = np.float(raw_input('enter new value for the correction parameter: \n'))
        elif  int(choice_user)==2:
            print('run another time with same correction value')
            texture = deepcopy(new_texture)
        elif int(choice_user)==4:
            print('run another time with new correction value')
            param_corr = np.float(raw_input('enter new value for the correction parameter: \n'))
            texture = deepcopy(new_texture)

