import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import numpy as np

PRD = os.environ['PRD']
CHECK = os.environ['CHECK']
if "DISPLAY" in os.environ:
    DISPLAY = os.environ['DISPLAY']
else:
    DISPLAY = ""
region_mapping_corr = float(os.environ['REGION_MAPPING_CORR'])
os.chdir(os.path.join(PRD, 'surface'))


def breadth_first_search(iposi, itrian, ilab):
    queue = [iposi]
    V = [iposi]
    while len(queue) > 0:
        iQ = queue.pop()
        iedges = list(itrian[np.argwhere(itrian==iQ)[:,0]].flatten())
        while len(iedges)>0:
            ineig = iedges.pop()
            if ineig not in V and texture[ineig]==ilab:
                V.append(ineig)
                queue.append(ineig)
    return V

def calculate_connected(texture, vert, trian):
    "find if the regions are connected components using Breadth-first seach"
    labels = np.unique(texture) 
    res = []
    for ilab in labels:
        ipos = np.nonzero(texture==ilab)
        ivert = vert[ipos]
        itrian=[]
        for itri in np.nonzero(texture==ilab)[0].tolist():
            itrian.extend(trian[np.nonzero(trian==itri)[0]])
        itrian = np.array(itrian).astype('int')
        V = breadth_first_search(ipos[0][0], itrian, ilab)
        res.append((ilab, ivert.shape[0]-len(V)))
    return res

def find_both_components(texture, vert, trian, ilab):
    " find the two subgraphs"
    ipos = np.nonzero(texture==ilab)
    ivert = vert[ipos]
    itrian=[]
    for itri in ipos[0].tolist():
        itrian.extend(trian[np.nonzero(trian==itri)[0]])
    itrian = np.array(itrian).astype('int')
    # first region
    V1 = breadth_first_search(ipos[0][0], itrian, ilab)
    # second region
    istart = 1
    while ipos[0][istart] in V1:
        istart += 1   
    V2 = breadth_first_search(ipos[0][istart], itrian, ilab)
    return (V1, V2)

def correct_sub_region(texture, trian, Vw):
    "correct the region mapping for the chosen component"
    new_texture = np.copy(texture)
    icount = 0
    while len(Vw)>0:
        iVw = Vw.pop()
        itrian = trian[np.nonzero(trian==iVw)[0]].flatten().astype('int').tolist()
        ir = list(filter(lambda x : new_texture[x] != new_texture[iVw], itrian))
        if len(ir)>0:
            new_texture[iVw] = new_texture[Counter(ir).most_common(1)[0][0]] 
        else:
            if icount<50: 
                Vw.insert(0, iVw)
                icount +=1
            else:
                # TODO: good error message
                print('error in correction')
                import pdb; pdb.set_trace()
    return new_texture

def check_region_mapping(texture, vert, trian, ilab):
    "drawing the region"
    ipos = np.nonzero(texture==ilab)
    itrian=[]
    for itri in ipos[0].tolist():
        itrian.extend(trian[np.nonzero(trian==itri)[0]])
    itrian = np.array(itrian).astype('int')
    bool_itrian = np.in1d(itrian, ipos[0]).reshape(itrian.shape[0], 3)
    itrian[np.nonzero(bool_itrian == False)] = 0
    citri = np.vstack([np.vstack([itrian[:,0], itrian[:,1]]).T, np.vstack([itrian[:,1],itrian[:,2]]).T, np.vstack([itrian[:,2],itrian[:,0]]).T])
    bcitri = (citri!=0).sum(1)
    valp = citri[bcitri==2]
    fig = plt.figure(figsize=(15, 15))
    fig.suptitle('region ' + str(int(ilab)))
    ax = fig.add_subplot(111, projection='3d') 
    plt.xlabel('x')
    plt.ylabel('y')
    for iv in np.arange(valp.shape[0]):
        ax.plot(vert[valp[iv], 0], vert[valp[iv], 1], vert[valp[iv], 2])
    # old function
    # xitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 0]
    # yitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 1]
    # zitrians = vert[np.hstack((itrian, itrian[:,0][:, newaxis])), 2]
    # fig = figure(figsize=(15, 15))
    # fig.suptitle('region ' + int(ilab))
    # ax = fig.add_subplot(111, projection='3d') 
    # for iv in range(xitrians.shape[0]): 
    #     ax.plot(xitrians[iv], yitrians[iv], zitrians[iv], alpha=0.4)
    # ax.scatter(vert[ipos, 0], vert[ipos, 1], vert[ipos, 2], c='b', s=45)
    show()

if __name__ == '__main__':

    rl = sys.argv[1]
    
    vert = np.loadtxt(rl+'_vertices_low.txt')
    trian = np.loadtxt(rl+'_triangles_low.txt')
    texture = np.loadtxt(rl + '_region_mapping_low.txt')

    res = np.array(calculate_connected(texture, vert, trian))
    wrong_labels = res[np.where(res[:,1]>0.),0][0]
    if len(wrong_labels)==0:
        print('evrything is fine, continuing')
    else:
        print("WARNING, some region have several components")
        for iwrong in wrong_labels:
            # TODO: handle more than two components
            (V1, V2) = find_both_components(texture, vert, trian, iwrong)
            if CHECK =="yes" and len(DISPLAY)>0:
                print("checking")
                check_region_mapping(texture, vert, trian, iwrong)
                i=0
                while True and i<10:
                    try:
                        choice_user = raw_input("""Do you want to get rid of region with:
                                            1) {0} nodes
                                            2) {1} nodes
                                            3) continue the pipeline anyway
                                            (answer: 1, 2, or 3)? \n""".format(len(V1), len(V2)))
                        print("you chose " + choice_user)

                        choice_user = int(choice_user)
                        if choice_user not in [1, 2, 3]:
                            raise ValueError
                        break 
                    except ValueError:
                        print('please choose 1, 2 or 3')
                        i += 1
                        continue
                else:
                    print('failure total, no check mode')
            else:
                print("no check, selecting automatically the smallest components")
                choice_user = np.argmin((len(V1), len(V2)))+1
            if choice_user==3:
                print("keep that correction")
                np.savetxt(rl + '_region_mapping_low.txt', new_texture)
            elif choice_user==1:
                texture = correct_sub_region(texture, trian, V1)
            elif  choice_user==2:
                texture =  correct_sub_region(texture, trian, V2)
            else: 
                print('failure of the choice')

        np.savetxt(rl + '_region_mapping_low.txt', texture)
