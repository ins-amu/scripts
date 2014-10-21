import os
import sys
PRD = os.environ['PRD']
CHECK = os.environ['CHECK']
DISPLAY = os.environ['DISPLAY']
region_mapping_corr = float(os.environ['region_mapping_corr'])
os.chdir(os.path.join(PRD, 'surface'))
from copy import deepcopy
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter

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
    new_texture = copy(texture)
    icount = 0
    while len(Vw)>0:
        iVw = Vw.pop()
        itrian = trian[np.nonzero(trian==iVw)[0]].flatten().tolist()
        ir = filter(lambda x : new_texture[x] != new_texture[iVw], itrian)
        if len(ir)>0:
            new_texture[iVw] = Counter(ir).most_common(1)[0][0] 
        else:
            if icount<10: 
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
    bool_itrian = np.in1d(itrian, ipos).reshape(itrian.shape[0], 3)
    itrian[np.nonzero(bool_itrian == False)] = 0
    citri = vstack([vstack([itrian[:,0], itrian[:,1]]).T, vstack([itrian[:,1],itrian[:,2]]).T, vstack([itrian[:,2],itrian[:,0]]).T])
    bcitri = (citri!=0).sum(1)
    valp = citri[bcitri==2]
    fig = figure(figsize=(15, 15))
    fig.suptitle('region ' + str(int(ilab)))
    ax = fig.add_subplot(111, projection='3d') 
    xlabel('x')
    ylabel('y')
    for iv in np.nditer(valp, flags=['external_loop'], order='C'):
        ax.plot(vert[iv, 0], vert[iv, 1], vert[iv, 2])

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
    vert = loadtxt(rl+'_vertices_low.txt')
    trian = loadtxt(rl+'_triangles_low.txt')
    texture = loadtxt(rl + '_region_mapping_low.txt')

    res = np.array(calculate_connected(texture, vert, trian))
    wrong_labels = res[np.where(res[:,1]>0.),0][0]
    if len(wrong_labels)==0:
        print 'evrything is fine, continuing'
    else:
        print "WARNING, some region have several components"
        for iwrong in wrong_labels:
            # TODO: handle more than two components
            (V1, V2) = find_both_components(texture, vert, trian, iwrong)
            if CHECK =="yes" and len(DISPLAY)>0:
                print "checking"
                check_region_mapping(texture, vert, trian, iwrong)
                choice_user = raw_input("""Do you want to get rid of region with:\n\
                                    1) {0} nodes\n\
                                    2) {1} nodes\n\
                                    3) continue the pipeline anyway\n\
                                    (answer: 1, 2, or 3)? \n"""
                                    .format(len(V1), len(V2)))
                print "you chose " + choice_user
            else:
                print "no check, selecting automatically the smallest components"
                choice_user=argmin((len(V1), len(V2)))+1
            if int(choice_user)==3:
                print "keep that correction"
                savetxt(rl + '_region_mapping_low.txt', new_texture)
            elif int(choice_user)==1:
                texture  =  correct_sub_region(texture, trian, V1)
            elif  int(choice_user)==2:
                texture =  correct_sub_region(texture, trian, V2)
            else: 
                print('please choose 1, 2, or 3')

        savetxt(rl + '_region_mapping_low.txt', texture)