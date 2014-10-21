from mpl_toolkits.mplot3d import Axes3D 

centres = loadtxt('centres_test.txt', usecols=(1,2,3))
weights = loadtxt('weights_test.txt')
lengths = loadtxt('tract_lengths_test.txt')
vertices = loadtxt('vertices.txt')

reg =12 
li = weights[reg]/np.max(weights[reg]) * 256.
figure()
subplot(221)
scatter(centres[reg,0], centres[reg,1], s=2000., marker=(5,2))
scatter(centres[:,0], centres[:,1], c=cm.jet(li), s=50.)
scatter(vertices[::10,0], vertices[::10,1], marker='+', alpha=.3)
subplot(222)
scatter(centres[reg,0], centres[reg,2], s=2000., marker=(5,2))
scatter(centres[:,0], centres[:,2], c=cm.jet(li), s=50.)
scatter(vertices[::10,0], vertices[::10,2], marker='+', alpha=.3)
subplot(223)
scatter(centres[reg,1], centres[reg,2], s=2000., marker=(5,2))
scatter(centres[:,1], centres[:,2], c=cm.jet(li), s=50.)
scatter(vertices[::10,1], vertices[::10,2], marker='+', alpha=.3)
ax = subplot(224, projection='3d')
ax.scatter(centres[:,0], centres[:,1], centres[:,2], c=cm.hot(li), s=50.)
ax.scatter(centres[reg,0], centres[reg,1], centres[reg,2], c=cm.hot(li), s=2000., marker=(5,2))
ax.scatter(vertices[::10,0], vertices[::10,1], vertices[::10,2], marker='+', alpha=.3)
