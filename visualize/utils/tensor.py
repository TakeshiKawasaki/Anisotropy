import numpy as np

def tensor(nx, ny, nz, scale, ratio):
	# input
	# orientation of the major axis: nx, ny, nz
	# scaling of the eigenvectors: e1, e2, e1 > e2

	# output 
	# tensor components: XX XY XZ YX YY YZ ZX ZY ZZ
	# [XX XY XZ]
	# [YX YY YZ]
	# [ZX ZY ZZ]

	# determine eigen vectors
	XZ = nx
	YZ = ny
	ZZ = nz

	#print(np.power(XZ,2) + np.power(YZ, 2) + np.power(ZZ,2))

	# you need a scale factor to ensure the eigenvalue is 1 for this eigen vector. Otherwise the cross product won't be orthogonal   
	d = 1/(np.sqrt(np.power(XZ,2)+ np.power(ZZ,2)))
	XX = d*ZZ
	YX = 0
	ZX = -d*XZ

	#print(np.power(XX,2) + np.power(YX, 2) + np.power(ZX,2))

	XY = -d*XZ*YZ
	YY = d*(np.power(XZ,2) + np.power(ZZ,2))
	ZY = -d*YZ*ZZ

	#print(np.power(XY,2) + np.power(YY, 2) + np.power(ZY,2))

	#print(d)

	#print(XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ)
	


	# scale eigenvectors using eigenvalues
	XX = scale*XX
	XY = scale*XY
	XZ = ratio*scale*XZ
	YX = scale*YX
	YY = scale*YY
	YZ = ratio*scale*YZ
	ZX = scale*ZX
	ZY = scale*ZY
	ZZ = ratio*scale*ZZ

	#print(XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ) # this is the correct one

	return XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ

if __name__ == '__main__':
	nx = 0.111660
	ny = -0.007014
	nz = 0.993722
	scale = 1
	ratio = 2
	tensor(nx, ny, nz, scale, ratio)



