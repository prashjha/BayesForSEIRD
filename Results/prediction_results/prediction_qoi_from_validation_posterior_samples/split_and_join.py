import sys
import numpy as np 

def split():

	qoi = np.load('qoi.npy')

	n = 5
	dn = int(qoi.shape[0]/n)

	qoi_n = [qoi[i*dn:(i+1)*dn] for i in range(n)]

	# write 
	for i in range(n):
		np.save('qoi_'+str(i)+'.npy', qoi_n[i])

def join():

	n = 5
	qoi = np.load('qoi_0.npy')
	for i in range(n-1):
		qoi_n = np.load('qoi_'+str(i+1)+'.npy')
		qoi = np.append(qoi, qoi_n, axis=0)

	# write
	np.save('qoi.npy')

run_type = int(sys.argv[1])
if run_type == 0:
	split()
else:
	join()
