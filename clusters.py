import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster import KMeans
import csv
# import alphashape
#from descartes import PolygonPatch
# synthetic classification dataset
from numpy import where
from sklearn.datasets import make_classification
from sklearn.cluster import DBSCAN
# birch clustering
from numpy import unique
from numpy import where
from sklearn.metrics import pairwise_distances
from sklearn.datasets import make_classification
from sklearn.cluster import Birch
from sklearn.cluster import MeanShift
import glob


import scipy
import scipy.cluster.vq
import scipy.spatial.distance
dst = scipy.spatial.distance.euclidean

def readcsv(filename):
    data = pd.read_csv(filename) #Please add four spaces here before this line
    return(np.array(data))

#yourArray = readcsv('positionsBreak.csv')#ShortPositions200D5nvalue0.csv') # ('positions.csv')

#yourArray = []
#yourArray = readcsv('CiL only/PositionsCiLOnly100D15nvalue9.csv')
#yourArray = readcsv('CoA CiL/PositionsCoACiL50D12nvalue5.csv')


# yourArray = readcsv('CoA CiL/PositionsCoACiL200D3nvalue11.csv')
# #
# #
# X = yourArray; #[:,[1,2]]

#files = glob.glob('trialNEWPositionsCoACiLeps200D6nvalue*.csv')

## Version 1


# plt.figure(figsize=(12, 3))
# for k in range(1,6):
#     kmeans = KMeans(n_clusters=k)
#     a = kmeans.fit_predict(X)
#     plt.subplot(1,5,k)
#     plt.scatter(X[:, 0], X[:, 1], c=a)
#     plt.xlabel('k='+str(k),fontsize = 16)
#
# #plt.tight_layout()

# plt.show()


def compute_inertia(a, y):
	W = [np.mean(pairwise_distances(y[a == c, :])) for c in np.unique(a)]
	return np.mean(W)

def bounding_box(X):
    xmin, xmax = min(X,key=lambda a:a[0])[0], max(X,key=lambda a:a[0])[0]
    ymin, ymax = min(X,key=lambda a:a[1])[1], max(X,key=lambda a:a[1])[1]
    return (xmin,xmax), (ymin,ymax)


def compute_gap(clustering, data, k_max=5, n_references=5):
	(xmin, xmax), (ymin, ymax) = bounding_box(X)
	# Create B reference datasets
	B = 5

	for i in range(B):
		Xb = []
		for n in range(len(X)):
			Xb.append([np.random.uniform(xmin, xmax),
					   np.random.uniform(ymin, ymax)])
		Xb = np.array(Xb)
	# if len(data.shape) == 1:
	# 	data = data.reshape(-1, 1)
	# reference = np.random.uniform( (xmin,xmax) ,*data.shape)
	reference = Xb
	#print("referce")
	#print(reference)

	reference_inertia = []
	for k in range(1, k_max + 1):
		local_inertia = []
		for _ in range(n_references):
			clustering.n_clusters = k
			assignments = clustering.fit_predict(reference)
			local_inertia.append(compute_inertia(assignments, reference))
		reference_inertia.append(np.mean(local_inertia))

	ondata_inertia = []
	for k in range(1, k_max + 1):
		clustering.n_clusters = k
		assignments = clustering.fit_predict(data)
		ondata_inertia.append(compute_inertia(assignments, data))
	#print("ref")
	#print(reference_inertia)
	#print("data")
	#print(ondata_inertia)
	gap = np.log(reference_inertia) - np.log(ondata_inertia)
	return gap, np.log(reference_inertia), np.log(ondata_inertia)

k_max = 5
maximum = []

Ds = [1, 3, 6, 9, 12, 15]
for d in Ds:
	files = glob.glob('CiL only/PositionsCiLOnly100D%dnvalue*.csv' %d)
#files = glob.glob('CoA CiL/PositionsCoACiL200D3nvalue*.csv')

	print(files)
	### Works well
	for file in files:
		X = readcsv(file)
		gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)

		maxind = np.argmax(gap)
		#print('Max index', maxind)
		# check if maximum is sufficiently significant
		if gap[maxind] < 0.1:
			maximum.append(0)
		else:
			maximum.append(maxind)

	#
	# 	maximum.append(maxind)
	#
	print('Max index', maximum)
	print(len(maximum))
	## gap statistic


# # plot weights and gap
# plt.plot(range(1, k_max + 1), reference_inertia,
# 		 '--o', label='reference')
# plt.plot(range(1, k_max + 1), ondata_inertia,
# 		 '-o', label='data')
# plt.xlabel('k',fontsize = 16)
# plt.ylabel('log($W_k$)',fontsize = 16)
# plt.xticks(fontsize=16)
# plt.xticks(np.arange(0, 7, step=1))
# plt.yticks(fontsize=16)
# plt.grid(True)
# plt.legend(fontsize=16)
# plt.show()
# #
# #
# # plot gap
# plt.plot(range(1, k_max+1), gap, '-o')
# plt.ylabel('gap',fontsize = 16)
# plt.xlabel('k', fontsize = 16)
# # plt.ylim((-6.6,-5.5))
# plt.xticks(fontsize=16)
# plt.xticks(np.arange(0, 7, step=1))
# plt.yticks(fontsize=16)
# plt.grid(True)
# plt.show()
#
#
#





















#print(X)
# plt.scatter(X[:,0], X[:,1])
# plt.show()

# wcss = []
# for i in range(1, 11):
#     kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
#     kmeans.fit(X)
#     wcss.append(kmeans.inertia_)
# plt.plot(range(1, 11), wcss)
# plt.title('Elbow Method')
# plt.xlabel('Number of clusters')
# plt.ylabel('WCSS')
# plt.show()
#
# kmeans = KMeans(n_clusters=2, init='k-means++', max_iter=300, n_init=10, random_state=0)
# pred_y = kmeans.fit_predict(X)
# plt.scatter(X[:,0], X[:,1])
# plt.scatter(kmeans.cluster_centers_[:, 0], kmeans.cluster_centers_[:, 1], s=300, c='red')
# plt.show()

# alpha_shape = alphashape.alphashape(X,0.01)
#
#
# fig, ax = plt.subplots()
# ax.scatter(*zip(*X))
# ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
# plt.show()


# BIRCH
# # define the model
# model = Birch(threshold=0.01, n_clusters=3)
# # fit the model
# model.fit(X)
# # assign a cluster to each example
# yhat = model.predict(X)
# # retrieve unique clusters
# clusters = unique(yhat)
# # create scatter plot for samples from each cluster
# for cluster in clusters:
# 	# get row indexes for samples with this cluster
# 	row_ix = where(yhat == cluster)
# 	# create scatter of these samples
# 	plt.scatter(X[row_ix, 0], X[row_ix, 1])
# # show the plot
# plt.show()

# # DBSCAN
# # define the model
# model = DBSCAN(eps=0.30, min_samples=9)
# # fit model and predict clusters
# yhat = model.fit_predict(X)
# # retrieve unique clusters
# clusters = unique(yhat)
# # create scatter plot for samples from each cluster
# for cluster in clusters:
# 	# get row indexes for samples with this cluster
# 	row_ix = where(yhat == cluster)
# 	# create scatter of these samples
# 	plt.scatter(X[row_ix, 0], X[row_ix, 1])
# # show the plot
# plt.show()


# # define the model
# model = MeanShift()
# # fit model and predict clusters
# yhat = model.fit_predict(X)
# # retrieve unique clusters
# clusters = unique(yhat)
# # create scatter plot for samples from each cluster
# print(clusters)
# for cluster in clusters:
# 	# get row indexes for samples with this cluster
# 	row_ix = where(yhat == cluster)
# 	# create scatter of these samples
# 	plt.scatter(X[row_ix, 0], X[row_ix, 1])
# # show the plot
# plt.show()





