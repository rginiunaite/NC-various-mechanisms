import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster import KMeans
import csv
import alphashape
from descartes import PolygonPatch
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


def readcsv(filename):
    data = pd.read_csv(filename) #Please add four spaces here before this line
    return(np.array(data))

yourArray = readcsv('positionsHalfbreak.csv')

# with open('positions.csv') as csvfile:
#     readCSV = (csv.reader(csvfile, delimiter=','))
#     arrayfromCSV = np.array(readCSV)
#     # for row in arrayfromCSV:
#     #     print(row[1],row[2],)
print(yourArray)

X = yourArray[:,[1,2]]


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




## gap statistic


plt.figure(figsize=(12, 3))
for k in range(1,6):
    kmeans = KMeans(n_clusters=k)
    a = kmeans.fit_predict(X)
    plt.subplot(1,5,k)
    plt.scatter(X[:, 0], X[:, 1], c=a)
    plt.xlabel('k='+str(k))
plt.tight_layout()
plt.show()


def compute_inertia(a, y):
	W = [np.mean(pairwise_distances(y[a == c, :])) for c in np.unique(a)]
	return np.mean(W)


def compute_gap(clustering, data, k_max=5, n_references=5):
	if len(data.shape) == 1:
		data = data.reshape(-1, 1)
	reference = np.random.rand(*data.shape)
	print("referce")
	print(reference)
	reference = reference 
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
	print("ref")
	print(reference_inertia)
	print("data")
	print(ondata_inertia)
	gap = np.log(reference_inertia) - np.log(ondata_inertia)
	return gap, np.log(reference_inertia), np.log(ondata_inertia)


k_max = 5
gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)

plt.plot(range(1, k_max + 1), reference_inertia,
		 '-o', label='reference')
plt.plot(range(1, k_max + 1), ondata_inertia,
		 '-o', label='data')
plt.xlabel('k')
plt.ylabel('log(inertia)')
plt.show()


# plot gap
plt.plot(range(1, k_max+1), gap, '-o')
plt.ylabel('gap')
plt.xlabel('k')
plt.ylim((-6.6,-5.5))

plt.show()