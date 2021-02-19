## I will use gap statistic to find the number of clusters, there are some other techniques at the end of thsi file but I did not use them
## I used this code https://glowingpython.blogspot.com/2019/01/a-visual-introduction-to-gap-statistics.html

## better, not yet implement but I should! https://datasciencelab.wordpress.com/tag/gap-statistic/

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.datasets.samples_generator import make_blobs
from sklearn.cluster import KMeans
import glob
import csv
# import alphashape
# from descartes import PolygonPatch
# synthetic classification dataset
# from numpy import where
# from sklearn.datasets import make_classification
# from sklearn.cluster import DBSCAN
# # # birch clustering
# from numpy import unique
# from numpy import where
from sklearn.metrics import pairwise_distances
# from sklearn.datasets import make_classification
# from sklearn.cluster import Birch
# from sklearn.cluster import MeanShift
# import glob


import scipy
import scipy.cluster.vq
import scipy.spatial.distance

dst = scipy.spatial.distance.euclidean

varywhat = 0;  # 0 - D and eps, 1 - beta and eps, 2 - beta and D


def readcsv(filename):
    data = pd.read_csv(filename)  # Please add four spaces here before this line
    return (np.array(data))


#yourArray = readcsv('CoA CiL biased vary D/PositionsLEADONLYDCoACiLeps200D6nvalue2.csv')
#yourArray = readcsv('CoA CiL biased vary D/PositionsAllbiasDCoACiLeps200D6nvalue2.csv')
# yourArray = readcsv('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps75D15nvalue1.csv')
#yourArray = readcsv('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickeps0D9nvalue1.csv')
# yourArray = readcsv('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps56D15nvalue50.csv')
# yourArray = readcsv('../Attr Rep LEAD ONLY Chick Data Files/AttrRepLEADONLYVARYDPositionsChickeps19D6nvalue0.csv')
# yourArray = readcsv('../Attr Rep ALLBIASED Chick Data Files/AttrRepALLBIASEDVARYDPositionsChickeps75D3nvalue0.csv')

# yourArray = readcsv('../CHICK DATA FINAL POSITIONS/RepOnlyALLBIASEDBiasedPositionsChickeps0D10beta4nvalue1.csv')
# X = yourArray;


#yourArray = readcsv('CoA CiL/PositionsCoACiL50D12nvalue5.csv')
yourArray = readcsv('positionsBreak.csv') # in thesis
#yourArray = readcsv('positionsVerycont.csv') # in thesis
#yourArray = readcsv('positions.csv')


X = yourArray[:,[1,2]] # for positions break and very cont


# plt.scatter(X[:, 0], X[:, 1])
# plt.show()



# # Version 1

#If I want to see the cells
plt.figure(figsize=(15, 10))  # before 12,3
for k in range(1, 6):
    kmeans = KMeans(n_clusters=k)
    a = kmeans.fit_predict(X)
    plt.subplot(5, 1, k)
    plt.scatter(X[:, 0], X[:, 1], s=10, c=a)
    plt.axis('scaled')
     # set the xlim to left, right
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    #plt.xlabel('k=' + str(k), fontsize=16)
    # if k == 1:
    #     plt.ylabel('Domain width, $\mu$m', fontsize=16)
    # plt.xlabel('Domain length, $\mu$m\nk=' + str(k), fontsize=16)

# plt.tight_layout()

plt.show()

def compute_inertia(a, y):
	W = [np.mean(pairwise_distances(y[a == c, :])) for c in np.unique(a)]
	return np.mean(W)


def Wk(mu, clusters):
    K = len(mu)
    return sum([np.linalg.norm(mu[i]-c)**2/(2*len(c)) \
               for i in range(K) for c in clusters[i]])



def bounding_box(X):
    xmin, xmax = min(X, key=lambda a: a[0])[0], max(X, key=lambda a: a[0])[0]
    ymin, ymax = min(X, key=lambda a: a[1])[1], max(X, key=lambda a: a[1])[1]
    return (xmin, xmax), (ymin, ymax)



## new version
#
# def gap_statistic(X):
#     (xmin,xmax), (ymin,ymax) = bounding_box(X)
#         # Dispersion for real distribution
#     ks = range(1,10)
#     Wks = zeros(len(ks))
#     Wkbs = zeros(len(ks))
#     sk = zeros(len(ks))
#     for indk, k in enumerate(ks):
#         mu, clusters = find_centers(X,k)
#         Wks[indk] = np.log(Wk(mu, clusters))
#             # Create B reference datasets
#         B = 10
#         BWkbs = zeros(B)
#         for i in range(B):
#             Xb = []
#             for n in range(len(X)):
#                 Xb.append([random.uniform(xmin,xmax),random.uniform(ymin,ymax)])
#             Xb = np.array(Xb)
#             mu, clusters = find_centers(Xb,k)
#             BWkbs[i] = np.log(Wk(mu, clusters))
#         Wkbs[indk] = sum(BWkbs)/B
#         sk[indk] = np.sqrt(sum((BWkbs-Wkbs[indk])**2)/B)
#     sk = sk*np.sqrt(1+1/B)
#     return(ks, Wks, Wkbs, sk)
#
#
# ## EXAMPLE
#
# ks, logWks, logWkbs, sk = gap_statistic(X)


 ## OLD VERSION
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
    # print("referce")
    # print(reference)

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

    # print("ref")
    # print(reference_inertia)
    # print("data")
    # print(ondata_inertia)
    gap = np.log(reference_inertia) - np.log(ondata_inertia)

    return gap, np.log(reference_inertia), np.log(ondata_inertia)


# For gap

k_max = 5
#
        ## gap statistic
maximum = []
gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)  # compute gap

maxind = np.argmax(gap)
                # print('Max index', maxind)
                # check if maximum is sufficiently significant
if gap[maxind] < 0.1:
    maximum.append(0)
else:
    maximum.append(maxind)

#plot weights and gap
plt.plot(range(1, k_max + 1), reference_inertia,
		 '--o', label='reference')
plt.plot(range(1, k_max + 1), ondata_inertia,
		 '-o', label='data')
plt.xlabel('k',fontsize = 16)
plt.ylabel('log($W_k$)',fontsize = 16)
plt.xticks(fontsize=16)
plt.xticks(np.arange(0, 7, step=1))
plt.yticks(fontsize=16)
plt.grid(True)
plt.legend(fontsize=16)
plt.show()
#
#
# plot gap
plt.plot(range(1, k_max+1), gap, '-o')
plt.ylabel('gap',fontsize = 16)
plt.xlabel('k', fontsize = 16)
# plt.ylim((-6.6,-5.5))
plt.xticks(fontsize=16)
plt.xticks(np.arange(0, 7, step=1))
plt.yticks(fontsize=16)
plt.grid(True)
plt.show()



# Ds = [1, 3, 6, 9, 12, 15]# values of D
Ds = [1, 4, 7, 10, 13]
# Ds = [1, 2, 3, 4, 5] ## BETA!!!
# run through all the files
epsilons = [0, 19, 38, 56, 75]
betas = [1, 4, 7, 10, 13]

# vary D and eps
#
if varywhat == 0:
    for eps in epsilons:
        for d in Ds:
            # files = glob.glob('CiL only/PositionsCiLOnly1D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL/NEWPositionsCoACiLeps250D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL BIASED all/PositionsBetaCoACiLeps100beta%dnvalue*.csv' % d)
            # files = glob.glob('Coa CiL Beta Lead Only/PositionsLEADONLYBetaCoACiLeps200beta%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsAllbiasDCoACiLeps250D%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsLEADONLYDCoACiLeps250D%dnvalue*.csv' % d)

            # files = glob.glob('Rep Only - New/PositionsRepOnlyeps19D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepeps94D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepLEADONLYVARYbetaeps%dbeta%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYbetaeps%dbeta%dnvalue*.csv' % (eps, d))

            ## Files with 100 simulations
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # # for plotting gap
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickepseps%dD%dnvalue10.csv' % (eps, d))
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue50.csv' % (eps, d))
            # files = glob.glob('../Attr Rep LEAD ONLY Chick Data Files/AttrRepLEADONLYVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Attr Rep ALLBIASED Chick Data Files/AttrRepALLBIASEDVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # Rep only all biased
            files = glob.glob('../CHICK DATA FINAL POSITIONS2/RepOnlyALLBIASEDBiasedPositionsChickeps%dD%dbeta4nvalue*.csv'%(eps,d))
            # Attr Rep all biased
            # files = glob.glob('../CHICK DATA FINAL POSITIONS2/AttrRepALLBIASEDBiasedPositionsChickeps%dD%dbeta4nvalue*.csv'%(eps,d))
            # Attr Rep biased leaders
            # files = glob.glob(
            #     '../CHICK DATA FINAL POSITIONS2/AttrRepLeadersBiasedPositionsChickeps%dD%dbeta4nvalue*.csv' % (eps, d))
            # # Attr Rep
            # files = glob.glob('../CHICK DATA FINAL POSITIONS2/AttrRepPositionsChickeps%dD%dbeta0nvalue*.csv'%(eps,d))
            # Rep Only
            # files = glob.glob('../CHICK DATA FINAL POSITIONS2/RepOnlyPositionsChickeps%dD%dbeta0nvalue*.csv'%(eps,d))

            # print(files)
            maximum = []
            ### Works well
            for file in files:
                X = readcsv(file)
                gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)  # compute gap

                maxind = np.argmax(gap)
                # print('Max index', maxind)
                # check if maximum is sufficiently significant
                if gap[maxind] < 0.1:
                    maximum.append(0)
                else:
                    maximum.append(maxind)

            #
            # 	maximum.append(maxind)
            #
            print('Epsilon ', eps, 'D', d)
            # print('Max index', maximum)
            # print('Just to check that I have 100', len(maximum))
            print('Number of non-zero entries', np.count_nonzero(maximum))




## end vary D and eps


# vary beta and eps

if varywhat == 1:

    for eps in epsilons:
        for b in betas:
            # files = glob.glob('CiL only/PositionsCiLOnly1D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL/NEWPositionsCoACiLeps250D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL BIASED all/PositionsBetaCoACiLeps100beta%dnvalue*.csv' % d)
            # files = glob.glob('Coa CiL Beta Lead Only/PositionsLEADONLYBetaCoACiLeps200beta%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsAllbiasDCoACiLeps250D%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsLEADONLYDCoACiLeps250D%dnvalue*.csv' % d)

            # files = glob.glob('Rep Only - New/PositionsRepOnlyeps19D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepeps94D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepLEADONLYVARYbetaeps%dbeta%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYbetaeps%dbeta%dnvalue*.csv' % (eps, d))

            ## Files with 100 simulations
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # # for plotting gap
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickepseps%dD%dnvalue10.csv' % (eps, d))
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue50.csv' % (eps, d))
            # files = glob.glob('../Attr Rep LEAD ONLY Chick Data Files/AttrRepLEADONLYVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Attr Rep ALLBIASED Chick Data Files/AttrRepALLBIASEDVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # Rep Only All biased
            # files = glob.glob('../CHICK DATA FINAL POSITIONS2/RepOnlyALLBIASEDBiasedPositionsChickeps%dD7beta%dnvalue*.csv'%(eps,b))
            # Attr Rep all biased
            # files = glob.glob(
            #	'../CHICK DATA FINAL POSITIONS/AttrRepALLBIASEDBiasedPositionsChickeps%dD7beta%dnvalue*.csv' % (eps, b))
            # Attr Rep biased leaders
            files = glob.glob(
                '../CHICK DATA FINAL POSITIONS/AttrRepLeadersBiasedPositionsChickeps%dD7beta%dnvalue*.csv' % (eps, b))

            # print(files)
            maximum = []
            ### Works well
            for file in files:
                X = readcsv(file)
                gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)  # compute gap

                maxind = np.argmax(gap)
                # print('Max index', maxind)
                # check if maximum is sufficiently significant
                if gap[maxind] < 0.1:
                    maximum.append(0)
                else:
                    maximum.append(maxind)

            #
            # 	maximum.append(maxind)
            #
            print('Epsilon ', eps, 'beta', b)
            # print('Max index', maximum)
            # print('Just to check that I have 100', len(maximum))
            print('Number of non-zero entries', np.count_nonzero(maximum))
        ## gap statistic

# plot weights and gap
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

## end vary beta and eps


## vary D and beta

if varywhat == 2:
    for D in Ds:
        for b in betas:
            # files = glob.glob('CiL only/PositionsCiLOnly1D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL/NEWPositionsCoACiLeps250D%dnvalue*.csv' %d)
            # files = glob.glob('CoA CiL BIASED all/PositionsBetaCoACiLeps100beta%dnvalue*.csv' % d)
            # files = glob.glob('Coa CiL Beta Lead Only/PositionsLEADONLYBetaCoACiLeps200beta%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsAllbiasDCoACiLeps250D%dnvalue*.csv' % d)
            # files = glob.glob('CoA CiL biased vary D/PositionsLEADONLYDCoACiLeps250D%dnvalue*.csv' % d)

            # files = glob.glob('Rep Only - New/PositionsRepOnlyeps19D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepeps94D%dnvalue*.csv' % d)
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepLEADONLYVARYbetaeps%dbeta%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/NEWPositionsAttrRepLEADONLYVARYDeps%dD%dnvalue*.csv' % (eps,d))
            # files = glob.glob('Attr and Rep - New/PositionsAttrRepALLBIASEDVARYbetaeps%dbeta%dnvalue*.csv' % (eps, d))

            ## Files with 100 simulations
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # # for plotting gap
            # files = glob.glob('../Rep Only Chick Data Files/RepOnlyVARYDPositionsChickepseps%dD%dnvalue10.csv' % (eps, d))
            # files = glob.glob('../Attr Rep Chick Data Files/AttrRepVARYDPositionsChickeps%dD%dnvalue50.csv' % (eps, d))
            # files = glob.glob('../Attr Rep LEAD ONLY Chick Data Files/AttrRepLEADONLYVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))
            # files = glob.glob('../Attr Rep ALLBIASED Chick Data Files/AttrRepALLBIASEDVARYDPositionsChickeps%dD%dnvalue*.csv' % (eps, d))

            # Rep Only all biased
            # files = glob.glob('../CHICK DATA FINAL POSITIONS/RepOnlyALLBIASEDBiasedPositionsChickeps38D%dbeta%dnvalue*.csv'%(D,b))
            # Attr Rep all biased
            # files = glob.glob(
            #	'../CHICK DATA FINAL POSITIONS/AttrRepALLBIASEDBiasedPositionsChickeps38D%dbeta%dnvalue*.csv' % (D, b))
            # Attr Rep biased leaders
            files = glob.glob('../CHICK DATA FINAL POSITIONS/AttrRepLeadersBiasedPositionsChickeps38D%dbeta%dnvalue*.csv' % (D, b))
            # print(files)
            maximum = []
            ### Works well
            for file in files:
                X = readcsv(file)
                gap, reference_inertia, ondata_inertia = compute_gap(KMeans(), X, k_max)  # compute gap

                maxind = np.argmax(gap)
                # print('Max index', maxind)
                # check if maximum is sufficiently significant
                if gap[maxind] < 0.1:
                    maximum.append(0)
                else:
                    maximum.append(maxind)

            #
            # 	maximum.append(maxind)
            #
            print('D ', D, 'beta ', b)
            # print('Max index', maximum)
            # print('Just to check that I have 100', len(maximum))
            print('Number of non-zero entries', np.count_nonzero(maximum))
        ## gap statistic

# plot weights and gap
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

# print(X)
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
