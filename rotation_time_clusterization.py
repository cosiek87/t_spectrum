# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 20:47:53 2021

@author: PC_Shark
"""
import numpy as np
import matplotlib.pyplot as pl
from sklearn.cluster import KMeans


liczba_obrotow = 59

pl.rcParams["font.family"] = "Arial"
pl.rcParams['figure.dpi'] = 300
pl.rcParams["legend.fancybox"] = False
pl.rcParams["legend.framealpha"] = 1

czas_kmeans = np.loadtxt("czas_do_klastra_heart_24.09.txt")
kroki_kmeans = np.array([1 for i in range(len(czas_kmeans))])
data = np.array([czas_kmeans,kroki_kmeans]).transpose()


kmeans = KMeans(init="k-means++",n_clusters=liczba_obrotow, n_init=20,max_iter=300,
                algorithm='elkan',tol=1e-6)
kmeans.fit(data)


pl.scatter(czas_kmeans,kroki_kmeans)
pl.plot(kmeans.cluster_centers_.transpose()[0],kmeans.cluster_centers_.transpose()[1],linewidth=0,c='r',marker='.')

kroki_clustered = kmeans.cluster_centers_.transpose()[0]

print(len(kroki_clustered),len(kmeans.labels_))

start_stop_times = np.array([[np.min(czas_kmeans[kmeans.labels_==i]),np.max(czas_kmeans[kmeans.labels_==i])] for i in set(kmeans.labels_)]).transpose()

times = np.hstack((start_stop_times[0],start_stop_times[1]))


print(np.sort(times)/1e12)
np.savetxt("czas_do_analizy_heart_24.09.txt",np.sort(times))
times = np.insert(times,0,0)
np.savetxt("os_czasu_heart_24.09.txt",np.sort(times)/1E12)
