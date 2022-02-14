from sklearn.cluster import KMeans
from sklearn import datasets
import matplotlib.pyplot as plt

iris = datasets.load_iris()

r = 10
plt.plot(range(1,11) ,[KMeans(n_clusters=k, random_state=0).fit(iris['data']).inertia_ for k in range(1, 11)])
plt.xticks(range(1,11))
plt.xlabel("k")
plt.ylabel("inertia")
plt.annotate('Elbow point', xy=(3, 80), xytext=(4.5, 250), arrowprops=dict(facecolor='black', shrink=0.05))
plt.plot(3, 80, 'o', ms=r*2, mec='orange', mfc='none', mew=2)

plt.savefig("elbow.png")
