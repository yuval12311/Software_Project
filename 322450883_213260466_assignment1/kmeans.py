import math
import sys
epsilon = 0.001

def input_error():
    print("Invalid Input!")
    exit(1)

def error():
    print("An Error Has Occurred")
    exit(1)

def read(f):
    return [[float(s) for s in line.split(",")] for line in f]


def closest_centroid(x, mu):
    min_i = 0
    min_dist = sqr_dist(x, mu[0])
    for i in range(len(mu)):
        if min_dist > sqr_dist(x, mu[i]):
            min_i = i
            min_dist = sqr_dist(x, mu[i])
    return min_i



def sqr_dist(x, y):
    return sum((x[i] - y[i])**2 for i in range(len(x)))

def add(x, y):
    return [x[i] + y[i] for i in range(len(x))]


def mult(a, x):
    return [a * x[i] for i in range(len(x))]


def add_to_cluster(cluster_sums, i, x):
    for j in range(len(x)):
        cluster_sums[i][j] += x[j]




def kmeans(data, k, maxiter):
    mu = data[:k]
    finished = False
    iter = 0
    while not finished:
        cluster_sizes = [0 for _ in mu]
        cluster_sums = [[0 for _ in mu[0]] for _ in mu]
        for x in data:
            j = closest_centroid(x, mu)
            cluster_sizes[j] += 1
            add_to_cluster(cluster_sums, j, x)
            

        iter += 1
        finished = True
        for i in range(len(mu)):
            newmui = mult(1/cluster_sizes[i], cluster_sums[i])
            finished &= sqr_dist(mu[i], newmui) < epsilon**2
            mu[i] = newmui
        finished |= iter >= maxiter
    return mu


def write(f, mu):
    for x in mu:
        f.write(','.join(f"{xi:.4f}" for xi in x) + "\n")



if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5 : input_error()
    maxiter = 200
    has_maxiter = 0 if len(sys.argv) == 4 else 1
    if not sys.argv[1].isdigit(): input_error()
    k = int(sys.argv[1])
    if k <= 1: input_error()
    try:
        f = open(sys.argv[2 + has_maxiter], 'r')
    except:
        error()
    maxiter = 200
    if has_maxiter and not sys.argv[2].isdigit(): input_error()
    elif has_maxiter: maxiter = int(sys.argv[2])
    if maxiter <= 0: input_error()
    data = read(f)
    if k >= len(data): input_error()
    f.close()
    mu = kmeans(data, k, maxiter)
    try:
        f = open(sys.argv[3 + has_maxiter], 'w')
    except:
        error()
    write(f, mu)
    f.close()

