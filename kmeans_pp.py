import pandas as pd
import numpy as np
import sys
import kmeans



def read(f1, f2):
    df1 = pd.read_csv(f1, header=None)
    df2 = pd.read_csv(f2, header=None)
    merged = pd.merge(df1, df2, on=0, sort=True)
    merged = merged.drop(0, axis=1)
    return merged.to_numpy()


def kmeans_pp(k, data):
    n = data.shape[0]
    mu = []
    np.random.seed(0)
    mu.append(np.random.choice(n))
    while len(mu) < k:
        D = np.apply_along_axis(lambda x : min(np.inner(x-data[i], x-data[i]) for i in mu), 1, data)
        p = (1/sum(D))*D
        mu.append(np.random.choice(n, p=p))
    return mu



def input_error():
    print("Invalid Input!")
    exit(1)

def error():
    print("An Error Has Occurred")
    exit(1)

if __name__ == "__main__":
    if len(sys.argv) < 5 or len(sys.argv) > 6 : input_error()
    maxiter = 300
    has_maxiter = 0 if len(sys.argv) == 5 else 1
    if not sys.argv[1].isdigit(): input_error()
    k = int(sys.argv[1])
    if k <= 1: input_error()
    eps = float(sys.argv[2+ has_maxiter])
    if eps < 0: input_error()
    
    f1 = sys.argv[3 + has_maxiter]
    f2 = sys.argv[4 + has_maxiter]

    if has_maxiter and not sys.argv[2].isdigit(): input_error()
    elif has_maxiter: maxiter = int(sys.argv[2])
    if maxiter <= 0: input_error()

    data = read(f1, f2)
    mu = kmeans_pp(k, data)
    print(','.join(f"{x}" for x in mu))
    data = data.tolist()
    mu = [data[i] for i in mu]
    mu = kmeans.fit(k, eps, maxiter, mu, data)
    for x in mu:
        print(','.join(f"{xi:.4f}" for xi in x))
    print()


    
