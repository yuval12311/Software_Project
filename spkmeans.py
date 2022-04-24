import spkmeans
import sys
import pandas as pd
import numpy as np

def extract_data(filename):
    return pd.read_csv(filename, header=None).values.tolist()

def print_matrix(mat):
    for row in mat:
        print(','.join(f"{x:.4f}" for x in row))
def print_matrix_jacobi(mat):
    for row in mat:
        print(','.join(f"{x:.4f}" if abs(x)>=0.00005 else "0.0000" for x in row))

def input_error():
    print("Invalid Input!")
    exit(1)

def error():
    print("An Error Has Occurred")
    exit(1)

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

def eigengap(eigenvals):
     return max(range(1, len(eigenvals)//2), key=lambda i: eigenvals[i]-eigenvals[i-1])

def normalize(row):
    norm = sqrt(sum(map(row, lambda x: x**2)))
    for i, x in enumerate(row):
        row[x] = 0 if norm == 0 else x/norm

def spk(data, k):
    normalized = spkmeans.lnorm(data)
    eigenvecs, eigenvals = spkmeans.jacobi(normalized)
    eigen_dataframe = pd.DataFrame(eigenvecs, columns = eigenvals)

    eigen_dataframe.sort_index(axis=1, inplace=True)

    if k == 0: 
        k = eigengap(eigen_dataframe.columns)
        if k == 1: error()
    u = eigen_dataframe.iloc[:,:k].to_numpy()
    u = np.apply_along_axis(lambda row: row if (norm := np.linalg.norm(row)) == 0 else row/norm, 1, u)
    mu_indexes = kmeans_pp(k, u)
    print(','.join(f"{x}" for x in mu_indexes))
    u = u.tolist()
    mu = [u[i] for i in mu_indexes]
    mu = spkmeans.fit(k, 0, 300, mu, u)
    for x in mu:
        print(','.join(f"{xi:.4f}" for xi in x))



if __name__ == "__main__":
    np.random.seed(0)
    if len(sys.argv) != 4 : input_error()
    if not sys.argv[1].isdigit(): input_error()
    k = int(sys.argv[1])
    if k < 0 or k == 1: input_error()
    if sys.argv[2] not in ["spk", "wam", "ddg", "lnorm", "jacobi"]: input_error()
    data = extract_data(sys.argv[3])

    if sys.argv[2] == "wam":
        print_matrix(spkmeans.wam(data))
    if sys.argv[2] == "ddg":
        print_matrix(spkmeans.ddg(data))
    if sys.argv[2] == "lnorm":
        print_matrix(spkmeans.lnorm(data))
    if sys.argv[2] == "jacobi":
        if len(data) != len(data[0]): input_error()
        eigen = spkmeans.jacobi(data)
        print(",".join(f"{x:.4f}" if abs(x)>=0.00005 else "0.0000" for x in eigen[1]))
        print_matrix_jacobi(eigen[0])
        
    if sys.argv[2] == "spk":
        spk(data, k)
        


