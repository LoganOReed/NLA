import numpy as np
import scipy as sp
from oct2py import octave
import matplotlib.pyplot as plt

def cg_bounds(N, k, tol=1e-14, M=None):
    """Generates values for testing specific bounds for convergence of CG."""
    a = k
    rho = (np.sqrt(np.float128(k)) - np.float128(1)) / (np.sqrt(np.float128(k)) + np.float128(1))
    sigma = np.linspace(1,a,N)
    A = np.diag(sigma)
    b = np.ones(N, dtype=np.float128) / np.sqrt(N, dtype=np.float128) 
    x_true = np.ones(N, dtype=np.float128) / (np.sqrt(N, dtype=np.float128) * sigma)
    e0 = np.sqrt(x_true.T @ A @ x_true)

    enorms = []

    def calculate_em(xk):
        """use xk at each iteration to get ek."""
        ek = (np.sqrt((x_true - xk).T @ A @ (x_true - xk))) / e0
        enorms.append(ek)

    xm, _ = sp.sparse.linalg.cg(A,b,x0=np.zeros(N), rtol=tol, M=M, callback=calculate_em)
    # rho_ms = np.array([rho ** i for i in range(N)])
    rho_ms = np.power(rho, np.arange(len(enorms), dtype=np.float128))

    bounds = enorms <= rho_ms
    print("Comparison of vectors:")
    for i, (v1, v2, comp) in enumerate(zip(enorms, rho_ms, bounds)):
        print(f"Iteration {i}: {v1} <= {v2} -> {comp}")

    return enorms, rho_ms
    
if __name__ == "__main__":
    cond = [2, 10, 50, 100, 1000]
    color = ["blue", "red", "green", "yellow", "black"]
    plt.figure(figsize=(8, 6))
    for i in range(len(cond)):
        em, rho = cg_bounds(50, cond[i])
        iter = np.arange(len(em))
        plt.plot(iter, rho, label=f'rho^m (k={cond[i]})', linestyle=':', color=color[i])  # Dotted line for rho_ms
        plt.plot(iter, em, label=f'||e_m||/||e_0|| (k={cond[i]})', linestyle='-', color=color[i])  # Solid line for enorms

    # plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.title('Comparison of enorms and rho_ms')
    plt.legend()
    plt.savefig(f"figures/{i}.png")
