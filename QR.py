import numpy as np


A = np.random.randint(1, 100,(3,3))
print("original matrix:")
print(A)

def householder(A):
    m,n=np.shape(A)
    Q = np.identity(m)
    R = np.copy(A)
    for i in range(m - 1):
        x = R[i:, i]
        e = np.zeros_like(x)
        e[0] = np.linalg.norm(x)
        u = x - e
        v = u / np.linalg.norm(u)
        H1 = np.identity(m)
        H1[i:, i:] -= 2.0 * np.outer(v, v)
        R = np.dot(H1, R)
        Q = np.dot(Q, H1)
    return (Q, R)


Q,R=householder(A)
print("found matrix Q:")
print(Q)
print("found matrix R:")
print(R)

q,r=np.linalg.qr(A)
print("check with np.linalg.qr matrix Q:")

print(q)

print("check with np.linalg.qr matrix R:")
print(r)





B=Q@R
for i in range(0,1000):
    B=R@Q
    Q,R=householder(B)

for i in range(0,3):
    print("eigenvalues ",i+1,B[i,i])




