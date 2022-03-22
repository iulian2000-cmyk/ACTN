import Crypto.Util.number
import math
import time


def gcd(a, b):
    if b == 0:
        return a
    else:
        return gcd(b, a % b)


def getPublicKey(fi, n):
    element = Crypto.Util.number.getRandomRange(1, fi)
    while gcd(element, fi) != 1:
        element = Crypto.Util.number.getRandomRange(1, fi)
        # print(gcd(element, fi))
    return n, element


def getPrivateKey(e, fi, n):
    d = pow(e, -1, fi) % fi
    # print('d',d)
    return n, d


print("\n PROBLEM I. RSA multi-prime \n ")
p = Crypto.Util.number.getPrime(512, randfunc=Crypto.Random.get_random_bytes)
q = Crypto.Util.number.getPrime(512, randfunc=Crypto.Random.get_random_bytes)
r = Crypto.Util.number.getPrime(512, randfunc=Crypto.Random.get_random_bytes)
n = p * q * r
fi = (p - 1) * (q - 1) * (r - 1)
print(f'p={p}')
print(f'r={r}')
print(f'q={q}')
print(f'n={n}')
print(f'fi(n)={fi}')

e = getPublicKey(fi, n)[1]
d = getPrivateKey(e, fi, n)[1]
m = 50
c = pow(m, e, n)
#
print(f'e={e}')
print(f'gcd(e,fi)={math.gcd(e, fi)}')
print(f'd={d}')
print(f'ed mod fi = {(e * d) % fi}')
print(f'm={m}')
print(f'c={c}')

print("\n Decryption using regular modular exponential algorithm \n")
start = time.time()
plaintext = pow(c, d, n)
end = time.time()
print(f'plaintext={plaintext}')
print('Time :', (end - start))
time_exponential = end-start

print("\n Decryption using Garner CRT \n")

# b_p = x_p , b_q = x_q , b_r = x_r

start = time.time()
b_p = pow(c % p, d % (p - 1), p)  # x_p
b_q = pow(c % q, d % (q - 1), q)  # x_q
b_r = pow(c % r, d % (r - 1), r)  # x_r

x_1 = b_p
alpha = ((b_q - x_1) * pow(p, -1, b_r)) % b_r
x_2 = x_1 + alpha * p
alpha = ((b_q - x_2) * pow(p * r, -1, b_q)) % b_q
x_3 = x_2 + alpha * p * q

# m_gauss = p * q * r
# c_p = m_gauss // p
# c_r = m_gauss // r
# c_q = m_gauss // q
# x_p = (pow(c_p % p, -1, p) * b_p) % p
# x_q = (pow(c_q % q, -1, q) * b_q) % q
# x_r = (pow(c_r % r, -1, r) * b_r) % r
# plaintext = (x_p * c_p + x_r * c_r + x_q * c_q) % m_gauss

plaintext = x_3
end = time.time()
print(f'plaintext={plaintext}')
print('Time :', (end - start))
time_garner = end - start

if time_garner == time_exponential :
    print("Same speed")
else:
    if time_garner < time_exponential:
        print("Garner algorithm is faster")
    else:
        print("Modular exponential is faster")


print("\n PROBLEM II. RSA multi-power \n ")
p = Crypto.Util.number.getPrime(512, randfunc=Crypto.Random.get_random_bytes)
q = Crypto.Util.number.getPrime(512, randfunc=Crypto.Random.get_random_bytes)
n = p * p * q
fi = p * (p - 1) * (q - 1)
print(f'p={p}')
print(f'r={r}')
print(f'q={q}')
print(f'n={n}')
print(f'fi(n)={fi}')

e = getPublicKey(fi, n)[1]
d = getPrivateKey(e, fi, n)[1]
m = 50
c = pow(m, e, n)
#
print(f'e={e}')
print(f'gcd(e,fi)={math.gcd(e, fi)}')
print(f'd={d}')
print(f'ed mod fi = {(e * d) % fi}')
print(f'm={m}')
print(f'c={c}')

print("\n Decryption using regular modular exponential algorithm \n")
start = time.time()
plaintext = pow(c, d, n)
end = time.time()
print(f'plaintext={plaintext}')
print('Time :', (end - start))
time_exponential = end-start
print("\n Decryption using Garner CRT \n")

start = time.time()
x_p2 = pow(c, d, n) % pow(p, 2)
x_q = pow(c, d, n) % q

x_0 = pow(c % p, d % (p - 1), p)
alpha = (c - pow(x_0, e, p * p)) // p

x_1 = (alpha * pow(e * (pow(x_0, e - 1, p * p)) % p, -1, p)) % p

x_p_square = p * x_1 + x_0

if x_p_square == x_p2:
    print("Hensel Lemma translated from modulo p to modulo p^2 ")

# m_gauss = p * p * q
# c_q = m_gauss // q
# c_p_square = m_gauss // (p * p)

b_q = x_q
b_p_square = x_p_square
# print(gcd(c_q,q),'-',gcd(c_p_square,p*p))

# x_sol_q = (pow(c_q % q, -1, q) * b_q) % q
# x_sol_p_square = (pow(c_p_square % (p * p), -1, p * p) * b_p_square) % (p * p)

# plaintext = (x_sol_q * c_q + x_sol_p_square * c_p_square) % m_gauss

x_1 = b_q
alpha = ((b_p_square - x_1) * pow(q, -1, b_p_square)) % b_p_square
x_2 = x_1 + alpha * q
plaintext = x_2


end = time.time()
time_garner = end - start
print(f'plaintext={plaintext}')
print(f'Time : {(end - start)}')

if time_garner == time_exponential :
    print("Same speed")
else:
    if time_garner < time_exponential:
        print("Garner algorithm is faster")
    else:
        print("Modular exponential is faster")

