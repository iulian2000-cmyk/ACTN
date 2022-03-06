import itertools
import Crypto.Util.number
import sys
import random
import math
import collections
import numpy
import array
import time


def P_x(x, m, p):
    res = 0
    for i in range(0, len(m)):
        res = res + (m[i] * math.pow(x, i + 1)) % p
    return res % p


def find_lcm(num1, num2):
    if num1 > num2:
        num = num1
        den = num2
    else:
        num = num2
        den = num1
    rem = num % den
    while rem != 0:
        num = den
        den = rem
        rem = num % den
    gcd = den
    lcm = int(int(num1 * num2) / int(gcd))
    return lcm


def poly_horner(A, x):
    if len(A) == 2:
        return x * (A[0] + x * A[1])
    else:
        p = A[-1]
        i = len(A) - 2
        while i >= 0:
            p = p * x + A[i]
            i -= 1
        return p


big_p = Crypto.Util.number.getPrime(161, randfunc=Crypto.Random.get_random_bytes)
print("FIRST TEST ")

print("Encoding")
print("Test instances simple encoding ")
k = 3
p = 11
s = 1

m = [7, 2]

n = k + 2 * s
i = 1
y = []
while i <= n:
    y.append(P_x(i, m, p))
    i = i + 1

print(f"y={y}")
print(f"n={n}")
print(f"k={k}")
print(f"m={m}")
print(f"p={p}")
print(f"s={s}")

print("Testing using Horner ")

k = 3
p = 11
s = 1
m = [7, 2]
n = k + 2 * s
i = 1
y = []
while i <= n:
    y.append(poly_horner(m, i) % p)
    i = i + 1
print(f"y={y}")


def calculate_fc_A(subset, z, p):
    start = time.time()
    result = 0
    for i in subset:
        produs = 1
        for j in subset:
            if j != i:
                produs = (produs * (j / (j - i))) % p
        result = result + (z[i - 1] * produs) % p
    end = time.time()
    return round(result % p), (end - start)


def calculate_fc_B(subset, z, p):
    s = time.time()
    subset_A = set(subset)
    products = []
    for i in subset_A:
        products.append((numpy.prod(list(subset_A - {i})) / (
                numpy.prod(list(subset_A - {i})) - i * sum(list(subset_A - {i})) + i * i)) % p)

    result = 0
    for i in subset_A:
        result = result + (products[list(subset_A).index(i)] * z[i - 1]) % p

    b = time.time()
    return round(result % p), (b - s)


def inverse_modular(x, p):
    result = -1
    for val in range(0, p):
        if (x * val) % p == 1:
            result = val
    return result


def calculate_fc_C(subset, z, p):
    a = time.time()
    sum_elements = 0
    vec_denominators = []
    vec_numerators = []
    for element in subset:
        prod_denominator = 1
        prod_numerator = 1
        for el_j in subset:
            if el_j != element:
                prod_numerator = prod_numerator * el_j
                prod_denominator = prod_denominator * (el_j - element)
        prod_numerator = (prod_numerator * z[element - 1]) % p
        vec_numerators.append(prod_numerator)
        vec_denominators.append(prod_denominator)

    lcm = find_lcm(vec_denominators[0], vec_denominators[1])
    for i in range(2, len(vec_denominators)):
        lcm = find_lcm(lcm, vec_denominators[i])

    aux = lcm

    for element in range(len(vec_denominators)):
        vec_numerators[element] = (vec_numerators[element] * (aux / vec_denominators[element])) % p
        vec_denominators[element] = aux
        sum_elements = sum_elements + vec_denominators[element]

    if sum_elements % aux == 0:
        b = time.time()
        return (sum_elements / aux) % p, (b - a)
    else:
        if inverse_modular(aux, p) != -1:
            sum_elements = sum_elements * inverse_modular(aux, p)
            b = time.time()
            return sum_elements, (b - a)


def multiply_polynom(A, B, m, n, p):
    prod = [0] * (m + n - 1);
    for i in range(m):
        for j in range(n):
            prod[i + j] += (A[i] * B[j]) % p
    return prod


def add_polynom(A, B, m, n, p):
    size = max(m, n);
    sum = [0 for i in range(size)]
    for i in range(0, m, 1):
        sum[i] = A[i]
    for i in range(n):
        sum[i] = (sum[i] + B[i]) % p
    return sum


def find_message(z, subset, p):
    m = []
    ret = []
    for i in subset:
        coeficients_x = []
        free_terms = []
        for j in subset:
            if j != i:
                coeficients_x.append(int(1 / (i - j)) % p)
                free_terms.append(int(j / (i - j)) % p)
        polynomials = []
        for index in range(0, len(coeficients_x)):
            polynomials.append([free_terms[index], coeficients_x[index]])
        polynom_z = []
        temp_polynom = [1, 0]
        for polynom in polynomials:
            temp_polynom = multiply_polynom(temp_polynom, polynom, len(temp_polynom), len(polynom), p)
        polynom_z = temp_polynom
        for index in range(0, len(polynom_z)):
            polynom_z[index] = (polynom_z[index] * z[i - 1]) % p
        m.append(polynom_z)

    temp_polynom = [0, 0]
    for polynom in m:
        temp_polynom = add_polynom(temp_polynom, polynom, len(temp_polynom), len(polynom), p)

    return temp_polynom


print()
print("Decoding")
z = []
z = y
random_index = random.randint(0, len(z) - 1)
print("Error at index :" + str(random_index))
old_val = z[random_index]
while z[random_index] == old_val:
    z[random_index] = random.randint(0, p - 1)

print(f"z={z}")
index_set = [random_index]
A = set([i + 1 for i in range(0, n)])
subsets = itertools.combinations(A - set(index_set), k)
fc_values_A = []
time_A = []
fc_values_B = []
time_B = []
fc_values_C = []
time_C = []
for set_input in subsets:
    fc_values_A.append(calculate_fc_A(set_input, z, p)[0])
    time_A.append(calculate_fc_A(set_input, z, p)[1])

    fc_values_B.append(calculate_fc_B(set_input, z, p)[0])
    time_B.append(calculate_fc_B(set_input, z, p)[1])

    fc_values_C.append(calculate_fc_C(set_input, z, p)[0])
    time_C.append(calculate_fc_C(set_input, z, p)[1])

    if calculate_fc_A(set_input, z, p)[0] == 0 or calculate_fc_B(set_input, z, p)[0] == 0:
        print(f"m={find_message(z, set_input, p)}")

print(f"Free coeficients using {k * (k - 1)} inversions : {fc_values_A}", f" TIME : {sum(time_A)}")
print(f"Free coeficients using {k} inversions : {fc_values_B}", f" TIME : {sum(time_B)}")
print(f"Using one inversion :{fc_values_C}", f"TIME : {sum(time_C)}")

print()
print()
print("SECOND TEST ")

second_message = []
for i in range(0, 2):
    second_message.append(Crypto.Util.number.getRandomRange(0, p, randfunc=None))

y = []
while i <= n:
    y.append(poly_horner(second_message, i) % p)
    i = i + 1

print("Encoding")
print(f"n={n}")
print(f"k={k}")
print(f"m={second_message}")
print(f"p={big_p}")
print(f"s={s}")
print(f"y={y}")
print()

print("Decoding")
z = []
z = y
random_index = random.randint(0, len(z) - 1)
print("Error at index :" + str(random_index))
old_val = z[random_index]
while z[random_index] == old_val:
    z[random_index] = Crypto.Util.number.getRandomRange(0, p)

print(f"z={z}")
index_set = [random_index]
A = set([i + 1 for i in range(0, n)])
subsets = itertools.combinations(A - set(index_set), k)
fc_values_A = []
time_A = []
fc_values_B = []
time_B = []
fc_values_C = []
time_C = []
for set_input in subsets:
    fc_values_A.append(calculate_fc_A(set_input, z, p)[0])
    time_A.append(calculate_fc_A(set_input, z, p)[1])

    fc_values_B.append(calculate_fc_B(set_input, z, p)[0])
    time_B.append(calculate_fc_B(set_input, z, p)[1])

    fc_values_C.append(calculate_fc_C(set_input, z, p)[0])
    time_C.append(calculate_fc_C(set_input, z, p)[1])

    if calculate_fc_A(set_input, z, p)[0] == 0 or calculate_fc_B(set_input, z, p)[0] == 0:
        print(f"m={find_message(z, set_input, p)}")

print(f"Free coeficients using {k * (k - 1)} inversions : {fc_values_A}", f" TIME : {sum(time_A)}")
print(f"Free coeficients using {k} inversions : {fc_values_B}", f" TIME : {sum(time_B)}")
print(f"Using one inversion :{fc_values_C}", f"TIME : {sum(time_C)}")
