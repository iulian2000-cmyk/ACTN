import Crypto.Util.number
import math
import time


def binaryToDecimal(binary):
    decimal = 0
    for digit in binary:
        decimal = decimal * 2 + int(digit)
    return decimal


def jacobi_symbol(a, n):
    if n <= 0:
        raise ValueError("'n' must be a positive integer.")
    if n % 2 == 0:
        raise ValueError("'n' must be odd.")
    a %= n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a /= 2
            n_mod_8 = n % 8
            if n_mod_8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    if n == 1:
        return result
    else:
        return 0


def solovay_strassen_prime_test(n, t):
    for i in range(1, t):
        a = Crypto.Util.number.getRandomRange(2, n - 2)
        r = pow(a, (n - 1) // 2, n)
        if r != 1 and r != n - 1:
            return False
        s = jacobi_symbol(a, n)
        if r != s % n:
            return False

    return True


print("Testing Jacobi Symbol algorithm AFCS:")
for i in [5, 7, 11, 13, 17, 19, 23, 29, 31, 33, 25, 21, 27, 39]:
    print(f'instance:{i}:result:{jacobi_symbol(i, 19)}')

print()
print("Testing Solovay Strassen algorithm:")
for i in [5, 7, 11, 13, 17, 19, 23, 29, 31, 33, 25, 21, 27, 39]:
    print(f'instance:{i}:result:{solovay_strassen_prime_test(i, 123)}')


def mersenne_reduction(a, M, n):
    length_M = M.bit_length()
    length_a = a.bit_length()
    if length_a < 2 * length_M:
        binary_a = '0' * (2 * length_M - length_a) + str(bin(a)).replace("0b", "")
        a0 = int(binary_a[:length_M], 2)
        a1 = int(binary_a[length_M:], 2)
        return (a1 + a0) % M
    elif length_a == 2 * length_M:
        binary_a = str(bin(a)).replace("0b", "")
        a0 = int(binary_a[:length_M], 2)
        a1 = int(binary_a[length_M:], 2)
        return (a1 + a0) % M


def lucas_lehmer_test(p):
    s = 4
    M = pow(2, p) - 1
    start = time.time()
    for i in range(p - 2):
        s = (s ** 2 - 2) % M
    if s == 0:
        print(f'{M} - is prime ')
    else:
        print(f'{M} - is composite')
    end = time.time()
    print(f'Normal reduction time :{end - start}')

    s = 4
    start = time.time()
    for i in range(p - 2):
        s = mersenne_reduction(s ** 2 - 2, M, p)

    if s == 0:
        print(f'{M} - is prime ')
    else:
        print(f'{M} - is composite')
    end = time.time()
    print(f'Dedicated reduction time :{end - start}')


print()
print("Testing Lucas-Lehmer algorithm :")
for i in [5, 7, 11, 13, 17, 19, 23, 29, 31, 33, 25, 21, 27, 39]:
    lucas_lehmer_test(i)
