# Based on:
# https://github.com/MRaczuk/padic, ver. 0.2.4
# mraczuk@gmail.com
# OSI Approved MIT License

from __future__ import annotations
from typing import Callable
from numpy import base_repr
from numpy.polynomial import Polynomial
from functools import lru_cache
# import sys

# sys.setrecursionlimit(2000)

class Qp:
    # Numbers are compared modulo p**PRECISION. Doesn't affect precision of computations.
    PRECISION: int = 32

    # Default precision for integer to Qp conversion.
    # This also may affect precision of some arithmetic operations.
    INTEGER_PRECISION: int = 64

    # Default value of p used in some functions for convenience so that one doesn't have to
    # keep inputting same prime over and over again. Can be set by set_prime function. Ignored
    # if None.
    DEFAULT_PRIME: int | None = None

    # Determines how many digits before period (on left) of a number will be displayed by default.
    # If set to None displays all of them. Should be non-negative.
    DISPLAY_PRECISION: int | None = None

    # Represents p-adic number as an interval s*(p^v) + O(p^N)
    def __init__(self, N: int, v: int, s: int, p: int | None = None) -> None:
        if p is None:
            p = Qp.DEFAULT_PRIME
        # Assumes p is prime
        # Assumes s _|_ p or s == 0
        if v >= N:
            s = 0
        if s == 0:
            v = N
        self.N: int = N
        self.v: int = v
        self.s: int = s % p ** (N - v)
        self.p: int = p

    def __abs__(self) -> int | float:
        return self.p ** (-self.v)

    def __eq__(self, other: Qp | int) -> bool:
        diff = self - other
        return diff.v >= Qp.PRECISION or diff.v == diff.N

    def __add__(self, other: Qp | int | float) -> Qp:
        if isinstance(other, float) and other == 0.0:
            return self
        if isinstance(other, Qp) and self.p == other.p:
            v_diff = self.v - other.v
            if v_diff > 0:
                return Qp(min(self.N, other.N), other.v, (self.p ** v_diff) * self.s + other.s, self.p)
            if v_diff < 0:
                return Qp(min(self.N, other.N), self.v, self.s + (self.p ** (-v_diff)) * other.s, self.p)
            return Qp.from_int(self.s + other.s, self.p, min(self.N, other.N), self.v)
        if isinstance(other, int):
            return self + Qp.from_int(other, self.p, self.N)
        else:
            raise RuntimeError(f"Can't add {str(self)} to {str(other)}")

    # This may look weird. That's a bypass to make class work with evaluation of numpy polynomials
    def __radd__(self, other: int | float) -> Qp:
        if isinstance(other, float) and other == 0.0:
            return self
        return self + other

    def __neg__(self) -> Qp:
        return Qp(self.N, self.v, -self.s, self.p)

    def __sub__(self, other: Qp | int) -> Qp:
        return self + (-other)

    def __rsub__(self, other: int) -> Qp:
        return -(self - other)

    def __mul__(self, other: Qp | int | float) -> Qp:
        if isinstance(other, float) and other == 1.0:
            return self
        if isinstance(other, Qp) and self.p == other.p:
            return Qp(min(self.v + other.N, other.v + self.N), self.v + other.v, self.s * other.s, self.p)
        if isinstance(other, int):
            return self * Qp.from_int(other, self.p, self.N + Qp.val(other, self.p) - Qp.val(self))
        else:
            raise RuntimeError(f"Can't multiply {self} with {other}")

    # This may look weird. That's a bypass to make class work with evaluation of numpy polynomials
    def __rmul__(self, other: int | float) -> Qp:
        if isinstance(other, float) and other == 1.0:
            return self
        return self * other

    def __truediv__(self, other: Qp | int) -> Qp:
        if isinstance(other, Qp) and self.p == other.p:
            N = min(self.v + other.N - 2 * other.v, self.N - other.v)
            v = self.v - other.v
            s = self.s * pow(other.s, -1, self.p ** (N - v))
            return Qp(N, v, s, self.p)
        if isinstance(other, int):
            return self / Qp.from_int(other, self.p, self.N + Qp.val(other, self.p) - Qp.val(self))
        else:
            raise RuntimeError(f"Can't divide {self} by {other}")

    def __rtruediv__(self, other: int) -> Qp:
        return Qp.from_int(other, self.p, self.N + Qp.val(other, self.p) - Qp.val(self)) / self

    def __mod__(self, other: Qp | int) -> Qp:
        if Qp.val(self) >= Qp.val(other, self.p):
            return Qp(max(Qp.INTEGER_PRECISION, self.N), max(Qp.INTEGER_PRECISION, self.N), 0, self.p)
        return Qp.from_int(self.s % self.p ** (min(Qp.val(other, self.p) - Qp.val(self), self.N)),
                              self.p, max(Qp.INTEGER_PRECISION, self.N), Qp.val(self))
    def __rmod__(self, other: int) -> Qp:
        return Qp.from_int(other, self.p, self.v) % self

    def __floordiv__(self, other: Qp | int) -> Qp:
        return (self - self % other) / other

    def __rfloordiv__(self, other: int) -> Qp:
        return (other - other % self) / self


    def wv(self):
        print(f's={self.s}')

        # digits = []
        # s = self.s
        # while True:
        #     (s, r) = divmod(s, self.p)
        #     out = f'{r}' + ("_" if len(out) > 1 else '') + out
        #     if s == 0:
        #         break











    def __str__(self) -> str:
        if self.p > 31:
            out = ''
            s = self.s
            while True:
              (s, r) = divmod(s, self.p)
              out = f'{r}' + ("_" if len(out) > 1 else '') + out
              if s == 0:
                break
            return f'{out} + O({self.p}^{self.N})'
        out = base_repr(self.s, self.p)
        if self.v >= 0:
            if Qp.DISPLAY_PRECISION is None:
                return out + ''.join(['0'] * self.v) + f' + O({self.p}^{self.N})'
            num = out + ''.join(['0'] * self.v)
            # Space is a workaround space cuz now (08.2023) 3-year-old bug causes pycharm not to
            # print ... if at the start of the string...
            return ((" ..." if Qp.DISPLAY_PRECISION < len(num) else "") +
                    num[-Qp.DISPLAY_PRECISION:] + f' + O({self.p}^{self.N})')
        else:
            if Qp.DISPLAY_PRECISION is None:
                return out[:self.v] + '.' + out[self.v:] + f' + O({self.p}^{self.N})'
            num = out[:self.v]
            # Same comment as above.
            return ((" ..." if Qp.DISPLAY_PRECISION < len(num) else "") +
                    num[-Qp.DISPLAY_PRECISION:] + '.' + out[self.v:] + f' + O({self.p}^{self.N})')

    def __repr__(self) -> str:
        return str(self)

    def __format__(self, format_spec: str) -> str:
        if format_spec == '':
            return str(self)
        if format_spec[0] == '.':
            digits = int(format_spec[1:])
            prev, Qp.DISPLAY_PRECISION = Qp.DISPLAY_PRECISION, digits
            out = str(self)
            Qp.DISPLAY_PRECISION = prev
            return out
        if format_spec in ['exact', 'all', 'a', '.None']:
            prev, Qp.DISPLAY_PRECISION = Qp.DISPLAY_PRECISION, None
            out = str(self)
            Qp.DISPLAY_PRECISION = prev
            return out
        return str(self)

    # Currently works for integer powers only. This may change in the future.
    # Modulo argument is for now ignored.
    def __pow__(self, power: int, modulo=None) -> Qp:
        assert isinstance(power, int)
        if power == 0:
            return Qp.from_int(1, self.p, max(Qp.INTEGER_PRECISION, self.N))
        if power == 1:
            return Qp(self.N, self.v, self.s, self.p)
        if power == -1:
            return 1 / self
        return self ** (power // 2) * self ** (power // 2) * self ** (power % 2)

    def __lshift__(self, other: int) -> Qp:
        assert isinstance(other, int)
        return self.p ** other * self

    def __rshift__(self, other: int) -> Qp:
        assert isinstance(other, int)
        return self // self.p ** other

    def __hash__(self) -> int:
        return self.s

    # Warning! Center isn't necessarily an integer!
    def center(self) -> int | float:
        return self.s * (self.p ** self.v)

    @staticmethod
    # No default prime on purpose.
    def val(n: Qp | int, p: int | None = None) -> int:
        if isinstance(n, Qp) and (n.p == p or p is None):
            return n.v
        if isinstance(n, int) and p is not None:
            if n == 0:
                return Qp.INTEGER_PRECISION
            out = 0
            while n % p == 0:
                out += 1
                n //= p
            return out
        raise RuntimeError("Valuation undefined for " + str(n), type(n))

    @staticmethod
    def _digit_value(c: str) -> int:
        o = ord(c)
        if ord('0') <= o <= ord('9'):
            return o - ord('0')
        if ord('A') <= o <= ord('Z'):
            return o - ord('A') + 10
        raise RuntimeError("Couldn't assign digit value for: " + c)

    @staticmethod
    def from_string(string: str, p: int | None = None) -> Qp:
        if p is None:
            p = Qp.DEFAULT_PRIME
        v = 0
        N = len(string) + 1
        non_zero_occurred = False
        dot_occurred = False
        stop = 0
        for i in range(len(string) - 1, -1, -1):
            char = string[i]
            if char == '.' and not dot_occurred:
                v -= len(string) - 1 - i
                N = i + 1
                dot_occurred = True
                continue
            if not char.isalnum():
                raise RuntimeError("Cannot parse string: " + string + " to p-adic integer.")
            if not non_zero_occurred and char != '0':
                non_zero_occurred = True
                stop = i
                v += len(string) - 2 - i if dot_occurred else len(string) - 1 - i
            if Qp._digit_value(char) >= p:
                raise RuntimeError("Cannot parse string: " + string + " to p-adic integer.")
        s = 0
        for char in string[:stop + 1]:
            if char == '.':
                continue
            s *= p
            s += Qp._digit_value(char)
        if s == 0:
            return Qp(N, N, 0, p)
        return Qp(N, v, s, p)

    @staticmethod
    # Calculates p^{v_adj}a + O(p^N) with a not necessarily coprime with p nor equal to 0
    # assumes that p is correct prime number
    def from_int(a: int, p: int | None = None, N: int | None = None, v_adj: int = 0) -> Qp:
        if p is None:
            p = Qp.DEFAULT_PRIME
        if N is None:
            N = Qp.INTEGER_PRECISION
        if a == 0:
            return Qp(N, N, a, p)
        v = Qp.val(a, p)
        return Qp(N, v + v_adj, a // (p ** v), p)

    @staticmethod
    # Creates p-adic number as fraction a/b. Doesn't check for corectness of given arguments
    def from_frac(a: int, b: int, p: int | None = None, N: int = None) -> Qp:
        if p is None:
            p = Qp.DEFAULT_PRIME
        if N is None:
            N = Qp.INTEGER_PRECISION
        return Qp.from_int(a, p, N) / Qp.from_int(b, p, N)


def gcd(a: int, b: int) -> int:
    if b == 0:
        return a
    return gcd(b, a % b)


class Rational:
    def __init__(self, p, q):
        g = gcd(p, q)
        self.num = p // g
        self.den = q // g

    def __add__(self, other):
        return Rational(self.num * other.den + self.den * other.num, self.den * other.den)

    def __mul__(self, other):
        return Rational(self.num * other.num, self.den * other.den)

    def __truediv__(self, other):
        return Rational(self.num * other.den, self.den * other.num)

    def __mod__(self, other):
        return Rational(self.num % other, self.den % other)

    def __neg__(self):
        return Rational(-self.num, self.den)

    def __str__(self):
        return str(self.num) + "/" + str(self.den)


def series(a: Callable[[int], int | Rational | Qp], n: int, z='p') -> \
                                    Callable[[int | Rational], Rational] \
                                    | Callable[[int | Qp, int | None], Qp]:
    def u1(x: int | Rational):
        if isinstance(x, int):
            x = Rational(x, 1)
        out = Rational(0, 1)
        for k in range(n, -1, -1):
            out *= x
            out += a(k)
        return out

    def u2(x: int | Qp, p: int | None):
        if p is None:
            p = x.p
        out = Qp.from_int(0, p)
        for k in range(n, -1, -1):
            out *= x
            out += a(k)
        return out

    if z == 'r':
        return lambda x: u1(x)
    else:
        return lambda x, p=None: u2(x, p)


@lru_cache(maxsize=1000)
def factorial(n):
    return n * factorial(n - 1) if n else 1


@lru_cache(maxsize=1000)
def binomial_coeff(a: Qp | int, b: int, p: int | None = None) -> Qp:
    if p is None:
        p = a.p
    return Qp.from_int(1, p) if b == 0 else binomial_coeff(a, b - 1, p) / b * (a - b + 1)


# Convergent for x = 1 + O(p)
def log(x: int | Qp, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        p = x.p
    return -series(lambda n: Qp.from_frac(1, n, p) if n != 0 else Qp.from_int(0, p), N)(1 - x, p)


# Convergent for |x|_p < p^{-1/{p-1}}
def exp(x: int | Qp, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        p = x.p
    return series(lambda n: Qp.from_frac(1, factorial(n), p), N)(x, p)


# Convergent for |x|_p < p^{-1/{p-1}}
def sin(x: int | Qp, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        p = x.p
    return series(lambda n: Qp.from_frac(1, factorial(n), p) * (-1) ** ((n - 1) // 2) if n % 2 else 0, N)(x, p)


# Convergent for |x|_p < p^{-1/{p-1}}
def cos(x: int | Qp, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        p = x.p
    return series(lambda n: Qp.from_frac(1, factorial(n), p) * (-1) ** ((n + 1) // 2) if (n - 1) % 2 else 0, N)(x, p)


# Convergence radius dependent on p and a
def binomial(x: int | Qp, a: int | Qp, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        if isinstance(x, Qp):
            p = x.p
        p = a.p
    return series(lambda n: binomial_coeff(a, n, p), N)(x - 1, p)


# Finds approximate root of a polynomial
# or doesn't.
# Currently checks for roots in Z_p only.
def find_approx_root(poly: Polynomial, p: int | None = None, depth: int = 5) -> Qp:
    der = poly.deriv()

    def check_path(x, i):
        v1 = Qp.val(poly(x), p)
        v2 = Qp.val(der(x), p)
        return 1 if v1 > 2 * v2 else 2 if v1 <= i else 0

    if p is None:
        p = poly.coef[0].p
    previous_paths = []
    current_paths = [Qp.from_int(0, p)]
    ppow = 1
    for i in range(depth):
        previous_paths = current_paths
        current_paths = []
        for path in previous_paths:
            for k in range(p):
                current_paths.append(path + k * ppow)
        for path in current_paths:
            res = check_path(path, i)
            if res == 1:
                return path
            if res == 2:
                current_paths.remove(path)
        ppow *= p
    raise RuntimeError("Root not found!")


# Calculates polynomial root closest to a given approximate root
# In case no root is given returns arbitrary root or raises NoRoot exception.
# Doesn't verify correctness of given approximate root.
# Note that in fact GHL (Generalised Hensel Lemma) for polynomial roots is used
# Currently checks for roots in Z_p only
def hensel(poly: Polynomial, approx: Qp | int | None = None, p: int | None = None, N: int = 100) -> Qp:
    if p is None:
        p = approx.p
    if approx is None:
        approx = find_approx_root(poly)
    if isinstance(approx, int):
        approx = Qp.from_int(approx, p)
    der = poly.deriv()
    out = approx
    for _ in range(N):
        out -= poly(out) / der(out)
    return out

xx = Qp.from_int(7,3) / 13
print(xx)
xx.wv()
