"""Exact multiquadratic arithmetic without factorization or symbolic simplification.

Represent sums of (a+i b)*sqrt(n), a,b rational, n positive integers.
Two radicands merge exactly when their ratio is a rational square.
No BFB construction or numerical tolerance is used here.
"""
from fractions import Fraction as F
from math import isqrt

def rational_sqrt(x):
    x=F(x)
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b) if a*a==x.numerator and b*b==x.denominator else None

class Radical:
    __slots__=('terms',)
    def __init__(self,value=0):
        self.terms={}
        if isinstance(value,Radical):self.terms=value.terms.copy()
        elif value:self.terms[1]=(F(value),F(0))
    @classmethod
    def sqrt(cls,x):
        x=F(x)
        if x<0:raise ValueError('Negative real radicand')
        out=cls()
        if x:out._add(x.numerator*x.denominator,F(1,x.denominator),F(0))
        return out
    @classmethod
    def imaginary(cls,x=1):
        out=cls()
        if x:out.terms[1]=(F(0),F(x))
        return out
    def _add(self,n,a,b):
        if not a and not b:return
        root=isqrt(n)
        if root*root==n:n,a,b=1,a*root,b*root
        for k in list(self.terms):
            factor=rational_sqrt(F(n,k))
            if factor is not None:
                x,y=self.terms[k];x+=factor*a;y+=factor*b
                if x or y:self.terms[k]=(x,y)
                else:del self.terms[k]
                return
        self.terms[n]=(a,b)
    def __add__(self,other):
        other=other if isinstance(other,Radical) else Radical(other);out=Radical(self)
        for k,(a,b) in other.terms.items():out._add(k,a,b)
        return out
    __radd__=__add__
    def __neg__(self):
        out=Radical();out.terms={k:(-a,-b) for k,(a,b) in self.terms.items()};return out
    def __sub__(self,other):return self+-Radical(other)
    def __rsub__(self,other):return Radical(other)+-self
    def __mul__(self,other):
        other=other if isinstance(other,Radical) else Radical(other);out=Radical()
        for n,(a,b) in self.terms.items():
            for m,(c,d) in other.terms.items():out._add(n*m,a*c-b*d,a*d+b*c)
        return out
    __rmul__=__mul__
    def __truediv__(self,other):return self*(F(1)/F(other))
    def conjugate(self):
        out=Radical();out.terms={k:(a,-b) for k,(a,b) in self.terms.items()};return out
    def rational(self):
        if not self.terms:return F(0)
        if set(self.terms)=={1} and not self.terms[1][1]:return self.terms[1][0]
        raise ArithmeticError('Expression is not real rational')
    def __bool__(self):return bool(self.terms)
    def __eq__(self,other):return not self-Radical(other)
    def payload(self):return [{'radicand':str(k),'real_coefficient':str(a),'imag_coefficient':str(b)} for k,(a,b) in sorted(self.terms.items())]
    def text(self):
        if not self.terms:return '0'
        return ' + '.join(f'({a} + ({b})*I)*sqrt({k})' if b else f'({a})*sqrt({k})' for k,(a,b) in sorted(self.terms.items()))
