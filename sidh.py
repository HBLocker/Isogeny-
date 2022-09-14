
from io import StringIO
from random import randint

class Complex(object):
	def __init__(self, real, imag=0):
		self.re = int(real)
		self.im = int(imag)
        
	def __add__(self, other):
		return Complex(self.re + other.re, self.im + other.im)

	def __sub__(self, other):
		return Complex(self.re - other.re, self.im - other.im)

	def __mul__(self, other):
		ab0=self.re*other.re
		ab1=self.im*other.im
		c=(self.re+self.im)*(other.re+other.im)
		return Complex((ab0-ab1)%p, (c-ab0-ab1)%p)

	def __neg__(self):  
		return Complex(-self.re, -self.im)

	def __eq__(self, other):
		return self.re == other.re and self.im == other.im

	def __ne__(self, other):
		return not self.__eq__(other)

	def __str__(self):
		return '(%u, %u)' % (self.re %p, self.im %p)

	def __repr__(self):
		return 'Complex(' + str(self.re) +  ', ' + str(self.im) + ')'

	def __pow__(self, power): #only squares required
		return Complex(((self.re+self.im)*(self.re-self.im))%p, (2*self.im*self.re)%p)
		
	def __mod__(self, p):
		return Complex(self.re % p, self.im % p)	



def j_invariant(A,C):

    """ The invariant function used. 
    Ref
    https://github.com/microsoft/PQCrypto-SIDH/blob/75ed5b09bd06a19cdad7660bef10e820074f808f/src/ec_isogeny.c#L247

    Args:
        A (Int): Projective of curve
        C (int): Projective of curve

    Returns output J-invariant
    """

    jinv = A**2  
    t1 = C**2
    t0 = t1 + t1
    t0 = jinv - t0
    t0 = t0 - t1
    jinv = t0 - t1
    t1 = t1**2
    jinv = jinv * t1
    t0 = t0 + t0
    t0 = t0 + t0
    t1 = t0**2
    t0 = t0 * t1
    t0 = t0 + t0
    t0 = t0 + t0
    jinv = inv(jinv)
    jinv = t0 * jinv

    return jinv


def xDBLADD(XP, ZP, XQ, ZQ, xPQ):
    """
    Montgomery Addition
    Input: projective Montgomery points xP=XP/ZP and xQ=XQ/ZQ
      affine difference x(P-Q)=xPQ

    Args:
        XP (int): 
        ZP (int): 
        XQ (int): 
        ZQ (int):
        xPQ (int): 

    Returns:
      projective coordinates x(Q+P)=XQP/XZP
    """

    t0 = XP + ZP   
    t1 = XP - ZP
    XP = XQ - ZQ
    ZP = XQ + ZQ
    t0 = (XP * t0)%p
    t1 = (ZP * t1)%p
    ZP = t0 - t1
    XP = t0 + t1
    ZP = (ZP**2)%p
    XQP = (XP**2)%p
    ZQP = (xPQ * ZP)%p
	
    return XQP, ZQP  


