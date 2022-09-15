
from io import StringIO
from random import randint

#ref  Ref https://github.com/microsoft/PQCrypto-SIDH/blob/75ed5b09bd06a19cdad7660bef10e820074f808f/src/ec_isogeny.c


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


def XDBL(X,Z,A24,C26):

  """
   Doubling of a Montgomery point in
  projective coordinates (X:Z)

  Args:
    projective coordinates xP=X/Z
     curve constant A24/C24 = (A/C+2)/4

      X (int): _description_
      Z (int): _description_
      A24 (int): _description_
      C26 (int): _description_

  Returns:
      _type_: 4M+2S+4a
  """
  
  t0 = X - Z      #code from msr
  t1 = X + Z
  t0 = t0**2
  t1 = t1**2 
  Z2 = C24 * t0
  X2 = Z2 * t1
  t1 = t1 - t0
  t0 = A24 * t1
  Z2 = Z2 + t0
  Z2 = Z2 * t1
	
  return X2, Z2   #cost: 4M+2S+4a

def edDBL(Y,Z,AE,DE):
	t0=Y**2
	t1=Z**2
	t2=t0+t1
	t2=t2**2
	t0=t0**2
	t1=t1**2
	t2=t2-t0
	t2=t2-t1
	Y2=t2*AE  
	t2=t2*DE  
	t0=t0*DE  
	t1=t1*AE  
	t0=t0+t1  
	Z2=t0-t2  
	Y2=Y2-t0 

	return Y2,Z2

def edDBLe(XP,ZP,A,C,e):

  """
  Computes [2^e](X:Z) on Montgomery curve with
   projective constant via e repeated doublings

  Returns:
      _type_: projective Montgomery x-coordinates Q <- (2^e)*P.
  """
  
  C2=C+C
  AE=A+C2
  DE=A-C2

  YeP = XP-ZP
  ZeP = ZP+XP
	
  for i in range(0,e):
    YeP, ZeP = edDBL(YeP, ZeP, AE, DE)
	
  XeP=YeP+ZeP
  ZeP=ZeP-YeP	
  return XeP, ZeP	


def xDBLADD(XP,ZP,XQ,ZQ,xPQ,A24):

  """Simultaneous doubling and differential addition.
    affine difference x(P-Q)=xPQ and curve constant A24=(A+2)/4.	

  Returns:
  projective coordinates of x(2P)=X2P/Z2P
      
  """
  t0 = XP + ZP                 
  t1 = XP - ZP 
  X2P = t0**2
  t2 = XQ - ZQ
  XQP = XQ + ZQ
  t0 = t0 * t2
  Z2P = t1**2
  t1 = t1 * XQP
  t2 = X2P - Z2P
  X2P = X2P * Z2P
  XQP = A24 * t2
  ZQP = t0 - t1
  Z2P = XQP + Z2P
  XQP = t0 + t1
  Z2P = Z2P * t2
  ZQP = ZQP**2
  XQP = XQP**2
  ZQP = xPQ * ZQP
    
  return X2P, Z2P, XQP, ZQP


def xDBLe(XP,ZP,A,C,e):
  """ Computes [2^e](X:Z) on Montgomery curve with projective
   constant via e repeated doublings.

  
  Returns:
      projective Montgomery x-coordinates Q <- (2^e)*P.
  """
  A24num = C + C                      
  A24den = A24num + A24num
  A24num = A24num + A
	
  XeP = XP
  ZeP = ZP
	
  for i in range(0,e):
    XeP, ZeP = xDBL(XeP, ZeP, A24num, A24den)
		
  return XeP, ZeP	                         



def xDBLADD(XP,ZP,XQ,ZQ,xPQ,A24):

  """ 
    projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ,
    affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.

  Returns:
      _type_: projective coordinates of x(2P)=X2P/Z2P   and x(Q+P)=XQP/ZQP
  """
  t0 = XP + ZP                  
  t1 = XP - ZP 
  X2P = t0**2
  t2 = XQ - ZQ
  XQP = XQ + ZQ
  t0 = t0 * t2
  Z2P = t1**2
  t1 = t1 * XQP
  t2 = X2P - Z2P
  X2P = X2P * Z2P
  XQP = A24 * t2
  ZQP = t0 - t1
  Z2P = XQP + Z2P
  XQP = t0 + t1
  Z2P = Z2P * t2
  ZQP = ZQP**2
  XQP = XQP**2
  ZQP = xPQ * ZQP
  
  return X2P, Z2P, XQP, ZQP

def xDBLe(XP,ZP,A,C,e):
  """
  Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.

  Returns:
      projective Montgomery x-coordinates Q <- (2^e)*P.
  """

  A24num = C + C                      
  A24den = A24num + A24num
  A24num = A24num + A
	
  XeP = XP
  ZeP = ZP
	
  for i in range(0,e):

    XeP, ZeP = xDBL(XeP, ZeP, A24num, A24den)
		
  return XeP, ZeP	                           
	

def xDBLADD(XP,ZP,XQ,ZQ,xPQ,A24,C24):
  """
   Simultaneous doubling and differential addition.
 
  Returns:
      projective coordinates of x(2P)=X2P/Z2P   and x(Q+P)=XQP/ZQP
       function assumes A24=1, C24=2 fixed
  """
	
  t0 = XP + ZP   
  t1 = XP - ZP
  XP = XQ - ZQ
  ZP = XQ + ZQ
  t2 = (XP * t0)%p
  t3 = (ZP * t1)%p
  ZP = t2 - t3
  XP = t2 + t3
  ZP = (ZP**2)%p
  XQP = (XP**2)%p
  ZQP = (xPQ * ZP)%p

  t1=(t1**2)%p  
  t0=(t0**2)%p  
  t2=t0+t1      
  t2=(t2**2)%p  
  t0=(t0**2)%p  
  t1=(t1**2)%p  
  t2=t2-t1
  X2P=(t2-t0)%p 
  Z2P=(t0-t1)%p 

  return X2P, Z2P, XQP, ZQP     


def xTPL(X, Z, AE, DE):
  """
  triple point in Edwards-coordinates
  projective x-coordinates of xP=X/Z  curve constant A/C

  Returns:
      projective x-coordinates of x(3P)=X3/Z3
  """
	
  Ye = X-Z     
  Ze = Z+X
	
  Ye, Ze = edDBL(Ye, Ze, AE, DE) #Edwards doubling
	
  t0 = Ze + Ze   
  t1 = Ye + Ye
  XP = X - Z
  ZP = X + Z
  t0 = XP * t0
  t1 = ZP * t1
  ZP = t0 - t1
  XP = t0 + t1
  ZP = ZP**2
  X3 = XP**2
  Z3 = X * ZP
  X3 = X3 * Z     

  return X3, Z3     

def xTPLe(X, Z, A, C, e):
  """
   Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants
   A24plus = A+2C and A24minus = A-2C.
  Output: projective Montgomery x-coordinates Q <- (3^e)*
  """
  XeP = X
  ZeP = Z
	
  C2=C+C  #Edwards coefficients
  AE=A+C2
  DE=A-C2
	
  for i in range (0, e):
    XeP, ZeP = xTPL(XeP, ZeP, AE, DE)

  return XeP, ZeP	          

