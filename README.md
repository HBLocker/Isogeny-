# Isogeny-
Trying to learn about isogeny  crypto


## The Basics 

##### Warning!
This is me trying to exaplin something quite rather complicated in a minimal fashion so it may be lacking to some crpyto people. 


## What even is an Isogony?

An isogeny is essently a special type of morphism between eleptic curves. The main idea about this is that we want to be able to map one curve from another curve.

![Example](https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSoCxuebe4c5Ci6mULqAMcMxqgk0GlR2KQXFaENeVaACuoHZEWZaPFI3Ykjzj1BclDMyxI&usqp=CAU)

With the use of isogony based cryptography there is no single curve, but many curves, adding complexity and robustness to the use of the algoithm. 

- ALice calcuates her kenel of 
```
R=mPB+nQB
```
 - Next she computes her isogony with [Velu formulars](https://eprint.iacr.org/2011/430.pdf) as:
```
 ϕa
```
- She uses ϕa to start hers random walk and ends up with trio:
```
(Ea,ϕa(PB),ϕa(PB))
```

Now she sends this trio to Bob (I know, she actually sends 3 points)
Bob does the same, let say his secret isogeny is ϕb. 

Upon receipt of corresponding trio from Bob 

```
(Eb,ϕb(PA),ϕb(QA))
```

 - Alice uses m,n to compute new kernel R′=m∗ϕb(PA)+n∗ϕb(QA) and new isogeny
Then she starts from Eb and does the random walk again ending on:
```
Eab
```
- Bob proceeds mutatis mutandis and ends up on a probably different curve Eba
 Then the j-invariant then can agree on a shared secret  j(Eab)==j(Eba)








### Refferences 


sagemath. (NA). Isogeny computation does not finish in Sage. [Online]. sagemath. Available at: https://ask.sagemath.org/question/40675/isogeny-computation-does-not-finish-in-sage/ [Accessed 20 November 2022].

Tanja Lange. (2021). Isogeny-based crypto. [Online]. Youtube. Available at: https://www.youtube.com/watch?v=ki_9p6CAdRc&list=PL6hzlGxGIS1Cnx3XS7ZD4wjcTmHqOEpTS [Accessed 20 November 2022].


ANALOGUES OF VELU’S FORMULAS FOR ISOGENIES ON ´ ALTERNATE MODELS OF ELLIPTIC CUR. (NA). DUSTIN MOODY AND DANIEL SHUMOW. NA. NA(NA), p.NA. [Online]. Available at: https://eprint.iacr.org/2011/430.pdf [Accessed 20 November 2022].


Will update this more 



