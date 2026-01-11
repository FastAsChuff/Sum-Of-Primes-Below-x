# Sum-Of-Primes-Below-x
Uses my adaptaion of Legendre's prime counting method to calculate the sum of all primes less than or equal to input x &lt;= 10^19.

This method is easy to understand and significantly faster than prime sieving all the way up to x. 

Only prime sieving up to sqrt(x) is required, though there are benefits to going further.

A cache is used for small intermediate values. The cache size can be altered, though the benefits are minor.

See https://math.stackexchange.com/questions/1378286/find-the-sum-of-all-primes-smaller-than-a-big-number/5118835#5118835

For a description of the equations used, see the PNG file included. 

Since c doesn't allow for 'infinite' recursion, I coded a stack which makes it look more complicated than it really is.

The results show an equivalent prime sieving speed in billions of integers per second (Gx/s). It is around 1000 times faster than my best prime sieving program for large x.
