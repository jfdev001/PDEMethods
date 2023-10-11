"""
    legendre(x, n)
  
Compute Legendre polynomial of order `n` at point `x` using explicit
representation.

# References
[1] : https://en.wikipedia.org/wiki/Legendre_polynomials
"""
function legendre(x, n)
    @assert n >= 0
    polysum = 0
    for k = 0:n
        intermed = binomial(n, k)*binomial(n+k, k)*((x - 1)/2)^k
        polysum += intermed
    end 
    return polysum 
end

"""
    hermite(x, n)

Compute physicist's Hermite polynomial of order `n` at point `x` using 
explicit expression.

# References
[1] : https://en.wikipedia.org/wiki/Hermite_polynomials
"""
function hermite(x, n)
    polysum = 0
    upper_bound = Int(floor(n/2))
    for m = 0:upper_bound
        polysum += ((2*x)^(n-2*m)*(-1)^m)/(factorial(m)*factorial(n - 2*m))
    end 
    return factorial(n)*polysum
end 
