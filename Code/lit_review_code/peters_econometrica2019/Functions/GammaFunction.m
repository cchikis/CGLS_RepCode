


function g = GammaFunction(Age,x,tau)

g = x * (1 - exp( - (tau - x)*Age))./(tau - x * exp(-(tau - x)*Age));

end