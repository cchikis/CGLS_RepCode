% This function calculates the expected number of products by age


function lnn = lnNbyAge(Age, x, tau)

gamma = GammaFunction(Age,x,tau);
lnn = (1-gamma)/gamma * sum(log(1:200).*gamma.^(1:200));

end