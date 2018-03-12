function out = logn_pdf(x,mu,sigma)

out = ( x.*sigma*sqrt(2*pi) ).^(-1) .*...
        exp( -(log(x)-mu).^2 ./ (2*sigma^2) );
    
end