function I = ZetaMN(m,n,x,y)
    
    integrand = @(t) (2*log(sqrt(1-x*y*cos(t)) + sqrt(1-y)) - log(y) - log(1-x*cos(t))) ./ ((1-x*cos(t)).^n.*(1-x*y*cos(t)).^(m-0.5));
    I = integral(@(t) integrand(t), 0, 2*pi);
    I = I .* (1-y).^(m-0.5) / (2*pi);
    
end