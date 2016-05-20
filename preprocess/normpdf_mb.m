function y = normpdf_mb(x,mean,stddev)
    
    y= (1/(stddev*sqrt(2*pi))*exp(-(x-mean).^2/(2*stddev^2)));