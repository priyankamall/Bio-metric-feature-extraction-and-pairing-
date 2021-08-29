function n = normalise(im, reqmean, reqvar)
    if ~(nargin == 1 | nargin == 3)
       error('No of arguments must be 1 or 3');
    end    
    if nargin == 1   
	if ndims(im) == 3        
	    hsv = rgb2hsv(im);
	    v = hsv(:,:,3);
	    v = v - min(v(:));    
	    v = v/max(v(:));
	    hsv(:,:,3) = v;
	    n = hsv2rgb(hsv);
	else                     
	    if ~isa(im,'double'), im = double(im); end
	    n = im - min(im(:));
	    n = n/max(n(:));
	end
    else  	
	if ndims(im) == 3   
	    error('cannot normalise colour image to desired mean and variance');
	end
	if ~isa(im,'double'), im = double(im); end	
	im = im - mean(im(:));    
	im = im/std(im(:));   
	n = reqmean + im*sqrt(reqvar);
    end  
