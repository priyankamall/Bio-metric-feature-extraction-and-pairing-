function newim = ridgefilter(im, orient, freq, kx, ky, showfilter)
    if nargin == 5
        showfilter = 0;
    end  
    angleInc = 3;      
    im = double(im);
    [rows, cols] = size(im);
    newim = zeros(rows,cols);
    [validr,validc] = find(freq > 0);  
    ind = sub2ind([rows,cols], validr, validc);
    freq(ind) = round(freq(ind)*100)/100;
    unfreq = unique(freq(ind)); 
    freqindex = ones(100,1);
    for k = 1:length(unfreq)
        freqindex(round(unfreq(k)*100)) = k;
    end
    filter = cell(length(unfreq),180/angleInc);
    sze = zeros(length(unfreq),1);
    for k = 1:length(unfreq)
        sigmax = 1/unfreq(k)*kx;
        sigmay = 1/unfreq(k)*ky;
        sze(k) = round(3*max(sigmax,sigmay));
        [x,y] = meshgrid(-sze(k):sze(k));
        reffilter = exp(-(x.^2/sigmax^2 + y.^2/sigmay^2)/2)...
                .*cos(2*pi*unfreq(k)*x);
        for o = 1:180/angleInc
            filter{k,o} = imrotate(reffilter,-(o*angleInc+90),'bilinear','crop'); 
        end
    end
    maxsze = sze(1);    
    finalind = find(validr>maxsze & validr<rows-maxsze & ...
                    validc>maxsze & validc<cols-maxsze);
    maxorientindex = round(180/angleInc);
    orientindex = round(orient/pi*180/angleInc);
    i = find(orientindex < 1);   orientindex(i) = orientindex(i)+maxorientindex;
    i = find(orientindex > maxorientindex); 
    orientindex(i) = orientindex(i)-maxorientindex; 
    for k = 1:length(finalind)
        r = validr(finalind(k));
        c = validc(finalind(k));
        filterindex = freqindex(round(freq(r,c)*100));
                s = sze(filterindex);   
        newim(r,c) = sum(sum(im(r-s:r+s, c-s:c+s).*filter{filterindex,orientindex(r,c)}));
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [freq, medianfreq] = ridgefreq(im, mask, orient, blksze, windsze, ...
                                        minWaveLength, maxWaveLength)     
    
    [rows, cols] = size(im);
    freq = zeros(size(im));
    for r = 1:blksze:rows-blksze
        for c = 1:blksze:cols-blksze
          blkim = im(r:r+blksze-1, c:c+blksze-1);   
          blkor = orient(r:r+blksze-1, c:c+blksze-1);       
          freq(r:r+blksze-1,c:c+blksze-1) =  ...
              freqest(blkim, blkor, windsze, minWaveLength, maxWaveLength);
        end
    end
    freq = freq.*mask;
    medianfreq = median(freq(find(freq>0)));  

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function [orientim, reliability] = ...
             ridgeorient(im, gradientsigma, blocksigma, orientsmoothsigma)
    [rows,cols] = size(im);
    sze = fix(6*gradientsigma);   if ~mod(sze,2); sze = sze+1; end
    f = fspecial('gaussian', sze, gradientsigma); 
    [fx,fy] = gradient(f);                       
    
    Gx = filter2(fx, im); 
    Gy = filter2(fy, im); 
    Gxx = Gx.^2;       
    Gxy = Gx.*Gy;
    Gyy = Gy.^2;
    sze = fix(6*blocksigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, blocksigma);
    Gxx = filter2(f, Gxx); 
    Gxy = 2*filter2(f, Gxy);
    Gyy = filter2(f, Gyy);
    denom = sqrt(Gxy.^2 + (Gxx - Gyy).^2) + eps;
    sin2theta = Gxy./denom;            
    cos2theta = (Gxx-Gyy)./denom;
       
    sze = fix(6*orientsmoothsigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, orientsmoothsigma);    
    cos2theta = filter2(f, cos2theta); % Smoothed sine and cosine of
    sin2theta = filter2(f, sin2theta); % doubled angles
    
    orientim = pi/2 + atan2(sin2theta,cos2theta)/2;
    Imin = (Gyy+Gxx)/2 - (Gxx-Gyy).*cos2theta/2 - Gxy.*sin2theta/2;
    Imax = Gyy+Gxx - Imin;
    reliability = 1 - Imin./(Imax+.001);
    reliability = reliability.*(denom>.001);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  
  function [normim, mask, maskind] = ridgesegment(im, blksze, thresh)
    im = normalise(im,0,1);  
    fun = inline('std(x(:))*ones(size(x))');
    stddevim = blkproc(im, [blksze blksze], fun);
    mask = stddevim > thresh;
    maskind = find(mask);
    im = im - mean(im(maskind));
    normim = im/std(im(maskind));    

