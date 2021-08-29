function nfimg = smoothen_frequency_image(fimg,RLOW,RHIGH,diff_cycles)
    valid_nbrs  =   3; 
    [ht,wt]     =   size(fimg);
    nfimg       =   fimg;
    N           =   1;
    h           =   fspecial('gaussian',2*N+1);
    cycles      =   0;
    invalid_cnt = sum(sum(fimg<RLOW | fimg>RHIGH));
    while((invalid_cnt>0 &cycles < diff_cycles) | cycles < diff_cycles)
        fimg    =   [flipud(fimg(1:N,:));fimg;flipud(fimg(ht-N+1:ht,:))]; 
        fimg    =   [fliplr(fimg(:,1:N)),fimg,fliplr(fimg(:,wt-N+1:wt))]; 
        for i=N+1:ht+N
         for j = N+1:wt+N
                blk = fimg(i-N:i+N,j-N:j+N);
                msk = (blk>=RLOW & blk<=RHIGH);
                if(sum(sum(msk))>=valid_nbrs)
                    blk           =blk.*msk;
                    nfimg(i-N,j-N)=sum(sum(blk.*h))/sum(sum(h.*msk));
                else
                    nfimg(i-N,j-N)=-1;
                end;
         end;
        end;
        fimg        =   nfimg;
        invalid_cnt =   sum(sum(fimg<RLOW | fimg>RHIGH));
        cycles      =   cycles+1;
    end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  function noimg = smoothen_orientation_image(oimg)
    gx      =   cos(2*oimg);
    gy      =   sin(2*oimg);
    msk     =   fspecial('gaussian',5);
    gfx     =   imfilter(gx,msk,'symmetric','same');
    gfy     =   imfilter(gy,msk,'symmetric','same');
    noimg   =   atan2(gfy,gfx);
    noimg(noimg<0) = noimg(noimg<0)+2*pi;
    noimg   =   0.5*noimg;
