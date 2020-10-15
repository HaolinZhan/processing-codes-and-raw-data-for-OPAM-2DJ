function sp_dp = dephase_fun(sp_ft,phc0,phc1)
% function sp_dp = dephase_fun(sp_ft,phc0,phc1)
% written by Chen Li
% input:   sp_ft - complex data after fourier transform
%         phc0  - the zero order phase correction
%         phc1  - the first order phase correction
% output:  sp_dp - spectral data after phase correction
%
% global k
k=size(sp_ft,2);
x1=0:1:k-1;
th=(phc0+x1.*phc1/k)*pi/180;
spr=real(sp_ft).*cos(th)-imag(sp_ft).*sin(th);
spi=real(sp_ft).*sin(th)+imag(sp_ft).*cos(th);
sp_dp=spr+j*spi;

