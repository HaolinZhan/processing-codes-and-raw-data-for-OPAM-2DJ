%% Matlab codes used for processing all orthogonal-pattern 2DJ spectra to obtain pure absorptive peaks 
%% Authered by Haolin Zhan, Yuqing Huang, Xiamen univerisity
clc;
clear all;
close all;
%% enter the right storage path of original FID in your computer
FIDpath='C:\Users\86157\Desktop\Raw_data_for_OPAM_2DJ_n_propanol.fid';
%% used for multi-band 2DJ exeriments
% FIDpath2='';
% FIDpath3='';
% FIDpath4='';
%% setting parameters
fn2=4096*4;
fn1=128;
%% Read the some parameters that are required to process the data
Para=getpar(FIDpath,'sw','np','np1','np2','retdly','ni','sw1','at','tauPS','droppts','lsfid');
sw=Para.sw;
np=Para.np;
np1=Para.np1/2;
np2=Para.np2;
retdly=Para.retdly;
ni=Para.ni;
at=Para.at;
sw1=Para.sw1;
isupdatefid=1;
% % parameters for window functions
win=0.6;
win2=2;
 
chunk1=0;
droppts=Para.droppts;
lsfid=Para.lsfid;
 if(droppts>=0)
 Droppts=droppts;
end
 if(lsfid>0 && droppts<0)
   Droppts=lsfid+1; 
end
if(lsfid<=0 && droppts<0)
    Droppts=1;
end
 na=fix(sw/(sw1)+0.5);
na1=na/2;
 if(chunk1>0)
    na1=chunk1;
end
  if(chunk1==0)
    na1=na;
 end
npre=(ni-1)*na+na1;
sw1=1/(2*retdly+2*np1/sw);
sw2=sw;
f1=linspace(-sw1/2,sw1/2,fn1);
f=linspace(-sw/2,sw/2,fn2);
p=linspace(-sw/(2*500),sw/(2*500),fn2);
acount1=0;
for k=1:ni;
    %% Read the FID data
     [R,I]=load_fid(FIDpath(1:end-4),k);
     FIDdata=R+1i*I;
    
    %% data rearrangement shown in Fig.2a and Fig. 2b in the manuscript to construct J couplings
    np1real=round(np1+retdly*sw);
    np2real=min(np2/2,round(np/np1real/2));
    for m=1:4*np2real
        FID2Dodd(:,m)=FIDdata((2*(m-1)*np1real+1):((2*m-1)*np1real));%uniform   
    end
    FID3Dodd(k,:,:)=FID2Dodd;
end
%% data rearrangement shown in Fig. 2b and Fig. 2c to construct pure shifts
FID2dodd=zeros(ni,np1real);
Fid2dresume=zeros(4*np2real,npre);
Fidresume=zeros(1,npre);
for j=1:4*np2real
    FID2dodd=reshape(FID3Dodd(:,:,j),ni,np1real);
    FIDproc=FID2dodd;
    Fidresume(1:na1)=FIDproc(1,(Droppts+1):(na1+Droppts));
   for i=1:(ni-1)
     Fidresume(((i-1)*na+na1+1):(i*na+na1))=FIDproc(i+1,(Droppts+1):(na+Droppts)); 
   end
     Fid2dresume(j,1:npre)=Fidresume;
    signal_1D(1,:,:)=Fid2dresume;
end
 
%% Fourier encoded or not, only available for multi-band experiments
% signal_1D1=fft(signal_1D);
signal_1D1=signal_1D;
 
Fid2dresume1=reshape(signal_1D1(1,:,:),[4*np2real,npre]);

%% add windows
scale=size(Fid2dresume1)
t1=1:scale(2);
t2=1:scale(1);
ht1=exp(0*pi.*(t1-t1(round(length(t1)/2)))/(max(t1)-min(t1))).^2;  
ht2=cos(win*pi.*(t2-t2(round(length(t2)/2)))/(max(t2)-min(t2))).^2;
ht1=ht2.'*ht1;
signal_ifft_window1=Fid2dresume1.*ht1;
    for i=1:scale(1)
 signal_ifft_window1(i,:)=addWindow(signal_ifft_window1(i,:),10000,win2,'LB');
   end
 %% original 2DJ spectrum
signal_window_fft1=fftshift(fft2(signal_ifft_window1,fn1,fn2));
 figure(11); contour(p,f1,real(signal_window_fft1),15);
 figure(12); plot(p,real(sum(signal_window_fft1)));
 figure(13); contour(p,f1,real(signal_window_fft1),5);
 
 %% anti-2DJ spectrum
 Fid2dfftre1=flipud(signal_window_fft1);
 figure(21); contour(p,f1,real(Fid2dfftre1),15);
 figure(22); plot(p,real(sum(Fid2dfftre1)));
 figure(23); contour(p,f1,real(Fid2dfftre1),5);
 
 %% combining both spectra
 Fid2dfftcom=Fid2dfftre1+signal_window_fft1;
 figure(31); contour(p,f1,real(Fid2dfftcom),18);
 figure(32); plot(p,real(sum(Fid2dfftcom)));
 figure(33); contour(p,f1,real(Fid2dfftcom),5);
 
 %% phase corrections along F2 dimension
 aaa3=Fid2dfftcom(fn1/2,:);
 k=size(aaa3,2);
 phcf=[0 0];
    for ikk=1:fn1
       aph1dat=Fid2dfftcom(ikk,:);
       [aph1datph,phc0,phc1]=acme(aph1dat,phcf);
       Fid2dfftcomph_all1(ikk,:)=aph1datph;
       phc0all(ikk)=phc0;
       phc1all(ikk)=phc1;
    end
 figure(41); contour(p,f1,real(Fid2dfftcomph_all1),65);
 figure(42); plot(p,real(sum(Fid2dfftcomph_all1)));
 figure(43); contour(p,f1,real(Fid2dfftcomph_all1),5);
 
%% peak alignment required to multi-band experiments
%  A=circshift(Fid2dfftcomph_all1,[0,-88]);
%  figure(104);
%  plot(real(sum(A)))
%  figure(105);
%  contour(p,f1,real(A),64);
%% in case of N-band experiments, code need to repeat N times with adjusting the corresponding subscripts
