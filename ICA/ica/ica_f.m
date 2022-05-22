function [Y1,Y2] = ica_f(X,NFFT,FS,OVERLAP,N)
%
%	function [Y1,Y2] = ica_f(X,NFFT,FS,OVERLAP,N)
%	split X based on the freequency domain Blind Separation
%	X must has 2 column,
%	
%	Shiro Ikeda 15,July,1998 

MT    =   40; % time for taking the Moving Average [msec]

[LENGTH,dim] = size(X);

if dim==2,
  fprintf(1,'Blind Separation in Time-Frequency domain.\n');

  win = hamming(NFFT);

  fprintf(1,'\nCorrelation Matrices.\n');
  M=correlation(X,NFFT,FS,win,OVERLAP,N);
  fprintf(1,'Done.\n');
  
  fprintf(1,'\nCalculating Decorrelation Matrices.\n');
  V=decorrelation(M,NFFT,N);
  fprintf(1,'Done.\n');
  
  fprintf(1,'\nSolving Permutation.\n');
  [P,S] =  permutation(V,X(1:min(2*FS,LENGTH),:),NFFT,FS,win,OVERLAP,MT,0.95);
  fprintf(1,'Done.\n');
  
  fprintf(1,'\nBuilding Filters.\n');
  [BT1,BT2]=sepfilter(S,P,V,NFFT);
  fprintf(1,'Done.\n');
  
  fprintf(1,'\nMaking Output Signals.\n');
  Y1(:,1)=conv(BT1(1,:),X(:,1))+conv(BT1(2,:),X(:,2));
  Y1(:,2)=conv(BT1(3,:),X(:,1))+conv(BT1(4,:),X(:,2));
  Y2(:,1)=conv(BT2(1,:),X(:,1))+conv(BT2(2,:),X(:,2));
  Y2(:,2)=conv(BT2(3,:),X(:,1))+conv(BT2(4,:),X(:,2));
  fprintf(1,'Done.\n');
else 
  fprintf(1,'This program cannot handle this data.\n');
end
