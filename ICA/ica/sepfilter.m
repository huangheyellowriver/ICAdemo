function [BT1,BT2]=sepfilter(S,P,V,NFFT)

NF = NFFT/2+1;

for f=1:NF,
  f1=2*f-1; f2=2*f;
  if S(f)==1,
    TMPV=inv(V(f1:f2,:));
    if P(f)==1,
      TMPV=flipud(TMPV);
    end
    TMPINV=inv(TMPV);
    BF1=TMPINV(:,1)*TMPV(1,:); BF2=TMPINV(:,2)*TMPV(2,:);
  else
    if P(f)==0,
      BF1=eye(2); BF2=zeros(2);
    else
      BF1=zeros(2); BF2=eye(2);
    end
  end

  BF1_11(f)=BF1(1,1); BF1_12(f)=BF1(1,2);
  BF1_21(f)=BF1(2,1); BF1_22(f)=BF1(2,2);

  BF2_11(f)=BF2(1,1); BF2_12(f)=BF2(1,2);
  BF2_21(f)=BF2(2,1); BF2_22(f)=BF2(2,2);
end

BF1_11(NF+1:NFFT)=conj(fliplr(BF1_11(2:NF-1)));
BF1_12(NF+1:NFFT)=conj(fliplr(BF1_12(2:NF-1)));
BF1_21(NF+1:NFFT)=conj(fliplr(BF1_21(2:NF-1)));
BF1_22(NF+1:NFFT)=conj(fliplr(BF1_22(2:NF-1)));
  
BF2_11(NF+1:NFFT)=conj(fliplr(BF2_11(2:NF-1)));
BF2_12(NF+1:NFFT)=conj(fliplr(BF2_12(2:NF-1)));
BF2_21(NF+1:NFFT)=conj(fliplr(BF2_21(2:NF-1)));
BF2_22(NF+1:NFFT)=conj(fliplr(BF2_22(2:NF-1)));
  
BT1(1,:)=ifft(BF1_11); BT1(2,:)=ifft(BF1_12);
BT1(3,:)=ifft(BF1_21); BT1(4,:)=ifft(BF1_22);
BT2(1,:)=ifft(BF2_11); BT2(2,:)=ifft(BF2_12);
BT2(3,:)=ifft(BF2_21); BT2(4,:)=ifft(BF2_22);
  
BT1=real(BT1); BT2=real(BT2);
