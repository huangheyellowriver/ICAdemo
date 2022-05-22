function [M] = correlation(X,NFFT,FS,win,OVERLAP,N)
%
%	function [M] = correlation(X,FS,NFFT,OVERLAP,N)
%	making correlation matrices for each frequency
%	X must has 2 column,
%	
%	Shiro Ikeda 14,July,1998 

SPLEN = 1000;

[LENGTH,dim] = size(X);

FL = floor((LENGTH-NFFT)/(NFFT-OVERLAP))+1;
DL = (SPLEN-1)*(NFFT-OVERLAP)+NFFT;
NF = NFFT/2+1;

repeat=floor(FL/SPLEN);

MARGIN=zeros(2*NF,N);
te=OVERLAP; M=zeros(2*NF,2*N+2); AVE=zeros(2*NF,1);
for ITER=1:repeat,
  ts = te-OVERLAP+1; te = ts-1+DL;
  B1=specgram(X(ts:te,1),NFFT,FS,win,OVERLAP);
  B2=specgram(X(ts:te,2),NFFT,FS,win,OVERLAP);
  if ITER==1,
    STKH1=B1(:,1:N); STKH2=B2(:,1:N);
  end
  for f=1:NF,
    f1=2*f-1; f2=2*f; 
    AVE(f1) = AVE(f1)+B1(f,:)*ones(SPLEN,1);
    AVE(f2) = AVE(f2)+B2(f,:)*ones(SPLEN,1);
    TMP=[MARGIN(f1,:),B1(f,:);MARGIN(f2,:),B2(f,:)];
    MARGIN(f1:f2,:)=[B1(f,SPLEN-N+1:SPLEN);B2(f,SPLEN-N+1:SPLEN)];
    for tau=0:N,
      TMP1=TMP(:,N-tau+1:SPLEN+N-tau)*TMP(:,N+1:SPLEN+N)';
      M(f1:f2,2*tau+1:2*tau+2)=M(f1:f2,2*tau+1:2*tau+2)+TMP1;
    end
  end
end

ts = te-OVERLAP+1; te = LENGTH;
B1=specgram(X(ts:te,1),NFFT,FS,win,OVERLAP);
B2=specgram(X(ts:te,2),NFFT,FS,win,OVERLAP);
RL=size(B1,2);
if repeat==0,
  STKH1=B1(:,1:N); STKH2=B2(:,1:N);
end
STKT1=B1(:,RL-N+1:RL); STKT2=B2(:,RL-N+1:RL);
for f=1:NF,
  f1=2*f-1; f2=2*f; 
  AVE(f1) = AVE(f1)+B1(f,:)*ones(RL,1);
  AVE(f2) = AVE(f2)+B2(f,:)*ones(RL,1);
  TMP=[MARGIN(f1,:),B1(f,:);MARGIN(f2,:),B2(f,:)];
  for tau=0:N,
    TMP1=TMP(:,N-tau+1:RL+N-tau)*TMP(:,N+1:RL+N)';
    M(f1:f2,2*tau+1:2*tau+2)=M(f1:f2,2*tau+1:2*tau+2)+TMP1;
  end
end

AVE=AVE/FL;

for f=1:NF,
  f1=2*f-1; f2=2*f; 
  TMPA = AVE(f1:f2)*AVE(f1:f2)';
  M(f1:f2,1:2)=M(f1:f2,1:2)-TMPA*FL;
  M(f1:f2,1:2)=M(f1:f2,1:2)+M(f1:f2,1:2)';
  for tau=1:N,
    t1=2*tau+1; t2=2*tau+2;
    T1 = AVE(f1:f2)*FL-[STKT1(f,N-tau+1:N);STKT2(f,N-tau+1:N)]*ones(tau,1);
    T2 = AVE(f1:f2)*FL-[STKH1(f,1:tau);STKH2(f,1:tau)]*ones(tau,1);
    TMP1=T1*AVE(f1:f2)'; TMP2=AVE(f1:f2)*T2';
    M(f1:f2,t1:t2)=M(f1:f2,t1:t2)-TMP1-TMP2+TMPA*(FL-tau);
    M(f1:f2,t1:t2)=M(f1:f2,t1:t2)+M(f1:f2,t1:t2)';
  end
end

for tau=0:N,
  M(:,2*tau+1:2*tau+2)=M(:,2*tau+1:2*tau+2)/(2*(FL-tau));
end
