function [P,S] = permutation(V,X,NFFT,FS,win,OVERLAP,MT,thres)

NF = NFFT/2+1;
D  = ceil(MT*FS/(1000*(NFFT-OVERLAP)));
MA = ones(D,1)/D;
S  = ones(1,NF); 
P  = zeros(1,NF);
TH = min(thres,1);

B1 = specgram(X(:,1),NFFT,FS,win,OVERLAP);
B2 = specgram(X(:,2),NFFT,FS,win,OVERLAP);
L  = size(B1,2);

for f=1:NF,
  
  Z = inv(V(2*f-1:2*f,:))*[B1(f,:);B2(f,:)];

  tmp1 = V(2*f-1:2*f,:)*[Z(1,:);zeros(1,L)]; 
  tmp2 = V(2*f-1:2*f,:)*[zeros(1,L);Z(2,:)];
  
  C11(f,:) = tmp1(1,:); C12(f,:) = tmp1(2,:);
  C21(f,:) = tmp2(1,:); C22(f,:) = tmp2(2,:);
  
  t1 = conv(MA,(ones(1,2)*abs(tmp1)));
  t2 = conv(MA,(ones(1,2)*abs(tmp2)));

  t1=t1/norm(t1); t2=t2/norm(t2);
  
  tmp(f)=t1*t2';
  if tmp(f)>TH,
    S(f)=0;
    C11(f,:) = B1(f,:); C21(f,:) = zeros(1,L); 
    C12(f,:) = B2(f,:); C22(f,:) = zeros(1,L);
  end
end

[tmp,l]=sort(tmp);
f=l(1); k=f;
y1 = zeros(1,L); y2 = y1;

for f=l(2:NF),
  
  t1=conv(MA,(ones(1,2)*abs([C11(f,:);C12(f,:)])));
  t2=conv(MA,(ones(1,2)*abs([C21(f,:);C22(f,:)])));

  if norm(t1)~=0,
    t1=t1/norm(t1);
  end
  if norm(t2)~=0,
    t2=t2/norm(t2);
  end
  
  y1=ones(1,3)*[y1;abs(C11(k,:));abs(C12(k,:))];
  y2=ones(1,3)*[y2;abs(C21(k,:));abs(C22(k,:))];
  
  ny1 = conv(MA,y1);ny1 = ny1/norm(ny1);
  ny2 = conv(MA,y2);ny2 = ny2/norm(ny2);
  
  if (abs(t1*ny1')+abs(t2*ny2') < abs(t1*ny2')+abs(t2*ny1')),
    P(f)=1;
  end
  k=f;
end
