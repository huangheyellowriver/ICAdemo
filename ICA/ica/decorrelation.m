function [V] = decorrelation(M,NFFT,N)

for f=1:NFFT/2+1,
  SP=inv(sqrtm(M(2*f-1:2*f,1:2)));
  for tau=1:N,
    M(2*f-1:2*f,2*tau+1:2*tau+2)=SP*M(2*f-1:2*f,2*tau+1:2*tau+2)*SP';
  end
  V(2*f-1:2*f,:)=joint_diag(M(2*f-1:2*f,3:2*N+2),0.000001);
  V(2*f-1:2*f,:)=inv(SP)*V(2*f-1:2*f,:);
  % V(2*f-1:2*f,:)=inv(V(2*f-1:2*f,:))*SP;
end

