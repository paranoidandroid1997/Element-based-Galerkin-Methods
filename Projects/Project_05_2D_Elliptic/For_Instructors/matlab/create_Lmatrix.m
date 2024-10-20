%---------------------------------------------------------------------%
%This function computes the 2D Laplacian Matrix on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Lmatrix = create_Lmatrix(intma,jac,wnq,ksi_x,ksi_y,eta_x,...         
                   eta_y,psi,dpsi,iperiodic,npoin,nelem,ngl,nq)

%Initialize
Lmatrix=zeros(npoin,npoin);
inode=zeros(ngl,ngl);

for e=1:nelem
   
   %Store Coordinates
   for j=1:ngl
   for i=1:ngl
      ip=intma(i,j,e);
      inode(i,j)=iperiodic(ip);
   end %i
   end %j
   
   %Do LGL Integration
   for l=1:nq
   for k=1:nq
       wq=wnq(k)*wnq(l)*jac(k,l,e);
       e_x=ksi_x(k,l,e);
       e_y=ksi_y(k,l,e);
       n_x=eta_x(k,l,e);
       n_y=eta_y(k,l,e);
       
       %Loop through I points
       for j=1:ngl
       for i=1:ngl
           ip=inode(i,j);
           h_e=dpsi(i,k)*psi(j,l);
           h_n=psi(i,k)*dpsi(j,l);
           dhdx_i=h_e*e_x + h_n*n_x;
           dhdy_i=h_e*e_y + h_n*n_y;
           
           %Loop through J points
           for n=1:ngl
           for m=1:ngl
               jp=inode(m,n);
               h_e=dpsi(m,k)*psi(n,l);
               h_n=psi(m,k)*dpsi(n,l);
               dhdx_j=h_e*e_x + h_n*n_x;
               dhdy_j=h_e*e_y + h_n*n_y;
               Lmatrix(ip,jp)=Lmatrix(ip,jp) - wq*(dhdx_i*dhdx_j + dhdy_i*dhdy_j);
           end %m
           end %n
       end %i
       end%j
   end %k
   end %l
end %e

%Boundary Conditions
for i=1:npoin
    j=iperiodic(i);
    if (i ~= j) 
        Lmatrix(i,:)=0;
        Lmatrix(i,i)=1;
    end
end


      
