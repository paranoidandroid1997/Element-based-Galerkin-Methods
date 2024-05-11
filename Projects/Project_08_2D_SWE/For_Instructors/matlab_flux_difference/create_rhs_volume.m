%----------------------------------------------------------------------%
%This subroutine builds the volume integral contribution using the Weak Form DGM-SEM
%on Quadrilateral Elements.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs_volume(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,flux_method)

%global arrays
rhs=zeros(npoin,1);
flux=zeros(ngl,ngl,2);
flux_es=zeros(ngl,ngl,ngl,ngl,2);
node=zeros(ngl,ngl);

%Construct Volume Integral Contribution
for e=1:nelem

    %Loop through element DOF and store Node Pointer and Fluxes
    I=0;
    for j=1:ngl
        for i=1:ngl
            I=I+1;
            I=intma(i,j,e);
            node(i,j)=I;
            flux(i,j,1)=q(I)*u(I);
            flux(i,j,2)=q(I)*v(I);
        end
    end

    %Compute Flux Matrix
    for j=1:ngl
        for i=1:ngl
            for l=1:ngl
                for k=1:ngl
                    flux_es(i,j,k,l,1)=flux(k,l,1);
                    flux_es(i,j,k,l,2)=flux(k,l,2);
                end
            end
        end
    end

    %Loop Integration Points
    for j=1:ngl
        for i=1:ngl
            
            I=iperiodic(node(i,j));
            wq=wnq(i)*wnq(j)*jac(i,j,e);
            e_x=ksi_x(i,j,e);
            e_y=ksi_y(i,j,e);
            n_x=eta_x(i,j,e);
            n_y=eta_y(i,j,e);

            for k=1:ngl
                %Ksi Terms
                h_e=dpsi(k,i)*psi(j,j);
                f_e=flux(k,j,1);
                g_e=flux(k,j,2);
                flux_e=h_e*( f_e*e_x + g_e*e_y);
                %Eta Terms
                h_n=psi(i,i)*dpsi(k,j);
                f_n=flux(i,k,1);
                g_n=flux(i,k,2);
                flux_n=h_n*( f_n*n_x + g_n*n_y);
                %RHS Contribution
                rhs(I)=rhs(I) - wq*(flux_e + flux_n);                
            end %k
        end %i
    end %j
end %e