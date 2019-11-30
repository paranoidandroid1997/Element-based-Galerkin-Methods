%----------------------------------------------------------------------%
%This subroutine constructs the Side Information for a High Order 
%Spectal Element Quads
%Written by Francis X. Giraldo 
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%
% INPUT LIST: iside = face information to know which points are on a face
%                     and which elements they belong to
%             intma = element connectivity
%             nface = number of faces/edges in the grid
%             ngl = number of interpolation points in one direction in an element
%
% OUTPUT LIST: face = face information stating which local edge number the face 
%                      is on and which elements they belong to (left and right neighbors).
%                      are the metric terms needed to
%               imapl and imapr give the tensor-product (i,j) points of the
%               edge points.
%
%----------------------------------------------------------------------%
function [face,imapl,imapr]=create_face(iside,intma,nface,ngl)

%global arrays
face=zeros(nface,4);
imapl=zeros(4,2,ngl);
imapr=zeros(4,2,ngl);

%local arrays
inode=zeros(4,1);
jnode=zeros(4,1);

%Construct Boundary Pointer
inode(1)=1;
inode(2)=ngl;
inode(3)=ngl;
inode(4)=1;
jnode(1)=1;
jnode(2)=1;
jnode(3)=ngl;
jnode(4)=ngl;

%Construct IMAP arrays
for l=1:ngl

   %eta=-1
   imapl(1,1,l)=l;
   imapl(1,2,l)=1;
   imapr(1,1,l)=ngl+1-l;
   imapr(1,2,l)=1;

   %ksi=+1
   imapl(2,1,l)=ngl;
   imapl(2,2,l)=l;
   imapr(2,1,l)=ngl;
   imapr(2,2,l)=ngl+1-l;

   %eta=+1
   imapl(3,1,l)=ngl+1-l;
   imapl(3,2,l)=ngl;
   imapr(3,1,l)=l;
   imapr(3,2,l)=ngl;

   %ksi=-1
   imapl(4,1,l)=1;
   imapl(4,2,l)=ngl+1-l;
   imapr(4,1,l)=1;
   imapr(4,2,l)=l;
end %l

%loop thru the sides
for i=1:nface

   ip1=iside(i,1);
   ip2=iside(i,2);
   iel=iside(i,3);
   ier=iside(i,4);

   %check for position on Left Element
   for j=1:4
      j1=j;
      j2=j+1;
      if (j2 > 4) 
         j2=1;
      end %j2   
      jp1=intma(iel,inode(j1),jnode(j1));
      jp2=intma(iel,inode(j2),jnode(j2));

      if (ip1 == jp1 && ip2 == jp2) 
         face(i,1)=j;
         break; %leave J loop
      end %ip1
   end %j

   %check for position on Right Element
   if (ier > 0)
      for j=1:4
         j1=j;
         j2=j+1;
         if (j2 > 4) 
            j2=1;
         end %j2   
         jp1=intma(ier,inode(j1),jnode(j1));
         jp2=intma(ier,inode(j2),jnode(j2));

         if (ip1 == jp2 && ip2 == jp1) 
            face(i,2)=j;
            break;          %leave J loop
         end %ip1
      end %j
   end %ier
   
   %Store Elements into face
   face(i,3)=iel;
   face(i,4)=ier;
end %i



