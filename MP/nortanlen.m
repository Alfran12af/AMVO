function [nv,tv,pl] = nortanlen(N,x,z)
% Outputs
%   nv: normal vectors at each control point
%   tv: tangent vectors at each control point
%   pl: panel lengths

nv=zeros(N,2); 
tv=zeros(N,2); 
pl=zeros(N,1); 
nx=zeros(N,1); nz=zeros(N,1);
tx=zeros(N,1); tz=zeros(N,1);
 
for i=1:N
    pl(i)=sqrt((x(i+1)-x(i))^2+(z(i+1)-z(i))^2); % panel length
    nx(i)=(z(i)-z(i+1))/pl(i); % normal vector x components (sin alpha)
    nz(i)=(x(i+1) - x(i))/pl(i); % normal vector z components (cos alpha)
    nv(i,:)=[nx(i), nz(i)];
    
    tx(i)= nz(i); % tangent vector x components (cos alpha)
    tz(i)=-nx(i); % tangent vector z components (-sin alpha)
    tv(i,:)=[tx(i), tz(i)];

end
end