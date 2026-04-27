function div = diverg(u, v, h, N)

div = zeros(N+2, N+2);

for i=2:N+1
    for j=2:N+1
        up = u(i,j);
        vp = v(i,j);
        uw = u(i-1,j);
        vs = v(i,j-1);
        div = (up + vp - uw - vs)*h;
    end
end

div = halo_update(div);

end