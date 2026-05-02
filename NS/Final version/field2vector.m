function b = field2vector(div, N)

% b = zeros(N, N);
b = div(2:end-1, 2:end-1); %Treu el Halo de la matriu div
b = reshape(b, [N*N, 1]); %Converteix la Matriu resultant sense Halo en Vector.

end