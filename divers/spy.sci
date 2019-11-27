
function spy(A)

    [i,j] = find(A~=0)
    [N,M] = size(A)

    xsetech([0,0,1,1],[1,0,M+1,N])
    xrects([j;N-i+1;ones(i);ones(i)],ones(i));
    xrect(1,N,M,N);

endfunction
