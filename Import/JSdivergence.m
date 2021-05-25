function out=JSdivergence(X,Y)

KL=@(P,Q) sum(P.*log(P./Q));

M=(X+Y)/2;
out=(KL(X,M)+KL(Y,M))/2;

end
