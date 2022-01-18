function create_linear_plt(y,X,ind,xlabel_,ylabel_,fname)
figure()
plot(X,y,ind)
hold on
m = [ones(size(X)) X]\y;
plot(X,m(1)+m(2).*X)
ylabel(ylabel_)
xlabel(xlabel_)
p = r2p(fastcorr(X,y),length(y));
title(sprintf('m: %.02f, r: %.02f, p: %.05f',m(2),fastcorr(X,y),p))
savefig(fname)