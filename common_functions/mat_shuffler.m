function mat = mat_shuffler(mat,ord)
% Shuffles the entries of a matrix sequentially
if isfloat(mat)
    for ii = 1:size(ord,1)
        m1 = mat(ord(ii,1));
        m2 = mat(ord(ii,2));
        mat(ord(ii,2))= m1;
        mat(ord(ii,1))= m2;
    end
elseif iscell(mat)
    for ii = 1:size(ord,1)
        m1 = mat{ord(ii,1)};
        m2 = mat{ord(ii,2)};
        mat{ord(ii,2)}= m1;
        mat{ord(ii,1)}= m2;
    end
end