function unmasked = unmasker(masked,mask)
% Finds unmasked vector such as unmasked(mask)= masked. Pad zeros wherever
% needed.
unmasked = zeros(size(mask));
mask_id = find(mask);
for ii = 1:length(masked)
    unmasked(mask_id(ii))=masked(ii);
end
