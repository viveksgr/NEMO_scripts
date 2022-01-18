function reshaped = nii_reshaper(vector_,voxel_id,sizer_)
% Reshapes NX1 vector_ into a 3D array specified by sizer_ and indices
% specified by voxel_id. 
reshaped = zeros(sizer_);
for ii = 1:length(vector_)
    reshaped(voxel_id(1,ii),voxel_id(2,ii),voxel_id(3,ii))=vector_(ii);
end
end


