function [output_arr] = vcat(input_array,row_sw)
% Concatenates vectors in a cell array with unequal vector sizes. If row_sw
% = true (default), concatenation is done across rows. Else, columns.
% Vectors in cells of input_array are column

if nargin<2
    row_sw = true;
end

% Row_concatenation
if row_sw
    output_arr = cell(1,size(input_array,2));
    for ii=1:size(input_array,2) % column
        temp = [];
        for jj=1:size(input_array,1) % row
            temp = cat(1,temp,input_array{jj,ii});
        end
        output_arr{ii} = temp; 
    end
end
    
% Column concatenation
if ~row_sw
    output_arr = cell(size(input_array,1),1);
    for ii=1:size(input_array,1) % row
        temp = [];
        for jj=1:size(input_array,2) % column
            temp = cat(1,temp,input_array{ii,jj});
        end
        output_arr{ii} = temp; 
    end
end