function [m_out] = sweeper_transform(m_in)
% For each column in m_in, count number of entries outside a linearly
% decreasing threshold in m_out. Also, provide slopes from sigmoid fits.

% thresh = linspace(0.5,-0.5,100);
thresh = linspace(0.5,-0.5,100);
F = @(b,x)(1/(1+exp(-b(1)*(x-b(2)))));
x = (1:1:length(thresh))';
m_out = zeros(length(thresh),size(m_in,2));
for zz = 1:size(m_in,2)
    temp = m_in(:,zz);
    for yy = 1:length(thresh)
        m_out(yy,zz)=sum((temp>thresh(yy)))/length(temp);
%         [ahat] = lsqcurvefit(F,[1 1],x,m_out(:,zz)');
%         ms(zz) = ahat(1);
    end
end


