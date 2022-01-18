function odornames = chunk_modifier(chunk, flow_rate, total_flow)
% Flow rate: fraction of input that should be the odor
% Total flow: flow of air + flow of odor

if nargin<3
    total_flow = 0.18;
end

flow_odor = flow_rate*total_flow*ones(length(chunk)/2,1);
flow_odor = [flow_odor; total_flow-flow_odor];
flow_air = total_flow-flow_odor;

odornames = cell(length(chunk),4);
odornames(:,1:2)=chunk(:,1:2);
odornames(:,3) = arr2vec(chunk(:,3),chunk(:,4),1);
odornames(:,4) = arr2vec(flow_odor, flow_air,0); 
end

