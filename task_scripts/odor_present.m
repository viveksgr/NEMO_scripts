function odor_present(odorname,options)
% Odorname = 2X4 cell. Odorname(1,:) = cid, name, valve_id, valve_flow for
% odor. Odorname(2,:) for air.
% Countdown  = 3 seconds before odor perception. Cue present for 1 second.
clearpict(1)
clearpict(2)
wait(2000)

if options.run_odors
    daq = OlfConfigDaq;
    OlfFlowControl(daq, odorname{1,4}(1), odorname{1,4}(2));
    OlfOpenLines(odorname{1,3}, daq);
end
wait(options.delay-3000);

setforecolour(1,1,1);
settextstyle('Arial',40);
preparestring('3',2,0,0);
drawpict(2);
wait(1000);

% Countdown = 2 seconds before odor perception. Present for 1 second.
clearpict(2)
setforecolour(1,1,1);
settextstyle('Arial',40);
preparestring('2',2,0,0);
drawpict(2);
wait(500)
% wait((2000+options.sniffDur)-options.delay);

% Set odor to clean
if options.run_odors
    OlfFlowControl(daq,odorname{end,4}(1),odorname{end,4}(2));
    OlfOpenLines(odorname{end,3},daq,16);
end
% wait(options.delay-(options.sniffDur+1))
wait(500)

% Countdown = 1 seconds.
clearpict(2)
setforecolour(1,1,1);
settextstyle('Arial',40);
preparestring('1',2,0,0);
drawpict(2);
wait(1000);

% Present sniff cue: "+" changes color
clearpict(2)
setforecolour(1,0,0);
settextstyle('Arial',40);
preparestring('+',2,0,0);
drawpict(2);
wait(options.sniffDur)


% Remove fixation
clearpict(2)
end

