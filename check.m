% Checking sound 
% Semin Cheon. M.D/Ph.D integrated course. MDAIL(KAIST). chsm0338@kaist.ac.kr
% 2021.11.22 (Y.M.D)
%
% !!  Sound Toolbox must be installed to use this function !!
%
% Input: trigger for check sound
% 
% Output: PPIPIK!
%%
function check(turnon)
WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);

if turnon == 1
pause(1)
playblocking(Audio);
elseif turnon == 0

else 
    
while true
    turnon = input('input = '); 
    if turnon == 1  
       break
    elseif turnon == 0
        break
    else
    end
end

if turnon == 0
    return
else
pause(1)
playblocking 
end

end

end