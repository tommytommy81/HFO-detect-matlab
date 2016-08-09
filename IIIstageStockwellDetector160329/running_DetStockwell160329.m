clear HFOobj results
HFOobj = input;
HFOobj.EVENTS = []
 
for ch = 1:size(data.x_bip,1)
    ch
 
    disp( ['STAGE I ch: ' num2str(ch)])
    HFOobj.ch = ch;
    % compute stage I and II
    [HFOobj]  = UMC_IIIstage_st1_2(data.x_bip(ch,:), HFOobj);
end
[HFOobj, results] = UMC_IIIstage_st3(data.x_bip, HFOobj);

 