function R = MoserShuffle(root)

    range = [30/1000 (root.b_ts(end)-root.b_ts(1))-30];
    

    for i = 1:100
        delta = range(1) + (range(2)-range(1))*rand(length(root.cel_x{1}),1);
        ts = mod(root.spike(root.cel(1), root.cel(2)).ts + delta, root.b_ts(end));
        Spk(i,1) = CMBHOME.Spike('ts',ts, 'vid_ts', root.b_ts);
    end
    
    root.spike = Spk;
    root = root.AlignSpike2Session;
    
    for i = 1:length(root.spike)
        root.cel = [i 1];
        R(i) = InstFR(root);
    end
    
end