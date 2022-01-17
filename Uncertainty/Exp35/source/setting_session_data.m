sessions = struct();

sessions.initial_cl_bar = 3;
sessions.initial_cl_wind = 1;
sessions.cue_combination= 4;
sessions.final_cl_bar = 6;
sessions.final_cl_wind = 5;
sessions.grating = [9,10];
sessions.starfield = [7,8];

save('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\low_reliability\20220114_60D05_7f_fly3\analysis\sessions_info.mat','sessions')