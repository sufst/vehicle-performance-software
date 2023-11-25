%% Set up the Lap Sim
veh = load('01 Vehicles/StagX.mat');
veh.mass = 270;
track = load('00 Tracks/FSUK_Sprint_MP.mat');
options.BOutputLog = 0;

ClA_vals = 0.0:0.1:5.0;
CdA_vals = 0.25:0.1:2.0;

sim_names = cell(length(ClA_vals),length(CdA_vals));
lap_times_skid = nan(size(sim_names));
energy = lap_times_skid;
vMax = energy;
sim_res = sim_names;
cnt = 0;

%% Main loop
fprintf('----- Running %d Sims -----\n', numel(lap_times_skid))
for cd = 1:length(CdA_vals)
    veh.CdA = CdA_vals(cd);
    
    for cl = 1:length(ClA_vals)
        cnt = cnt+1;
        fprintf('Sim %02d... ',cnt);
        veh.ClA = ClA_vals(cl);
        
        sim_names{cl,cd} = sprintf('%dNm_%dkg',veh.CdA,veh.ClA);
        
        opts.PlotGGV = 0;
        opts.BOutputLog = 0;
        try
            res = runSkidpad(veh);
            lap_times_skid(cl,cd) = res.time;
            fprintf('LapTime : %05.2f \n',lap_times_skid(cl,cd));
        catch err
            keyboard
            fprintf('FAILED\n')
        end
    end
end