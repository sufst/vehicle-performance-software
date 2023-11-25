[dmesh,lmesh] = meshgrid(CdA_vals, ClA_vals);

%% Plot 1: Lift vs Drag Lap Time Sensitivity
effmesh = lmesh./dmesh; % efficiency

figure(4)

ax(1) = axes;
contourf(lmesh,dmesh,lap_times,15,'linecolor','none','ShowText','off')

c = colorbar;
colormap(flipud(parula));
c.Label.String = 'Lap Time (s)';
xlabel('Lift Coefficient')
ylabel('Drag Coefficient')
title('iso-Efficiency lines for Drag vs Lift')
grid on
grid minor

ax(2) = axes;
contour(lmesh,dmesh,effmesh,0.5:0.5:6.0,'k-','ShowText','on')
ax(2).Color = 'none';
ax(2).Position = ax(1).Position;
ax(2).XAxis.Visible = 'off';
ax(2).YAxis.Visible = 'off';



%% Plot 2: Lift vs Drag Points Sensitivity
points = nan(size(lap_times));

for i = 1:size(lap_times,1)
    for j = 1:size(lap_times,2)
        points(i,j) = calcSprintScore(lap_times(i,j));
    end
end

figure(5)

ax(1) = axes;
contourf(lmesh,dmesh,points,15,'linecolor','none','ShowText','off')
c = colorbar;
colormap(parula);
c.Label.String = 'Sprint Points';
xlabel('Lift Coefficient')
ylabel('Drag Coefficient')
title('iso-Points lines for Drag vs Lift')
grid on
grid minor

ax(2) = axes;
contour(lmesh,dmesh,points,80:0.5:100,'k-','ShowText','on')
ax(2).Color = 'none';
ax(2).Position = ax(1).Position;
ax(2).XAxis.Visible = 'off';
ax(2).YAxis.Visible = 'off';


%% Plot 3: Lift vs Drag Energy Consumption Sensitivity over Endurance
figure(6)

ax(1) = axes;
contourf(lmesh,dmesh,energy,15,'linecolor','none','ShowText','off')
c = colorbar;
colormap(flipud(parula));
c.Label.String = 'Energy Consumption per Lap (MJ)';
xlabel('Lift Coefficient')
ylabel('Drag Coefficient')
title('iso-MJ lines for Drag vs Lift')
grid on
grid minor

ax(2) = axes;
contour(lmesh,dmesh,energy,13:0.1:15.6,'k-','ShowText','on')
ax(2).Color = 'none';
ax(2).Position = ax(1).Position;
ax(2).XAxis.Visible = 'off';
ax(2).YAxis.Visible = 'off';