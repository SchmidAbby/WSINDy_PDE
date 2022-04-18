function displacement_data = integrate_data(velocity_data, time)
    
    [x_pts, t_pts] = size(velocity_data);

    % initialize matrix for displacement data
    displacement_data = zeros(x_pts, t_pts);

    for i=1:x_pts
        displacement_data(i, :) = cumtrapz(velocity_data(i,:), time);
    end    

%     % movie 
%     vidfile = VideoWriter('2022-04-14-movie-time1to4540.mp4','MPEG-4');
%     open(vidfile);
%     for i=1:4540 
%         figure(1)
%         %p1 = plot(1:x_pts, velocity_data(:,i));
%         %hold on
%         p2 = plot(1:x_pts, displacement_data(:,i));
%         ylim([-10 10])
%         %ylabel('Velocity (mm/s)')
%         %legend('Velocity', 'Displacement')
%     
%         hold off
%         set(gca,'XColor', 'none')
%         frame = getframe(gcf);
%         writeVideo(vidfile, frame);
%     end
%     close(vidfile);

end