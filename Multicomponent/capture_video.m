function capture_video(number_of_timesteps,nodes,elements,phi,tau)

video_filename = 'test';
%% open figure and video file (define framerate, qualty, etc.)
fig = figure('position',[50 100 1200 680]);
n_now=now; % to prevent overwriting existing files
vidObj = VideoWriter([video_filename,'_',num2str(n_now)]);
vidObj.Quality = 100;
vidObj.FrameRate = 24;
open(vidObj);

%% drawing plots and capturing their frames

for j = 1:number_of_timesteps
    
    trisurf(elements, nodes(:,1), nodes(:,2),0*nodes(:,2),phi(:,j));
        view([0 0 1]);
        shading interp
        title(['t = ',num2str(tau*j)])
        axis('equal');
        colorbar
        colormap("jet")
        pbaspect([1,1,1])
    % drawing and capturing the current frame
    drawnow
    writeVideo(vidObj, getframe(fig));
end

%% closing files
close(fig); % in some instances it is worthwhile to close the figure after each frame
close(vidObj);