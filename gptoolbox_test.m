t = tsurf(F,V);
set(gcf,'COlor',0.94*[1 1 1]);
teal = [144 216 196]/255;
pink = [254 194 194]/255;
bg_color = pink;
fg_color = teal;
for pass = 1:10
  switch pass
  case 1
    % blank run
    axis([-209.4       119.38      -181.24       262.67      -247.28 247.38]);
  case 2
    axis equal;
    axis([-209.4       119.38      -181.24       262.67      -247.28 247.38]);
    axis vis3d;
  case 3
    t.EdgeColor = 'none';
  case 4
    set(t,fphong,'FaceVertexCData',repmat(fg_color,size(V,1),1));
  case 5
    set(t,fsoft);
  case 6
    l = light('Position',[0.2 -0.2 1]);
  case 7
    set(gca,'Visible','off');
  case 8
    set(gcf,'Color',bg_color);
  case 9
    s = add_shadow(t,l,'Color',bg_color*0.8,'BackgroundColor',bg_color,'Fade','infinite');
  case 10
    apply_ambient_occlusion(t,'AddLights',false,'SoftLighting',false);
  end

  vidObj = VideoWriter(sprintf('nefertiti-%02d.mp4',pass),'MPEG-4');
  vidObj.Quality = 100;
  vidObj.open;
  thetas = linspace(30,-30,450);
  for theta = thetas(1:end-1)
    view(theta,30);
    drawnow;
    vidObj.writeVideo(getframe(gcf));
  end
  vidObj.close;

end