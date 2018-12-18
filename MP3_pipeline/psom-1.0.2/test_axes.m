cm = [0 1 1 0 0;1 0 0 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 0 0];
cm2 = [0 1 1 0 1;1 0 1 1 1;1 0 0 0 0;0 0 0 0 1;1 0 1 1 0];

h(1) = figure;
view(biograph(cm));
fig(1)=gcf;       ax1=gca;   % Figure and axes handles
c1 = get(fig(1), 'Children');
copyobj(c1,h(1));            % Copying the object to h(1)

h(2) = figure;
view(biograph(cm2));
fig(2)=gcf;       ax2=gca;   % Figure and axes handles
c2 = get(fig(2), 'Children'); 
copyobj(c2,h(2));            % Copying the object to h(2)

figure;
% Creating and getting handles of subplot axes
% Retrieving properties of children and copying to new parents and hiding the axes lines
s1 = subplot(1,2,1);    bg1 = get(ax1,'children');  copyobj(bg1,s1);    axis off;
s2 = subplot(1,2,2);    bg2 = get(ax2,'children');  copyobj(bg2,s2);    axis off;

close(fig);    close(h);    % Closing previously opened figures