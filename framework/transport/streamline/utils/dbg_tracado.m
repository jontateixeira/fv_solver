% Scritp para plotar as linhas - on debugging...

global elem coord 

Sl = reshape([S, repmat({[nan, nan]}, [numel(S),1])]',[], 1);
Ss = vertcat(Sl{:});
% plot(Ss(:,1), Ss(:,2), 'r-')

patch('Faces',elem(:,1:4),'Vertices',[coord(:,1),coord(:,2)],'FaceColor',...
'white','EdgeColor',[1,1,1]*0.75);
hold on;
plot(Ss(:,1), Ss(:,2), 'r-');