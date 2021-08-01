


%% Permeability Field Selection
% 1 - Monofasico do Mestrado do Fernando - Ativar o termo Fonte no Assembly
% 2 - Bifasico Altamente Heterogeneo e Isotropico
% 3 - Bifasico Aleatorio Chueh
permFlag  = 5;

if permFlag == 1
    %global fontEach
    %%gerando termo fonte
    %eleTag = (elem(:,end) == 2);
    %fontEach = zeros(size(elem,1),1);
    %fontEach(eleTag) = elemarea(eleTag)*(1/0.0625);
    %% gerando funcao de campo de permeabilidade
    ek = 5 * 10^-2;
    permFF = @(x,y)[(y.^2 + ek.*(x.^2)) , -(1-ek).*x.*y, -(1-ek).*x.*y, (x.^2 + ek.*(y.^2))];
    %% alterando campos de permeabilidade
    kmap = [ [1:size(elem,1)]' permFF(centelem(:,1),centelem(:,2))];
    elem(:,end) = kmap(:,1);
elseif permFlag == 2
    %% gerando funcao de campo de permeabilidade
    ek = 5 * 10^-2;
    perISO = @(x,y) ek .^( 1*sqrt(4).*cos(6.*pi.*x).*cos(6.*pi.*y) );
   % perISO = @(x,y) ek .^( 0.5*sqrt(4).*sin(2.*pi.*(x-0.25)).*cos(2.*pi.*(y)) );
    weights = perISO(centelem(:,1),centelem(:,2)); 

   % flag = (centelem(:,1) <= 0.25 & centelem(:,2) <= 0.25) | (centelem(:,1) >= 0.75 & centelem(:,2) >= 0.75);
    %weights(flag) = 1;
    
    kmap = [ [1:size(elem,1)]'  weights zeros(size(elem,1),1) zeros(size(elem,1),1) weights];
    elem(:,end) = [1:size(elem,1)]';
    postprocessorName(kmap(:,2),0,superFolder, 'PermeHet')
elseif permFlag == 3
    alet = load ('chue.mat');
    alet = struct2cell(alet.chue);
    alet = cell2mat(alet);
    normP = @(x,y)(x.^2 + y.^2).^.5;
   % phiIso = @(x,y,numb) exp(-  ((1/0.05).*( normP((x - alet(numb,2)),((y - alet(numb,3))) ) ) ).^2 );
   
    
    phiIso = @(x,y,numb) exp(1).^(- ((normP((x - alet(numb,2)),((y - alet(numb,3)))) ./0.05).^2)  );

    
    acum = zeros(size(elem,1),1);
    for ii = 1:40
       acum = acum +  phiIso(centelem(:,1),centelem(:,2), ii*ones(size(elem,1),1));       
    end
    flagT = (acum <= 0.01);
    acum(flagT) = 0.01;
    
    flagT = (acum > 4);
    acum(flagT) = 4;
    weights = acum;
    %weights = 4*ones(size(elem,1),1) - weights;
    
    kmap = [ [1:size(elem,1)]'  weights zeros(size(elem,1),1) zeros(size(elem,1),1) weights];
    elem(:,end) = [1:size(elem,1)]';
    postprocessorOP(kmap,1,  superFolder, 'PermHet')
%    postprocessorName(kmap(:,2),0,superFolder, 'PermeHet')

elseif permFlag == 4
    alet = load('chuehartur.mat', 'perm_points')
    %alet = struct2cell(alet.chue);
    alet = alet.perm_points;
    normP = @(x,y)(x.^2 + y.^2).^.5;
   % phiIso = @(x,y,numb) exp(-  ((1/0.05).*( normP((x - alet(numb,2)),((y - alet(numb,3))) ) ) ).^2 );
   
    
    phiIso = @(x,y,numb) exp(1).^(- ((normP((x - alet(numb,2)),((y - alet(numb,3)))) ./0.05).^2)  );

    
    acum = zeros(size(elem,1),1);
    for ii = 1:40
       acum = acum +  phiIso(centelem(:,1),centelem(:,2), ii*ones(size(elem,1),1));       
    end
    flagT = (acum <= 0.01);
    acum(flagT) = 0.01;
    
    flagT = (acum > 4);
    acum(flagT) = 4;
    weights = acum;
    %weights = 4*ones(size(elem,1),1) - weights;
    
    kmap = [ [1:size(elem,1)]'  weights zeros(size(elem,1),1) zeros(size(elem,1),1) weights];
    pormap = pormap .* ones(size(elem,1),1);
    elem(:,end) = [1:size(elem,1)]';
    postprocessorName(kmap(:,2),0,superFolder, 'PermeHet')

elseif permFlag == 5
    %% gerando funcao de campo de permeabilidade
    ang = pi/2;
    R1 = [cos(ang) -sin(ang); sin(ang) cos(ang)]; 
    R2 = [cos(ang) sin(ang); -sin(ang) cos(ang)]; 
    KK = [100 0 ; 0 .01];
    Ka1 = R1*KK*R2;
    
    tag_k1 = centelem(:,1) <= 0.5;
    kmap = zeros(size(elem,1),5);
    kmap(:,1) = 1:size(elem,1);
    kmap(tag_k1,2) = KK(1,1);
    kmap(tag_k1,3) = KK(1,2);
    kmap(tag_k1,4) = KK(2,1);
    kmap(tag_k1,5) = KK(2,2);
    ek = 1 * 10^-3;
    li = 10^-3;
    kmap(~tag_k1,2) = (centelem(~tag_k1,2) + li) + ek*(centelem(~tag_k1,1) + li);
    kmap(~tag_k1,3) = -(1-ek)* ((centelem(~tag_k1,1) + li) .*(centelem(~tag_k1,2) + li) );
    kmap(~tag_k1,4) = -(1-ek)* ((centelem(~tag_k1,1) + li) .*(centelem(~tag_k1,2) + li) );
    kmap(~tag_k1,5) = (centelem(~tag_k1,1) + li) + ek*(centelem(~tag_k1,2) + li);
    
    
%     perISO = @(x,y) ek .^( 1*sqrt(4).*cos(6.*pi.*x).*cos(6.*pi.*y) );
%    % perISO = @(x,y) ek .^( 0.5*sqrt(4).*sin(2.*pi.*(x-0.25)).*cos(2.*pi.*(y)) );
%     weights = perISO(centelem(:,1),centelem(:,2)); 
% 
%    % flag = (centelem(:,1) <= 0.25 & centelem(:,2) <= 0.25) | (centelem(:,1) >= 0.75 & centelem(:,2) >= 0.75);
%     %weights(flag) = 1;
%     
%     kmap = [ [1:size(elem,1)]'  weights zeros(size(elem,1),1) zeros(size(elem,1),1) weights];
    elem(:,end) = [1:size(elem,1)]';
%     postprocessorName(kmap(:,2),0,superFolder, 'PermeHet')    
    
end



