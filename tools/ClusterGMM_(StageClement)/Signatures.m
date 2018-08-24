
%Si on ne veut afficher qu'une seule signature :
% b = bar([ProbVolume1.ProbVolume.'; nan(size(ProbVolume1.ProbVolume.'))],'Stacked');
% set(gca, 'xtick',1:size(ProbVolume1.ProbVolume,1));

function Signatures(Informations,Statistiques, Couleurs)

% La syntaxe classique de la creation d'histogramme ne permet pas de
% n'afficher qu'un seul histogramme. Pour faire cela, il faut la modifier
% legerement. On doit donc traiter les deux cas (affichage d'un ou de
% plusieurs histogrammes) separement.

% Cas ou le clustering ne concerne qu'un seul couple Patient/Timepoint
% (tres rare ...). On ne doit afficher qu'un histogramme (ie qu'une signature).
% PAS ENCORE TESTE 
if size(Informations(1).Sign,2) == 1
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    b = bar([Statistiques(1).ProbVolume; nan(size(Statistiques(1).ProbVolume))],'Stacked');
    set(gca, 'xtick',1:size(Statistiques(1).ProbVolume,1));
    
    for ii = 1:length(Statistiques(1).ProbVolume)
        b(1,ii).FaceColor = Couleurs(ii,:);
    end
    
    MoyClust = zeros(size(Statistiques(1).MoyCartesVolume));
    for i = 1:length(Statistiques)
        MoyClust = MoyClust + Statistiques(i).MoyCartesVolume/length(Statistiques);
    end
    
    
    colonnes = cell(1,length(Informations.Cartes));
    lignes = cell(1,length(Statistiques(1).ProbVolume)*2);
    for i = 1:length(Informations.Cartes)
        colonnes{1,i} = char(Informations.Cartes(i));
    end
    
    for i = 1:size(Statistiques(1).MoyCartesVolume,1)/2
        colergen = @(color,text) ['<html><table border=0 width=20 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
        lignes{1,2*i-1} = strcat('Cluster_',num2str(i),'_Mean');
        lignes{1,2*i} = strcat('Cluster_',num2str(i),'_SD');
    end
    
    
    f1 = gcf;
    t1 = uitable(f1,'Data',MoyClust,'ColumnName',colonnes,'RowName',lignes);%,'FontS',8);
    
    
    t1.Position(1) = 300;
    t1.Position(2) = 300;
    t1.Position(3) = t1.Extent(3);
    t1.Position(4) = t1.Extent(4);
    
    
    subplot(2,2,4)
    radarPlot(MoyClust.',Informations)
    
    
    
else            % Cas courant, avec plusieurs couples Patient/Timepoint
    
    %On va afficher des tableaux, on definit les titres de leurs lignes et
    %colonnes.
    colonnes = cell(1,length(Informations.Cartes));
    lignes = cell(1,length(Statistiques(1).ProbVolume)*2);
    for i = 1:length(Informations.Cartes)
        colonnes{1,i} = char(Informations.Cartes(i));
    end
    
    for i = 1:length(Statistiques(1).ProbVolume)
        %colergen = @(color,text) ['<html><table border=0 width=20 bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
        lignes{1,2*i-1} = strcat('Cluster_',num2str(i),'_Mean');
        lignes{1,2*i} = strcat('Cluster_',num2str(i),'_SD');
    end
    
   
    
    ProbasClass = [];
    IndClass = [];
    NomsGroupes = [];
    UniGroup = unique(Informations.Sign(3,:));
    TP = [];
    % On classe les donnees par groupe pathologique
    for k = 1:length(UniGroup)
        for i = 1:length(Statistiques)
            if strcmp(char(Informations.Sign(3,i)),char(UniGroup(k)))
                ProbasClass = [ProbasClass; Statistiques(i).ProbVolume];
                IndClass = [IndClass Informations.Sign(1,i)];
                TP = [TP Informations.Sign(2,i)];
            end
        end
        Longueur(k) = sum(Informations.Sign(3,:) == UniGroup(k));
        NomsGroupes = strcat(char(NomsGroupes), char(UniGroup(k)),', ');
    end
    
    % Permet d'agrandir automatiquement la figure creee.
    figure('units','normalized','outerposition',[0 0 1 1])
    
    subplot(2,2,1)
    %Histogramme ( = Signature) de chaque couple Patient/Timepoint. Les
    %histogrammes sont classes par groupes pathologiques.
    b = bar(ProbasClass,'stacked');  
    
    % On attribue � chaque partie des histogrammes la couleur du cluster
    % qu'elle represente.
    for ii = 1:length(Statistiques(1).ProbVolume)
        b(1,ii).FaceColor = Couleurs(ii,:);
    end
    
    title(strcat('Signatures ',num2str(length(Statistiques(1).ProbVolume)), ' Clusters, ',num2str(length(colonnes)), ' Cartes'));
    ylabel('Pourcents')
    xlabel(NomsGroupes)
    
    
    hold on
    long = 0;
    % On trace des lignes verticales entre les differents groupes
    % pathologiques afin de bien les differencier.
    for k = 1:length(UniGroup)-1
        long = long + Longueur(k);
        line([long+0.5 long+0.5],[0 120],'LineWidth',4)
    end
    hold off
    

colonnesgroups = cell(1,length(UniGroup));
ProbStat = zeros(length(Statistiques(1).ProbVolume),length(UniGroup));
EcTy = zeros(length(Statistiques(1).ProbVolume),length(UniGroup));
long = 1;
% Calcul des statistiques (moyenne et ecart type) des signatures moyennes.
for k = 1:length(UniGroup)
    for j = 1:length(Statistiques(1).ProbVolume)
        ProbStat(j,k) = nanmean(ProbasClass(long:(long+Longueur(k)-1),j));
        EcTy(j,k) = nanstd(ProbasClass(long:(long+Longueur(k)-1),j));
    end
    long = long + Longueur(k);
    colonnesgroups{1,k} = strcat(char(UniGroup(k)), ' (%)');
end

MoyProb = zeros(size(Statistiques(1).MoyCartesVolume,1),length(UniGroup));
for i = 1:size(Statistiques(1).MoyCartesVolume,1)/2
    MoyProb(2*i-1,:) = ProbStat(i,:);
    MoyProb(2*i,:) = EcTy(i,:);
end
    
% Si il n'y a qu'un seul groupe (assez courant), rebelote, il faut modifier la syntaxe afin
% de n'afficher qu'un seul histogramme
if length(UniGroup) == 1
    subplot(2,2,2)
    bmoy = bar([ProbStat.'; nan(size(ProbStat'))],'Stacked');
    set(gca, 'xtick',1:size(ProbStat,1));
    
    
else         %On revient au cas classique � plusieurs groupes donc plusieurs histogrammes
    subplot(2,2,2)
    bmoy = bar(ProbStat.','stacked');
end

% On fait correspondre les couleurs des histogrammes et des clusters.
for ii = 1:length(Statistiques(1).ProbVolume)
    bmoy(ii).FaceColor = Couleurs(ii,:);
end
title(strcat('Signatures moyennes ',num2str(length(Statistiques(1).ProbVolume)), ' Clusters, ',num2str(length(colonnes)), ' Cartes'))
ylabel('Pourcents')
xlabel(NomsGroupes)


% Lignes a priori inutiles avec la correction de l'erreur de calcul des
% moyennes mais au cas ou ...
ConvertStat = struct2cell(Statistiques);
MoyClust = nanmean(cell2mat(ConvertStat(3,1,:)),3);


    
%Lignes qui permettent de corriger les erreurs de calcul de moyenne et
%d'ecart type. L'erreur etait de calculer la moyenne globale a partir des
%moyennes de chaque patient, alors que ces patients ne comportent pas le
%meme nombre de pixels et biaisent donc le poids de leur moyennes. J'ai
%donc rajoute dans la fonction AnalyseClusterGMM le calcul global de la
%moyenne et des ecarts types. C'est ces valeurs que l'on affiche desormais.
for i = 1:size(Statistiques(1).MoyCartesVolume,1)/2
    MoyClust(2*i-1,:) = Statistiques(1).MoyGlobal(2*i-1,:);
    MoyClust(2*i,:) = Statistiques(1).Ecart_Type_Global(2*i,:);
end

    
    % Ici on choisit la position initiale des tableaux. Il peut etre utile
    % de savoir qu'on peut modifier les positions des objets, leurs titres,
    % leurs legendes, etc... directement sur la figure en activant l'option
    % "edit plot" dans le menu "Tools".
    
    f1 = gcf;
    t = uitable(f1,'Data',MoyClust,'ColumnName',colonnes,'RowName',lignes);%,'FontS',8);
    
    t.Position(1) = 250;
    t.Position(2) = 275;
    t.Position(3) = t.Extent(3);
    t.Position(4) = t.Extent(4);
    
    
    t2 = uitable(f1,'Data',MoyProb,'ColumnName',colonnesgroups,'RowName',lignes);%,'FontS',8);
    
    t2.Position(1) = 250;
    t2.Position(2) = 0;
    t2.Position(3) = t2.Extent(3);
    t2.Position(4) = t2.Extent(4);
    
    NomsPatients =  cell(length(IndClass),1);
    for k = 1:length(IndClass)
        NomsPatients{k,1} = char(IndClass(k));
        NomsPatients{k,2} = char(TP(k));
    end
    col = cell(1,2);
    col{1,1} = 'Patient';
    col{1,2} = 'TimePoint';
    
    
    t3 = uitable(f1,'Data',NomsPatients,'ColumnName',col);
    
    
    t3.Position(3) = t3.Extent(3);
    t3.Position(4) = t3.Extent(4);
    t3.Position(1) = 0;
    t3.Position(2) = 980-t3.Position(4);
  
    subplot(2,2,4)
    
    %Puisque le graphique araignee n'affiche que les moyennes, il est
    %inutile de conserver les ecarts types.
    for i = 1:size(MoyClust,1)/2
        MoyClust2(i,:) = MoyClust(2*i-1,:);
    end
    
    
    % La fonction radarplot (recuperee sur le file exchange) parmet de
    % tracer des graphiques araignee. Je l'ai legerement modifie pour
    % qu'elle affiche les bonnes couleurs et les bons titres.
    radarPlot(MoyClust2.',Informations, Couleurs)

   
   
    
end
end
