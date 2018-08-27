
function[MoyCartesVolume, ProbVolume, Ecart_Type_Global, Sign, MoyGlobal] = AnalyseClusterGMM(data_in_table)


% MoyCartesVolume : moyenne et ecart type des valeurs des cartes par
%       cluster sur chaque volume, donc sur l'ensemble des tranches de
%       chaque couple Patient/Timepoint.
% ProbVolume : proportion de chaque cluster sur chaque volume, donc sur
%       l'ensemble des tranches de chaque couple Patient/Timepoint
% Ecart_Type_Global : ecart types des valeurs des cartes par cluster sur
%       toutes les donnees utilisees pour le clustering.
% Sign : Classement des couples Patient/Timepoint avec le groupe
%       d'appartenance. On peut ainsi relier chaque valeur des statistiques
%       � son patient.
% MoyGlobal : moyenne des valeurs des cartes par cluster sur toutes les
%       donnees utilisees pour le clustering.
data_in_table.Group = categorical(data_in_table.Group);
data_in_table.Patient = categorical(data_in_table.Patient);
data_in_table.Tp = categorical(data_in_table.Tp);

Pat = unique(data_in_table.Patient);
Noms = [];
Timp = [];
Group = [];

%Classement des patients et timepoints. On associe en plus le groupe de
%chaque patient et on obtient la matrice Sign
for i = 1: length(Pat)
    data_in_table_Pat = data_in_table(data_in_table.Patient == Pat(i),:);
    Timep = unique(data_in_table_Pat.Tp);
    for j = 1: length(Timep)
        Noms = [Noms Pat(i)];
        Timp = [Timp Timep(j)];
        Group = [Group unique(data_in_table_Pat.Group)];
    end
end
Sign = [Noms ; Timp ; Group];

%Une cellule par couple Patient/Timepoint
MoyCartesTranches = cell(1,size(Sign,2));
SommeTranches = cell(1,size(Sign,2));
ProbTranches = cell(1,size(Sign,2));
MoyCartesVolume = cell(1,size(Sign,2));
SommeVolume = cell(1,size(Sign,2));
ProbVolume = cell(1,size(Sign,2));

for h = 1:size(Sign,2)
    data_in_table_temp = data_in_table(data_in_table.Patient == Sign(1,h),:);  %On ne selectionne que les donnees correspondant au bon patient
    data_in_table_tmp = data_in_table_temp(data_in_table_temp.Tp == Sign(2,h),:);  %On ne selectionne que les donnees correspondant au bon timepoint
    %Z = unique(data_in_table_tmp.Coord_Z); %On recupere les valeurs de la coordonnee Z, c'est ainsi que l'on va separer les tranches.
    %MoyCartesTranches{h} = zeros(max(data_in_table.Cluster)*2,size(data_in_table_tmp,2)-8,length(Z));
    %SommeTranches{h} = zeros(max(data_in_table.Cluster),length(Z));
    %ProbTranches{h} = zeros(max(data_in_table.Cluster),length(Z));
    
%     for i=1:length(Z)                    %Tranches
%         for c = 1:max(data_in_table.Cluster)         % Clusters
%             for p=8:size(data_in_table_tmp,2)-1                   %Cartes
%                 if sum(double((data_in_table_tmp.Coord_Z == Z(i) & data_in_table_tmp.Cluster == c)))== 0
%                     %Si le cluster n'est pas present sur la tranche, on
%                     %fixe la valeur des statistiques a NaN, pour ne pas
%                     %perturber les prochains calculs de moyenne.
%                     MoyCartesTranches{h}(2*c-1,p-7,i) = NaN;
%                     MoyCartesTranches{h}(2*c,p-7,i) = NaN;
%                 else
%                     %Moyenne des valeurs des cartes des pixels de la
%                     %tranche de coordonnee Z = Z(i) et appartenant au
%                     %cluster c.
%                     MoyCartesTranches{h}(2*c-1,p-7,i) = nanmean(double(data_in_table_tmp(data_in_table_tmp.Coord_Z == Z(i) & data_in_table_tmp.Cluster == c,p)));  
%                     %Ecart type des valeurs des cartes des pixels de la
%                     %tranche de coordonnee Z = Z(i) et appartenant au
%                     %cluster c.
%                     MoyCartesTranches{h}(2*c,p-7,i) = nanstd(double(data_in_table_tmp(data_in_table_tmp.Coord_Z == Z(i) & data_in_table_tmp.Cluster == c,p)));
%                 end
%             end
%             if mean(double(data_in_table_tmp.Coord_Z == Z(i) & data_in_table_tmp.Cluster == c)) == 0
%                 SommeTranches{h}(c,i) = 0;
%             else
%                 %Somme de chaque cluster sur chacune des tranches concernees.
%                 SommeTranches{h}(c,i) = nansum(double(data_in_table_tmp.Cluster(data_in_table_tmp.Coord_Z == Z(i) & data_in_table_tmp.Cluster == c)))/c;   % Probabilites
%             end
%         end
%         %Proportion de chaque cluster sur chacune des tranches concernees.
%         ProbTranches{h}(:,i) = SommeTranches{h}(:,i)/nansum(SommeTranches{h}(:,i))*100;
%     end
%     
%     
%     %%Stats sur le volume (=Toutes les tranches)
%     
%     
%     SommeVolume{h} = nansum(SommeTranches{h},2);
%     
%     %Proportion de chaque cluster sur toutes les tranches de chaque couple
%     %Patient/Timepoint
%     ProbVolume{h} = SommeVolume{h}/nansum(SommeVolume{h})*100;
%     
    for c = 1:max(data_in_table.Cluster)
        for p=4:size(data_in_table_tmp,2)-1
            %Moyenne des valeurs des cartes des pixels appartenant au
            %cluster c sur toutes les tranches de chaque couple
            %Patient/Timepoint
            MoyCartesVolume{h}(2*c-1,p-3) = mean(table2array(data_in_table_tmp(data_in_table_tmp.Cluster ==c,p)));
            %Ecart type des valeurs des cartes des pixels appartenant au
            %cluster c sur toutes les tranches de chaque couple
            %Patient/Timepoint
            MoyCartesVolume{h}(2*c,p-3) = std(table2array(data_in_table_tmp(data_in_table_tmp.Cluster == c,p)));
            
            
            
            
            
            
            
            
            
            
            if mean(data_in_table_tmp.Cluster == c) == 0
                SommeVolume{h}(c) = 0;
            else
                %Somme de chaque cluster sur chacune des tranches concernees.
                SommeVolume{h}(c) = nansum(double(data_in_table_tmp.Cluster(data_in_table_tmp.Cluster == c)))/c;   % Probabilites
            end

        %Proportion de chaque cluster sur chacune des tranches concernees.
        ProbVolume{h} = SommeVolume{h}/nansum(SommeVolume{h})*100;

            
            
            
            
            
            
            
            
            
            
        end
    end
    
end

%On calcule les ecart types des valeurs des cartes par cluster sur toutes
%les donnees utilisees pour le clustering (= sur la database initiale :
%data_in_table)
Ecart_Type_Global = zeros(max(data_in_table.Cluster)*2,size(data_in_table_tmp,2)-4);
for c = 1:max(data_in_table.Cluster)
    for p=4:size(data_in_table_tmp,2)-1 
        Ecart_Type_Global(2*c,p-3) = nanstd(table2array(data_in_table(data_in_table.Cluster == c,p)));
    end
end

%On calcule les moyennes des valeurs des cartes par cluster sur toutes
%les donnees utilisees pour le clustering (= sur la database initiale :
%data_in_table)
MoyGlobal = zeros(max(data_in_table.Cluster)*2,size(data_in_table_tmp,2)-4);
for c = 1:max(data_in_table.Cluster)
    for p=4:size(data_in_table_tmp,2)-1 
        MoyGlobal(2*c-1,p-3) = nanmean(table2array(data_in_table(data_in_table.Cluster == c,p)));
    end
end



end