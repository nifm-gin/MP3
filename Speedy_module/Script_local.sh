#!/bin/bash

# Command to run the standard cluster profil for this script :
# oar=1 core=8 walltime=48 cluster=mistis repetitionNumber=8 referenceModelTime=new referenceDelineationTime=new atypicalModelTime=new testDataTime=new ./Script_v2016_02_08.sh

# Commande pour lancer en local
# oar=0 core=1 repetitionNumber=1 referenceModelTime=new referenceDelineationTime=new atypicalModelTime=new testDataTime=new K=9 quantileOfInterestOnWeights=0.01 KSignature=5 nbCores=4 ./Script_v2016_02_08.sh
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#   Script du protocole de traitement des images IRM via                       #
#   un modele de melange de lois, en particulier via                           #
#   l'algorithme EM pour les melanges de distributions multiple scaled         #
#                                                                              #
#   Ajustement du script pour les donnees Hommes de Benjamin 2016-06           #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Dates d'execution du script et des points de sauvegarde :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Date d'execution du script
# . identifiant unique pour sauvegarder les resultats du protocole (un format de
#   date n'est pas necessaire, une chaine de caractere suffit)
# . valeur par defaut : l'heure actuelle du systeme
if [ -z "$launchTime" ]; then
    launchTime=$(date +'%Y_%m_%d_-_%H_%M_%S')
fi

# Identifiant du modele de reference a utiliser
# . valeur par defaut : 2016_03_25_-_09_36_33
# . valeur particuliere : new -> creation d'un nouveau model de reference
#   ayant pour identifiant '$launchTime'
if [ -z "$referenceModelTime" ]; then
    referenceModelTime="2016_06_04_-_09_22_33"
fi

# Identifiant de la delineation a utiliser
# . valeur par defaut : 2016_03_25_-_09_36_33
# . valeur particuliere : new -> creation d'une nouvelle delineation
#   ayant pour identifiant '$launchTime'
if [ -z "$referenceDelineationTime" ]; then
    referenceDelineationTime="2016_03_25_-_09_36_33"
fi

# Identifiant du modele atypique a utiliser
# . valeur par defaut : 2016_03_25_-_09_36_33
# . valeur particuliere : new -> creation d'un nouveau model atypique
#   ayant pour identifiant '$launchTime'
if [ -z "$atypicalModelTime" ]; then
    atypicalModelTime="2016_03_25_-_09_36_33"
fi

# Identifiant du test a effectuer
# . valeur par defaut : 2016_03_25_-_09_36_33
# . valeur particuliere : new -> creation d'un nouveau test
#   ayant pour identifiant '$launchTime'
if [ -z "$testDataTime" ]; then
    testDataTime="2016_06_04_-_09_22_33"
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Protocole a executer :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Chemin absolu du dossier contenant les fichiers a executer
# . valeur par defaut : le dossier courant
if [ -z "$filesDirectory" ]; then
    # filesDirectory="/Users/Lemasson/Desktop/Code/"
    filesDirectory=$(pwd)
fi

# Chemin relatif du dossier contenant les donnees :
# . valeur par defaut : le dossier "Data"
if [ -z "$dataDirectory" ]; then
    # dataDirectory=$(pwd)/Data
    dataDirectory="/services/scratch/mistis/arnaud/Code_R/Humans_data"
fi

# Chemin absolu du dossier pour sauvegarder les resultats :
# . valeur par defaut : le scratch du cluster a l'Inria
#   "/services/scratch/mistis/arnaud/Code_R"
if [ -z "$saveDirectory" ]; then
   # saveDirectory=$(pwd)
   saveDirectory="/services/scratch/mistis/arnaud/Code_R"
fi

# Fichier R a executer
# . valeur par defaut : la derniere version du protocole
#   "00_Protocol_v2016_02_08.R"
if [ -z "$file" ]; then
    file="00_Protocol_v2016_02_08.R"
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Parametres du cluster OAR de l'INRIA :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Utilisation du cluster OAR
# . valeur par defaut : 1 -> utilisation de OAR
# . valeur specifique : 0 -> lancement en local du code R
if [ -z "$oar" ]; then
    oar=0
fi

# Nombre de noeuds a reserver dans OAR
# . valeur par defaut : 1
# . le code R actuel ne peut pas gerer plus de 1 noeud
if [ -z "$nodes" ]; then
    nodes=1   
fi

# Nombre de coeurs a reserver par noeud dans OAR
# . valeur par defaut : 1
if [ -z "$core" ]; then
    core=1
fi

# Duree de reservation dans OAR
# . valeur par defaut : 6:00:00 -> 6h
# . il faut prevoir 30h pour l'ajustement complet avec 20 classes et MMSD
if [ -z "$walltime" ]; then
    walltime=6:00:00
fi

# Nom du cluster sur lequel faire la reservation
# . valeur par defaut : ("mistis" "SIC") -> repartition de la reservation sur
#   les clusters mistis et SIC
# . pour specifier N clusters, il faut utiliser la syntaxe suivante
#   ("cluster1" "cluster2" ... "clusterN")
if [ -z "$cluster" ]; then
    cluster=("mistis" "SIC")
fi

# Modification de la variable 'cluster' pour quelle soit reconnue par OAR
if [ "${#cluster[*]}" -gt 1 ]; then
    temp="cluster='${cluster[0]}'"
    for i in $(seq 1 1 $(( ${#cluster[*]} - 1 )) );
    do
    temp="$temp OR cluster='${cluster[$i]}'"
    done
else
    temp="cluster='$cluster'"
fi
cluster=$temp
unset temp

# Utilisation du type besteffort pour le job a soumettre (ie priorite minimale)
# . valeur par defaut : 
if [ -z "$besteffort" ]; then
    besteffort=0
fi



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Parametres du mode parallele :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Execution parallele :
if [ -z "$tryParallel" ]; then
    tryParallel=TRUE
fi

# Detection automatique du nombre de coeurs :
if [ -z "$detectCores" ]; then
    detectCores=FALSE
fi

# Nombre de coeurs a utiliser :
if [ -z "$nbCores" ]; then
    nbCores=$((nodes*core))
    # nbCores=4
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Parametres du mode d'apprentissage :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Noms des parametres IRM :
if [ -z "$nomPar" ]; then
    #nomPar=rresliceFLAIRSENSEone,rresliceToneW_SE_POSTCLEARone
    #nomPar=Ttwomap,DCErehaus,rCBV,resliceDWI-SENSE
    nomPar=Ttwomap,DCErehaus,rCBV
fi

# Plages des valeurs admissibles :
if [ -z "$valPhysio" ]; then
 #ADC #AUC#AUCmax #AUCrehaus #AUCttp   #T1map#T2map#T2starPre#T2starPost  #T2starCorr3D  #deltaR2   #deltaR2star#R2prim#CBV #CBF#SO2 #VSI #Density#CMRO2 #rresliceFLAIRSENSEone#rresliceToneW_SE_POSTCLEARone #DCEAUC#Ttwomap #DCErehaus#rCBV#reslicedDWISENSE #Coreg_Rtwoprim#DCEMax#DCEttp#MTT#Tzero#TTP,#rCBF
    valPhysio=0,4e3,-2e5,1.5e6,2e4,2e5,0,150,0,300,0,1.5e4,0,1e3,0,3e2,0,1.5e2,0,3e2,-1,1,-1,1,-5e1,5e1,-5e1,5e1,-4e2,8e2,-3e2,3e2,0,1e2,-2e2,3e3,-5e1,1e2,0,5.0912e4,0,1.4692e5,0,1e8,0,5e4,0,1e3,0,1e2,0,1e7,0,50,0,1e8,0,1e8,0,50,0,1000,0,1000,0,1000          

    # valPhysio=0,4e3,  #ADC
    # -2e5,1.5e6      #AUC
    # ,2e4,2e5         #AUCmax
    # ,0,150           #AUCrehaus
    # ,0,300           #AUCttp
    # ,0,1.5e4         #T1map
    # ,0,1e3           #T2map
    # ,0,3e2           #T2starPre
    # ,0,1.5e2         #T2starPost
    # ,0,3e2           #T2starCorr3D
    # ,-1,1            #deltaR2
    # ,-1,1            #deltaR2star
    # ,-5e1,5e1        #R2prim
    # ,-5e1,5e1        #CBV
    # ,-4e2,8e2        #CBF
    # ,-3e2,3e2        #SO2
    # ,0,1e2           #VSI
    # ,-2e2,3e3        #Density
    # ,-5e1,1e2        #CMRO2
    # ,0,5.0912e4      #rresliceFLAIRSENSEone
    # ,0,1.4692e5      #rresliceToneW_SE_POSTCLEARone
    # ,0,1e8           #DCEAUC
    # ,0,5e4           #Ttwomap
    # ,0,1e3           #DCErehaus
    # ,0,1e2           #rCBV
    # ,0,1e7           #reslicedDWISENSE
    # ,0,50            #Coreg_Rtwoprim
    # ,0,1e8           #DCEMax
    # ,0,1e8           #DCEttp
    # ,0,50            #MTT
    # ,0,1000          #Tzero
    # ,0,1000          #TTP
    # ,0,1000          #rCBF
fi

# Noms des plages de valeurs admissibles :
if [ -z "$valPhysioNames" ]; then
    valPhysioNames=ADC,AUC,AUCmax,AUCrehaus,AUCttp,T1map,T2map,T2starPre,T2starPost,T2starCorr3D,deltaR2,deltaR2star,R2prim,CBV,CBF,SO2,VSI,Density,CMRO2,rresliceFLAIRSENSEone,rresliceToneW_SE_POSTCLEARone,DCEAUC,Ttwomap,DCErehaus,rCBV,reslicedDWISENSE,Coreg_Rtwoprim,DCEMax,DCEttp,MTT,Tzero,TTP,rCBF
fi

# Type du modele EM a utiliser :
if [ -z "$modelEm" ]; then
    modelEm="mmsd"
fi

# Type de la methode pour choisir le nombre de clusters a utiliser :
if [ -z "$clusterNumberChoice" ]; then
    clusterNumberChoice="sh"
fi

# Nombre de points consecutifs a prendre en compte pour l'heuristique de pente :
if [ -z "$slopeHeuristicPointNumber" ]; then
    slopeHeuristicPointNumber=4
fi

# Nombre d'iterations pour "mmsd.estimate" :
if [ -z "$iterNumber" ]; then
    iterNumber=1e3
fi

# Methode d'optimisation pour les matrices orthogonales :
if [ -z "$orthogonalMethod" ]; then
    orthogonalMethod=S
    # orthogonalMethod=FG
fi

# Nombre d'iterations pour la methode d'optimisation sur les matrices orthogonales :
if [ -z "$orthogonalIterations" ]; then
    orthogonalIterations=1e2
fi

# Nombre de repetitionNumber pour la classification :
if [ -z "$repetitionNumber" ]; then
    repetitionNumber=2
fi

# Nombre de groupes :
if [ -z "$K" ]; then
    K=20
fi

# Vecteur des nombres de groupes a tester :
if [ -z "$KVect" ]; then
    KVect=NULL
fi

# Doit-on enlever les classes vides :
if [ -z "$isEmptyClusterDeleted" ]; then
    isEmptyClusterDeleted=FALSE
fi

# Doit-on ne garder que les groupes majoritaires :
if [ -z "$keepOnlyMajorGroups" ]; then
    keepOnlyMajorGroups=FALSE
fi

# Combien de groupes majoritaires garder :
if [ -z "$nbMajorGroups" ]; then
    nbMajorGroups=-1
fi

# Est-ce qu'on enleve les donnees de reference apres la partie de delimitation automatique :
if [ -z "$removeReferenceDataAfterAutomaticDelimitation" ]; then
    removeReferenceDataAfterAutomaticDelimitation=false
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Parametres de construction des signatures :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Noms des ROI a utiliser pour l'analyse discriminante :
if [ -z "$voiDA" ]; then
    voiDA=(Brain)
fi

# Methodes de segmentation des ROI :
if [ -z "$voiSegmentationType" ]; then
    # voiSegmentationType="classes"
    voiSegmentationType="threshold"
fi

# Seuil sur l'image des poids pour la construction des signatures :
if [ -z "$thresholdOfInterestOnWeights" ]; then
    thresholdOfInterestOnWeights=-1
fi

# Seuil d'exclusion des patients atypiques par rapport au pourcentage de voxels non atypiques :
if [ -z "$thresholdOfAtypicalOnInlierBrain" ]; then
    thresholdOfAtypicalOnInlierBrain=0
fi

# Quantile des poids sains a considerer pour la construction des signatures :
if [ -z "$quantileOfInterestOnWeights" ]; then
    quantileOfInterestOnWeights=1e-2
fi

# Utilisation des regles d'experts pour l'exclusion de donnees :
if [ -z "$useExpertRules" ]; then
    useExpertRules=FALSE
fi

# Nombre de groupes pour la segmentation des ROI :
if [ -z "$KRoi" ]; then
    KRoi=2
fi

# Vecteur des nombres de groupes a tester pour la segmentation des ROI :
if [ -z "$KRoiVect" ]; then
    KRoiVect=NULL
fi

# Nombre de groupes pour la caracterisation des ROI :
if [ -z "$KSignature" ]; then
    KSignature=20
fi

# Vecteur des nombres de groupes a tester pour la caracterisation des ROI  :
if [ -z "$KSignatureVect" ]; then
    KSignatureVect=NULL
fi

# Fonctions d'analyse discriminante :
if [ -z "$discriminantAnalysis" ]; then
    discriminantAnalysis=("lda","hdda","mclustda")
fi

# Models pour l'analyse lda :
if [ -z "$ldaMethode" ]; then
    ldaMethode=("moment","mle","mve","t")
fi

# Models pour l'analyse HDDA :
if [ -z "$hddaModels" ]; then
    hddaModels=("ABQD","AjBQD","ABQkD","AkBQkD","AkjBQkD","ABkQkD","AkBkQkD","AkjBkQkD","ABQkDk","AkBQkDk","AkjBQkDk","ABkQkDk","AkBkQkDk","AkjBkQkDk")
fi

# Models pour l'analyse MclustDA :
if [ -z "$mclustdaModels" ]; then
    mclustdaModels=("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV")
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Parametres graphiques :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sauvegarde des images individuelles :
if [ -z "$saveFigureDirectory" ]; then
    saveFigureDirectory=TRUE
fi

# Fonction pour generer les couleurs des classes :
if [ -z "$colorFunction" ]; then
    colorFunction=rainbow
fi

# Couleur de fond des graphiques :
if [ -z "$bgPlot" ]; then
    bgPlot=black
fi

# Couleur de fond du texte :
if [ -z "$bgPlotText" ]; then
    bgPlotText=white
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Fichiers a importer :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Fichiers pour l'apprentissage :
if [ -z "$referenceType" ]; then
    referenceType="AVC,BG,GBM,Leucemie,Menin,Meta"
    # referenceType="GBM"
fi

if [ -z "$referenceTypeRoi" ]; then
    referenceTypeRoi="Contro"
fi

if [ -z "$listFileName" ]; then
    listFileName="Humans_AVC_allMaps_ROIContro,Humans_BG_allMaps_ROIContro,Humans_GBM_allMaps_ROIContro,Humans_Leucemie_allMaps_ROIContro,Humans_Menin_allMaps_ROIContro,Humans_Meta_allMaps_ROIContro"
    # listFileName="Humans_GBM_allMaps_ROIContro"

    # listFileName="Humans_GBMetBG_allMaps_ROITumorflair"
    # listFileName="Humans_allPatients_allMaps_ROIContro"
    # listFileName="Rats_Sain_-_Brain_-_Lemasson_2016"
    # ,Rats_Sain_-_Brain_-_Lemasson_2016
    # ,Rats_Sain_-_Tumor_-_Lemasson_2016
fi

if [ -z "$listPatientToRemove" ]; then
    listPatientToRemove="NULL"
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Fichiers pour la prediction :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if [ -z "$atypicalType" ]; then
    # atypicalType="Sain,9L,C6,F98,RG2"
    atypicalType="AVC,BG,GBM,Leucemie,Menin,Meta"
    # atypicalType="GBM"
fi

if [ -z "$atypicalTypeRoi" ]; then
    atypicalTypeRoi="Brain"
fi

if [ -z "$reconListFileName" ]; then
reconListFileName="Humans_AVC_allMaps_ROIBrain,Humans_BG_allMaps_ROIBrain,Humans_GBM_allMaps_ROIBrain,Humans_Leucemie_allMaps_ROIBrain,Humans_Menin_allMaps_ROIBrain,Humans_Meta_allMaps_ROIBrain"
    # reconListFileName="Humans_GBM_allMaps_ROIBrain"
    #reconListFileName="Humans_allPatients_allMaps_ROIBrain,Humans_allPatients_allMaps_ROIContro,Humans_allPatients_allMaps_ROITumorflair,Humans_allPatients_allMaps_ROITumorgado"
    #reconListFileName="Rats_Sain_-_Brain_-_Lemasson_2016,Rats_9L_-_Brain_-_Lemasson_2016,Rats_C6_-_Brain_-_Lemasson_2016,Rats_F98_-_Brain_-_Lemasson_2016,Rats_RG2_-_Brain_-_Lemasson_2016"
    # ,Rats_Sain_-_Brain_-_Lemasson_2016
    # ,Rats_Sain_-_Tumor_-_Lemasson_2016
    # ,Rats_9L_-_Brain_-_Lemasson_2016
    # ,Rats_9L_-_Tumor_-_Lemasson_2016
    # ,Rats_C6_-_Brain_-_Lemasson_2016
    # ,Rats_C6_-_Tumor_-_Lemasson_2016
    # ,Rats_F98_-_Brain_-_Lemasson_2016
    # ,Rats_F98_-_Tumor_-_Lemasson_2016
    # ,Rats_RG2_-_Brain_-_Lemasson_2016
    # ,Rats_RG2_-_Tumor_-_Lemasson_2016
fi

if [ -z "$reconListPatientToRemove" ]; then
    reconListPatientToRemove="NULL"
fi

# Roi pour la comparaison :
if [ -z "$comparisonAtypicalTypeRoi" ]; then
comparisonAtypicalTypeRoi="tumor-flair"
fi

# Fichiers pour la comparaison :
if [ -z "$comparisonListFileName" ]; then
    comparisonListFileName="Humans_AVC_allMaps_ROITumorflair,Humans_BG_allMaps_ROITumorflair,Humans_GBM_allMaps_ROITumorflair,Humans_Leucemie_allMaps_ROITumorflair,Humans_Menin_allMaps_ROITumorflair,Humans_Meta_allMaps_ROITumorflair"
#comparisonListFileName="Humans_GBM_allMaps_ROITumorflair"

    # ,Rats_Sain_-_Tumor_-_Lemasson_2016
    # ,Rats_9L_-_Tumor_-_Lemasson_2016
    # ,Rats_C6_-_Tumor_-_Lemasson_2016
    # ,Rats_F98_-_Tumor_-_Lemasson_2016
    # ,Rats_RG2_-_Tumor_-_Lemasson_2016
fi

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Fichiers de test :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if [ -z "$testType" ]; then
    testType="AVC,BG,GBM,Leucemie,Menin,Meta"
        #testType="Sain,9L,C6,F98,RG2,Blind_test_9L,Blind_test_C6,Blind_test_F98,Blind_test_Sain,Outlier_9L,Outlier_C6"
fi

if [ -z "$testTypeRoi" ]; then
    testTypeRoi="Brain"
fi

if [ -z "$testListFileName" ]; then
testListFileName="Humans_AVC_allMaps_ROIBrain,Humans_BG_allMaps_ROIBrain,Humans_GBM_allMaps_ROIBrain,Humans_Leucemie_allMaps_ROIBrain,Humans_Menin_allMaps_ROIBrain,Humans_Meta_allMaps_ROIBrain"
    
    #testListFileName="Rats_Sain_-_Brain_-_Lemasson_2016,Rats_9L_-_Brain_-_Lemasson_2016,Rats_C6_-_Brain_-_Lemasson_2016,Rats_F98_-_Brain_-_Lemasson_2016,Rats_RG2_-_Brain_-_Lemasson_2016,Rats_Blind_test_9L_-_Brain_-_Lemasson_2016,Rats_Blind_test_C6_-_Brain_-_Lemasson_2016,Rats_Blind_test_F98_-_Brain_-_Lemasson_2016,Rats_Blind_test_Sain_-_Brain_-_Lemasson_2016,Rats_Outlier_9L_-_Brain_-_Lemasson_2016,Rats_Outlier_C6_-_Brain_-_Lemasson_2016"
    # ,Rats_Sain_-_Brain_-_Lemasson_2016
    # ,Rats_9L_-_Brain_-_Lemasson_2016
    # ,Rats_C6_-_Brain_-_Lemasson_2016
    # ,Rats_F98_-_Brain_-_Lemasson_2016
    # ,Rats_RG2_-_Brain_-_Lemasson_2016
    # ,Rats_Blind_test_9L_-_Brain_-_Lemasson_2016
    # ,Rats_Blind_test_C6_-_Brain_-_Lemasson_2016
    # ,Rats_Blind_test_Sain_-_Brain_-_Lemasson_2016
    # ,Rats_Outlier_9L_-_Brain_-_Lemasson_2016
    # ,Rats_Outlier_C6_-_Brain_-_Lemasson_2016
fi

if [ -z "$testListPatientToRemove" ]; then
    testListPatientToRemove="NULL"
fi

if [ -z "$testComparisonListFileName" ]; then
    testComparisonListFileName="Humans_AVC_allMaps_ROITumorflair,Humans_BG_allMaps_ROITumorflair,Humans_GBM_allMaps_ROITumorflair,Humans_Leucemie_allMaps_ROITumorflair,Humans_Menin_allMaps_ROITumorflair,Humans_Meta_allMaps_ROITumorflair"

    #testComparisonListFileName="Rats_Sain_-_Tumor_-_Lemasson_2016,Rats_9L_-_Tumor_-_Lemasson_2016,Rats_C6_-_Tumor_-_Lemasson_2016,Rats_F98_-_Tumor_-_Lemasson_2016,Rats_RG2_-_Tumor_-_Lemasson_2016,Rats_Blind_test_9L_-_Tumor_-_Lemasson_2016,Rats_Blind_test_C6_-_Tumor_-_Lemasson_2016,Rats_Blind_test_F98_-_Tumor_-_Lemasson_2016,Rats_Blind_test_Sain_-_Tumor_-_Lemasson_2016,Rats_Outlier_9L_-_Tumor_-_Lemasson_2016,Rats_Outlier_C6_-_Tumor_-_Lemasson_2016"
    # ,Rats_Sain_-_Tumor_-_Lemasson_2016
    # ,Rats_9L_-_Tumor_-_Lemasson_2016
    # ,Rats_C6_-_Tumor_-_Lemasson_2016
    # ,Rats_F98_-_Tumor_-_Lemasson_2016
    # ,Rats_RG2_-_Tumor_-_Lemasson_2016
    # ,Rats_Blind_test_9L_-_Tumor_-_Lemasson_2016
    # ,Rats_Blind_test_C6_-_Tumor_-_Lemasson_2016
    # ,Rats_Blind_test_Sain_-_Tumor_-_Lemasson_2016
    # ,Rats_Outlier_9L_-_Tumor_-_Lemasson_2016
    # ,Rats_Outlier_C6_-_Tumor_-_Lemasson_2016
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Sauvegarde des calculs :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sauvegarde du fichier RData :
if [ -z "$saveRData" ]; then
    saveRData=TRUE
fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Lancement du protocole :
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
commandToExecute="R --vanilla --file=$file --args --oar=$oar --launchTime=$launchTime --referenceModelTime=$referenceModelTime --referenceDelineationTime=$referenceDelineationTime --atypicalModelTime=$atypicalModelTime --testDataTime=$testDataTime --filesDirectory=$filesDirectory --dataDirectory=$dataDirectory --saveDirectory=$saveDirectory --tryParallel=$tryParallel --detectCores=$detectCores --nbCores=$nbCores --nomPar=$nomPar --valPhysio=$valPhysio --valPhysioNames=$valPhysioNames --modelEm=$modelEm --clusterNumberChoice=$clusterNumberChoice --slopeHeuristicPointNumber=$slopeHeuristicPointNumber --iterNumber=$iterNumber --orthogonalMethod=$orthogonalMethod --orthogonalIterations=$orthogonalIterations --repetitionNumber=$repetitionNumber --K=$K --KVect=$KVect --isEmptyClusterDeleted=$isEmptyClusterDeleted --keepOnlyMajorGroups=$keepOnlyMajorGroups --nbMajorGroups=$nbMajorGroups --removeReferenceDataAfterAutomaticDelimitation=$removeReferenceDataAfterAutomaticDelimitation --voiDA=$voiDA --voiSegmentationType=$voiSegmentationType --thresholdOfInterestOnWeights=$thresholdOfInterestOnWeights --thresholdOfAtypicalOnInlierBrain=$thresholdOfAtypicalOnInlierBrain --quantileOfInterestOnWeights=$quantileOfInterestOnWeights --useExpertRules=$useExpertRules --KRoi=$KRoi --KRoiVect=$KRoiVect --KSignature=$KSignature --KSignatureVect=$KSignatureVect --discriminantAnalysis=$discriminantAnalysis --ldaMethode=$ldaMethode --hddaModels=$hddaModels --mclustdaModels=$mclustdaModels --saveFigureDirectory=$saveFigureDirectory --colorFunction=$colorFunction --bgPlot=$bgPlot --bgPlotText=$bgPlotText --referenceType=$referenceType --referenceTypeRoi=$referenceTypeRoi --listFileName=$listFileName --listPatientToRemove=$listPatientToRemove --atypicalType=$atypicalType --atypicalTypeRoi=$atypicalTypeRoi --reconListFileName=$reconListFileName --reconListPatientToRemove=$reconListPatientToRemove --comparisonAtypicalTypeRoi=$comparisonAtypicalTypeRoi --comparisonListFileName=$comparisonListFileName --testType=$testType --testTypeRoi=$testTypeRoi --testListFileName=$testListFileName --testListPatientToRemove=$testListPatientToRemove --testComparisonListFileName=$testComparisonListFileName --saveRData=$saveRData"

if [ "$oar" == 0 ]; then
    $commandToExecute
elif [ $besteffort == 1 ]; then
    oarsub -t besteffort -t idempotent -l /nodes=$nodes/core=$core,walltime=$walltime -p "$cluster" --directory=$filesDirectory --name="$launchTime" --stdout="OAR/OAR.%jobid%.$launchTime.stdout" --stderr="OAR/OAR.%jobid%.$launchTime.stderr" --notify="[END,ERROR]exec:MMSD_-_Script_notification_0.2.sh $launchTime" "$commandToExecute"
else
    oarsub -l /nodes=$nodes/core=$core,walltime=$walltime -p "$cluster" --directory=$filesDirectory --name="$launchTime" --stdout="OAR/OAR.%jobid%.$launchTime.stdout" --stderr="OAR/OAR.%jobid%.$launchTime.stderr" --notify="[END,ERROR]exec:MMSD_-_Script_notification_0.2.sh $launchTime" "$commandToExecute"
fi

echo "Launch_time="$launchTime;
