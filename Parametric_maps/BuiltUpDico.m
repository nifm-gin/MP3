function Dic = BuiltUpDico(DicoPrePath,DicoPostPath,tout)

% Built the ratio signal Post/Pre and resampled the dictionary to the
% time tout.
% Input: - DicoPrePath: path to the dictionary corresponding to the signal
%           before injection
%        - DicoPostPath: path to the dictionary corresponding to the signal
%        after injection
%        - tout: echo times (in s) for the resampling.
%
% Output: - Dic: Structure containing the Dicionary of the ratio (Dic.dat)
%        and the corresponding input values used to generate the data (Dic.par)
%
%  Nicolas Pannetier, Fev 28th 2013, UCSF

%% ComputeDico
% load
load(DicoPrePath,'par');
Dic.par     = par;
clear par
load(DicoPrePath,'Dico');
Pre = bsxfun(@rdivide,Dico,Dico(:,1));
clear Dico
load(DicoPostPath,'Dico');
Post = bsxfun(@rdivide,Dico,Dico(:,1));
clear Dico

load(DicoPrePath,'Seq');

Dic.dat = abs(Post)./abs(Pre);
tin = Seq.Tacq;
Dic.dat = spline(tin(:),Dic.dat,tout(:));
