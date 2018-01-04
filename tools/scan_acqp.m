function res=scan_acqp(strref,texte,numeric)
lmax=length(texte);
temp1=findstr(texte,strref);
temp1=temp1+length(strref)-1;
temp2=temp1;
while ( ~(strcmp(texte(temp2),'#')) & ~(strcmp(texte(temp2),'$')) & (temp2 < lmax))
    temp2=temp2+1;
end
if temp2 < lmax
    switch numeric
        case 1
            temp4=texte(temp1+1:temp2-1);
            openpar=findstr(temp4,'(');
            if ~isnan(openpar)
            closepar=findstr(temp4,')');
            temp3=strrep(temp4(openpar+1:closepar-1),',',' ');
            elsize=sscanf(temp3,'%f');
            else
                elsize=1;
                closepar=0;
                temp3=temp4;
            end
            res=sscanf(temp4(closepar+1:end),'%f');
            if numel(elsize)>1
                elsize=fliplr(elsize(:)');
                res=reshape(res,elsize)';
            end
%             temp3=temp4(openpar:closepar);
%             temp3=strrep(texte(temp1+1:temp2-1) , '(' , ' ');
%             temp3=strrep(temp3 , ')' , ' ');
%             res=sscanf(temp3,'%f');
%             if numel(res)>1 %On supprime le nombre d'éléments si celui-ci existe
%                 res=res(2:end);
%             end
        case 0
            temp4=texte(temp1+1:temp2-2);
            openpar=findstr(temp4,'(');
            if ~isnan(openpar)
                closepar=findstr(temp4,')')+2;
            else
                closepar=1;
            end
            res=temp4(closepar:end);
        case 2
            temp3=texte(temp1+1:temp2-1);
            temp1=1;
            while ~(strcmp(temp3(temp1),'<'))
                temp1=temp1+1;
            end
            res=temp3(temp1+1:end-2);
    end
else
    res=NaN; %Not a Number
end
return

