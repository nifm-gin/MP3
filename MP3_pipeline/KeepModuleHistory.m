function J = KeepModuleHistory(J, Struct, Name_Mod)

if ~isfield(J, 'Bricks')
    J.Bricks = struct();
else
    FieldName = fieldnames(J.Bricks);
    Ext = '';
    Ind = 2;
    while any(strcmp([Name_Mod, Ext], FieldName))
        Ext = ['_', num2str(Ind)];
        Ind = Ind+1;
    end
    Name_Mod = [Name_Mod, Ext];
end
J.Bricks.(Name_Mod) = Struct;

end