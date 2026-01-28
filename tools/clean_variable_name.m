function var_name = clean_variable_name(var_name, option)


var_name= strrep(var_name,'-','');
var_name= strrep(var_name,':','');
var_name= strrep(var_name,'/','_');
var_name= strrep(var_name,' ','_');
var_name= strrep(var_name,'.','');
var_name= strrep(var_name,'+','');
var_name= strrep(var_name,'*','star');

switch option
    case 1 % also remove numerical value
        var_name=strrep(var_name,'0','zero');
        var_name=strrep(var_name,'1','one');
        var_name=strrep(var_name,'2','two');
        var_name=strrep(var_name,'3','three');
        var_name=strrep(var_name,'4','four');
        var_name=strrep(var_name,'5','five');
        var_name=strrep(var_name,'6','six');
        var_name=strrep(var_name,'7','seven');
        var_name=strrep(var_name,'8','eight');
        var_name=strrep(var_name,'9','ten');
    otherwise
        
end

