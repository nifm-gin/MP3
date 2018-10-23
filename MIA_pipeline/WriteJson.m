function WriteJson(Json,jsonfile)
%WriteJson Write the data contained in Json in the jsonfile file.
JMod = jsonencode(Json);
fidmod = fopen(jsonfile, 'w');
fwrite(fidmod, JMod, 'uint8');
fclose(fidmod);
end

