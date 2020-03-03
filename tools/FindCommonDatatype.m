function CommonDatatype = FindCommonDatatype(Types)
%FINDCOMMONDATATYPE Find the highest Datatype among all the ones in a cell
%vector
%   The function does : 
%     if any(double) 
%            return double
%     elseif any(single) 
%            return single
%     elseif only int 
%            return int(max)
%     elseif only uint 
%            return uint(max)
%     elseif uint and int 
%             if max <= 32
%                 return int(max*2)
%             else
%                 return single
%             end
%     else
%             warning
%             return double
%     end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Documentation about data types %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double : 64-bit float : 1 bit for the sign, 11 for the exponent (included 1 for its sign), 52 for the fraction
%         The range for double is:
%             -1.79769e+308 (=(2-2^52)*2^1023) to -2.22507e-308 and
%             2.22507e-308 to  1.79769e+308
%         Type eps(x) to display the positive distance from abs(x) to the next larger in
%             magnitude floating point number of the same precision as x 
%             
%             
% single : 32-bit float : 1 bit for the sign, 8 for the exponent (included 1 for its sign), 23 for the fraction
%         The range for single is:
%             -3.40282e+38 (=(2-2^23)*2^127) to -1.17549e-38 and
%             1.17549e-38 to  3.40282e+38
%         Type eps(x) to display the positive distance from abs(x) to the next larger in
%             magnitude floating point number of the same precision as x
%             
%             
% uint8 :  unsigned integer in the range 0 <--> 255
% uint16 : unsigned integer in the range 0 <--> 65 535
% uint32 : unsigned integer in the range 0 <--> 4 294 967 295
% uint64 : unsigned integer in the range 0 <--> 18 446 744 073 709 551 615
% 
% 
% int8 :  signed integer in the range -128 <--> 127
% int16 : signed integer in the range -32 768 <--> 32 767
% int32 : signed integer in the range -2 147 483 648 <--> 2 147 483 647
% int64 : signed integer in the range -9 223 372 036 854 775 808 <--> 9 223 372 036 854 775 807


if any(strcmp(Types, 'double'))
    CommonDatatype = 'double';
elseif any(strcmp(Types, 'single'))
    CommonDatatype = 'single';
elseif all(contains(Types, 'int') & ~contains(Types, 'uint')) % if there is only int
    BitNumber = zeros(1,length(Types));
    for i=1:length(Types)
        Dt = Types{i};
       BitNumber(i) = str2double(Dt(4:end));
    end
    CommonDatatype = ['int',num2str(max(BitNumber))];
elseif all(contains(Types, 'uint')) % if there is only uint
    BitNumber = zeros(1,length(Types));
    for i=1:length(Types)
        Dt = Types{i};
       BitNumber(i) = str2double(Dt(5:end));
    end
    CommonDatatype = ['int',num2str(max(BitNumber))];
elseif any(contains(Types, 'int') & ~contains(Types, 'uint')) && any(contains(Types, 'uint')) % if there is a mix of int and uint
    BitNumber = zeros(1,length(Types));
    IsInt = contains(Types, 'int') & ~contains(Types, 'uint');
    %IsUint = contains(Types, 'uint');
    for i=1:length(Types)
        Dt = Types{i};
        if IsInt(i)
            BitNumber(i) = str2double(Dt(4:end));
        else
            BitNumber(i) = str2double(Dt(5:end));
        end
    end
    if max(BitNumber) <= 32
        CommonDatatype = ['int',num2str(max(BitNumber*2))];
    else
        CommonDatatype = 'single';
    end
else
    warning('Datatype non supported, ''double'' returned');
    CommonDatatype = 'double';
end

end
