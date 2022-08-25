function [fW,fL,zone]=filename_reader(files,i,c)
%This function outputs the values for transistors size contained in the filename
%It requires to be called in a loop in order to adress each file individually
%Each file must conatin in order:  Width and then Length.
switch c
    case 'out' % output curves
        if length(files(i).name) < 55
            % getting transistors sizes registered in the filename
            auxread = textscan(files(i).name,'%s%s%s%s%s*');
            fWL = split(auxread{1,5},'_');
            zone = auxread{1,4};
            fW = str2double(fWL{1,1}); %filename width
            fL = str2double(fWL{2,1}); %filename length
        else
            auxread = textscan(files(i).name,'%s%s[%s;%s%s%s');
            fWL = split(auxread{1,3},'_');
            fW = str2double(fWL{1,1}); %filename width
            fL = str2double(fWL{2,1}); %filename length
            zone = fWL{5,1}; % filename zone
        end
    case 'lin' %linear curves
        auxread = textscan(files(i).name,'%s%s[%s');
        fWL = split(auxread{1,3},'_');
        fW = str2double(fWL{1,1}); %filename width
        fL = str2double(fWL{2,1}); %filename length
         z = split(fWL{5,1},'(');
         zone = z{1,1};
    case 'sat' %saturation curves
        auxread = textscan(files(i).name,'%s%s[%s');
        fWL = split(auxread{1,3},'_');
        fW = str2double(fWL{1,1}); %filename width
        fL = str2double(fWL{2,1}); %filename length
        z = split(fWL{5,1},'(');
        zone = z{1,1};
    case 'cov'
        auxread = textscan(files(i).name,'%s[%s;%s%s%s');
        fWL = split(auxread{1,2},'_');
        fW = str2double(fWL{1,1});
        fL = str2double(fWL{2,1});
        zone = fWL{5,1};
    case 'cox'
        auxread = textscan(files(i).name,'%s');
        aux = split(auxread{1,1},'_');
        fW = aux{2,1}
%       aux = split(auxread{1,1},'_');
%       fL = aux(2,1);
        
    
end
end
