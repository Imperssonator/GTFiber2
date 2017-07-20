function cell2csv_cols(fileName, cellArray, separator, excelYear, decimal)
% Writes cell array content into a *.csv file.
%
% CELL2CSV(fileName, cellArray, separator, excelYear, decimal)
%
% fileName     = Name of the file to save. [ i.e. 'text.csv' ]
% cellArray    = Name of the Cell Array where the data is in
% separator    = sign separating the values (default = ';')
% excelYear    = depending on the Excel version, the cells are put into
%                quotes before they are written to the file. The separator
%                is set to semicolon (;)
% decimal      = defines the decimal separator (default = '.')
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
% added the choice of decimal separator, 11/2010, S.Fiedler

% Checking für optional Variables
if ~exist('separator', 'var')
    separator = ',';
end

if ~exist('excelYear', 'var')
    excelYear = 1997;
end

if ~exist('decimal', 'var')
    decimal = '.';
end

% Setting separator for newer excelYears
if excelYear > 2000
    separator = ';';
end

% Write file
datei = fopen(fileName, 'w');

cellSizes = zeros(size(cellArray));
for i = 1:size(cellArray,1)
    for j = 1:size(cellArray,2)
        if isnumeric(cellArray{i,j})
            cellSizes(i,j) = length(cellArray{i,j});
        else
            cellSizes(i,j) = 1;
        end
    end
end

maxLens = max(cellSizes,[],1); % Each column now knows how wide it must be
% disp(maxLens)
totalLen = sum(maxLens,2);
% disp(totalLen)


for z=1:size(cellArray, 1)
    rowCount = 1;
    for s=1:size(cellArray, 2)
        
        var = eval(['cellArray{z,s}']);
%         disp(var)
        numEntries = 1;
        % If zero, then empty cell
        if size(var, 1) == 0
            var = '';
            numEntries = 1;
        end
        % If numeric -> String
        if isnumeric(var)
            numEntries = length(var);
            for ii = 1:length(var)
                varii = num2str(var(ii));
                disp(varii)
                % Conversion of decimal separator (4 Europe & South America)
                % http://commons.wikimedia.org/wiki/File:DecimalSeparator.svg
                if decimal ~= '.'
                    varii = strrep(varii, '.', decimal);
                end
%                 disp(varii)
                % If newer version of Excel -> Quotes 4 Strings
                if excelYear > 2000
                    varii = ['"' varii '"'];
                end
                
                % OUTPUT value
                fprintf(datei, '%s', varii);
                
                % OUTPUT separator
                if rowCount ~= totalLen
                    fprintf(datei, separator);
                    rowCount=rowCount+1;
                end
            end
        else
            % If logical -> 'true' or 'false'
            if islogical(var)
                if var == 1
                    var = 'TRUE';
                else
                    var = 'FALSE';
                end
                numEntries = 1;
            end
            % If newer version of Excel -> Quotes 4 Strings
            if excelYear > 2000
                var = ['"' var '"'];
            end
            
            % OUTPUT value
            fprintf(datei, '%s', var);
            
            % OUTPUT separator
            if rowCount ~= totalLen
                fprintf(datei, separator);
                rowCount = rowCount+1;
            end
        end
        for jj = 1:maxLens(s)-numEntries
            if rowCount ~= totalLen
                fprintf(datei,separator);
                rowCount = rowCount+1;
            end
        end
        
    end
    if z ~= size(cellArray, 1) % prevent a empty line at EOF
        % OUTPUT newline
        fprintf(datei, '\n');
    end
end
% Closing file
fclose(datei);
% END

end