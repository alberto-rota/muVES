function x4_bifurcation_QUARTO = importPtsFile(filename, dataLines)
    %IMPORTFILE Import data from a text file
    %  X4_BIFURCATION_QUARTO = IMPORTFILE(FILENAME) reads data from text
    %  file FILENAME for the default selection.  Returns the data as a cell
    %  array.
    %
    %  X4_BIFURCATION_QUARTO = IMPORTFILE(FILE, DATALINES) reads data for
    %  the specified row interval(s) of text file FILENAME. Specify
    %  DATALINES as a positive scalar integer or a N-by-2 array of positive
    %  scalar integers for dis-contiguous row intervals.
    %
    %  Example:
    %  x4_bifurcation_QUARTO = importfile("/MATLAB Drive/pts/4_bifurcation_QUARTO.pts", [1, Inf]);
    %
    %  See also READTABLE.
    %
    % Auto-generated by MATLAB on 21-May-2020 22:34:36
    
    %% Input handling
    
    % If dataLines is not specified, define defaults
    if nargin < 2
        dataLines = [1, Inf];
    end
    
    %% Setup the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 5);
    
    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = "\t";
    
    % Specify column names and types
    opts.PreserveVariableNames = true;
    opts.VariableNames = ["BEGIN_LIST", "VarName2", "VarName3", "VarName4", "VarName5"];
    opts.VariableTypes = ["char", "double", "double", "double", "char"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, ["BEGIN_LIST", "VarName5"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["BEGIN_LIST", "VarName5"], "EmptyFieldRule", "auto");
    
    % Import the data
    x4_bifurcation_QUARTO = readtable(filename, opts);
    
    %% Convert to output type
    x4_bifurcation_QUARTO = table2cell(x4_bifurcation_QUARTO);
    numIdx = cellfun(@(x) ~isnan(str2double(x)), x4_bifurcation_QUARTO);
    x4_bifurcation_QUARTO(numIdx) = cellfun(@(x) {str2double(x)}, x4_bifurcation_QUARTO(numIdx));
end