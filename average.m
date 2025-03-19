function average()
    % Step 1: Generate the complete date range
    % Initialize arrays to hold all dates
    daysIn2020 = 366 - 129 + 1;  % 2020 was a leap year
    daysIn2021 = 365;
    daysIn2022 = 365;
    totalDays = daysIn2020 + daysIn2021 + daysIn2022;
    years = zeros(totalDays, 1);
    doys = zeros(totalDays, 1);

    % Fill in 2020 dates
    idx = 1;
    for doy = 129:366
        years(idx) = 2020;
        doys(idx) = doy;
        idx = idx + 1;
    end

    % Fill in 2021 dates
    for doy = 1:365
        years(idx) = 2021;
        doys(idx) = doy;
        idx = idx + 1;
    end

    % Fill in 2022 dates
    for doy = 1:365
        years(idx) = 2022;
        doys(idx) = doy;
        idx = idx + 1;
    end

    % Create a matrix with all dates
    allDates = [years, doys];
    numDays = size(allDates, 1);

    % Step 2: Define output file name first to exclude it from processing
    dataFolder = 'C:\Users\Lenovo 15\Desktop\METDATANASApower';
    outputFileName = 'nepal_daily_averages_complete.csv';
    outputFilePath = fullfile(dataFolder, outputFileName);
    
    % Get list of all files in the directory
    allFiles = dir(fullfile(dataFolder, '*.*'));
    
    % Filter to only include data files and exclude the output file
    % Filter files to include only those with specified extensions
% Preallocate fileList to improve performance
fileList = repmat(struct('name', '', 'folder', '', 'date', '', 'bytes', 0, 'isdir', false, 'datenum', 0), 1, length(allFiles));
validFileCount = 0; % Counter for valid files

for i = 1:length(allFiles)
    % Skip directories, current directory (.), parent directory (..), and files with unwanted extensions
    if allFiles(i).isdir || ...
       strcmpi(allFiles(i).name, '.') || ...
       strcmpi(allFiles(i).name, '..') || ...
       ~any(strcmpi(getFileExtension(allFiles(i).name), {'csv', 'xls', 'xlsx'}))
        continue; % Skip this file if it doesn't meet the criteria
    end
    
    % Append the valid file to the fileList array
    validFileCount = validFileCount + 1; % Increment valid file counter
    fileList(validFileCount) = allFiles(i); % Assign the valid file to the preallocated array
end

% Trim fileList to the actual number of valid files
fileList = fileList(1:validFileCount);

% Check if any valid files were found
if isempty(fileList)
    error('No data files found in the specified folder. Please ensure the files have .csv, .xls, or .xlsx extensions.');
end    
fprintf('Found %d data files to process.\n', length(fileList));

    % Step 3: Initialize arrays to store daily averages
    avgT2M = nan(numDays, 1);
    avgRH2M = nan(numDays, 1);
    avgWS2M = nan(numDays, 1);
    avgQV2M = nan(numDays, 1); % In case QV2M is present
    validCounts = zeros(numDays, 1);  % Track number of valid readings for each day
    
    % Arrays to store sums for each day
    t2mSums = zeros(numDays, 1);
    rh2mSums = zeros(numDays, 1);
    ws2mSums = zeros(numDays, 1);
    qv2mSums = zeros(numDays, 1);
    
    % Step 4: Process each file
    for fileIdx = 1:length(fileList)
        % Show progress
        fprintf('Processing file %d of %d: %s\n', fileIdx, length(fileList), fileList(fileIdx).name);
        
        % Construct full file path
        filePath = fullfile(dataFolder, fileList(fileIdx).name);
        
        try
            % Try to read the file
            try
                % First attempt - CSV format
                data = readtable(filePath, 'FileType', 'text', 'Delimiter', ',');
            catch
                try
                    % Second attempt - Excel format
                    data = readtable(filePath);
                catch
                    warning('Failed to read file %s as CSV or Excel. Skipping.', fileList(fileIdx).name);
                    continue;
                end
            end
            
            % Get the variable names and normalize them to uppercase
            varNames = data.Properties.VariableNames;
            upperVarNames = upper(varNames);
            
            % Check if required columns exist
            requiredCols = {'YEAR', 'DOY', 'T2M', 'RH2M', 'WS2M'};
            colIndices = zeros(1, length(requiredCols));
            
            % Find indices of required columns
            for i = 1:length(requiredCols)
                colIdx = find(strcmp(requiredCols{i}, upperVarNames));
                if isempty(colIdx)
                    warning('File %s is missing required column %s. Skipping file.', ...
                            fileList(fileIdx).name, requiredCols{i});
                    colIndices = [];
                    break;
                end
                colIndices(i) = colIdx;
            end
            
            % Skip this file if any required column is missing
            if isempty(colIndices)
                continue;
            end
            
            % Find index of QV2M column if it exists
            qvColIdx = find(strcmp('QV2M', upperVarNames));
            hasQV = ~isempty(qvColIdx);
            
            % Extract and convert data to numeric if needed
            yearCol = varNames{colIndices(1)};
            doyCol = varNames{colIndices(2)};
            t2mCol = varNames{colIndices(3)};
            rh2mCol = varNames{colIndices(4)};
            ws2mCol = varNames{colIndices(5)};
            
            if hasQV
                qv2mCol = varNames{qvColIdx};
            end
            
            % Ensure numeric data
            yearData = ensureNumeric(data.(yearCol));
            doyData = ensureNumeric(data.(doyCol));
            t2mData = ensureNumeric(data.(t2mCol));
            rh2mData = ensureNumeric(data.(rh2mCol));
            ws2mData = ensureNumeric(data.(ws2mCol));
            
            if hasQV
                qv2mData = ensureNumeric(data.(qv2mCol));
            end
            
            % Process each day's data
            for dayIdx = 1:numDays
                currentYear = allDates(dayIdx, 1);
                currentDOY = allDates(dayIdx, 2);
                
                % Find rows for this day
                dayRows = (yearData == currentYear) & (doyData == currentDOY);
                
                if any(dayRows)
                    % Calculate means for this location and day
                    t2mVal = mean(t2mData(dayRows), 'omitnan');
                    rh2mVal = mean(rh2mData(dayRows), 'omitnan');
                    ws2mVal = mean(ws2mData(dayRows), 'omitnan');
                    
                    % Add to sums if values are valid
                    if ~isnan(t2mVal) && ~isnan(rh2mVal) && ~isnan(ws2mVal)
                        t2mSums(dayIdx) = t2mSums(dayIdx) + t2mVal;
                        rh2mSums(dayIdx) = rh2mSums(dayIdx) + rh2mVal;
                        ws2mSums(dayIdx) = ws2mSums(dayIdx) + ws2mVal;
                        
                        if hasQV
                            qv2mVal = mean(qv2mData(dayRows), 'omitnan');
                            if ~isnan(qv2mVal)
                                qv2mSums(dayIdx) = qv2mSums(dayIdx) + qv2mVal;
                            end
                        end
                        
                        validCounts(dayIdx) = validCounts(dayIdx) + 1;
                    end
                end
            end
            
        catch e
            warning('Error processing file %s: %s', fileList(fileIdx).name, e.message);
            continue;
        end
    end
    
    % Step 5: Calculate averages from the sums
    for dayIdx = 1:numDays
        if validCounts(dayIdx) > 0
            avgT2M(dayIdx) = t2mSums(dayIdx) / validCounts(dayIdx);
            avgRH2M(dayIdx) = rh2mSums(dayIdx) / validCounts(dayIdx);
            avgWS2M(dayIdx) = ws2mSums(dayIdx) / validCounts(dayIdx);
            
            if any(qv2mSums > 0)  % If any QV2M data was found
                avgQV2M(dayIdx) = qv2mSums(dayIdx) / validCounts(dayIdx);
            end
        end
    end

    % Step 6: Create output table
    if any(qv2mSums > 0)  % If QV2M data was found
        results = array2table([allDates, avgT2M, avgRH2M, avgQV2M, avgWS2M, validCounts], ...
            'VariableNames', {'YEAR', 'DOY', 'avgT2M', 'avgRH2M', 'avgQV2M', 'avgWS2M', 'LocationCount'});
    else
        results = array2table([allDates, avgT2M, avgRH2M, avgWS2M, validCounts], ...
            'VariableNames', {'YEAR', 'DOY', 'avgT2M', 'avgRH2M', 'avgWS2M', 'LocationCount'});
    end

    % Step 7: Export results
    writetable(results, outputFilePath);
    fprintf('Results successfully exported to %s\n', outputFilePath);

    % Display the first few rows to verify
    fprintf('\nFirst few rows of the results:\n');
    disp(head(results));

    % Display summary statistics
    fprintf('\nSummary Statistics:\n');
    fprintf('Total days analyzed: %d\n', numDays);
    fprintf('Days with valid data: %d\n', sum(validCounts > 0));
    
    if sum(validCounts > 0) > 0
        fprintf('Average number of locations per day: %.1f\n', mean(validCounts(validCounts > 0)));
        fprintf('Maximum locations for a single day: %d\n', max(validCounts));
        fprintf('Minimum locations for a single day with data: %d\n', min(validCounts(validCounts > 0)));
    end
end

% Helper function to get file extension
function ext = getFileExtension(filename)
    [~, ~, ext] = fileparts(filename);
    if ~isempty(ext) && ext(1) == '.'
        ext = ext(2:end);  % Remove the leading dot
    end
end

% Helper function to ensure data is numeric
function numData = ensureNumeric(data)
    if isnumeric(data)
        numData = data;
    else
        % Try to convert to numeric, replacing non-convertible values with NaN
        numData = str2double(data);
        
        % Handle cell arrays of strings
        if all(isnan(numData)) && iscell(data)
            numTemp = zeros(size(data));
            for i = 1:numel(data)
                if isnumeric(data{i})
                    numTemp(i) = data{i};
                elseif ischar(data{i}) || isstring(data{i})
                    val = str2double(data{i});
                    if isnan(val)
                        numTemp(i) = NaN;
                    else
                        numTemp(i) = val;
                    end
                else
                    numTemp(i) = NaN;
                end
            end
            numData = numTemp;
        end
    end
end