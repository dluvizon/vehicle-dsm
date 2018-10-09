% KLT_READ_FEATURE_TABLE  Read a feature table written by KLTWriteFeatureTable.
% 
%    [X, Y, VAL] = KLT_READ_FEATURE_TABLE(FILENAME) reads the feature table
%    in FILENAME into the two-dimensional arrays X, Y, and VAL.  The
%    dimensions of X, Y, and VAL are M by N, where M is the number of
%    features and N is the number of image frames.  
%
%    The file may be ASCII or binary.
function [x, y, val] = klt_read_featuretable(filename)

fp = fopen(filename);
if (fp == -1)  error('Unable to open file'), end
a = fread(fp, 6, 'uchar');
fclose(fp);

if all(a' == 'KLTFT1') == 1, 
    [x, y, val] = klt_read_featuretable_binary(filename);
else
    [x, y, val] = klt_read_featuretable_ascii(filename);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binary version
function [x, y, val] = klt_read_featuretable_binary(filename)

fp = fopen(filename);
if (fp == -1)  error('Unable to open file'), end
a = fread(fp, 6, 'uchar');
if all(a' == 'KLTFT1') ~= 1, warning('File is not a binary KLT feature table'), end
a = fread(fp, 2, 'int');
nframes = a(1);
nfeatures = a(2);
x = zeros(nfeatures, nframes);
y = zeros(nfeatures, nframes);
val = zeros(nfeatures, nframes);
for j = 1 : nfeatures,
    for i = 1 : nframes,
        x(j, i) = fread(fp, 1, 'float32');
        y(j, i) = fread(fp, 1, 'float32');
        val(j, i) = fread(fp, 1, 'int');
    end
end
fclose(fp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASCII version

function [x, y, val] = klt_read_featuretable_ascii(filename)

str1 = '!!! Warning:  This is a KLT data file.  Do not modify below this line !!!';
preamble = 1;  % whether in preamble of file (in which comments are allowed)
counter = 0;   % current line number of file
ln_dim = 5;    % line number of file in which dimensions are specified
ln_data = 10;  % first line number of file containing data

fp = fopen(filename);
if (fp == -1)  error('Unable to open file'), end
while (~feof(fp)),
    tline = fgets(fp);
    if preamble & strncmp(tline, str1, size(str1,2)) == 1,
        preamble = 0;
    elseif ~preamble,
        if counter == ln_dim,
            a = sscanf(tline, 'nFrames = %d, nFeatures = %d');
            nframes = a(1);
            nfeatures = a(2);
            x = zeros(nfeatures, nframes);
            y = zeros(nfeatures, nframes);
            val = zeros(nfeatures, nframes);
            
        elseif counter >= ln_data,
        
            % parse line
            feat_number = sscanf(tline, '%d |') + 1;
            if feat_number ~= counter - ln_data + 1, warning('Unexpected feature number'), end
            index = strfind(tline, '|');
            tline = tline( index+2 : end );
            a = sscanf(tline, '(%f,%f)=%d ');
            if size(a,1) ~= nframes*3 | size(a,2) ~= 1, warning('Unexpected number of values in row'), end
            a = reshape(a, 3, nframes);
            x(feat_number, :)   = a(1, :);
            y(feat_number, :)   = a(2, :);
            val(feat_number, :) = a(3, :);
%            disp(tline)
        end
        
        counter = counter + 1;
    end
end
if counter - ln_data ~= nfeatures, warning('Unexpected number of rows'), end
fclose(fp);