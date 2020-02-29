function      [data, result]= readtext(fname, delimiter, comment, quotes, options)

%  Usage:     [data, result]= readtext(fname, delimiter, comment, quotes, options)
% 
% Whatever text file you give it, readtext returns an array of the contents (or send me a 
%   bug report). Matlab can't read variable length lines or variable type values with the standard 
%   library. readtext can read any text file. Any string (or even regexp) can be delimiting, 
%   default is a comma. Everything after (and including) a comment character, until the line end, 
%   is ignored. Quote characters may also be given, everything between them is treated as one item. 
%   There are options to control what will be converted to numbers and how empty items are saved. 
% 
% If you find any errors, please let me know: peder at axensten dot se
% 
% fname:      the file to be read.
% 
% delimiter:  (default: ',') any string. May be a regexp, but this is a bit slow on large files. 
% 
% comment:    (default: '') zero or one character. Anything after (and including) this character, 
%   until the end of the line, will be ignored. 
% 
% quotes:     (default: '') zero, one (opening quote equals closing), or two characters (opening 
%   and closing quote) to be treated as paired braces. Everything between the quotes will be 
%   treated as one item. The quotes will remain. Quotes may be nested.
% 
% options:    (default: '') may contain (concatenate combined options): 
% - 'textual': no numeric conversion ('data' is a cell array of strings only), 
% - 'numeric': everything is converted to a number or NaN ('data' is a numeric array, empty items 
%   are converted to NaNs unless 'empty2zero' is given), 
% - 'empty2zero': an empty field is saved as zero, and 
% - 'empty2NaN': an empty field is saved as NaN. 
% - 'usewaitbar': call waitbar to report progress. If you find the wait bar annoying, get 'waitbar 
%   alternative' at http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=11398
% 
% data:       A cell array containing the read text, divided into cells by delimiter and line 
%   endings. 'data' will be empty if the file is not found, could not be opened, or is empty. 
%   With the option 'numeric', 'data' will be a numeric array, with 'textual', 'data' will be a 
%   cell array of strings only, and otherwise it will be a mixed cell array. For Matlab < version 7, 
%   returned strings may contain leading white-space.
% 
% result:     a structure:
% .min: minimum number of columns found in a line.
% .max: number of columns in 'data', before removing empty columns.
% .rows: number of rows in 'data', before removing empty rows. 
% .numberMask: true, if numeric conversion ('NaN' converted to NaN counts).
% .number: number of numeric conversions ('NaN' converted to NaN counts).
% .emptyMask: true, if empty item in file.
% .empty: number of empty items in file.
% .stringMask: true, if non-number and non-empty.
% .string: number of non-number, non-empty items.
% .quote: number of quotes. 
% 
% INSPIRATION: loadcell.m (id 1965). The goal of readtext is to be at least as flexible (you be 
%   the judge) and quicker. On my test file (see below), readtext is about 3--4 times 
%   as quick, maybe even more on large files. In readtext you may use a regexp as 
%   delimiter and it can ignore comments in the text file. 
% 
% SPEED:      Reading a 1MB file (150000 items!) with 'numeric' takes about 100 seconds on a 
%   fairly slow system. Time scales approximately linearly with input file size. 
% - Conversion from string to numeric is slow (I can't do anything about this), but using the 
%   option 'textual' is a lot quicker (above case takes 12 seconds).
% - Using a regexp delimiter is slower (during initializing), it adds 250 seconds! 
% 
% EXAMPLE:    [a,b]= readtext('txtfile', '[,\t]', '#', '"', 'numeric-empty2zero')
% This will load the file 'txtfile' into variable a, treating any of tab or comma as
%   delimiters. Everything from and including # to the next newline will be ignored. 
%   Everything between two double quotes will be treated as a string. Everything will 
%   be converted to numbers and a numeric array returned. Non-numeric items will become 
%   NaNs and empty items are converted to zero. 
% 
% Copyright (C) Peder Axensten (peder at axensten dot se), 2006.

% HISTORY:
% Version 1.0, 2006-05-03.
% Version 1.1, 2006-05-07:
% - Made 'options' case independent. 
% - Removed result.errmess -- now use error(errmess) instead. 
% - Removed result.nan -- it was equivalent to result.string, so check this instead.
% - Added some rows', 'result' fields: 'numberMask', 'emptyMask', and 'stringMask' 
%   (see 'result:', above).
% - A few small bug fixes.
% Version 1.2, 2006-06-06:
% - Now works in Matlab 6.5.1 (R13SP1) (maybe 6.5 too), versions <6.5 will NOT work.
% Version 1.3, 2006-06-20:
% - Better error report when file open fails. 
% - Somewhat quicker. 
% - Recommends 'waitbar alternative'. Ok with Matlab orig. waitbar too, of course. 
% Version 1.4, 2006-07-14:
% - Close waitbar instead of deleting it, and some other minor waitbar compatibility fixes. 
% Version 1.5, 2006-08-13:
% - No more (La)TeX formatting of file names. 
% - Prefixed waitbar messages with '(readtext)'. 
% Version 1.6, 2006-10-02:
% - Better removal of comments. Could leave an empty first row before. 
% - Added a 'usewaitbar' option. 
% - Now removes empty last columns and rows. 
% 
% TO DO:
% - Write a better text macro. 
% - Add result.quoteMask.
% - Add 'removeemptycolumns' and 'removeemptyrows' options. 

% KEYWORDS:     import, read, load, text, delimited, cell, numeric, array, flexible
% 
% REQUIREMENTS: Works in Matlab 6.5.1 (R13SP1) (probably 6.5 too), versions <6.5 will NOT work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Read (or set to default) the input arguments.
	if((nargin < 1) || ~ischar(fname) || isempty(fname))		% Is there a file name?
		error('First argument must be a file name!'); 
	end
	if(nargin < 2), delimiter=	',';				end			% Default delimiter value.
	if(nargin < 3), comment=	'';					end			% Default comment value.
	if(nargin < 4), quotes=		'';					end			% Default quotes value.
	if(nargin < 5), options=	[];					end			% Default options value.
	
	options=		lower(options);
	op_waitbar=		~isempty(strfind(options, 'usewaitbar'));	% Do waitbar calls. 
	op_numeric=		~isempty(strfind(options, 'numeric'));		% Read as numerical. 
	op_textual=		~isempty(strfind(options, 'textual')) && ~op_numeric;	% Read as textual. 
	op_empty=		[];											% Ignore empties, ...
	if(~isempty(strfind(options, 'empty2zero')))
		op_empty=		0;										% ... or replace by zero ...
	elseif(op_numeric || ~isempty(strfind(options, 'empty2nan')))
		op_empty=		NaN;									% ... or replace by NaN.
	end
	if(op_textual), op_empty= num2str(op_empty);	end			% Textual 'empty'.
	if(~ischar(comment) || (length(comment) > 1))
		error('Argument ''comment'' must be a string of maximum one character.');
	end
	if(~ischar(quotes) || (length(quotes) > 2))
		error('Argument ''quotes'' must be a string of maximum two characters.');
	end
	
	% Set the default return values.
	result.min=		Inf;
	result.max=		0;
	result.quote=	0;
	
	% Read the file.
	[fid, errmess]=	fopen(fname, 'r');							% Open the file.
	if(fid < 0), error(['Trying to open ' fname ': ' errmess]); end
	text=			fread(fid, 'uchar=>char')';					% Read the file.
	fclose(fid);												% Close the file.
	
	if(op_waitbar)
		th= waitbar(0, '(readtext) Initialising...');			% Show waitbar.
		thch=			findall(th, '-property', 'Interpreter');
		set(thch, 'Interpreter', 'none');						% No (La)TeX) formatting. 
	end
	
	% Clean up the text.
	eol=			char(10);
	text=			strrep(text, [char(13) char(10)], eol);		% Replace Windows-style eol.
	text=			strrep(text, char(13), eol);				% Replace MacClassic-style eol.
	if(~isempty(comment))										% Remove comments.
		text=	regexprep(text, ['^' comment '[^\n]*\n'], '');	% Remove commented lines. 
		text=	regexprep(text, [comment '[^\n]*'], '');		% Remove commented line endings. 
	end
	if(text(end) ~= eol), text= [text eol];				end		% End string with eol, if none.
	
	% Find column and row dividers.
	delimiter=		strrep(delimiter, '\t', char( 9));			% Convert to one char, quicker?
	delimiter=		strrep(delimiter, '\n', char(10));
	delimiter=		strrep(delimiter, '\r', char(10));
	delimiter=		strrep(delimiter, '\f', char(12));
	if(1 == length(delimiter))									% Find column dividers quickly.
		delimS=		find((text == delimiter) | (text == eol));
		delimE=		delimS;
	elseif(isempty(regexp(delimiter, '[\+\*\?\|\[^$<>]', 'once')))	% Find them rather quickly.
		delimS=		strfind(text, delimiter);
		eols=		find(text == eol);
		delimE=		union(eols, delimS + length(delimiter) - 1);
		delimS=		union(eols, delimS);
	else														% Find them with regexp.
		[delimS, delimE]=	regexp(text, [delimiter '|' eol]);
	end
	divRow=			[0, find(text == eol), length(text)];		% Find row dividers+last.
	
	% Keep quoted text together.
	if(~isempty(quotes))										% Should we look for quotes?
		if((length(quotes) == 1) || (quotes(1) == quotes(2)))	% Opening char == ending.
			exclE=			find(text == quotes(1));
			exclS=			exclE(1:2:end);
			exclE=			exclE(2:2:end);
		else													% Opening char ~= closing.
			exclS=			find(text == quotes(1));
			exclE=			find(text == quotes(2));
		end
		if((length(exclS) ~= length(exclE)) || (sum(exclS > exclE) > 0))
			if(op_waitbar), close(th); 	end						% Close waitbar or it'll linger.
			error('Opening and closing quotes don''t match in file %s.', fname); 
		end
		if(~isempty(exclS))										% We do have quoted text.
			if(op_waitbar), waitbar(0, th, '(readtext) Doing quotes...'); end	% Inform user.
			r=		1;
			rEnd=	length(exclS);
			n=		1;
			nEnd=	length(delimS);
			result.quote=	rEnd;
			while((n < nEnd) && (r < rEnd)) % "Remove" delimiters and newlines within quyotes.
				while((r <= rEnd) && (delimS(n) > exclE(r))), r= r+1;	end
				while((n <= nEnd) && (delimS(n) < exclS(r))), n= n+1;	end
				while((n <= nEnd) && (delimS(n) >= exclS(r)) && (delimS(n) <= exclE(r)))
					delimS(n)=	0;
					n=			n+1;
				end
				if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd); end	% Update waitbar.
			end
			if(op_waitbar), waitbar(1);	end;
			delimE=	delimE(delimS > 0);
			delimS=	delimS(delimS > 0);
		end
	end
	delimS=		delimS-1;										% Last char before delimiter.
	delimE=		[1 delimE(1:end-1)+1];							% First char after delimiter.
	
	% Do the stuff: convert text to cell (and maybe numeric) array.
	if(op_waitbar), waitbar(0, th, sprintf('(readtext) Reading ''%s''...', fname));	end
	r=				1;
	c=				1;											% Presize data to optimise speed.
	data=			cell(length(divRow), ceil(length(delimS)/(length(divRow)-1)));
	nums=			zeros(size(data));							% Presize nums to optimise speed.
	nEnd=			length(delimS);								% Prepare for a waitbar.
	for n=1:nEnd
		temp=			text(delimE(n):delimS(n));
		data{r, c}= 	temp;									% Textual item.
		if(~op_textual), nums(r, c)= str2double(temp);	end		% Quicker(!) AND better waitbar.
		if(text(delimS(n)+1) == eol)							% Next row.
			result.min=		min(result.min, c);					% Find shortest row.
			result.max=		max(result.max, c);					% Find longest row.
			r=				r+1;
			c=				0;
		end
		c=				c+1;
		if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd);	end	% Update waitbar.
	end
	if(op_waitbar), waitbar(1);	end
	
	% Clean up the conversion and do the result statistics.
	if(op_waitbar), waitbar(0, th, '(readtext) Cleaning up...');	end	% Inform user.
	data=				data(1:(r-1), 1:result.max);			% In case we started off to big.
	if(~op_textual), nums= nums(1:(r-1), 1:result.max);	end		% In case we started off to big.
	if(exist('strtrim', 'builtin')), data= strtrim(data);		% Not in Matlab 6.5...
	else							 data= deblank(data);		
	end
	while(all(cellfun('isempty', data(end, :))))				% Remove empty last lines. 
		data=	data(1:end-1, :); 
		nums=	nums(1:end-1, :); 
		r=		r-1;
	end 
	while(all(cellfun('isempty', data(:, end))))				% Remove empty last columns. 
		data=	data(:, 1:end-1); 
		nums=	nums(:, 1:end-1); 
		c=		c-1;
	end 
	result.rows=		r-1;
	empties=			cellfun('isempty', data);				% Find empty items.
	result.emptyMask=	empties;
	if(op_textual)
		result.numberMask=	repmat(false, size(data));			% No numbers, all strings.
		result.stringMask=	~empties;							% No numbers, all strings.
		data(empties)=		{op_empty};							% Set correct empty value.
	else
		result.numberMask=	~(isnan(nums) & ~strcmp(data, 'NaN'));	% What converted well.
		if(op_numeric)
			nums(empties)=		op_empty;						% Set correct empty value.
			data=				nums;							% Return the numeric array.
			result.stringMask=	~(empties | result.numberMask);	% Didn't convert well: so strs.
		else
			data(result.numberMask)= num2cell(nums(result.numberMask));	% Copy back numerics.
			data(empties)=		{op_empty};						% Set correct empty value.
			result.stringMask=	cellfun('isclass', data, 'char');	% Well, the strings.
		end
	end
	result.empty=		sum(result.emptyMask(:));				% Count empties.
	result.numberMask=	result.numberMask & ~result.emptyMask;	% Empties don't count.
	result.number=		sum(result.numberMask(:));				% Count numbers.
	result.stringMask=	result.stringMask & ~result.emptyMask;	% Empties don't count.
	result.string=		sum(result.stringMask(:));				% Count strings.
	
 	if(op_waitbar), close(th);	end								% Removing the waitbar. 
end

