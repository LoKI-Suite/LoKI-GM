function LXCat2LoKI(LXCatFileName)
%LXCat2LoKI function to help in the process of adapting any LXCat file to the LoKI code 
%   This function help the end user in the process of adapting any LXCat txt file to the LoKI code. In particular, the
%   function creates a new LXCat file that includes in the first comment line of each electron impact cross section the
%   proper description (in LoKI format) of the reactants and products of the collision. 

  % definition of regular expressions to parse LXCat file
  LXCatRegExp1 = 'PROCESS: (?<reactants>.+?)(?<direction>->|<->) (?<products>.+?), (?<type>\w+)';

  % open LXCat file
  fileIDin = fopen(['Input' filesep LXCatFileName], 'r');
  if fileIDin<0
    error(' Unable to open LXCat file: %s\n', ['Input' filesep LXCatFileName]);
  end
  
  % create output file
  LXCatModFileName = [LXCatFileName(1:end-4) '_modified.txt' ];
  fileIDout = fopen(['Input' filesep LXCatModFileName], 'wt');
  if fileIDout<0
    error(' Unable to open modified LXCat file: %s\n', LXCatModFileName);
  end
  
  % parse LXCat input file and create modified LXCat output file
  while ~feof(fileIDin)
    % read line from original LXCat file
    line = fgetl(fileIDin);
    % copy line from original to modified LXCat file
    fprintf(fileIDout, '%s\n', line);
    % analise if a cross section has been found
    description = regexp(line, LXCatRegExp1, 'names', 'once');
    if ~isempty(description)
      % ask the user for the correct description of the collision
      userDescrition = inputdlg(line);
      % read and copy "PARAMETER:" line from original LXCat file into the modified one
      line = fgetl(fileIDin);
      fprintf(fileIDout, '%s\n', line);
      % write comment line with user prescribed description 
      fprintf(fileIDout, 'COMMENT: [%s, %s]\n', userDescrition{1}, description.type);
    end
  end
  
  % close LXCat files, booth original and modified (output) ones
  fclose(fileIDin);
  fclose(fileIDout);

end