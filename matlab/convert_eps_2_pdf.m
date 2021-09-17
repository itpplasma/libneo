% Convert an eps figure to pdf using epstopdf.
%
% input:
% ------
% pr_name: string, name+path of the eps file.
%
% output:
% -------
% none
%
% side effects:
% -------------
% Creates a new file.
% Path variable might be changed during execution, but is reset at the end.
function convert_eps_2_pdf(pr_name)
  path = getenv('PATH');
  if strfind(computer,'MAC')
    setenv('PATH',[path,':/usr/local/bin:/usr/texbin']);
    system(['epstopdf ',pr_name]);
  else
    system(['/usr/bin/epstopdf ',pr_name]);
  end
  setenv('PATH',path)
end
