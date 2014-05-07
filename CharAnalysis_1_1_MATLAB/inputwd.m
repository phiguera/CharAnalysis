function r = inputwd(prompt, default)

% input prompt with default value
% Usage: inputwd('prompt', defaultvalue)

if nargin ~= 2, 
   error('Usage: inputwd(prompt, defaultvalue)');
end

rep = input(char(prompt));

if isempty(rep),
   r = default;
else
   r = rep;

end

