function req_ok = check_req(req)
%REQ_OK(v,req) check if the current version of MATLAB has the requirements
%specified in 'req'. These must be correctly named (case sensitive). 'req'
%is a cell array of string/character vectors
%
%Example: 
%req = {'Image Processing Toolbox','Bioinformatics Toolbox'};
%req_ok = check(req);
%'req_ok' is 1 if the current MATLAB installation has all the
%requirements, 0 otherwise
v = ver;
req_ok = zeros(1,length(req));
for j = 1:length(req)
    for i = 1:length(v)
        req_ok(j) = strcmp(v(i).Name,req(j)) || req_ok(j);
    end
end
req_ok = all(req_ok);
end