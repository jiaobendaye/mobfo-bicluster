function resid = residu(data)
%Find residues         
  trandata = data';  
  vecdata = trandata(:);
  bicRow = ones(size(data,1),1);
  bicCol = ones(size(data,2),1);
  resid = mprintres(size(data,1), size(data,2), vecdata, bicRow, bicCol);
end
