function ztype = atypeToZtype(atype)
if (atype > 99)
   ztype = round(atype/100);
else
   ztype = atype;
end
end

