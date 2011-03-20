function err = skipline(fid)

size = fread(fid,1,'integer*4');
fread(fid,size,'char*1');
size2 = fread(fid,1,'integer*4');
err = size-size2;