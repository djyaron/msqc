function [frag, dir] = tmpfrag(cfg, tplpath)
  dir = tempname;
  if (mkdir(dir) ~= 1), error('Could not create temporary directory.'), end;
  copyfile(tplpath, [dir filesep cfg.template '.tpl']);
  frag = Fragment(dir, cfg);
end