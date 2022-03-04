function [songs] = read_songs(file_name)

fid = fopen(file_name, 'r');

if fid == -1
      error('Error: check that path has been set to find appropriate file.')
end

songs = fread(fid, 'float32');

fclose(fid);