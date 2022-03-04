function write_songs(filename, song)

%write_songs(filename, song) 

fid = fopen(filename, 'w');

fwrite(fid, song, 'float32');

fclose(fid); 