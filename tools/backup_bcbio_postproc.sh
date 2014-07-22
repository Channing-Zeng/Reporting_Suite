find . -regex ".*${1}" | while read found_file;
 do 
    new_fpath=${found_file}_bk;
    echo "${found_file} -> ${new_fpath}";
    mv ${found_file} ${new_fpath};
 done