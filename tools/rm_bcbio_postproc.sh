find . -regex ".*${1}" | while read found_file;
 do
    echo "${found_file};
    rm -r ${found_file};
 done