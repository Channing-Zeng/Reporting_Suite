#!/bin/sh
while read -r -a cols; do
    (
		if [[ ${cols} == \#\#* ]] ;
		then
			echo "${cols[*]}"
		else  
	        IFS=$'\t'
		    echo "${cols[*]}"
		fi
    )
done <$1
