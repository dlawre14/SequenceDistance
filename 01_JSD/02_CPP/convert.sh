awk '{key=$0; getline; print $0 "\t" substr(key,2);}' $1
