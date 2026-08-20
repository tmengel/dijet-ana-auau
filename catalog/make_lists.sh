#/bin/bash

python3 catalog.py 
for j in 10 20 30 ; do 
    python3 make_list.py -d trees -s "jet${j}_hijing_scaled" -n 100
    find "$(pwd)/jet${j}_hijing_scaled" -type f -name "*.list" | sort -u > jet${j}_hijing_scaled_master.list
done