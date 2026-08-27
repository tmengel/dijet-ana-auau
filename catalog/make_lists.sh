#/bin/bash

python3 catalog.py 
# for j in 10 20 30 ; do 
#     python3 make_list.py -d trees -s "jet${j}_hijing_unsclaed" -n 100
#     find "$(pwd)/jet${j}_hijing_unsclaed" -type f -name "*.list" | sort -u > jet${j}_hijing_unsclaed_master.list
# done
for j in 10 20 30 ; do 
    python3 make_list.py -d trees -s "jet${j}_hijing_scaled" -n 100 -t 08_24_2026
    find "$(pwd)/jet${j}_hijing_scaled" -type f -name "*.list" | sort -u > jet${j}_hijing_scaled_master.list
done

# python3 make_list.py -d trees -s "mb_hijing_unsclaed" -n 100
# find "$(pwd)/mb_hijing_unsclaed" -type f -name "*.list" | sort -u > mb_hijing_unsclaed_master.list
# python3 make_list.py -d trees -s "mb_hijing_scaled" -n 100 
# find "$(pwd)/mb_hijing_scaled" -type f -name "*.list" | sort -u > mb_hijing_scaled_master.list