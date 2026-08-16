#/bin/bash

python3 catalog.py 
# python3 make_list.py -d overlay -s hijing_sub1
# for js in 30 ; do
#     for sub in 0 1 ; do
#         for jetR in 3 4 ; do
#             jetSample="jet${js}_r0${jetR}_sub${sub}"
#             rm -f ${jetSample}/*.list
#             echo "Making list for ${jetSample}"
#             python3 make_list.py -d overlay -s ${jetSample}
#             find "$(pwd)/${jetSample}" -name "*.list" | sort -u > "${jetSample}.list"
#         done
#     done
# done