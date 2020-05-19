#!/bin/bash
el=$(head -n2 RDFDAT | tail -n 1 | awk '{print $2+1}')
tail -n $(wc -l RDFDAT | awk '{print ($1-2)}') RDFDAT | split --lines=$el
paste x*| head -n1|awk '{printf "#           # ";for (i=1;i<=NF;i+=2) printf "%2s %2s        # ", $i,$(i+1);printf "\n"}' > rdfhead
paste x*| awk  '{printf "%s  ", $1 ;for (i=2;i<=NF;i+=2) printf "%s  ", $i; printf "\n" }' > rdf_temp.dat
cat rdfhead rdf_temp.dat > rdf_alla.dat 
sed "2d" rdf_alla.dat > rdf_all.dat
rm rdfhead rdf_temp.dat rdf_alla.dat
python3 /home/oadicks/bin/DL_tool_DDCT/qimaker/make_Dr.py
python3 /home/oadicks/bin/DL_tool_DDCT/qimaker/make_Qi.py

#python3 ~/bin/grtotmaker/make_rdf.py
#num_files=0
#for i in x*; do
#num_files+=1
#for i in x*; do
#a=$(head -n 1 $i)
#echo ${a/ /}

#done
#echo $num_files
