#!/bin/bash

pids=$(ps -u $USER -o pid -o comm | grep R | awk '{print $1}')
for pid in $pids
do
  kill -9 $pid
done

fs=(progress stderr stdout)
for ((i=0; i<${#fs[@]}; i++));
do
  if [ -f ${fs[i]}.txt ]
  then
    rm ${fs[i]}.txt 
  fi
done

ds=(fig rds JournalFigures)
for ((i=0; i<${#ds[@]}; i++));
do
  rm -rf ../${ds[i]}
done

ds=(slides manuscript)
es=(aux bbl blg log nav out pdf snm synctex.gz toc vrb)

for ((i=0; i<${#ds[@]}; i++));
do
  for((j=0; j<${#es[@]}; ++j));
  do
    if [ -f ../${ds[i]}/*.${es[j]} ]
    then
      rm ../${ds[i]}/*.${es[j]}
    fi
  done
done