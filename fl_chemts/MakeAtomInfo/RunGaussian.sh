#!/bin/bash

i=0
thread=$(nproc --all)
for com in work/*.com; do
  echo $com
  (g16 $com) &
  i=$((i+1))
  if [[ $i -ge $thread ]];then
    wait
    i=0
  fi
done
