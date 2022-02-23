#!/bin/usr/env bash

array=($(seq 10 10 40))

echo "array is ${array[*]}"

echo "array size is ${#array[@]}"

for ((i=0 ; i<${#array[@]} ; i++))
do
  echo "Number: ${array[$i]}"
done
