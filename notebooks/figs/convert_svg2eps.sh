#!/bin/bash
for file in *.svg; do 
   inkscape "$file" -E "${file%svg}eps"; 
done

