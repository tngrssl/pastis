#!/bin/bash
for file in *.svg; do 
   inkscape "$file" --export-width 2048 --export-png "${file%svg}png"; 
done

