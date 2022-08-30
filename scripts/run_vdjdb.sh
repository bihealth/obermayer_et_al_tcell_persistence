#!/bin/bash

vdjmatch match -S human -R TRA ../data/VDJdb/input/*TRA.txt ../data/VDJdb/TRA/ 
vdjmatch match -S human -R TRB ../data/VDJdb/input/*TRB.txt ../data/VDJdb/TRB/ 
