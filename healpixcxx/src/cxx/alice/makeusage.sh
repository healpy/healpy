#! /bin/bash

awk 'BEGIN {print "cout <<"} // {print"\""$0"\\n\""} END{print ";"}' alice_usage.txt > alice_usage.h
