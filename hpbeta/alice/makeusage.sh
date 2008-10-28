#! /bin/bash

yes 'cout << "' | head -`wc -l alice_usage.txt | awk '{print $1}'` > foo.txt
yes '" << endl;' | head -`wc -l alice_usage.txt | awk '{print $1}'` > foo2.txt
paste --delimiters=  foo.txt alice_usage.txt foo2.txt > alice_usage.h
