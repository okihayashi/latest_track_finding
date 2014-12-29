#!/bin/sh

for i in {4..10}
do
    echo "---This is $i th run---"
    ./TrackFinding $i
done

