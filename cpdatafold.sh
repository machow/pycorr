#!/bin/bash
FILES="$@"
for f in $(cat "$@")
do
	cp -r /mnt/cd/hasson/janice/Pieman/piesky/subjects/$f pieNDiv/subjects
done
