#!/usr/bin/env bash

split_into_chunks () {
	IN_DIR="$1"
	OUT_DIR="$2"
	LINES=$3

	split -l $LINES -d <(cat $IN_DIR/*) $OUT_DIR/chunk_

	for file in $OUT_DIR/*
	do
    mv "$file" "$file.txt"
	done

}

DATA_IN_DIR="$1"
DATA_OUT_DIR="$2"
ANNOTATION_IN_DIR="$3"
ANNOTATION_OUT_DIR="$4"
CHUNKSIZE=$5

split_into_chunks $DATA_IN_DIR $DATA_OUT_DIR $CHUNKSIZE
split_into_chunks $ANNOTATION_IN_DIR $ANNOTATION_OUT_DIR $CHUNKSIZE