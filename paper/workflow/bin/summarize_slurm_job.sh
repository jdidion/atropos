#!/bin/bash
jobId=`cat $1`
sacct -P -o ALL -j $jobId
