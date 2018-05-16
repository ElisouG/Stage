#!/bin/bash


date=`date`
messageCommit="$1 Mise a jour du : $date"
git add *
git commit -m "'$messageCommit'"
git push