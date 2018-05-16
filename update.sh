#!/bin/bash
git pull;
git add *;
message = date ;
git commit -m 'Mise a jour du $message';
git push;
