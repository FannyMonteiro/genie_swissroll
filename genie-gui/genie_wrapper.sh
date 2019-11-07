#! /bin/bash

if [ -f ./genie_active ]; then rm ./genie_active; sleep 2; fi

./genie.exe > ./genie.out 2>&1 &

touch ./genie_active

while [ -f ./genie_active ]; do \
  if [ -f ./genie_pause ]; then kill -STOP $!; rm ./genie_pause; fi; \
  if [ -f ./genie_continue ]; then kill -CONT $!; rm ./genie_continue; fi; \
sleep 1; done

kill -HUP $!

