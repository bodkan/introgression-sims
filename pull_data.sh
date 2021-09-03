#!/usr/bin/env bash

rm -rf model/ results/

scp -r taco:/science/willerslev/users-shared/science-snm-willerslev-krd114/alba/introgression-sims/{model,results} .
