#!/bin/bash
for j in 10 20 30 ; do
    condor_submit match.job -a batch_name=jet${j}_hijing_scaled
done