#!/bin/bash
set -eu
for i in gd10 yn56;do
	Rscript dge.R ${i}
done
