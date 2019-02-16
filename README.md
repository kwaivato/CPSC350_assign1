# CPSC350_assign1

1) Header
Vidal Arroyo
2253413
arroy118@mail.chapman.edu
CPSC 350-01
Assignment 1

2) Brief Overview
i. main.cpp: this file acts as the main function for the program.
ii. menu.cpp/menu.h: these files create a public interface that works with the user to enter multiple DNA-read files to be analyzed.
iii. logger.cpp/logger.h: these files create a class that can open DNA-read files and write both statistical parameters (along with new DNA reads) to a results file titled 'VidalArroyo.txt'.
iv: geneCalculator.cpp/geneCalculator.h: these files create a 'gene calculator' that can compute statistics concerning the lengths and composition of DNA reads from a DNA read file.
v: genePool.cpp/genePool.h: these files create a 'gene pool' that outputs DNA strings based on the statistical parameters of a previous DNA-read file.
vi: boxMuller.cpp/boxMuller.h: these files create a box-muller transform program, which generates a normal random variable D from a previously known mean and standard deviation.

3) Inline Comments
To run the programs in Docker using g++, type the following commands:

g++ *.cpp -o assign1
assign1 {fileName.txt}
