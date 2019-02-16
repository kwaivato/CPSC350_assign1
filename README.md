# CPSC350_assign1

(Header) Name: Vidal Arroyo, ID: 2253413, Email: arroy118@mail.chapman.edu, Class: CPSC 350-01, Assignment: Assignment 1

(Brief Overview) i. main.cpp: this file acts as the main function for the program; ii. menu.cpp/menu.h: these files create a public interface that works with the user to enter multiple DNA-read files to be analyzed; iii. logger.cpp/logger.h: these files create a class that can open DNA-read files and write both statistical parameters (along with new DNA reads) to a results file titled 'VidalArroyo.txt'; iv: geneCalculator.cpp/geneCalculator.h: these files create a 'gene calculator' that can compute statistics concerning the lengths and composition of DNA reads from a DNA read file; v: genePool.cpp/genePool.h: these files create a 'gene pool' that outputs DNA strings based on the statistical parameters of a previous DNA-read file; vi: boxMuller.cpp/boxMuller.h: these files create a box-muller transform program, which generates a normal random variable D from a previously known mean and standard deviation.

(Inline Comments) To run the programs in Docker using g++, type the following commands: 1) g++ *.cpp -o assign1 2) ./assign1 {fileName.txt}

(Other Notes) I used inspiration from the following sources for the following topics: 1) Converting font to all lower/upper case - http://www.cplusplus.com/forum/general/28440/ ; 2) IO Stream - example on pg. 502 of the Absolute C++ Book (5th Edition); 3) Use of the pi constant - https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c ; 4) File Appending - https://stackoverflow.com/questions/2393345/how-to-append-text-to-a-text-file-in-c ; 5) Checking if a file exists - https://stackoverflow.com/questions/25225948/how-to-check-if-a-file-exists-in-c-with-fstreamopen?rq=1
