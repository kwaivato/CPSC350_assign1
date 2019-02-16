// menu.h
// A declared menu interface used to facilitate a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "logger.h"

using namespace std;

// Allows a user to continuously enter DNA-read files and have their results outputted to a results file. Must use
// with the logger, geneCalculator, genePool, and boxMuller classes.
class Menu
{
  public:
    Menu();
    ~Menu();
    // The Run method is the user-friendly interface used for a DNA file reader interface.
    void Run(string file_name); //this method is going to start a menu and close it once the user is done based on the user's choice
  private:
    string m_choice;
    string m_file_name;
    Logger m_logger;
};
