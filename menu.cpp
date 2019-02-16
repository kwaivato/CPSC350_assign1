// menu.cpp
// An implemented menu interface used to facilitate a DNA reader program.
// @author Vidal M. Arroyo
// @author arroy118@mail.chapman.edu
// @version 1.0

#include "menu.h"

using namespace std;

// The m_choice variable "Y" is defaulted to ensure that we don't have junk hanging out in the string.
Menu::Menu()
{
  m_choice = "Y";
  Logger m_logger;
}

Menu::~Menu()
{

}

// This run method will act as the public interface for the user. Comments will be appended where deemed necessary.
void Menu::Run(string file_name)
{
  cout << "Welcome to MyGenePool: the only computer program that can compute statistics for your gene pools with style! :-)" << endl;
  cout << "My name is Gene, and I am your virtual assistant for MyGenePool." << endl;
  cout << "It looks like you want to analyze '" << file_name << "'. I will do that now." << endl;
  m_file_name = file_name;
  Logger m_logger;
  m_logger.setGeneCalculator();
  m_logger.setGenePool();
  // The Out method is used here since this is the first file.
  m_logger.Out(m_file_name);
  cout << "Do you want to continue? Type 'Y' to continue and anything else to exit: ";
  // A cin statement must be used for extra files since they can no longer be entered via the command line.
  cin >> m_choice;
  // This while loop is run as long as the uder wants to enter files to be analyzed
  while (m_choice == "Y")
  {
    cout << "What is the name of your .txt file that you would like MyGenePool to analyze?: ";
    cin >> m_file_name;
    Logger m_logger;
    m_logger.setGeneCalculator();
    m_logger.setGenePool();
    m_logger.AppendOut(m_file_name);
    cout << "Do you want to continue? Type 'Y' to continue and anything else to exit: ";
    cin >> m_choice;
  }
  cout << "Thank you for using MyGenePool, the only computer program that can compute statistics for your gene pools with style! :-)" << endl;
}
