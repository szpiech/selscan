#ifndef __PARAM_T_H__
#define __PARAM_T_H__

#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <cctype>

using namespace std;

const string ARG_HELP = "--help";

class param_t
{
public:

    bool addFlag(string flag, bool value, string label, string description);
    bool addFlag(string flag, double value, string label, string description);
    bool addFlag(string flag, int value, string label, string description);
    bool addFlag(string flag, char value, string label, string description);
    bool addFlag(string flag, string value, string label, string description);
    bool addFlag(string flag, const char value[], string label, string description);

    bool addListFlag(string flag, string value, string label, string description);
    bool addListFlag(string flag, const char value[], string label, string description);
    bool addListFlag(string flag, int value, string label, string description);


    void printHelp();

    bool parseCommandLine(int argc, char *argv[]);

    bool getBoolFlag(string flag);
    double getDoubleFlag(string flag);
    int getIntFlag(string flag);
    char getCharFlag(string flag);
    string getStringFlag(string flag);

    vector<string> getStringListFlag(string flag);
    vector<int> getIntListFlag(string flag);

    void setPreamble(string str);

    param_t();


private:

    map<string, bool> argb;
    map<string, double> argd;
    map<string, int> argi;
    map<string, char> argch;
    map<string, string> args;

    map<string, vector< string > > listargs;
    map<string, vector< int > > listargi;    

    map<string, string> help;
    map<string, bool> isSet;

    map<string, string> labels;

    bool goodDouble(string str);
    bool goodInt(string str);
    bool goodChar(string str);

    bool flagExists(string flag);

    string preamble;
};

#endif
