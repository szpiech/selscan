#include "param_t.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

bool param_t::addFlag(string flag, bool value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        if (value) buffer = "true";
        else buffer = "false";
        argb[flag] = value;
        help[flag] = "<bool>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, double value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%.2f", value);
        buffer = charBuffer;
        argd[flag] = value;
        help[flag] = "<double>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%d", value);
        buffer = charBuffer;
        argi[flag] = value;
        help[flag] = "<int>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, char value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%c", value);
        buffer = charBuffer;
        argch[flag] = value;
        help[flag] = "<char>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        args[flag] = value;
        help[flag] = "<string>: " + description + "\n\tDefault: " + value;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addFlag(string flag, const char value[], string label, string description)
{
    return this->addFlag(flag, string(value), label, description);
}

bool param_t::addListFlag(string flag, int value, string label, string description)
{
    if (!flagExists(flag))
    {
        string buffer;
        char charBuffer[100];
        sprintf(charBuffer, "%d", value);
        buffer = charBuffer;
        listargi[flag].push_back(value);
        help[flag] = "<int1> ... <intN>: " + description + "\n\tDefault: " + buffer;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool  param_t::addListFlag(string flag, string value, string label, string description)
{
    if (!flagExists(flag))
    {
        listargs[flag].push_back(value);
        help[flag] = "<string1> ... <stringN>: " + description + "\n\tDefault: " + value;
        labels[flag] = label;
    }
    else
    {
        cerr << "ERROR: " << flag << " already exists.\n";
        throw 0;
    }

    return true;
}

bool param_t::addListFlag(string flag, const char value[], string label, string description)
{
    return this->addListFlag(flag, string(value), label, description);
}

void param_t::printHelp()
{
    map<string, string>::iterator it;

    cerr << preamble << endl;

    cerr << "----------Command Line Arguments----------\n\n";

    for (it = help.begin(); it != help.end(); it++)
    {
        if (labels[(*it).first].compare("SILENT") != 0)
        {
            cerr << (*it).first << " " << (*it).second << "\n\n";
        }
    }

    return;
}

bool param_t::goodDouble(string str)
{
    string::iterator it;
    //int dashCount = 0;
    int decimalCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '.' && *it != '-') return 0;
        if (*it == '.') decimalCount++;
        if (*it == '-' && it != str.begin()) return 0;
        if (/*dashCount > 1 || */decimalCount > 1) return 0;
    }
    return 1;
}

bool param_t::goodInt(string str)
{
    string::iterator it;
    //int dashCount = 0;
    for (it = str.begin(); it != str.end(); it++)
    {
        if (!isdigit(*it) && *it != '-') return 0;
        if (*it == '-' && it != str.begin()) return 0;
        //if (dashCount > 1) return 0;
    }
    return 1;
}

bool param_t::goodChar(string str)
{
    if (str.length() > 1) return 0;
    return 1;
}

bool param_t::parseCommandLine(int argc, char *argv[])
{
    int badFlags = 0;

    for (int i = 1; i < argc; i++)
    {
        if (isSet.count(argv[i]) > 0)
        {
            cerr << "ERROR: Duplicate " << argv[i] << " found.\n";
            badFlags++;
            break;
        }
        else if (argb.count(argv[i]) > 0)
        {
            argb[argv[i]] = !argb[argv[i]];
            isSet[argv[i]] = true;
        }
        else if (argi.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodInt(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                badFlags++;
                break;
            }
            else
            {
                argi[argv[i]] = atoi(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargi.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargi[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (goodInt(string(argv[i + 1]))) //make sure the next value is OK
                    {
                        listargi[argv[flagIndex]].push_back(atoi(argv[i + 1]));
                        i++;
                    }
                    //if it is a bad int...
                    else if (!goodInt(string(argv[i + 1])) && !flagExists(string(argv[i + 1]))) 
                    {
                        cerr << "ERROR: " << argv[i + 1] << " is not a valid integer.\n";
                        badFlags++;
                        break;
                    }
                    else //if the next value is another flag..
                    {
                        if (listargi[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argd.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodDouble(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid double.\n";
                badFlags++;
                break;
            }
            else
            {
                argd[argv[i]] = atof(argv[i + 1]);
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (args.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                args[argv[i]] = argv[i + 1];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else if (listargs.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else
            {
                listargs[argv[i]].clear();//clear the default value
                int flagIndex = i;//remember where the flag is in argv
                while (i + 1 < argc) //go until the end of the list
                {
                    if (argv[i + 1][0] != '-') //make sure the next value isn't another flag
                    {
                        listargs[argv[flagIndex]].push_back(string(argv[i + 1]));
                        i++;
                    }
                    else //if the next value is another flag...
                    {
                        if (listargs[argv[flagIndex]].size() == 0)
                        {
                            cerr << "ERROR: No arguments found for " << argv[flagIndex] << ".\n";
                            badFlags++;
                        }
                        break;
                    }
                }
            }
        }
        else if (argch.count(argv[i]) > 0)
        {
            if (i + 1 >= argc)
            {
                cerr << "ERROR: No argument found for " << argv[i] << ".\n";
                badFlags++;
                break;
            }
            else if (!goodChar(string(argv[i + 1])))
            {
                cerr << "ERROR: " << argv[i + 1] << " is not a valid character.\n";
                badFlags++;
                break;
            }
            else
            {
                argch[argv[i]] = argv[i + 1][0];
                isSet[argv[i]] = true;
                i++;
            }
        }
        else //if (argv[i][0] == '-')
        {
            cerr << "ERROR: command line flag " << argv[i] << " not recognized.\n";
            badFlags++;
        }
    }

    if (getBoolFlag(ARG_HELP))
    {
        this->printHelp();
        throw 0;
    }

    if (badFlags) throw 0;

    return 0;
}

bool param_t::flagExists(string flag)
{
    return (help.count(flag) > 0);
}

param_t::param_t()
{
    this->addFlag(ARG_HELP, false, "__help", "Prints this help dialog.");
}

bool param_t::getBoolFlag(string flag)
{
    if (argb.count(flag) > 0) return argb[flag];

    cerr << "ERROR: There are no bool flags named " << flag << "\n";
    throw 0;
}

double param_t::getDoubleFlag(string flag)
{
    if (argd.count(flag) > 0) return argd[flag];

    cerr << "ERROR: There are no double flags named " << flag << "\n";
    throw 0;
}

int param_t::getIntFlag(string flag)
{
    if (argi.count(flag) > 0) return argi[flag];

    cerr << "ERROR: There are no int flags named " << flag << "\n";
    throw 0;
}

char param_t::getCharFlag(string flag)
{
    if (argch.count(flag) > 0) return argch[flag];

    cerr << "ERROR: There are no char flags named " << flag << "\n";
    throw 0;
}

string param_t::getStringFlag(string flag)
{
    if (args.count(flag) > 0) return args[flag];

    cerr << "ERROR: There are no string flags named " << flag << "\n";
    throw 0;
}

vector<string> param_t::getStringListFlag(string flag)
{
    if (listargs.count(flag) > 0) return listargs[flag];

    cerr << "ERROR: There are no string list flags named " << flag << "\n";
    throw 0;
}

vector<int> param_t::getIntListFlag(string flag)
{
    if (listargi.count(flag) > 0) return listargi[flag];

    cerr << "ERROR: There are no int list flags named " << flag << "\n";
    throw 0;
}

void param_t::setPreamble(string str)
{
    preamble = str;
    return;
}
