#ifndef COMMAND_H
#define COMMAND_H

#include <string>
#include <map>

using namespace std;

class Command
{
public:
    Command();
    Command(const string &s);
    map<string,string> arguments;
    string command;
    string object_name;

};

#endif // COMMAND_H
