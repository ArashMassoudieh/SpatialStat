#ifndef COMMAND_H
#define COMMAND_H

#include <string>
#include <map>
#include <vector>

using namespace std;

enum class command_type {creator, modifier};

struct command_parameters
{
    string Object;
    string Output;
    command_type CommandType;
};

class Command
{
public:
    Command();
    Command(const string &s);
    map<string,string> arguments;
    string command;
    string object_name;
    static vector<char> deliminators;
    command_type CommadType;
    static map<string,command_parameters> Commands_Structure;

};

#endif // COMMAND_H
