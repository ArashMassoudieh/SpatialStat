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
    Command operator = (const Command &C);
    Command(const Command &C);
    ~Command();
    map<string,string> arguments;
    string command;
    string object_name;
    static vector<char> deliminators;
    static map<string,command_parameters> Command_Structures;
    bool Command_Structures_Initialized = false;
    void Initialize_Command_Structure();

};

#endif // COMMAND_H
