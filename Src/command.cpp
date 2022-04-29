#include "command.h"
#include "vector"
#include "Utilities.h"

vector<char> Command::deliminators = vector<char>({'=', '(', ')', ','});
map<string,command_parameters> Command::Command_Structures = map<string,command_parameters>();

Command::Command()
{
    if (!Command_Structures_Initialized)
        Initialize_Command_Structure();
}

Command::Command(const string &cmdline)
{
    if (!Command_Structures_Initialized)
        Initialize_Command_Structure();


    vector<string> s = aquiutils::split(cmdline,deliminators);
    if (s.size()>=2)
    {
        command = s[1];
        if (Command_Structures.count(command)==1)
        {   if (Command_Structures[command].CommandType == command_type::creator)
            {   object_name = s[0];
                for (unsigned int i=2; i<s.size(); i++)
                {
                    vector<string> arg = aquiutils::split(s[i],':');
                    if (arg.size()==2)
                        arguments[arg[0]] = arg[1];
                }
            }
            else if (Command_Structures[command].CommandType == command_type::modifier)
            {
                object_name = s[0];
                for (unsigned int i=2; i<s.size(); i++)
                {
                    vector<string> arg = aquiutils::split(s[i],':');
                    if (arg.size()==2)
                        arguments[arg[0]] = arg[1];
                }
            }
        }
    }
}

void Command::Initialize_Command_Structure()
{
    command_parameters cmd;
    Command_Structures["CreateGrid"].CommandType = command_type::creator; Command_Structures["CreateGrid"].Object = "Grid"; Command_Structures["CreateGrid"].Output = "Grid";
    Command_Structures["CreateDistribution"].CommandType = command_type::creator; Command_Structures["CreateDistribution"].Object = "Distribution"; Command_Structures["CreateDistribution"].Output = "Distribution";
    Command_Structures_Initialized = true;
}

