#include "command.h"
#include "vector"
#include "Utilities.h"

vector<char> Command::deliminators = vector<char>({'=', '(', ')', ',','*'});
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

Command Command::operator = (const Command &C)
{
    arguments = C.arguments;
    command = C.command;
    object_name = C.object_name;
    return *this;
}


Command::~Command(void)
{

}

Command::Command(const Command &C)
{
    arguments = C.arguments;
    command = C.command;
    object_name = C.object_name;
}

void Command::Initialize_Command_Structure()
{
    command_parameters cmd;
    Command_Structures["CreateGrid"].CommandType = command_type::creator; Command_Structures["CreateGrid"].Object = "Grid"; Command_Structures["CreateGrid"].Output = "Grid";
    Command_Structures["CreateDistribution"].CommandType = command_type::creator; Command_Structures["CreateDistribution"].Object = "Distribution"; Command_Structures["CreateDistribution"].Output = "Distribution";
    Command_Structures["WriteToFile"].CommandType = command_type::modifier; Command_Structures["WriteToFile"].Object = "Distribution"; Command_Structures["WriteToFile"].Output = "";
    Command_Structures["SetInverseCumulative"].CommandType = command_type::modifier; Command_Structures["SetInverseCumulative"].Object = "Distribution"; Command_Structures["WriteToFile"].Output = "";
    Command_Structures["WriteInverseCumulativeToFile"].CommandType = command_type::modifier; Command_Structures["WriteInverseCumulativeToFile"].Object = "Distribution"; Command_Structures["WriteInverseCumulativeToFile"].Output = "";
    Command_Structures["AssignKField"].CommandType = command_type::modifier; Command_Structures["AssignKField"].Object = "Grid"; Command_Structures["AssignKField"].Output = "";
    Command_Structures["WriteKFieldToVTP"].CommandType = command_type::modifier; Command_Structures["WriteKFieldToVTP"].Object = "Grid"; Command_Structures["WriteKFieldToVTP"].Output = "";
    Command_Structures["RenormalizeKField"].CommandType = command_type::modifier; Command_Structures["RenormalizeKField"].Object = "Grid"; Command_Structures["RenormalizeKField"].Output = "";
    Command_Structures["SolveHydro"].CommandType = command_type::modifier; Command_Structures["SolveHydro"].Object = "Grid"; Command_Structures["SolveHydro"].Output = "";
    Command_Structures["WriteHydroSolutionToVTP"].CommandType = command_type::modifier; Command_Structures["WriteHydroSolutionToVTP"].Object = "Grid"; Command_Structures["WriteHydroSolutionToVTP"].Output = "";
    Command_Structures["SolveTransport"].CommandType = command_type::modifier; Command_Structures["SolveTransport"].Object = "Grid"; Command_Structures["SolveTransport"].Output = "";
    Command_Structures["WriteConcentrationToVTP"].CommandType = command_type::modifier; Command_Structures["WriteConcentrationToVTP"].Object = "Grid"; Command_Structures["WriteConcentrationToVTP"].Output = "";


    Command_Structures_Initialized = true;
}

