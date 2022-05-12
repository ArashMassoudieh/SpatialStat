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

    vector<string> argument_list;
    string the_rest;
    if (aquiutils::contains(cmdline,"("))
    {
        argument_list = aquiutils::split(aquiutils::split(aquiutils::split(cmdline,'(')[1],')')[0],',');

        for (unsigned int i=0; i<argument_list.size(); i++)
        {
            vector<string> arg = aquiutils::split(argument_list[i],'=');
            if (arg.size()==2)
                arguments[arg[0]] = arg[1];
        }
        the_rest = aquiutils::split(cmdline,'(')[0];
    }
    string rhs;
    if (aquiutils::contains(the_rest,"="))
    {
        object_name=aquiutils::split(the_rest,'=')[0];
        rhs = aquiutils::split(the_rest,'=')[1];
    }
    else
    {
        rhs = the_rest;
    }
    if (aquiutils::contains(rhs,"."))
    {
        object_name = aquiutils::split(rhs,'.')[0];
        command = aquiutils::split(rhs,'.')[1];
    }
    else
    {
        command = aquiutils::split(rhs,'.')[0];
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
    Command_Structures["GetConcentrationBTCAtX"].CommandType = command_type::creator; Command_Structures["GetConcentrationBTCAtX"].Object = "Grid"; Command_Structures["GetConcentrationBTCAtX"].Output = "timeseries";


    Command_Structures_Initialized = true;
}

