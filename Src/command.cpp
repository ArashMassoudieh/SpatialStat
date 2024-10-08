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
        output_name=aquiutils::split(the_rest,'=')[0];
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
    if (object_name=="") object_name = output_name;

}

Command Command::operator = (const Command &C)
{
    output_name = C.output_name;
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
    output_name = C.output_name;
    arguments = C.arguments;
    command = C.command;
    object_name = C.object_name;
}

void Command::Initialize_Command_Structure()
{
    Command_Structures["CreateGrid"].Object = object_type::grid; Command_Structures["CreateGrid"].Output = object_type::grid;
    Command_Structures["CreateDistribution"].Object = object_type::distribution; Command_Structures["CreateDistribution"].Output = object_type::distribution;
    Command_Structures["WriteDistToFile"].Object = object_type::distribution; Command_Structures["WriteDistToFile"].Output = object_type::none;
    Command_Structures["WriteTimeSeriesToFile"].Object = object_type::distribution; Command_Structures["WriteTimeSeriesToFile"].Output = object_type::none;
    Command_Structures["SetInverseCumulative"].Object = object_type::distribution; Command_Structures["WriteToFile"].Output = object_type::none;
    Command_Structures["WriteInverseCumulativeToFile"].Object = object_type::distribution; Command_Structures["WriteInverseCumulativeToFile"].Output = object_type::none;
    Command_Structures["AssignKField"].Object = object_type::grid; Command_Structures["AssignKField"].Output = object_type::none;
    Command_Structures["WriteKFieldToVTP"].Object = object_type::grid; Command_Structures["WriteKFieldToVTP"].Output = object_type::none;
    Command_Structures["WritePathwayToVTP"].Object = object_type::pathwayset; Command_Structures["WritePathwayToVTP"].Output = object_type::none;
    Command_Structures["RenormalizeKField"].Object = object_type::grid; Command_Structures["RenormalizeKField"].Output = object_type::none;
    Command_Structures["SolveHydro"].Object = object_type::grid; Command_Structures["SolveHydro"].Output = object_type::none;
    Command_Structures["WriteHydroSolutionToVTP"].Object = object_type::grid; Command_Structures["WriteHydroSolutionToVTP"].Output = object_type::none;
    Command_Structures["SolveTransport"].Object = object_type::grid; Command_Structures["SolveTransport"].Output = object_type::none;
    Command_Structures["WriteConcentrationToVTP"].Object = object_type::grid; Command_Structures["WriteConcentrationToVTP"].Output = object_type::none;
    Command_Structures["GetConcentrationBTCAtX"].Object = object_type::grid; Command_Structures["GetConcentrationBTCAtX"].Output = object_type::timeseries;
    Command_Structures["CreateTrajectories"].Object = object_type::grid; Command_Structures["CreateTrajectories"].Output = object_type::pathwayset;
    Command_Structures["Uniformize"].Object = object_type::pathwayset; Command_Structures["Uniformize"].Output = object_type::none;
    Command_Structures["WriteFieldToVTP"].Object = object_type::grid; Command_Structures["WriteFieldToVTP"].Output = object_type::none;
    Command_Structures["GetMarginalVelocityDistribution"].Object = object_type::grid; Command_Structures["GetMarginalVelocityDistribution"].Output = object_type::timeseries;
    Command_Structures["WriteTimeSeriesToFile"].Object = object_type::timeseries; Command_Structures["WriteTimeSeriesToFile"].Output = object_type::none;
    Command_Structures["GetVelocityDistributionAtXSection"].Object = object_type::timeseries; Command_Structures["GetVelocityDistributionAtXSection"].Output = object_type::none;
    Command_Structures["AssignKField_2nd_order"].Object = object_type::grid; Command_Structures["AssignKField_2nd_order"].Output = object_type::grid;
    Command_Structures["GetMarginalDistribution"].Object = object_type::grid; Command_Structures["GetMarginalDistribution"].Output = object_type::distribution;

    Command_Structures_Initialized = true;
}

