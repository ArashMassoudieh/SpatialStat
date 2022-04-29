#include "environment.h"
#include "interface.h"
#include "grid.h"
#include "Distribution.h"

Environment::Environment()
{

}


bool Environment::Execute(const Command &cmd)
{
    if (Grid::HasCommand(cmd.command))
    {
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::creator)
        {
            Objects[cmd.object_name]=new Grid();
            dynamic_cast<Grid*>(Objects[cmd.object_name])->CreateGrid(cmd.arguments);
            return true;
        }
    }
    if (CDistribution::HasCommand(cmd.command))
    {
        if (cmd.Command_Structures[cmd.command].CommandType == command_type::creator)
        {
            Objects[cmd.object_name]=new Grid();
            dynamic_cast<Grid*>(Objects[cmd.object_name])->CreateGrid(cmd.arguments);
            return true;
        }
    }
}
