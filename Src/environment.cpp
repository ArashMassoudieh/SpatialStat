#include "environment.h"
#include "interface.h"
#include "grid.h"

Environment::Environment()
{

}


bool Environment::Execute(const Command &cmd)
{
    if (Grid::HasCommand(cmd.command))
    {
        if (cmd.CommadType == command_type::creator)
        {

        }
    }
}
