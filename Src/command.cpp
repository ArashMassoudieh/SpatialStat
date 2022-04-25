#include "command.h"
#include "vector"
#include "Utilities.h"

Command::Command()
{

}

Command::Command(const string &cmdline)
{
    vector<string> s = aquiutils::split(cmdline,' ');
    if (s.size()>=2)
    {
        command = s[0];
        object_name = s[1];
        for (unsigned int i=2; i<s.size(); i++)
        {
            vector<string> arg = aquiutils::split(s[i],'=');
            if (arg.size()==2)
                arguments[arg[0]] = arg[1];
        }
    }
}

